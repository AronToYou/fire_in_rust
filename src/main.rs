use indicatif::ProgressIterator;
use std::ops::{Add, Sub, Mul};

// ---------------- Utility Functions ----------------
#[inline] fn clamp(v: f32, lb: f32, ub: f32) -> f32 { v.max(lb).min(ub) }

fn new_sim(p: Params) -> Sim<impl Fn((f32, f32)) -> (f32, f32)> {
    let maxx = p.nx as f32 - 1.001;
    let maxy = p.ny as f32 - 1.001;
    let clamp = move |(x, y): (f32, f32)| (x.clamp(0.0, maxx), y.clamp(0.0, maxy));
    Sim::new(p, clamp)
}

#[derive(Clone, Copy)]
struct P(f32, f32);

impl Add for P {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        P(self.0 + rhs.0, self.1 + rhs.1)
    }
}

impl Sub for P {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        P(self.0 - rhs.0, self.1 - rhs.1)
    }
}

impl Mul<f32> for P {
    type Output = Self;
    fn mul(self, rhs: f32) -> Self {
        P(self.0 * rhs, self.1 * rhs)
    }
}

// ---------------- Simulation Parameters ----------------
#[derive(Clone, Copy)]
struct Params {
    // Grid Parameters //
    nx: usize, ny: usize,  // grid lengths [number of cells]
    h: f32,     // grid cell length     [m]
    dt: f32,    // simulation time step [s]

    // Physical Constants //
    visc: f32,      // kinematic viscosity [m^2/s] (for velocity diffusion)
    k_diff: f32,    // diffusion constant  [m^2/s] (for scalar diffusion)
    k_buoy: f32,    // buoyancy constant   [N/K]   (for buoyancy force)
    temp_air: f32,  // ambient air temperature [K]
    temp_max: f32,  // peak flame temperature  [K]
    k_cool: f32,    // cooling constant        [K/s]
    vconf: f32,     // vorticity confinement   [1/s]
    k_react: f32,   // reaction rate      [m/s]
    d_fuel: f32,    // denisty of fuel    [kg/m^3]
    d_hgas: f32,    // density of hot gas [kg/m^3]

    // scaling factor for rendering //
    render_scale: f32
}

// ---------------- Simulation State ----------------
struct Sim<C> where C: Fn((f32, f32)) -> (f32, f32) {
    p: Params,  // Simulation parameters (defined above)
    clamp: C,

    ux: Vec<f32>, uy: Vec<f32>,  // velocity field (x, y)
    p_h: Vec<f32>,    // pressure for hot gas
    div_h: Vec<f32>,  // divergence (hot gas)

    phi: Vec<f32>,   // level set (+pos in fuel region, -neg outside, 0 at boundary)
    temp: Vec<f32>,  // temperature field (hot gas domain)
    rt: Vec<f32>,    // reaction-time tracker (1 at fuel; decreases after crossing)
    dns: Vec<f32>,   // smoke density (simple)

    tmp: Vec<f32>, tmp2: Vec<f32>  // temporary intermediate fields for calculations
}

/// Which field from `Sim` to print
enum Field { 
    Ux, Uy,
    Ph, Divh,
    Phi, Temp, Rt, Dns
}

impl<C> Sim<C> where C: Fn((f32, f32)) -> (f32, f32) {
    fn new(p: Params, clamp: C) -> Self {
        let n = p.nx*p.ny;
        let mut s = Self {
            p, clamp,
            ux: vec![0.0; n], uy: vec![0.0; n],
            p_h: vec![0.0; n],
            div_h: vec![0.0; n],

            phi: vec![1.0; n],
            temp: vec![p.temp_air; n],
            rt: vec![0.0; n],
            dns: vec![0.0; n],

            tmp: vec![0.0; n], tmp2: vec![0.0; n]
        };
        s.init_fuel_inlet();
        s
    }

    /// Initialize fuel inlet at bottom of domain
    fn init_fuel_inlet(&mut self) {
        let (nx, ny) = (self.p.nx, self.p.ny);
        let mut idx: usize = 0;
        for y in 0..ny {
            for _ in 0..nx {
                if y < 10 {
                    self.uy[idx] = 1.5;
                    self.phi[idx] = 5.0;
                } else {
                    self.phi[idx] = -5.0;
                }
                idx += 1;
            }
        }
    }

    fn step(&mut self) {
        // (1) Update level set field //
        self.update_levelset();

        // (2) Semi-Lagrangian advection of velocity fields //
        self.semi_lagrangian_advect();
    }

    // ----------------- (1) Thin-flame Level Set Propagation -----------------
    /// Updates the level set using upwind one-sided differencing to estimate spatial derivatives
    fn update_levelset(&mut self) {
        let (nx, ny, k_react) = (self.p.nx, self.p.ny, self.p.k_react);
        let mut idx = nx;
        for _ in 1..ny-1 {
            idx += 1;
            for _ in 1..nx-1 {
                // (1.3) (unscaled) Central differencing for normed gradient (∇φ/|∇φ|) //
                let gx = self.phi[idx+1] - self.phi[idx-1];  // gradient x-component
                let gy = self.phi[idx+nx] - self.phi[idx-nx];  // gradient y-component
                let norm = (gx*gx + gy*gy).sqrt().max(1e-8);  // gradient norm
                
                // (1.2) Velocity of implicit surface (where φ==0) //
                let wx = self.ux[idx] + k_react*gx/norm;
                let wy = self.uy[idx] + k_react*gy/norm;

                // (unscaled) Upwind one-sided differencing //
                let ddx = if wx > 0.0 {
                    self.phi[idx] - self.phi[idx-1]
                } else {
                    self.phi[idx+1] - self.phi[idx]
                };
                let ddy = if wy > 0.0 {
                    self.phi[idx] - self.phi[idx-nx]
                } else {
                    self.phi[idx+nx] - self.phi[idx]
                };

                // (1.4) (scaled) Application of time derivative //
                self.tmp[idx] = self.phi[idx] - (wx*ddx + wy*ddy)*(self.p.dt/self.p.h);

                idx += 1;
            }
            idx += 1;
        }
        std::mem::swap(&mut self.phi, &mut self.tmp);
    }

    // ------------- (2) Semi-Lagrangian advection of velocity fields -------------
    /// Runge-Kutta 2-stage backtrace, bilinear velocity sampling, clamped at boundaries
    fn semi_lagrangian_advect(&mut self) {
        let (ux, uy) = (&self.ux, &self.uy);
        let (nx, ny, dt, h) = (self.p.nx, self.p.ny, self.p.dt, self.p.h);
        for y in 0..ny {
            for x in 0..nx {
                let idx = x + y*nx;
                let P(x1, y1) = P(x as f32, y as f32);

                // backtrace half-step to midpoint //
                let P(x0, y0) = P(x1, y1) - P(ux[idx], uy[idx])*(0.5*dt/h);
                // midpoint velocity //
                let ux0 = self.sample_bilin(ux, x0, y0);
                let uy0 = self.sample_bilin(uy, x0, y0);

                // backtrace full step //
                let P(x0, y0) = P(x1, y1) - P(ux0, uy0)*(dt/h);
                // final velocity //
                self.tmp[idx] = self.sample_bilin(ux, x0, y0);
                self.tmp2[idx] = self.sample_bilin(uy, x0, y0);
            }
        }
        std::mem::swap(&mut self.ux, &mut self.tmp);
        std::mem::swap(&mut self.uy, &mut self.tmp2);
    }

    // ----------------- Utility Functions -----------------
    /// Bilinear sampling of scalar field, clamped at boundaries
    fn sample_bilin(&self, f: &Vec<f32>, x: f32, y: f32) -> f32 {
        let (nx, ny) = (self.p.nx, self.p.ny);
        let x = clamp(x, 0.0, (nx as f32) - 1.001);
        let y = clamp(y, 0.0, (ny as f32) - 1.001);  
        let x0 = x.floor() as usize; let y0 = y.floor() as usize;
        let x1 = (x0+1).min(nx-1);   let y1 = (y0+1).min(ny-1);
        let tx = x - x0 as f32;      let ty = y - y0 as f32;
        
        let f00 = f[x0 + y0*nx];       let f01 = f[x0 + y1*nx];
        let f10 = f[x1 + y0*nx];       let f11 = f[x1 + y1*nx];
        let a = f00*(1.0-tx) + f10*tx; let b = f01*(1.0-tx) + f11*tx;
        a*(1.0-ty) + b*ty
    }

    /// Print a downsampled version of a field to the console
    fn print_field(&self, which: Field) {
        let (label, field): (&str, &[f32]) = match which {
            Field::Ux => ("Velocity x Component", &self.ux),
            Field::Uy => ("Velocity y Component", &self.uy),
            Field::Ph => ("Hot Gas Pressure", &self.p_h),
            Field::Divh => ("Divergence of Hot Gas Pressure", &self.div_h),
            Field::Phi => ("Level Set", &self.phi),
            Field::Temp => ("Temperature", &self.temp),
            Field::Rt => ("Reaction Parameter", &self.rt),
            Field::Dns => ("Smoke Density", &self.dns)
        };
        let (nx, ny) = (self.p.nx, self.p.ny);
        let (stride_x, stride_y) = (nx/30, ny/30);
        println!("{label} field ({nx}x{ny}), sampled with stride: ({stride_x}, {stride_y})");
        for y in (0..ny).rev().step_by(stride_y) {
            for x in (0..nx).step_by(stride_x) {
                print!("{:4.1} ", field[x + y*nx]);
            }
            println!();
        }
    }
}

// ---------------- Main ----------------
fn main() {
    let mut sim = new_sim(Params {
        nx: 240, ny: 320,
        h: 1.0,
        dt: 0.5,

        visc: 0.0005,
        k_diff: 0.0001,
        k_buoy: 0.05,
        temp_air: 300.0,
        temp_max: 2200.0,
        k_cool: 0.002,
        vconf: 2.0,
        k_react: 0.7,
        d_fuel: 1.0,
        d_hgas: 0.1,

        render_scale: 1.0
    });

    // Prints the initialized level set field
    sim.print_field(Field::Phi);
    
    let steps = 300;
    println!("Simulating {steps} steps...");
    for _ in (0..steps).progress() {
        sim.step();
    }
    sim.print_field(Field::Phi);
    println!("Done!");
}
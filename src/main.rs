use indicatif::ProgressIterator;
use std::ops::{Add, Sub, Mul};

// ---------------- Utility Functions ----------------
fn new_sim<const NX: usize, const NY: usize>(p: Params) -> Sim<NX, NY, impl Fn((f32, f32)) -> (f32, f32)> {
    let maxx = (NX as f32) - 1.001;
    let maxy = (NY as f32) - 1.001;
    let clamp_xy = move |(x, y): (f32, f32)| (x.clamp(0.0, maxx), y.clamp(0.0, maxy));
    Sim::<NX, NY, _>::new(p, clamp_xy)
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
struct Sim<const NX: usize, const NY: usize, C> where C: Fn((f32, f32)) -> (f32, f32) {
    p: Params,  // Simulation parameters (defined above)
    clamp_xy: C,  // clamping function for coordinates

    ux: [[f32; NY]; NX], uy: [[f32; NY]; NX],  // velocity field (x, y)
    p_h: [[f32; NY]; NX],    // pressure for hot gas
    div_h: [[f32; NY]; NX],  // divergence (hot gas)

    phi: [[f32; NY]; NX],   // level set (+pos in fuel region, -neg outside, 0 at boundary)
    temp: [[f32; NY]; NX],  // temperature field (hot gas domain)
    rt: [[f32; NY]; NX],    // reaction-time tracker (1 at fuel; decreases after crossing)
    dns: [[f32; NY]; NX],   // smoke density (simple)

    tmp: [[f32; NY]; NX], tmp2: [[f32; NY]; NX]  // temporary intermediate fields for calculations
}

/// Which field from `Sim` to print
enum Field { 
    U,
    Ph, Divh,
    Phi, Temp, Rt, Dns
}

impl<const NX: usize, const NY: usize, C> Sim<NX, NY, C> where C: Fn((f32, f32)) -> (f32, f32) {
    fn new(p: Params, clamp_xy: C) -> Self {
        let mut s = Self {
            p, clamp_xy,
            ux: [[0.0; NY]; NX], uy: [[0.0; NY]; NX],
            p_h: [[0.0; NY]; NX],
            div_h: [[0.0; NY]; NX],

            phi: [[1.0; NY]; NX],
            temp: [[p.temp_air; NY]; NX],
            rt: [[0.0; NY]; NX],
            dns: [[0.0; NY]; NX],

            tmp: [[0.0; NY]; NX], tmp2: [[0.0; NY]; NX]
        };
        s.init_fuel_inlet();
        s
    }

    /// Initialize fuel inlet at bottom of domain
    fn init_fuel_inlet(&mut self) {
        for x in 0..NX {
            for y in 0..NY {
                if y < 10 {
                    self.uy[x][y] = 1.5;
                    self.phi[x][y] = 5.0;
                } else {
                    self.phi[x][y] = -5.0;
                }
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
        let k_react = self.p.k_react;
        for x in 1..NX-1 {
            for y in 1..NY-1 {
                // (1.3) (unscaled) Central differencing for normed gradient (∇φ/|∇φ|) //
                let gx = self.phi[x+1][y] - self.phi[x-1][y];  // gradient x-component
                let gy = self.phi[x][y+1] - self.phi[x][y-1];  // gradient y-component
                let norm = (gx*gx + gy*gy).sqrt().max(1e-8);  // gradient norm
                
                // (1.2) Velocity of implicit surface (where φ==0) //
                let wx = self.ux[x][y] + k_react*gx/norm;
                let wy = self.uy[x][y] + k_react*gy/norm;

                // (unscaled) Upwind one-sided differencing //
                let ddx = if wx > 0.0 {
                    self.phi[x][y] - self.phi[x-1][y]
                } else {
                    self.phi[x+1][y] - self.phi[x][y]
                };
                let ddy = if wy > 0.0 {
                    self.phi[x][y] - self.phi[x][y-1]
                } else {
                    self.phi[x][y+1] - self.phi[x][y]
                };

                // (1.4) (scaled) Application of time derivative //
                self.tmp[x][y] = self.phi[x][y] - (wx*ddx + wy*ddy)*(self.p.dt/self.p.h);
            }
        }
        std::mem::swap(&mut self.phi, &mut self.tmp);
    }

    // ------------- (2) Semi-Lagrangian advection of velocity fields -------------
    /// Runge-Kutta 2-stage backtrace, bilinear velocity sampling, clamped at boundaries
    fn semi_lagrangian_advect(&mut self) {
        let (ux, uy) = (&self.ux, &self.uy);
        let (dt, h) = (self.p.dt, self.p.h);
        for x in 0..NX {
            for y in 0..NY {
                let P(x1, y1) = P(x as f32, y as f32);

                // backtrace half-step to midpoint //
                let P(x0, y0) = P(x1, y1) - P(ux[x][y], uy[x][y])*(0.5*dt/h);
                // midpoint velocity //
                let ux0 = self.sample_bilin(ux, x0, y0);
                let uy0 = self.sample_bilin(uy, x0, y0);

                // backtrace full step //
                let P(x0, y0) = P(x1, y1) - P(ux0, uy0)*(dt/h);
                // final velocity //
                self.tmp[x][y] = self.sample_bilin(ux, x0, y0);
                self.tmp2[x][y] = self.sample_bilin(uy, x0, y0);
            }
        }
        std::mem::swap(&mut self.ux, &mut self.tmp);
        std::mem::swap(&mut self.uy, &mut self.tmp2);
    }

    // ----------------- Utility Functions -----------------
    /// Bilinear sampling of scalar field, clamped at boundaries
    fn sample_bilin(&self, f: &[[f32; NY]; NX], x: &f32, y: &f32) -> f32 {
        let (x, y) = (self.clamp_xy)((x, y));
        let x0 = x.floor() as usize; let y0 = y.floor() as usize;
        let x1 = (x0+1).min(NX-1);   let y1 = (y0+1).min(NY-1);
        let tx = x - x0 as f32;      let ty = y - y0 as f32;
        
        let f00 = f[x0][y0];       let f01 = f[x0][y1];
        let f10 = f[x1][y0];       let f11 = f[x1][y1];  // store x1 in x0
        let a = f00*(1.0-tx) + f10*tx; let b = f01*(1.0-tx) + f11*tx;  // store tx in x1/y1 in x0/
        a*(1.0-ty) + b*ty  // store ty in tx in x1/y1 in x0/
    }

    /// Print a downsampled version of a field to the console
    fn print_field(&self, which: Field) {
        let (label, field): (&str, &[[f32; NY]; NX]) = match which {
            Field::U => ("Velocity (ux, uy)", &self.ux),
            Field::Ph => ("Hot Gas Pressure", &self.p_h),
            Field::Divh => ("Divergence of Hot Gas Pressure", &self.div_h),
            Field::Phi => ("Level Set", &self.phi),
            Field::Temp => ("Temperature", &self.temp),
            Field::Rt => ("Reaction Parameter", &self.rt),
            Field::Dns => ("Smoke Density", &self.dns)
        };
        let (stride_x, stride_y) = (NX/30, NY/30);
        println!("{label} field ({NX}x{NY}), sampled with stride: ({stride_x}, {stride_y})");
        for y in (0..NY).rev().step_by(stride_y) {
            for x in (0..NX).step_by(stride_x) {
                print!("{:4.1} ", field[x][y]);
            }
            println!();
        }
    }
}

// ---------------- Main ----------------
fn main() {
    let mut sim = new_sim::<240, 320>(Params {
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
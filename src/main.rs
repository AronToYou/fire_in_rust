use crate::grid::{P, GridDisp, Linterp};
use indicatif::ProgressIterator;

mod grid;

const MIN_NORM: f32 = 1e-8;

// ---------------- Utility Functions ----------------
fn new_sim<const NX: usize, const NY: usize>(p: Params) -> Sim<NX, NY, impl Fn((f32, f32)) -> (f32, f32)> {
    let maxx = (NX as f32) - 1.001;
    let maxy = (NY as f32) - 1.001;
    let clamp_xy = move |(x, y): (f32, f32)| (x.clamp(0.0, maxx), y.clamp(0.0, maxy));
    Sim::<NX, NY, _>::new(p, clamp_xy)
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

    u: Box<[[P<f32>; NY]; NX]>,        // velocity field (x, y)
    p_h: Box<[[f32; NY]; NX]>,    // pressure for hot gas
    div_h: Box<[[f32; NY]; NX]>,  // divergence (hot gas)

    phi: Box<[[f32; NY]; NX]>,       // level set (+pos in fuel region, -neg outside, 0 at boundary)
    temp_gas: Box<[[f32; NY]; NX]>,  // temperature field (hot gas domain)
    rt: Box<[[f32; NY]; NX]>,        // reaction-time tracker (1 at fuel; decreases after crossing)
    dns: Box<[[f32; NY]; NX]>,       // smoke density (simple)

    tmp: Box<[[f32; NY]; NX]>, tmp2: Box<[[P<f32>; NY]; NX]>  // temporary intermediate fields for calculations
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
            u: Box::new([[P(0.0, 0.0); NY]; NX]),
            p_h: Box::new([[0.0; NY]; NX]),
            div_h: Box::new([[0.0; NY]; NX]),

            phi: Box::new([[1.0; NY]; NX]),
            temp_gas: Box::new([[p.temp_air; NY]; NX]),
            rt: Box::new([[0.0; NY]; NX]),
            dns: Box::new([[0.0; NY]; NX]),

            tmp: Box::new([[0.0; NY]; NX]), tmp2: Box::new([[P(0.0, 0.0); NY]; NX])
        };
        s.init_fuel_inlet();
        s
    }

    /// Initialize fuel inlet at bottom of domain
    fn init_fuel_inlet(&mut self) {
        for x in 0..NX {
            for y in 0..NY {
                if y < 10 {
                    self.u[x][y] = P(0.0, 1.5);
                    self.phi[x][y] = 5.0;
                    self.temp_gas[x][y] = self.p.temp_air + 10.0;
                } else {
                    self.phi[x][y] = -5.0;
                }
            }
        }
    }

    fn step(&mut self) {
        // (1) Update level set field //
        self.update_levelset();

        // (2) Add Forces {Bouyancy, vorticity confinement} //
        self.add_forces();

        // (3) Semi-Lagrangian advection of velocity fields //
        self.semi_lagrangian_advect();
    }

    // ----------------- (1) Thin-flame Level Set Propagation -----------------
    /// Updates the level set using upwind one-sided differencing to estimate spatial derivatives
    fn update_levelset(&mut self) {
        let k_react = self.p.k_react;
        self.tmp.copy_from_slice(&*self.phi);  // TODO
        for x in 1..NX-1 {
            for y in 1..NY-1 {
                // (1.3) (unscaled) Central differencing for normed gradient (∇φ/|∇φ|) //
                let gx = self.phi[x+1][y] - self.phi[x-1][y];  // gradient x-component
                let gy = self.phi[x][y+1] - self.phi[x][y-1];  // gradient y-component
                let norm = (gx*gx + gy*gy).sqrt().max(MIN_NORM);  // gradient norm
                
                // (1.2) Velocity of implicit surface (where φ==0) //
                let P(wx, wy) = self.u[x][y] + P(gx, gy)*(k_react/norm);
                
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
        std::mem::swap(&mut *self.phi, &mut *self.tmp);
    }

    fn add_forces(&mut self) {
        let force = &mut *self.tmp2;
        let (k_buoy, temp_air) = (self.p.k_buoy, self.p.temp_air);
        for x in 0..NX {
            for y in 0..NY {
                // (2.1) Buoyancy force α(T - T_air)ŷ //
                force[x][y] = P(0.0, k_buoy*(self.temp_gas[x][y] - temp_air));

                // (2.2.1) (unscaled) Vorticity ω //
                let P(_, dv_dx) = match x {
                    0 =>             self.u[x+1][y],
                    _ if x < NX-1 => self.u[x+1][y] - self.u[x-1][y],
                    _ =>                            - self.u[x-1][y],
                };
                let P(du_dy, _) = match y {
                    0 =>             self.u[x][y+1],
                    _ if y < NY-1 => self.u[x][y+1] - self.u[x][y-1],
                    _ =>                            - self.u[x][y-1],
                };
                self.tmp[x][y] = dv_dx - du_dy;
            }
        }
        let (vconf, h, dt) = (self.p.vconf, self.p.h, self.p.dt);
        for x in 1..NX-1 {
            for y in 1..NY-1 {
                // (2.2.2) (unscaled) Central differencing for normed gradient N = (∇|ω|/|∇|ω||) //
                let gx = self.tmp[x+1][y].abs() - self.tmp[x-1][y].abs();  // gradient x-component
                let gy = self.tmp[x][y+1].abs() - self.tmp[x][y-1].abs();  // gradient y-component
                let norm = (gx*gx + gy*gy).sqrt().max(MIN_NORM);  // gradient norm

                // (2.2.3) (scaled) Force of vorticity confinement εh(N x ω) //
                force[x][y] += P(-gy, gx)*self.tmp[x][y]*(vconf*h/norm);

                // (2.3) Add force //
                self.u[x][y] += force[x][y]*dt;
            }
        }
    }

    // ------------- (2) Semi-Lagrangian advection of velocity fields -------------
    /// Runge-Kutta 2-stage backtrace, bilinear velocity sampling, clamped at boundaries
    fn semi_lagrangian_advect(&mut self) {
        let (u, phi) = (&*self.u, &*self.phi);
        let (dt, h) = (self.p.dt, self.p.h);
        for x in 0..NX {
            for y in 0..NY {
                let P(x1, y1) = P(x as f32, y as f32);

                // Proceed accordingly if whether implicit surface is crossed //
                if self.sample_bilin(phi, (x1, y1)) > 0.0 {  // already in fuel region initially
                    // backtrace half-step to midpoint //
                    let P(x0, y0) = P(x1, y1) - u[x][y]*(0.5*dt/h);
                    let u0 = self.sample_bilin(u, (x0, y0));  // midpoint velocity

                    // backtrace full step //
                    let P(x0, y0) = P(x1, y1) - u0*(dt/h);
                    self.tmp2[x][y] = self.sample_bilin(u, (x0, y0));  // final velocity

                } else {
                    // backtrace half-step to midpoint //
                    let P(x0, y0) = P(x1, y1) - u[x][y]*(0.5*dt/h);
                    let phi_0 = self.sample_bilin(phi, (x0, y0));
                    let u0 = if phi_0 > 0.0 {  // if boundary crossed...
                        self.sample_ghost_velocity(P(x0, y0))  // ...sample ghost velocity
                    } else {
                        self.sample_bilin(u, (x0, y0))
                    };  // midpoint velocity

                    // backtrace full step //
                    let P(x0, y0) = P(x1, y1) - u0*(dt/h);
                    let phi_0 = self.sample_bilin(phi, (x0, y0));
                    self.tmp2[x][y] = if phi_0 > 0.0 {  // if boundary crossed...
                        self.sample_ghost_velocity(P(x0, y0))  // ...sample ghost velocity
                    } else {
                        self.sample_bilin(u, (x0, y0))
                    };  // final velocity
                }
            }
        }
        std::mem::swap(&mut *self.u, &mut *self.tmp2);
    }

    fn sample_ghost_velocity(&self, p: P<f32>) -> P<f32> {
        let P(x, y) = p;
        let (phi, u) = (&*self.phi, &*self.u);
        let (d_h, d_f, s) = (self.p.d_hgas, self.p.d_fuel, self.p.k_react);

            // (unscaled) Central differencing for normed gradient at non-integer coordinate
            let nx = self.sample_bilin(phi, (x+1.0, y)) - self.sample_bilin(phi, (x-1.0, y));
            let ny = self.sample_bilin(phi, (x, y+1.0)) - self.sample_bilin(phi, (x, y-1.0));
            let norm = (nx*nx + ny*ny).sqrt().max(MIN_NORM);
            let n = P(nx, ny)*(1.0/norm);  // implicit surface unit normal vector
            
            let u_ghost = self.sample_bilin(u, (x, y)) + n*(d_f/d_h - 1.0)*s;
            
            u_ghost
    }

    // ----------------- Utility Functions -----------------
    /// Bilinear sampling of scalar field, clamped at boundaries
    fn sample_bilin<T: Linterp>(&self, f: &[[T; NY]; NX], p: (f32, f32)) -> T {
        let (x, y) = (self.clamp_xy)(p);
        let P(x0, y0) = P(x, y).floor();
        let P(tx, ty) = P(x, y) - P(x0 as f32, y0 as f32);

        let f00 = f[x0][y0];       let f01 = f[x0][y0+1];
        let f10 = f[x0+1][y0];       let f11 = f[x0+1][y0+1];  // store x1 in x0
        let a = f00*(1.0-tx) + f10*tx; let b = f01*(1.0-tx) + f11*tx;  // store tx in x1/y1 in x0/
        a*(1.0-ty) + b*ty  // store ty in tx in x1/y1 in x0/
    }

    /// Print a downsampled version of a field to the console
    fn print_field(&self, which: Field) {
        let (label, field): (&str, &dyn GridDisp) = match which {
            Field::U => ("Velocity (ux, uy)", &*self.u),
            Field::Ph => ("Hot Gas Pressure", &*self.p_h),
            Field::Divh => ("Divergence of Hot Gas Pressure", &*self.div_h),
            Field::Phi => ("Level Set", &*self.phi),
            Field::Temp => ("Temperature", &*self.temp_gas),
            Field::Rt => ("Reaction Parameter", &*self.rt),
            Field::Dns => ("Smoke Density", &*self.dns)
        };
        let (stride_x, stride_y) = ((NX/30).max(1), (NY/30).max(1));
        println!("{label} field ({NX}x{NY}), sampled with stride: ({stride_x}, {stride_y})");
        for y in (0..NY).rev().step_by(stride_y) {
            for x in (0..NX).step_by(stride_x) {
                print!("{:4.1} ", field.at(x, y));
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
    sim.print_field(Field::U);
    
    let steps = 100;
    println!("Simulating {steps} steps...");
    for _ in (0..steps).progress() {
        sim.step();
    }
    sim.print_field(Field::Phi);
    sim.print_field(Field::U);
    println!("Done!");
}
use indicatif::ProgressIterator;
use std::time::Duration;
use std::thread::sleep;

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
struct Sim {
    p: Params,  // Simulation parameters (defined above)
    
    uf: Vec<f32>, vf: Vec<f32>,  // velocity field fuel    (x, y)
    uh: Vec<f32>, vh: Vec<f32>,  // velocity field hot gas (x, y)
    p_h: Vec<f32>,    // pressure for hot gas
    div_h: Vec<f32>,  // divergence (hot gas)

    phi: Vec<f32>,   // level set (+pos in fuel region, -neg outside, 0 at boundary)
    temp: Vec<f32>,  // temperature field (hot gas domain)
    rt: Vec<f32>,    // reaction-time tracker (1 at fuel; decreases after crossing)
    dns: Vec<f32>,   // smoke density (simple)

    tmp: Vec<f32>  // temporary intermediate field for calculations
}

impl Sim {
    fn new(p: Params) -> Self {
        let n = p.nx*p.ny;
        let mut s = Self {
            p,
            uf: vec![0.0; n], vf: vec![0.0; n],
            uh: vec![0.0; n], vh: vec![0.0; n],
            p_h: vec![0.0; n],
            div_h: vec![0.0; n],

            phi: vec![1.0; n],
            temp: vec![p.temp_air; n],
            rt: vec![0.0; n],
            dns: vec![0.0; n],

            tmp: vec![0.0; n]
        };
        s.init_fuel_inlet();
        s
    }

    fn init_fuel_inlet(&mut self) {
        let (nx, ny) = (self.p.nx, self.p.ny);
        let mut idx: usize = 0;
        for y in 0..ny {
            for _ in 0..nx {
                if y < 10 {
                    self.vf[idx] = 1.5;
                    self.phi[idx] = 5.0;
                } else {
                    self.phi[idx] = -5.0;
                }
                idx += 1;
            }
        }
    }
}

// ---------------- Main ----------------
fn main() {
    let sim = Sim::new(Params {
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

    
    let steps = 300;
    let ten_ms = Duration::from_millis(10);
    println!("Simulating {steps} steps...");
    for _ in (0..steps).progress() {
        sleep(ten_ms);
    }
    println!("Done!");
}
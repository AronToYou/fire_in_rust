# Outline
Incompressible Navier-Stokes
$$\pdv{\vec{u}}{t} = -(\vec{u}\cdot\nabla)\vec{u} - \frac{1}{\rho}\nabla p + \nu\nabla^2\vec{u} + \vec{f}$$
Velocity updated for incompressible flow (Stam's 4-step loop)
1. add force
2. advect
3. diffuse
4. project
## Stable Fluids
## A) **Thin flame** tracked by **level set** $\phi$[^1]
1. Implicit surface where ($\phi=0$) defines the interface between
	- **Fuel** ($\phi>0$) blue-core reaction zone
	- **Hot Gas** ($\phi<0$)
2. velocity of surface
$$w = \vec{u}_{fuel} + S\vec{n}$$
3. Implicit surface unit normal vector (central differencing)
$$\vec{n} = \nabla\phi/|\nabla\phi|$$
4. Time derivative of level set (upwind differencing for $\nabla\phi$)
$$\phi_t = -(\vec{u}_f + S\vec{n})\cdot\nabla\phi$$
- evolve **fuel** and **hot gas** velocity fields separately
## B) Velocity updated via Stam's 4-step loop[^3]
$$\vec{\text{w}}_0(\vec{x}) := \vec{u}(\vec{x},\; t)$$
$$\vec{\text{w}}_0 = \vec{u}_t$$
$$\vec{\text{w}}_0\rightarrow\text{add force}\rightarrow\vec{\text{w}}_1\rightarrow\text{advect}\rightarrow\vec{\text{w}}_2\rightarrow\text{diffuse}\rightarrow\vec{\text{w}}_3\rightarrow\text{project}\rightarrow\vec{\text{w}}_4$$
$$\vec{u}(\vec{x},\; t+\Delta t) = \vec{\text{w}}_4(\vec{x})$$
### 1. **Add Force** (before advection)[^2]
 1. Buoyancy
$$f_{buoy} = \alpha(T - T_{air})\hat{z}$$
2. Vorticity Confinement
	1. Vorticity vector
$$\vec{\omega} = \nabla\times\vec{u}$$
	2. Normalized vorticity location vector (central differencing)
$$\vec{N} = \nabla|\vec{\omega}|/|\nabla|\vec{\omega}||$$
	 3. Force of vorticity confinement
 $$f_{conf} = \varepsilon h(\vec{N}\times\vec{\omega})$$
 3. Add Force
$$\vec{\text{w}}_1 = \vec{u}_t + \triangle t\ \vec{f}$$
### 2. **Advection**[^3]
- Stable Semi-Lagrangian
	- `semi_lagrangian_advect`
$$\text{Advection of }\vec{\text{w}}_1(P_0)\text{ where } P_0 := (x_0, y_0)$$
1. backtrace half-step
$$P_{-1/2} = P_0 - \vec{\text{w}}_1(P_0)\frac{\Delta t}{2h}$$
2. sample velocity* 
$$\vec{\text{w}}_1(P_{-1/2}) \approx \text{bilinear\_sample}(\vec{\text{w}}_1,\; P_{-1/2})$$
3. backtrace full-step
$$P_{-1} = P_0 - \vec{\text{w}}_1(P_{-1/2})\frac{\Delta t}{h}$$
4. sample velocity*
$$\vec{\text{w}}_2(P_0) \approx \vec{\text{w}}_1(P_{-1}) \approx \text{bilinear\_sample}(\vec{\text{w}}_1,\; P_{-1})$$
#### **\*Expansion coupling** across front via **ghost normal velocity** 
- (mass conservation) $\nabla\cdot\vec{u}=0$
- i.e. when **hot gas** advection samples **fuel**, synthesize correct cross-interface value
	- $V_h^{ghost} = V_f + (\rho_f/\rho_h - 1)S$
		- (normal velocity of fuel) $V_f = \vec{u}_f\cdot\vec{n}$
	- $\vec{u}_h^{ghost} = V_h^{ghost}\vec{n} + \vec{u}_f - (\vec{u}_f\cdot\vec{n})\vec{n}$
	- $\vec{u}_h^{ghost} = (\rho_f/\rho_h - 1)S + \vec{u}_f$
### 3. **Diffusion**
- with implicit diffusion
	- `diffuse_velocity`
$$\pdv{\vec{\text{w}}_2}{t} = \nu\nabla^2\vec{\text{w}}_2$$
$$\vec{\text{w}}_3 = \vec{\text{w}}_2 + (\nu\Delta t)\nabla^2\vec{\text{w}}_2$$
### 4. **Projection**
- Onto divergent free fields
- and a Poisson solve for pressure
	- `project_hot`
$$\nabla^2q = \nabla\cdot\vec{\text{w}}_3$$
$$\vec{\text{w}}_4 = \vec{\text{w}}_3 - \nabla q$$
## **Reaction-time scalar Y** (advected + linear source)
- (advected and decreased by 1 per unit time)
- time since crossing blue core
#### **Temperature** (advected) 
- with radiative cooling 
	- proportional to $(T - T_{air})^4$
$T_t = -(\vec{u}\cdot\nabla)T - c_T\Big(\frac{T - T_{air}}{T_{max} - T_{air}}\Big)^4$
- $c_T$  :  cooling constant
	- `p.k_cool`
- By Fourier's law for an isotropic medium, the rate of flow of heat energy per unit area through a surface is proportional to the negative temperature gradient across it:
	- $\vec{q}_{heat flow}  = -k\nabla T(\vec{x}, t)$
## **Offline rendering** simple emission/absorption with blackbody colormap
- not full Monte-Carlo ray marcher, but good first pass)
- for full correctness, see stochastic ray marching of the RTE (stanford)
- Diffusion (via Jacobi or CG) not applies to `T`, `Y`, or `D`

| Variable | Description                      |
| -------- | -------------------------------- |
| N        | grid length (Vectors are N x N)  |
| h        | grid spacing                     |
| phi      | (grid-center) level set function |
| p        | (grid-center) pressure           |
| T        | (grid-center) temperature        |
| vf       | velocity of **fuel**             |
| vg       | velocity of **gas**              |
- (phi = 0) defines interface between **fuel** and **gas**
- velocity units of 'cells / time'
- cell-centered fields
- Conjugate-Gradient Poisson solver (SPD 5-point Laplacian)
- Ghost-velocity
	- (applied during semi-Lagrangian sampling of hot-gas field)
	- pressure jump conditions not enforced in the linear system




# Miscellaneous Notes
## Discretization
### Staggered Field Grids[^2]
#### Velocity Field
Defined at faces/edges of grid cells
$$\vec{u} = (u, v)$$
$$u_{i+1/2,\; j}\quad i\in\{0,\cdots, N_x\}\quad j\in\{1,\cdots,N_y\}$$
$$v_{i,\; j+1/2}\quad i\in\{1,\cdots, N_x\}\quad j\in\{0,\cdots,N_y\}$$
##### Cell-Centered Velocities
$$\bar{u}_{i,\; j} = (u_{i+1/2,\; j} + u_{i-1/2,\; j})/2 \quad \bar{v}_{i,\; j} = (v_{,\; j+1/2} + v_{i,\; j-1/2})/2$$
##### Vorticity
$$|\omega| = \omega^3_{i,\; j} = (\bar{v}_{i+1,\; j} - \bar{v}_{i-1,\; j} - \bar{u}_{i,\; j+1} + \bar{u}_{i,\; j-1})/2h$$
#### All Other Fields
Defined at center of grid cells
$$F_{i,\; j}\quad i\in\{1,\cdots, N_x\}\quad j\in\{1,\cdots,N_y\}$$


## Advection Equation
Advection for a conserved quantity described by a scalar field $\psi(t, x, y, z)$ that flows with velocity $\vec{u}$
$$\pdv{\psi}{t} + \nabla\cdot(\psi\vec{u}) = 0$$
$$\big(\nabla\cdot\vec{u} = 0\big) \implies \pdv{\psi}{t} + \vec{u}\cdot\nabla\psi = 0$$
## Advection-Diffusion-Dissipation Equation
For a generic scalar field $S$ that flows with velocity $\vec{u}$
$$\pdv{S}{t} = -(\vec{u}\cdot\nabla)S + k_S\nabla^2S - a_SS + source$$
- $k_S$  :  diffusion constant
- $a_S$  :  dissipation rate

[^1]: Nguyen, Duc & Fedkiw, Ronald & Jensen, Henrik. (2002). [Physically Based Modeling and Animation of Fire](http://dx.doi.org/10.1145/566570.566643). ACM Transactions on Graphics. 21. 10.1145/566570.566643. 
[^2]: Fedkiw, Ronald & Stam, Jos & Jensen, Henrik. (2001). [Visual Simulation of Smoke](http://dx.doi.org/10.1145/383259.383260). ACM SIGGRAPH2001, 2001.. 10.1145/383259.383260. 
[^3]: Stam, Jos. (2001). [Stable Fluids](http://dx.doi.org/10.1145/311535.311548). ACM SIGGRAPH 99. 1999. 10.1145/311535.311548. 

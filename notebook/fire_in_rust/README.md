# Outline
Incompressible Navier-Stokes
$$\pdv{\mathbf{u}}{t} = -(\mathbf{u}\cdot\nabla)\mathbf{u} - \frac{1}{\rho}\nabla p + \nu\nabla^2\mathbf{u} + \mathbf{f}$$
## Stable Fluids
## 1. **Thin flame** tracked by **level set** $\phi$
1. Implicit surface where ($\phi=0$) defines the interface between
	- **Fuel** ($\phi>0$) blue-core reaction zone
	- **Hot Gas** ($\phi<0$)
2. velocity of surface
$$w = u_{fuel} + S\cdot n$$
3. Implicit surface unit normal vector (central differencing)
$$\mathbf{n} = \nabla\phi/|\nabla\phi|$$
4. Time derivative of level set (upwind differencing for $\nabla\phi$)
$$\phi_t = -(\mathbf{u}_f + S\mathbf{n})\cdot\nabla\phi$$
- evolve **fuel** and **hot gas** velocity fields separately
### 2. **Add Force** (before advection)
#### Buoyancy
- $f_{buoy} = \alpha(T - T_{air})\hat{z}$
#### Vorticity Confinement
- $f_{conf} = \varepsilon h(N\times \omega)$ 
### 3. **Advection**
- for incompressible flow (Stam's 4-step loop)
- Stable Semi-Lagrangian
	- `semi_lagrangian_advect`
### 4. **Diffusion**
- with implicit diffusion
	- `diffuse_velocity`
### 4. **Projection**
- and a Poisson solve for pressure
	- `project_hot`
## **Expansion coupling** across front via **ghost normal velocity** 
- (mass conservation) $\nabla\cdot\mathbf{u}=0$
- i.e. when **hot gas** advection samples **fuel**, synthesize correct cross-interface value
	- $V_h^{ghost} = V_f + (\rho_f/\rho_h)S$
		- (normal velocity of fuel) $V_f = \mathbf{u}_f\cdot\mathbf{n}$
	- $u_h^{ghost} = V_h^{ghost}\mathbf{n} + \mathbf{u}_f - (\mathbf{u}_f\cdot\mathbf{n})\mathbf{n}$
## **Reaction-time scalar Y** (advected + linear source)
- (advected and decreased by 1 per unit time)
- time since crossing blue core
#### **Temperature** (advected) 
- with radiative cooling 
	- proportional to $(T - T_{air})^4$
$T_t = -(\mathbf{u}\cdot\nabla)T - c_T\Big(\frac{T - T_{air}}{T_{max} - T_{air}}\Big)^4$
- $c_T$  :  cooling constant
	- `p.k_cool`
- By Fourier's law for an isotropic medium, the rate of flow of heat energy per unit area through a surface is proportional to the negative temperature gradient across it:
	- $\mathbf{q}_{heat flow}  = -k\nabla T(\mathbf{x}, t)$
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
## Advection Equation
Advection for a conserved quantity described by a scalar field $\psi(t, x, y, z)$ that flows with velocity $\mathbf{u}$
$$\pdv{\psi}{t} + \nabla\cdot(\psi\mathbf{u}) = 0$$
$$\big(\nabla\cdot\mathbf{u} = 0\big) \implies \pdv{\psi}{t} + \mathbf{u}\cdot\nabla\psi = 0$$
## Advection-Diffusion-Dissipation Equation
For a generic scalar field $S$ that flows with velocity $\mathbf{u}$
$$\pdv{S}{t} = -(\mathbf{u}\cdot\nabla)S + k_S\nabla^2S - a_SS + source$$
- $k_S$  :  diffusion constant
- $a_S$  :  dissipation rate
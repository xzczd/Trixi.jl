
@doc raw"""
    LatticeBoltzmannEquation2D(; Ma, Re, collision_op=collision_bgk,
                               c=1, L=1, rho0=1, u0=nothing, nu=nothing)

The Lattice-Boltzmann equation
```math
\partial_t u_\alpha + v_{\alpha,1} \partial_1 u_\alpha + v_{\alpha,2} \partial_2 u_\alpha = 0
```
in two space dimensions for the D2Q9 scheme.

The characteristic Mach number and Reynolds numbers are specified as `Ma` and `Re`. By the
default, the collision operator `collision_op` is set to the BGK model. `c`, `L`, and `rho0`
specify the mean thermal molecular velocity, the characteristic length, and the reference density,
respectively. They can usually be left to the default values. If desired, instead of the Mach
number, one can set the macroscopic reference velocity `u0` directly (`Ma` needs to be set to
`nothing` in this case). Likewise, instead of the Reynolds number one can specify the kinematic
viscosity `nu` directly (in this case, `Re` needs to be set to `nothing`).


The nine discrete velocity directions of the D2Q9 scheme are sorted as follows:

```
  6  2  5
   ╲ │ ╱
    ╲│╱
  3──9──1
    ╱│╲
   ╱ │ ╲
  7  4  8
```

The corresponding opposite directions are:
* 1 ←→ 3
* 2 ←→ 4
* 3 ←→ 1
* 4 ←→ 2
* 5 ←→ 7
* 6 ←→ 8
* 7 ←→ 5
* 8 ←→ 6
* 9 ←→ 9

The main sources for the base implementation were
1. Misun Min, Taehun Lee, **A spectral-element discontinuous Galerkin lattice Boltzmann method for
   nearly incompressible flows**, J Comput Phys 230(1), 2011
   [doi:10.1016/j.jcp.2010.09.024](https://doi.org/10.1016/j.jcp.2010.09.024)
2. Karsten Golly, **Anwendung der Lattice-Boltzmann Discontinuous Galerkin Spectral Element Method
   (LB-DGSEM) auf laminare und turbulente nahezu inkompressible Strömungen im dreidimensionalen
   Raum**, Master Thesis, University of Cologne, 2018.
3. Dieter Hänel, **Molekulare Gasdynamik**, Springer-Verlag Berlin Heidelberg, 2004
   [doi:10.1007/3-540-35047-0](https://doi.org/10.1007/3-540-35047-0)
"""
struct LatticeBoltzmannEquation2D{RealT<:Real, CollisionOp} <: AbstractLatticeBoltzmannEquation{2, 9}
  c::RealT    # mean thermal molecular velocity
  c_s::RealT  # isothermal speed of sound
  rho0::RealT # macroscopic reference density

  Ma::RealT   # characteristic Mach number
  u0::RealT   # macroscopic reference velocity

  Re::RealT   # characteristic Reynolds number
  L::RealT    # reference length
  nu::RealT   # kinematic viscosity

  weights::SVector{9, RealT}  # weighting factors for the equilibrium distribution
  v_alpha1::SVector{9, RealT} # discrete molecular velocity components in x-direction
  v_alpha2::SVector{9, RealT} # discrete molecular velocity components in y-direction

  collision_op::CollisionOp   # collision operator for the collision kernel
end

function LatticeBoltzmannEquation2D(; Ma, Re, collision_op=collision_bgk,
                                    c=1, L=1, rho0=1, u0=nothing, nu=nothing)
  # Sanity check that exactly one of Ma, u0 is not `nothing`
  if isnothing(Ma) && isnothing(u0)
    error("Mach number `Ma` and reference speed `u0` may not both be `nothing`")
  elseif !isnothing(Ma) && !isnothing(u0)
    error("Mach number `Ma` and reference speed `u0` may not both be set")
  end

  # Sanity check that exactly one of Re, nu is not `nothing`
  if isnothing(Re) && isnothing(nu)
    error("Reynolds number `Re` and visocsity `nu` may not both be `nothing`")
  elseif !isnothing(Re) && !isnothing(nu)
    error("Reynolds number `Re` and visocsity `nu` may not both be set")
  end

  # Calculate isothermal speed of sound
  # The relation between the isothermal speed of sound `c_s` and the mean thermal molecular velocity
  # `c` depends on the used phase space discretization, and is valid for D2Q9 (and others). For
  # details, see, e.g., [3] in the docstring above.
  c_s = c / sqrt(3)

  # Calculate missing quantities
  if isnothing(Ma)
    Ma = u0 / c_s
  elseif isnothing(u0)
    u0 = Ma * c_s
  end
  if isnothing(Re)
    Re = u0 * L / nu
  elseif isnothing(nu)
    nu = u0 * L / Re
  end

  # Promote to common data type
  Ma, Re, c, L, rho0, u0, nu = promote(Ma, Re, c, L, rho0, u0, nu)

  # Source for weights and speeds: https://cims.nyu.edu/~billbao/report930.pdf
  weights  = @SVector [1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36, 4/9]
  v_alpha1 = @SVector [ c,   0,  -c,   0,   c,   -c,   -c,    c,    0 ]
  v_alpha2 = @SVector [ 0,   c,   0,  -c,   c,    c,   -c,   -c,    0 ]

  LatticeBoltzmannEquation2D(c, c_s, rho0, Ma, u0, Re, L, nu,
                             weights, v_alpha1, v_alpha2,
                             collision_op)
end


get_name(::LatticeBoltzmannEquation2D) = "LatticeBoltzmannEquation2D"
varnames_cons(equations::LatticeBoltzmannEquation2D) = ntuple(v -> "pdf"*string(v), nvariables(equations))
varnames_prim(::LatticeBoltzmannEquation2D) = @SVector ["rho", "v1", "v2", "p"]

# Set initial conditions at physical location `x` for time `t`
"""
    initial_condition_constant(x, t, equation::LatticeBoltzmannEquation2D)

A constant initial condition to test free-stream preservation.
"""
function initial_condition_constant(x, t, equation::LatticeBoltzmannEquation2D)
  @unpack u0 = equation
  rho = pi
  v1 = u0
  v2 = u0

  return equilibrium_distribution(rho, v1, v2, equation)
end


"""
    boundary_condition_wall_noslip(u_inner, orientation, direction, x, t,
                                   surface_flux_function,
                                   equation::LatticeBoltzmannEquation2D)

No-slip wall boundary condition using the bounce-back approach.
"""
function boundary_condition_wall_noslip(u_inner, orientation, direction, x, t,
                                        surface_flux_function,
                                        equation::LatticeBoltzmannEquation2D)
  # For LBM no-slip wall boundary conditions, we set the boundary state to
  # - the inner state for outgoing particle distribution functions
  # - the *opposite* inner state for all other particle distribution functions
  # See the list of (opposite) directions in the docstring of `LatticeBoltzmannEquation2D`.
  if direction == 1 # boundary in -x direction
    pdf1 = u_inner[3]
    pdf2 = u_inner[4]
    pdf3 = u_inner[3] # outgoing
    pdf4 = u_inner[2]
    pdf5 = u_inner[7]
    pdf6 = u_inner[6] # outgoing
    pdf7 = u_inner[7] # outgoing
    pdf8 = u_inner[6]
    pdf9 = u_inner[9]
  elseif direction == 2 # boundary in +x direction
    pdf1 = u_inner[1] # outgoing
    pdf2 = u_inner[4]
    pdf3 = u_inner[1]
    pdf4 = u_inner[2]
    pdf5 = u_inner[5] # outgoing
    pdf6 = u_inner[8]
    pdf7 = u_inner[5]
    pdf8 = u_inner[8] # outgoing
    pdf9 = u_inner[9]
  elseif direction == 3 # boundary in -y direction
    pdf1 = u_inner[3]
    pdf2 = u_inner[4]
    pdf3 = u_inner[1]
    pdf4 = u_inner[4] # outgoing
    pdf5 = u_inner[7]
    pdf6 = u_inner[8]
    pdf7 = u_inner[7] # outgoing
    pdf8 = u_inner[8] # outgoing
    pdf9 = u_inner[9]
  else # boundary in +y direction
    pdf1 = u_inner[3]
    pdf2 = u_inner[2] # outgoing
    pdf3 = u_inner[1]
    pdf4 = u_inner[2]
    pdf5 = u_inner[5] # outgoing
    pdf6 = u_inner[6] # outgoing
    pdf7 = u_inner[5]
    pdf8 = u_inner[6]
    pdf9 = u_inner[9]
  end
  u_boundary = SVector(pdf1, pdf2, pdf3, pdf4, pdf5, pdf6, pdf7, pdf8, pdf9)

  # Calculate boundary flux
  if direction in (2, 4) # u_inner is "left" of boundary, u_boundary is "right" of boundary
    flux = surface_flux_function(u_inner, u_boundary, orientation, equation)
  else # u_boundary is "left" of boundary, u_inner is "right" of boundary
    flux = surface_flux_function(u_boundary, u_inner, orientation, equation)
  end

  return flux
end


function boundary_condition_moving_wall_ypos(u_inner, orientation, direction, x, t,
                                             surface_flux_function,
                                             equation::LatticeBoltzmannEquation2D)
  @assert direction == 4 "moving wall assumed in +y direction"

  @unpack rho0, u0, weights, c_s = equation
  cs_squared = c_s^2

  pdf1 = u_inner[3] + 2 * weights[1] * rho0 * u0 / cs_squared
  pdf2 = u_inner[2] # outgoing
  pdf3 = u_inner[1] + 2 * weights[3] * rho0 * (-u0) / cs_squared
  pdf4 = u_inner[2]
  pdf5 = u_inner[5] # outgoing
  pdf6 = u_inner[6] # outgoing
  pdf7 = u_inner[5] + 2 * weights[7] * rho0 * (-u0) / cs_squared
  pdf8 = u_inner[6] + 2 * weights[8] * rho0 * u0 / cs_squared
  pdf9 = u_inner[9]

  u_boundary = SVector(pdf1, pdf2, pdf3, pdf4, pdf5, pdf6, pdf7, pdf8, pdf9)

  # Calculate boundary flux (u_inner is "left" of boundary, u_boundary is "right" of boundary)
  flux = surface_flux_function(u_inner, u_boundary, orientation, equation)

  return flux
end


"""
    initial_condition_lid_driven_cavity(x, t, equation::LatticeBoltzmannEquation2D)

Initial state for a lid-driven cavity flow setup. To be used in combination with
[`boundary_condition_lid_driven_cavity`](@ref) and [`boundary_condition_wall_noslip`](@ref).
"""
function initial_condition_lid_driven_cavity(x, t, equation::LatticeBoltzmannEquation2D)
  @unpack L, u0, nu = equation

  rho = 1
  v1 = 0
  v2 = 0

  return equilibrium_distribution(rho, v1, v2, equation)
end


"""
    boundary_condition_lid_driven_cavity(u_inner, orientation, direction, x, t,
                                         surface_flux_function,
                                         equation::LatticeBoltzmannEquation2D)
    
Boundary condition for a lid-driven cavity flow setup, where the top lid (+y boundary) is a moving
no-slip wall. To be used in combination with [`initial_condition_lid_driven_cavity`](@ref).
"""
function boundary_condition_lid_driven_cavity(u_inner, orientation, direction, x, t,
                                              surface_flux_function,
                                              equation::LatticeBoltzmannEquation2D)
  return boundary_condition_moving_wall_ypos(u_inner, orientation, direction, x, t,
                                             surface_flux_function, equation)
end


"""
    initial_condition_couette_unsteady(x, t, equation::LatticeBoltzmannEquation2D)

Initial state for an *unsteady* Couette flow setup, which is also the exact solution for the
incompressible Navier-Stokes equations. To be used in combination with
[`boundary_condition_couette`](@ref) and [`boundary_condition_wall_noslip`](@ref). In the limit,
this setup will converge to the state set in [`initial_condition_couette_steady`](@ref).
"""
function initial_condition_couette_unsteady(x, t, equation::LatticeBoltzmannEquation2D)
  @unpack L, u0, rho0, nu = equation

  x1, x2 = x
  v1 = u0*x2/L
  for m in 1:100
    lambda_m = m * pi / L
    v1 += 2 * u0 * (-1)^m/(lambda_m * L) * exp(-nu * lambda_m^2 * t) * sin(lambda_m * x2)
  end

  rho = 1
  v2 = 0

  return equilibrium_distribution(rho, v1, v2, equation)
end


"""
    initial_condition_couette_unsteady(x, t, equation::LatticeBoltzmannEquation2D)

Initial state for a *steady* Couette flow setup. To be used in combination with
[`boundary_condition_couette`](@ref) and [`boundary_condition_wall_noslip`](@ref).
"""
function initial_condition_couette_steady(x, t, equation::LatticeBoltzmannEquation2D)
  @unpack L, u0, rho0 = equation

  rho = rho0
  v1 = u0 * x[2] / L
  v2 = 0

  return equilibrium_distribution(rho, v1, v2, equation)
end


"""
    boundary_condition_couette(u_inner, orientation, direction, x, t,
                               surface_flux_function,
                               equation::LatticeBoltzmannEquation2D)

Moving *upper* wall boundary condition for a Couette flow setup. To be used in combination with
[`boundary_condition_wall_noslip`](@ref) for the lower wall and
[`boundary_condition_periodic`](@ref) for the lateral boundaries.
"""
function boundary_condition_couette(u_inner, orientation, direction, x, t,
                                    surface_flux_function,
                                    equation::LatticeBoltzmannEquation2D)
  return boundary_condition_moving_wall_ypos(u_inner, orientation, direction, x, t,
                                             surface_flux_function, equation)
end


# Pre-defined source terms should be implemented as
# function source_terms_WHATEVER(u, x, t, equation::LatticeBoltzmannEquation2D)


# Calculate 1D flux in for a single point
@inline function calcflux(u, orientation, equation::LatticeBoltzmannEquation2D)
  if orientation == 1
    v_alpha = equation.v_alpha1
  else
    v_alpha = equation.v_alpha2
  end
  return v_alpha .* u
end


function flux_lax_friedrichs(u_ll, u_rr, orientation, equation::LatticeBoltzmannEquation2D)
  if orientation == 1
    v_alpha = equation.v_alpha1
  else
    v_alpha = equation.v_alpha2
  end
  return 0.5 * ( v_alpha .* (u_ll + u_rr) - abs.(v_alpha) .* (u_rr - u_ll) )
end


"""
    density(u, equation::LatticeBoltzmannEquation2D)

Calculate the macroscopic density from the particle distribution functions `u`.
"""
@inline density(u, equation::LatticeBoltzmannEquation2D) = sum(u)


"""
    velocity(u, orientation, equation::LatticeBoltzmannEquation2D)

Calculate the macroscopic velocity for the given `orientation` (1 -> x, 2 -> y) from the
particle distribution functions `u`.
"""
@inline function velocity(u, orientation, equation::LatticeBoltzmannEquation2D)
  if orientation == 1
    v_alpha = equation.v_alpha1
  else
    v_alpha = equation.v_alpha2
  end
  
  return sum(v_alpha .* u)/density(u, equation)
end


"""
    velocity(u, equation::LatticeBoltzmannEquation2D)

Calculate the macroscopic velocity vector from the particle distribution functions `u`.
"""
@inline function velocity(u, equation::LatticeBoltzmannEquation2D)
  @unpack v_alpha1, v_alpha2 = equation
  rho = density(u, equation)
  
  return SVector(sum(v_alpha1 .* u)/rho, sum(v_alpha2 .* u)/rho)
end


"""
    pressure(u, equation::LatticeBoltzmannEquation2D)

Calculate the macroscopic pressure from the particle distribution functions `u`.
"""
pressure(u, equation::LatticeBoltzmannEquation2D) = density(u, equation) * equation.c_s^2


"""
    equilibrium_distribution(alpha, rho, v1, v2, equation::LatticeBoltzmannEquation2D)

Calculate the local equilibrium distribution for the distribution function with index `alpha` and
given the macroscopic state defined by `rho`, `v1`, `v2`.
"""
@inline function equilibrium_distribution(alpha, rho, v1, v2, equation::LatticeBoltzmannEquation2D)
  @unpack weights, c_s, v_alpha1, v_alpha2 = equation

  va_v = v_alpha1[alpha]*v1 + v_alpha2[alpha]*v2
  cs_squared = c_s^2
  v_squared = v1^2 + v2^2

  return weights[alpha] * rho * (1 + va_v/cs_squared
                                   + va_v^2/(2*cs_squared^2)
                                   - v_squared/(2*cs_squared))
end


function equilibrium_distribution(alpha, u, equation::LatticeBoltzmannEquation2D)
  rho = density(u, equation)
  v1, v2 = velocity(u, equation)

  return equilibrium_distribution(alpha, rho, v1, v2, equation)
end


@inline function equilibrium_distribution(rho, v1, v2, equation::LatticeBoltzmannEquation2D)
  return SVector(equilibrium_distribution(1, rho, v1, v2, equation),
                 equilibrium_distribution(2, rho, v1, v2, equation),
                 equilibrium_distribution(3, rho, v1, v2, equation),
                 equilibrium_distribution(4, rho, v1, v2, equation),
                 equilibrium_distribution(5, rho, v1, v2, equation),
                 equilibrium_distribution(6, rho, v1, v2, equation),
                 equilibrium_distribution(7, rho, v1, v2, equation),
                 equilibrium_distribution(8, rho, v1, v2, equation),
                 equilibrium_distribution(9, rho, v1, v2, equation))
end


function equilibrium_distribution(u, equation::LatticeBoltzmannEquation2D)
  rho = density(u, equation)
  v1, v2 = velocity(u, equation)

  return equilibrium_distribution(rho, v1, v2, equation)
end


"""
    collision_bgk(u, dt, equation::LatticeBoltzmannEquation2D)

Collision operator for the Bhatnagar, Gross, and Krook (BGK) model.
"""
@inline function collision_bgk(u, dt, equation::LatticeBoltzmannEquation2D)
  @unpack c_s, nu = equation
  tau = nu / (c_s^2 * dt)
  return -(u - equilibrium_distribution(u, equation))/(tau + 1/2)
end



@inline have_constant_speed(::LatticeBoltzmannEquation2D) = Val(true)

@inline function max_abs_speeds(equation::LatticeBoltzmannEquation2D)
  @unpack c = equation

  return c, c
end


# Convert conservative variables to primitive (macroscopic)
@inline function cons2prim(u, equation::LatticeBoltzmannEquation2D)
  return SVector(density(u, equation),
                 velocity(u, equation)...,
                 pressure(u, equation))
end

# Convert conservative variables to entropy variables
@inline cons2entropy(u, equation::LatticeBoltzmannEquation2D) = u

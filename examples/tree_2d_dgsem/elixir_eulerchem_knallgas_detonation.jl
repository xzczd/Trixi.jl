
using OrdinaryDiffEq
using Trixi
using KROME

###############################################################################
# semidiscretization of the compressible Euler equations

equations = CompressibleEulerMultichemistryEquations2D(gammas             = (1.4, 1.4, 1.4, 1.4, 1.4),
                                                       gas_constants      = (4.1242, 0.2598, 0.4, 0.4615, 0.2968),
                                                       heat_of_formations = (0.0, 0.0, -50.0, -100.0, 0.0))

initial_condition = initial_condition_knallgas_detonation

boundary_conditions = boundary_condition_knallgas_detonation

#boundary_condition = BoundaryConditionDirichlet(initial_condition)
#boundary_conditions = (x_neg=boundary_condition,
#                       x_pos=boundary_condition,
#                       y_neg=boundary_condition,
#                       y_pos=boundary_condition)

chemistry_term = chemistry_knallgas_detonation

surface_flux = flux_hllc
volume_flux  = flux_central
basis = LobattoLegendreBasis(3)
indicator_sc = IndicatorHennemannGassner(equations, basis,
                                         alpha_max = 1.0,
                                         alpha_min = 0.0001,
                                         alpha_smooth = true,
                                         variable=density_pressure)
volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                 volume_flux_dg=volume_flux,
                                                 volume_flux_fv=surface_flux)
solver = DGSEM(basis, surface_flux, volume_integral)

coordinates_min = ( 0.0, 0.0)
coordinates_max = ( 6.0, 6.0)   # eigentlich (6, 2)

#mesh = StructuredMesh((100, 100), coordinates_min, coordinates_max, periodicity=false)

mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level=8,
                n_cells_max=100_000,
                periodicity=false)


semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver, 
                                    boundary_conditions=boundary_conditions, source_terms=nothing,
                                    chemistry_terms=chemistry_term)


###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 0.1)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000

analysis_callback = AnalysisCallback(semi, interval=analysis_interval)

alive_callback = AliveCallback(analysis_interval=analysis_interval)

save_solution = SaveSolutionCallback(interval=10000,
                                     save_initial_solution=true,
                                     save_final_solution=true,
                                     solution_variables=cons2prim)

stepsize_callback = StepsizeCallback(cfl=0.1)

chemistry_callback = KROMEChemistryCallback()

callbacks = CallbackSet(summary_callback,
                        analysis_callback, 
                        alive_callback, 
                        save_solution,
                        stepsize_callback,
                        chemistry_callback)

###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition=false),   #stage_limiter!, step_limiter!, 
            dt=1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep=false, callback=callbacks, maxiters=1e5);
summary_callback() # print the timer summary

PLOTTING = false 
@testset "extract R modes" begin 

  DESC.desc_jl_use_gpu_if_available()

  surface = DESC.jl_fourierRZToroidalSurface(
    R_lmn = [10, 1],
    modes_R = [[0, 0], [1, 0]],
    Z_lmn = [0, -1],
    modes_Z = [[0, 0], [-1, 0]],
  )

  knots = LinRange(0, 1, 20)
  pressure_values = zeros(20)
  iota_values = 1 .+ 1.5 * (knots.^2)

  pressure = DESC.jl_profiles_spline_profile(
    values = pressure_values, knots = knots
  )

  iota = DESC.jl_profiles_spline_profile(
    values = iota_values, 
    knots = knots
  )

  rtol = 1.0E-6
  @test isapprox(pressure._knots[length(pressure._knots) - 1], 0.94736842)
  
  eq = DESC.jl_equilibrium(
    surface=surface,
    pressure=pressure,
    iota=iota,
    Psi=1.0,  # flux (in Webers) within the last closed flux surface
    NFP=1,  # number of field periods
    L=6,  # radial spectral resolution
    M=6,  # poloidal spectral resolution
    N=0,  # toroidal spectral resolution (it's a tokamak, so we don't need any toroidal modes)
    L_grid=12,  # real space radial resolution, slightly oversampled
    M_grid=9,  # real space poloidal resolution, slightly oversampled
    N_grid=0,  # real space toroidal resolution
    sym=true,  # explicitly enforce stellarator symmetry
  )

  if PLOTTING
    DESC.jl_plotting_plot_section(eq, "|F|", norm_F = true, log=true)
  end 

  optimizer = DESC.jl_create_optimizer("lsq-exact")
  constraints = (
        DESC.jl_objective_fix_boundary_r(), 
        DESC.jl_objective_fix_boundary_z(), 
        DESC.jl_objective_fix_pressure(), 
        DESC.jl_objective_fix_iota(), 
        DESC.jl_objective_fix_psi()
  )

  objectives = DESC.jl_objective_force_balance()
  obj = DESC.jl_objective_function(objectives)


  DESC.jl_solve_equilibrium(eq, 
    verbose=2, ftol=1e-8, objective=obj, optimizer=optimizer, constraints=constraints
  )

  if PLOTTING
    DESC.jl_plotting_plot_surfaces(
      eq
    )

    DESC.jl_plotting_plot_section(
      eq, 
      "|F|", 
      norm_F=true, log=true
    )
    DESC.jl_plotting_plot_1d(eq, "p");
  end 

  delta_p = zeros(Int, size(eq.p_l))
  p_values = 1000 * (1 .- eq.pressure._knots .^2)
  delta_p = p_values 

  eq1 = DESC.jl_perturb(
    eq, 
    Dict("p_l" => delta_p), 
    order = 2
  )

  # TODO: FIX
  # @test isapprox(eq1.pressure.params, 1000 * (1 .- eq.pressure._knots .^ 2))

  if PLOTTING 
    DESC.jl_plotting_plot_section(eq1, "|F|", norm_F=true, log=true);
  end 
  
  constraints = (
    DESC.jl_objective_fix_boundary_r(), 
    DESC.jl_objective_fix_boundary_z(), 
    DESC.jl_objective_fix_pressure(), 
    DESC.jl_objective_fix_iota(), 
    DESC.jl_objective_fix_psi()
  )

  objective = DESC.jl_objective_force_balance()
  obj = DESC.jl_objective_function(objectives=objectives)
  DESC.jl_solve_equilibrium(eq, verbose=2, ftol = 1e-4, optimizer=optimizer, constraints=constraints, objective=obj)
  DESC.jl_plotting_plot_section(eq1, "|F|", norm_F=true, log=true)


#   println("Pressure values at these knots:\n", pressure.params)
#   println("Rotational Transform rho knots:\n", iota._knots)
#   println("Rotational Transform  values at these knots:\n", iota.params)

end 

PLOTTING = false 
@testset "extract R modes" begin 

  DESC.desc_jl_use_gpu_if_available()

  surface = DESC.FourierRZToroidalSurface(
    R_lmn = [10, 1],
    modes_R = [[0, 0], [1, 0]],
    Z_lmn = [0, -1],
    modes_Z = [[0, 0], [-1, 0]],
  )

  knots = LinRange(0, 1, 20)
  pressure_values = zeros(20)
  iota_values = 1 .+ 1.5 * (knots.^2)

  pressure = DESC.SplineProfile(
    values = pressure_values, knots = knots
  )

  iota = DESC.SplineProfile(
    values = iota_values, 
    knots = knots
  )

  rtol = 1.0E-6
  @test isapprox(pressure._knots[length(pressure._knots) - 1], 0.94736842)
  
  eq = DESC.Equilibrium(
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

    DESC.plot_section(eq, "|F|", norm_F = true, log=true, output_plots = false)

  optimizer = DESC.Optimizer("lsq-exact")
  constraints = (
        DESC.FixBoundaryR(), 
        DESC.FixBoundaryZ(), 
        DESC.FixPressure(), 
        DESC.FixIota(), 
        DESC.FixPsi()
  )

  objectives = DESC.ForceBalance()
  obj = DESC.ObjectiveFunction(objectives)


  eq.solve(
    verbose=2, ftol=1e-8, objective=obj, optimizer=optimizer, constraints=constraints
  )

    DESC.plot_surfaces(
      eq, 
      output_plots = false
    )

    DESC.plot_section(
      eq, 
      "|F|", 
      norm_F=true, log=true, output_plots = false
    )
    DESC.plot_1d(eq, "p", output_plots = false);

  delta_p = zeros(Int, size(eq.p_l))
  p_values = 1000 * (1 .- eq.pressure._knots .^2)
  delta_p = p_values 

  eq1 = DESC.equilibrium_perturb(
    eq, 
    Dict("p_l" => delta_p), 
    order = 2
  )


    DESC.plot_section(eq1, "|F|", norm_F=true, log=true, output_plots = false);
  
  constraints = (
    DESC.FixBoundaryR(), 
    DESC.FixBoundaryZ(), 
    DESC.FixPressure(), 
    DESC.FixIota(), 
    DESC.FixPsi()
  )

  objective = DESC.ForceBalance()
  obj = DESC.ObjectiveFunction(objectives=objectives)
  eq.solve(verbose=2, ftol = 1e-4, optimizer=optimizer, constraints=constraints, objective=obj)
  DESC.plot_section(eq1, "|F|", norm_F=true, log=true, output_plots = false)


#   println("Pressure values at these knots:\n", pressure.params)
#   println("Rotational Transform rho knots:\n", iota._knots)
#   println("Rotational Transform  values at these knots:\n", iota.params)

end 

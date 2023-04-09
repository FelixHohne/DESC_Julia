FULL_RUN = false

@testset "saveEquilibrium" begin 
  eq = DESC.jl_equilibrium()
  DESC.jl_save_equilibrium(eq, "test_eq.hdf5")
end 

@testset "loadEquilibrium" begin 
  eq = DESC.jl_load_equilibrium("test_eq.hdf5", "hdf5")
end 


@testset "Construct Linear Objectives" begin 
  objective = DESC.jl_objective_fix_boundary_r()
end 

@testset "constructEqFamily" begin 
  py"""
  import desc 
  import desc.equilibrium 

  eq = desc.equilibrium.Equilibrium()
  eq_fam = desc.equilibrium.EquilibriaFamily(eq)
  """
  eq = DESC.jl_equilibrium()
  eq_fam = DESC.jl_equilibria_family(eq)
end 

@testset "construct Linear Grid" begin 

  grid = DESC.jl_linear_grid(rho=[0.6, 0.8, 1.0], sym=true)

end 


@testset "extract R modes" begin 
    surf = DESC.jl_fourierRZToroidalSurface(
      [1, 0.125, 0.1],
      [-0.125, -0.1],
      [[0, 0], [1, 0], [0, 1]],
      [[-1, 0], [0, -1]],
      4
    )

    eq = DESC.jl_equilibrium(
        M = 8, 
        N = 8, 
        Psi=0.04, 
        surface=surf
    )

end 


@testset "set eq for QAS_output.h5" begin 
  if FULL_RUN 
    py"""
    import numpy as np
    import desc
    from desc import set_device
    set_device('gpu')
    """
    
    surf = DESC.jl_fourierRZToroidalSurface(
      [1, 0.125, 0.1],
      [-0.125, -0.1],
      [[0, 0], [1, 0], [0, 1]],
      [[-1, 0], [0, -1]],
      4
    )

    eq = DESC.jl_equilibrium(
        M = 8, 
        N = 8, 
        Psi=0.04, 
        surface=surf
    )

    eq = last(DESC.jl_solve_continuation_automatic(eq, objective = "force", verbose=3, bdry_step=0.5))
    DESC.jl_save_equilibrium(eq, "QAS_output_continuation_automatic.hdf5")
  end 
end 


@testset "QAS_output.h5" begin 
    eq = DESC.jl_load_equilibrium("QAS_output_continuation_automatic.hdf5", "hdf5")
    println(eq)
    eq_fam = DESC.jl_equilibria_family(eq)
    grid = DESC.jl_linear_grid(M=eq.M, N=eq.N, NFP=eq.NFP, rho=[0.6, 0.8, 1.0], sym=true)

    for n in 1:eq.M
      print(n)
      println(" Optimizing boundary modes with M, N <= %d\n", n); 
      objective = DESC.jl_objective_function(
        (DESC.jl_objective_quasisymmetry_two_term(helicity = (1, eq.NFP), grid=grid, normalize=false),
        DESC.jl_objective_aspect_ratio(target=8, weight=1e1, normalize=false)), 
        verbose = 0
      )

      R_abs = abs.(eq.surface.R_basis.modes)
      # println("R_abs")
      # display(R_abs)
      max_Rabs = maximum(R_abs, dims=1)
      # println("max Rabs")
      # display(max_Rabs)
      # println("Selection")
      R_elem_mode = eq.surface.R_basis.modes[findall(>(n), max_Rabs), :]
      R_modes = vcat([0, 0, 0], R_elem_mode)

      Z_abs = abs.(eq.surface.Z_basis.modes)
      max_Rabs = maximum(Z_abs, dims=1)
      Z_modes = eq.surface.Z_basis.modes[findall(>(n), max_Rabs), :]


  end 



end 
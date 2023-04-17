using PyCall
FULL_RUN = false

@testset "extract R modes" begin 
  DESC.desc_jl_use_gpu_if_available()
  surf = DESC.jl_fourierRZToroidalSurface(
    R_lmn = [1, 0.125, 0.1],
    Z_lmn = [-0.125, -0.1],
    modes_R = [[0, 0], [1, 0], [0, 1]],
    modes_Z = [[-1, 0], [0, -1]],
    NFP = 4
  )

  eq = DESC.jl_equilibrium(
      M = 8, 
      N = 8, 
      Psi=0.04, 
      surface=surf
  )

end 

@testset "Construct boolean R modes" begin 
  DESC.desc_jl_use_gpu_if_available()
  constraints = (
      DESC.jl_objective_fix_boundary_r(), 
      DESC.jl_objective_fix_boundary_z(), 
    )
end 



@testset "saveEquilibrium" begin 
  DESC.desc_jl_use_gpu_if_available()
  eq = DESC.jl_equilibrium()
  DESC.jl_save_equilibrium(eq, "test_eq.hdf5")
end 

@testset "loadEquilibrium" begin 
  DESC.desc_jl_use_gpu_if_available()
  eq = DESC.jl_load_equilibrium("test_eq.hdf5", "hdf5")
end 


@testset "constructEqFamily" begin 
  DESC.desc_jl_use_gpu_if_available()
  py"""
  import desc 
  import desc.equilibrium 

  eq = desc.equilibrium.Equilibrium()
  eq_fam = desc.equilibrium.EquilibriaFamily(eq)
  """
  eq = DESC.jl_equilibrium()
  eq_fam = DESC.jl_equilibria_family(eq)

  eq2 = DESC.jl_equilibrium()
  DESC.jl_equilibrium_family_append(eq_fam, eq2)

end 

@testset "construct Linear Grid" begin 
  DESC.desc_jl_use_gpu_if_available()
  grid = DESC.jl_linear_grid(rho=[0.6, 0.8, 1.0], sym=true)

end 

@testset "set eq for QAS_output.h5" begin 
  DESC.desc_jl_use_gpu_if_available()
  if FULL_RUN 
    py"""
    import numpy as np
    import desc

    """
    
    surf = DESC.jl_fourierRZToroidalSurface(
      R_lmn = [1, 0.125, 0.1],
      Z_lmn = [-0.125, -0.1],
      modes_R = [[0, 0], [1, 0], [0, 1]],
      modes_Z = [[-1, 0], [0, -1]],
      NFP = 4
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
  if FULL_RUN 
    DESC.desc_jl_use_gpu_if_available()
    eq = DESC.jl_load_equilibrium("QAS_output_continuation_automatic.hdf5", "hdf5")
    println(eq)
    eq_fam = DESC.jl_equilibria_family(eq)
    grid = DESC.jl_linear_grid(M=eq.M, N=eq.N, NFP=eq.NFP, rho=[0.6, 0.8, 1.0], sym=true)

      for n in 1:eq.M
        objective = DESC.jl_objective_function(
          (DESC.jl_objective_quasisymmetry_two_term(helicity = (1, eq.NFP), grid=grid, normalize=false),
          DESC.jl_objective_aspect_ratio(target=8, weight=1e1, normalize=false)), 
          verbose = 0
        )

        R_abs = abs.(eq.surface.R_basis.modes)
        max_Rabs = maximum(R_abs, dims=2)
        # equivalent to np.squeeze. See https://stackoverflow.com/questions/52505760/dropping-singleton-dimensions-in-julia
        # ensures indexing with findall broadcasts correctly
        max_Rabs = dropdims(max_Rabs, dims = tuple(findall(size(max_Rabs) .== 1)...))
        R_elem_mode = eq.surface.R_basis.modes[findall(>(n), max_Rabs), :]
        R_modes = vcat([0 0 0], R_elem_mode)
      
        Z_abs = abs.(eq.surface.Z_basis.modes)
        max_Zabs = maximum(Z_abs, dims=2)
        max_Zabs = dropdims(max_Zabs, dims = tuple(findall(size(max_Zabs) .== 1)...))
        Z_modes = eq.surface.Z_basis.modes[findall(>(n), max_Zabs), :]

        constraints = (
          DESC.jl_objective_force_balance(), 
          DESC.jl_objective_fix_boundary_r(modes=R_modes), 
          DESC.jl_objective_fix_boundary_z(modes=Z_modes), 
          DESC.jl_objective_fix_pressure(), 
          DESC.jl_objective_fix_current(), 
          DESC.jl_objective_fix_psi()
        )

        println("Constructed DESC constraints")

        optimizer = DESC.jl_create_optimizer("lsq-exact")

        println("Defined optimizer")

        println("Beginning equilibrium optimization")
        eq_new, out = DESC.jl_optimize_equilibrium(
          last(eq_fam); 
          objective=objective, 
          constraints=constraints, 
          optimizer=optimizer, 
          maxiter=1, 
          verbose=3, 
          copy=true, 
          options = Dict(
            "initial_trust_radius" => 0.5,
            "perturb_options" => Dict("verbose" => 0),
            "solve_options" => Dict("verbose" => 0)
          )
        )
        DESC.jl_equilibrium_family_append(eq_fam, eq_new)
        
    end 

    DESC.jl_save_equilibrium_family(eq_fam, "qas_julia_test_results.hdf5")
  end 
end 

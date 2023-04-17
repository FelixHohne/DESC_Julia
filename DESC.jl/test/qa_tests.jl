@testset "set eq for QAS_output.h5" begin 

    py"""
    import numpy as np
    import desc
   
    """
    
    surf = DESC.jl_fourierRZToroidalSurface(
      R_lmn = [1, 0.166, 0.1],
      Z_lmn = [-0.166, -0.1],
      modes_R = [[0, 0], [1, 0], [0, 1]],
      modes_Z = [[-1, 0], [0, -1]],
      NFP = 2
    )

    eq = DESC.jl_equilibrium(M=8, N = 8, Psi=0.087, surface=surf)

    eq = last(DESC.jl_solve_continuation_automatic(
        eq, 
        objective = "force", 
        bdry_step = 0.5, 
        verbose = 3
    ))

    DESC.jl_save_equilibrium(eq, "Qa_continuation_automatic.hdf5")

end 

@testset "Qa_output.h5" begin 
    eq = DESC.jl_load_equilibrium("Qa_continuation_automatic.hdf5", "hdf5")
    println(eq)
    eq_fam = DESC.jl_equilibria_family(eq)
    grid = DESC.jl_linear_grid(M=eq.M, N=eq.N, NFP=eq.NFP, rho=[0.6, 0.8, 1.0], sym=true)
  
      for n in 1:eq.M

        objective = DESC.jl_objective_function(
          (DESC.jl_objective_quasisymmetry_two_term(helicity = (1, 0), grid=grid, normalize=false),
          DESC.jl_objective_aspect_ratio(target=0.42, weight=1e1, normalize=false), 
          DESC.jl_objective_rotational_transform(target = 0.42, weight = 10, normalize = false)), 
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
  
    DESC.jl_save_equilibrium_family(eq_fam, "qa_julia_test_results.hdf5")
  
  end 
@testset "set eq for QAS_output.h5" begin 

    py"""
    import numpy as np
    import desc
   
    """
    
    surf = DESC.FourierRZToroidalSurface(
      R_lmn = [1, 0.166, 0.1],
      Z_lmn = [-0.166, -0.1],
      modes_R = [[0, 0], [1, 0], [0, 1]],
      modes_Z = [[-1, 0], [0, -1]],
      NFP = 2
    )

    eq = DESC.Equilibrium(M=8, N = 8, Psi=0.087, surface=surf)

    eq = last(DESC.solve_continuation_automatic(
        eq, 
        objective = "force", 
        bdry_step = 0.5, 
        verbose = 3
    ))

    eq.save("Qa_continuation_automatic.hdf5")

end 

@testset "Qa_output.h5" begin 
    eq = DESC.equilibrium_load("Qa_continuation_automatic.hdf5", "hdf5")
    println(eq)
    eq_fam = DESC.EquilibriaFamily(eq)
    grid = DESC.LinearGrid(M=eq.M, N=eq.N, NFP=eq.NFP, rho=[0.6, 0.8, 1.0], sym=true)
  
      for n in 1:eq.M

        objective = DESC.ObjectiveFunction(
          (DESC.QuasisymmetryTwoTerm(helicity = (1, 0), grid=grid, normalize=false),
          DESC.AspectRatio(target=0.42, weight=1e1, normalize=false), 
          DESC.RotationalTransform(target = 0.42, weight = 10, normalize = false)), 
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
          DESC.ForceBalance(), 
          DESC.FixBoundaryR(modes=R_modes), 
          DESC.FixBoundaryZ(modes=Z_modes), 
          DESC.FixPressure(), 
          DESC.FixCurrent(), 
          DESC.FixPsi()
        )
  
        println("Constructed DESC constraints")
  
        optimizer = DESC.Optimizer("lsq-exact")
  
        println("Defined optimizer")
  
        println("Beginning equilibrium optimization")
        eq_new, out = DESC.equilibrium_optimize(
          eq, 
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
  
    eq_fam.save("qa_julia_test_results.hdf5")
  
  end 
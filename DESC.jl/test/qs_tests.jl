using PyCall
FULL_RUN = false

@testset "extract R modes" begin 
  DESC.desc_jl_use_gpu_if_available()
  surf = DESC.FourierRZToroidalSurface(
    R_lmn = [1, 0.125, 0.1],
    Z_lmn = [-0.125, -0.1],
    modes_R = [[0, 0], [1, 0], [0, 1]],
    modes_Z = [[-1, 0], [0, -1]],
    NFP = 4
  )

  eq = DESC.Equilibrium(
      M = 8, 
      N = 8, 
      Psi=0.04, 
      surface=surf
  )

end 

@testset "Construct boolean R modes" begin 
  DESC.desc_jl_use_gpu_if_available()
  constraints = (
      DESC.FixBoundaryR(), 
      DESC.FixBoundaryZ(), 
    )
end 



@testset "saveEquilibrium" begin 
  DESC.desc_jl_use_gpu_if_available()
  eq = DESC.Equilibrium()
  eq.save("test_eq.hdf5")
end 

@testset "loadEquilibrium" begin 
  DESC.desc_jl_use_gpu_if_available()
  eq = DESC.equilibrium_load("test_eq.hdf5", "hdf5")
end 


@testset "constructEqFamily" begin 
  DESC.desc_jl_use_gpu_if_available()
  py"""
  import desc 
  import desc.equilibrium 

  eq = desc.equilibrium.Equilibrium()
  eq_fam = desc.equilibrium.EquilibriaFamily(eq)
  """
  eq = DESC.Equilibrium()
  eq_fam = DESC.EquilibriaFamily(eq)

  eq2 = DESC.Equilibrium()
  DESC.equilibrium_family_append(eq_fam, eq2)

end 

@testset "construct Linear Grid" begin 
  DESC.desc_jl_use_gpu_if_available()
  grid = DESC.LinearGrid(rho=[0.6, 0.8, 1.0], sym=true)

end 

@testset "set eq for QAS_output.h5" begin 
  DESC.desc_jl_use_gpu_if_available()
  if FULL_RUN 
    py"""
    import numpy as np
    import desc

    """
    
    surf = DESC.FourierRZToroidalSurface(
      R_lmn = [1, 0.125, 0.1],
      Z_lmn = [-0.125, -0.1],
      modes_R = [[0, 0], [1, 0], [0, 1]],
      modes_Z = [[-1, 0], [0, -1]],
      NFP = 4
    )

    eq = DESC.Equilibrium(
        M = 8, 
        N = 8, 
        Psi=0.04, 
        surface=surf
    )

    eq = last(DESC.solve_continuation_automatic(eq, objective = "force", verbose=3, bdry_step=0.5))
    eq.save("QAS_output_continuation_automatic.hdf5")
  end 
end 


@testset "QAS_output.h5" begin 
  if FULL_RUN 
    DESC.desc_jl_use_gpu_if_available()
    eq = DESC.equilibrium_load("QAS_output_continuation_automatic.hdf5", "hdf5")
    println(eq)
    eq_fam = DESC.EquilibriaFamily(eq)
    grid = DESC.LinearGrid(M=eq.M, N=eq.N, NFP=eq.NFP, rho=[0.6, 0.8, 1.0], sym=true)

      for n in 1:eq.M
        objective = DESC.ObjectiveFunction(
          (DESC.QuasisymmetryTwoTerm(helicity = (1, eq.NFP), grid=grid, normalize=false),
          DESC.AspectRatio(target=8, weight=1e1, normalize=false)), 
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
        DESC.equilibrium_family_append(eq_fam, eq_new)
        
    end 

    eq_fam.save("qas_julia_test_results.hdf5")
  end 
end 

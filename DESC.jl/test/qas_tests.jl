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



@testset "QAS_output.h5" begin 
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
    println(eq)

    eq_fam = DESC.jl_equilibria_family(eq)







  end 
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
    eq_fam = DESC.jl_equilibria_family(eq)

    grid = LinearGrid(M=eq.M, N=eq.N, NFP=eq.NFP, rho=[0.6, 0.8, 1.0], sym=true)

    for n in 1:eq.M
      print(n)
      println(" Optimizing boundary modes with M, N <= %d\n", n); 
      objective = DESC.jl_objective_function(
        DESC.jl_objective_quasisymmetry_two_term(helicity = (1, eq.NFP), grid=grid, normalize=false),
        DESC.jl_objective_aspect_ratio(target=8, weight=1e1, normalize=false)
        verbose = 0
      )
    end 








  end 
@testset "Simple Continuation Test" begin 

eq = DESC.Equilibrium()
eq_fam = DESC.EquilibriaFamily(eq)

new_eq_fam = DESC.solve_continuation(
    eq_fam, 
    objective = "vacuum", 
    optimizer = "lsq-exact", 
    ftol = 1e-8, 
    gtol = 1e-4, 
    verbose = 3
)

surf = DESC.FourierRZToroidalSurface(
    R_lmn = [1, 0.125, 0.1],
    Z_lmn = [-0.125, -0.1],
    modes_R = [[0, 0], [1, 0], [0, 1]],
    modes_Z = [[-1, 0], [0, -1]],
    NFP = 4
  )

  new_eq = DESC.Equilibrium(
      M = 9, 
      N = 8, 
      Psi=0.04, 
      surface=surf
  )

new_eq_fam_v2 = DESC.solve_continuation_automatic(
    new_eq, 
    objective = "vacuum", 
    optimizer = "lsq-exact", 
    ftol = 1e-2, 
    nfev=2,
    gtol = 1e-2, 
    verbose = 3
)

end 

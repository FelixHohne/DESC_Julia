using Pkg
Pkg.activate(".") 
using DESC 

function code_to_time(use_gpu)
    if use_gpu
        DESC.desc_jl_use_gpu_if_available()
    end 
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

    DESC.jl_solve_continuation_automatic(eq, objective = "force", verbose=3, bdry_step=0.5)
end 


@time code_to_time(false)

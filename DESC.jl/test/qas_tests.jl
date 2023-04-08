@testset "QAS_output.h5" begin 
    py"""
    import numpy as np
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

  end 
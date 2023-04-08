using PyCall


function jl_linear_grid(;
    L = nothing, 
    M = nothing, 
    N = nothing, 
    NFP = 1, 
    sym = false, 
    axis = true, 
    endpoint = false, 
    rho = [1.0], 
    theta = [0.0], 
    zeta = [0.0]
)

    py"""
    import desc 
    import desc.grid
    import numpy as np  

    def create_linear_grid():
        return desc.grid.LinearGrid(
            L = $L, 
            M = $M, 
            N = $N, 
            NFP = $NFP, 
            sym = $sym, 
            axis = $axis, 
            endpoint = $endpoint, 
            rho = $rho, 
            theta = $theta, 
            zeta = $zeta
        )
    """
    output = py"create_linear_grid"()

end 



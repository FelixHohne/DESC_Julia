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

function jl_quadrature_grid(;
    L, 
    M, 
    N; 
    NFP = 1
)
    py"""
    import desc 
    import desc.grid
    import numpy as np  

    def create_quadrature_grid():
        return desc.grid.QuadratureGrid(
            L = $L, 
            M = $M, 
            N = $N, 
            NFP = $NFP
        )
    """
    output = py"create_quadrature_grid"()
end 

function jl_concentric_grid(;
    L, 
    M, 
    N; 
    NFP = 1, 
    sym=false, 
    axis=false, 
    node_pattern = "jacobi"
)
    py"""
    import desc 
    import desc.grid
    import numpy as np  

    def create_concentric_grid():
        return desc.grid.ConcentricGrid(
            L = $L, 
            M = $M, 
            N = $N, 
            NFP = $NFP, 
            sym=$sym, 
            axis=$axis, 
            node_pattern = $node_pattern
        )
    """
    output = py"create_concentric_grid"()
end 




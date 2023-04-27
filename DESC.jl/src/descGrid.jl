using PyCall


function Grid(
    nodes;
    sort = true
)
    py"""
    import desc 
    import desc.grid
    import numpy as np  

    if isinstance($nodes, np.ndarray):
        n_nodes = np.ascontiguousarray($nodes)
        assert n_nodes.flags['C_CONTIGUOUS']
    else:
        n_nodes = $nodes

    def create_grid():
        return desc.grid.Grid(
            n_nodes, 
            sort = $sort
        )
    """
    output = py"create_grid"()
end 

function LinearGrid(;
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

    if isinstance($rho, np.ndarray):
        n_rho = np.ascontiguousarray($rho)
        assert n_rho.flags['C_CONTIGUOUS']
    else:
        n_rho = $rho

    if isinstance($theta, np.ndarray):
        n_theta = np.ascontiguousarray($theta)
        assert n_theta.flags['C_CONTIGUOUS']
    else:
        n_theta = $theta

    if isinstance($zeta, np.ndarray):
        n_zeta = np.ascontiguousarray($zeta)
        assert n_zeta.flags['C_CONTIGUOUS']
    else:
        n_zeta = $zeta

    def create_linear_grid():
        return desc.grid.LinearGrid(
            L = $L, 
            M = $M, 
            N = $N, 
            NFP = $NFP, 
            sym = $sym, 
            axis = $axis, 
            endpoint = $endpoint, 
            rho = n_rho, 
            theta = n_theta, 
            zeta = n_zeta
        )
    """
    output = py"create_linear_grid"()

end 

function QuadratureGrid(
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

function ConcentricGrid(
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






function QuadratureGrid(
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

function grid_load(
    load_from;
    file_format = nothing
)
    py"""
    import desc 

    def load():
        return desc.grid.load(
           $load_from, 
           file_format = $file_format
        )
    """
    output = py"load"()
end 




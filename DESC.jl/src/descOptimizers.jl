using PyCall


#=
Functions implement desc.optimize API
=#

function Optimizer(method) 
    py"""
    def create_optimizer():
        return desc.optimize.Optimizer(
            $method, 
        )
    """
    output = py"create_optimizer"()
end 

function optimize_load(
    load_from, 
    file_format = nothing 
)
py"""
import desc
import desc.optimize 

def optimize_load():
    return desc.optimize.Optimizer.load(
        $load_from, 
        file_format = $file_format
    )
"""
result = py"optimize_load"()
end 


function optimize_optimize(
    optimizer, 
    eq, 
    objective;
    constraints = (), 
    ftol = nothing, 
    xtol = nothing, 
    gtol = nothing, 
    x_scale = "auto", 
    verbose = 1, 
    maxiter = nothing, 
    options = Dict())

    py"""
    import desc
    import desc.optimize

    if isinstance($x_scale, np.ndarray):
        new_x_scale = np.ascontiguousarray($x_scale)
        assert new_x_scale.flags['C_CONTIGUOUS']
    else:
        new_x_scale = $x_scale

    def execute_optimize_command():
        $optimizer.optimize(
            $eq, 
            $objective, 
            constraints = $constraints, 
            ftol = $ftol, 
            xtol = $xtol, 
            gtol = $gtol, 
            x_scale = new_x_scale, 
            verbose = $verbose, 
            maxiter = $maxiter, 
            options = $options
        )
    """
    optimize_result = py"execute_optimize_command"()
end 

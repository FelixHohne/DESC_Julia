using PyCall


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

# function jl_optimize_fmintr(
#     fun, 
#     x0, 
#     grad;
#     hess = "bfs", 
#     args = (), 
#     method = "dogleg", 
#     x_scale = 1, 
#     ftol = 1e-6, 
#     xtol = 1e-6, 
#     gtol = 1e-6, 
#     verbose = 1, 
#     maxiter = nothing, 
#     callback = nothing, 
#     options = Dict()
# )

# py"""
# import numpy as np
# import desc

# def fmintr():
#     return desc.optimize.fmintr(
#         fun = $fun, 
#         x0 = $x0, 
#         grad = $grad, 
#         hess = $hess, 
#         args = $args, 
#         method = $method, 
#         x_scale = $x_scale, 
#         ftol = $ftol, 
#         xtol = $xtol, 
#         gtol = $gtol, 
#         verbose = $verbose, 
#         maxiter = $maxiter, 
#         callback = $callback, 
#         options = $options 
#     )
# """
# result = py"fmintr"()

# end 

# function jl_optimize_lsqtr(
#     fun, 
#     x0, 
#     jac; 
#     args = (), 
#     x_scale = 1, 
#     ftol = 1e-6, 
#     xtol = 1e-6, 
#     gtol = 1e-6, 
#     verbose = 1, 
#     maxiter = nothing, 
#     tr_method = "svd", 
#     callback = nothing, 
#     options = Dict()
# )

# py"""
# import numpy as np
# import desc

# def lsqtr():
#     return desc.optimize.lsqtr(
#         fun = $fun, 
#         x0 = $x0, 
#         jac = $jac; 
#         args = $args, 
#         x_scale = $x_scale, 
#         ftol = $ftol, 
#         xtol = $xtol, 
#         gtol = $gtol, 
#         verbose = $verbose, 
#         maxiter = $maxiter, 
#         tr_method = $tr_method, 
#         callback = $callback, 
#         options = $options
#     )
# """
# result = py"lsqtr"()

# end 
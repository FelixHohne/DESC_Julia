using PyCall


function jl_create_optimizer(method) 
    py"""
    def create_optimizer():
        return desc.optimize.Optimizer(
            $method, 
        )
    """
    output = py"create_optimizer"()
end 


function jl_create_example_constraints() 

    py"""
    from desc.objectives import (
        get_fixed_boundary_constraints,
        ObjectiveFunction,
        FixBoundaryR,
        FixBoundaryZ,
        FixLambdaGauge,
        FixPressure,
        FixIota,
        FixPsi,
        ForceBalance,
    )
    def create_constraints():
        constraints = (
            FixBoundaryR(fixed_boundary=True),  # enforce fixed  LCFS for R
            FixBoundaryZ(fixed_boundary=True),  # enforce fixed  LCFS for R
            FixLambdaGauge(),  # Fix the gauge for Lambda (in stellarator symmetric cases, this sets lambda to zero at the magnetic axis)
            FixPressure(),  # enforce that the pressure profile stay fixed
            FixIota(),  # enforce that the rotational transform profile stay fixed
            FixPsi(),  # enforce that the enclosed toroidal stay fixed
        ) 
        return constraints  
    """
    constraints = py"create_constraints"() 
end 

function jl_optimize(
    optimizer, 
    equilibrium, 
    objective, 
    constraints = (), 
    ftol = nothing, 
    xtol = nothing, 
    gtol = nothing, 
    x_scale = "auto", 
    verbose = 1, 
    maxiter = nothing, 
    options = Dict())

    py"""
    from desc.equilibrium import Equilibrium
    def execute_optimize_command():
        $optimizer.optimize(
            $equilibrium, 
            $objective, 
            $constraints, 
            $ftol, 
            $xtol, 
            $gtol, 
            $x_scale, 
            $verbose, 
            $maxiter, 
            $options
        )
    """
    optimize_result = py"execute_optimize_command"()
end 


function jl_optimize_equilibrium(
    eq; 
    objective=nothing,
    constraints= nothing,
    optimizer="proximal-lsq-exact",
    ftol= nothing,
    xtol= nothing,
    gtol= nothing,
    maxiter=50,
    x_scale="auto",
    options=Dict(),
    verbose=1,
    copy=false
)

    py"""
    import numpy as np
    import desc
   
    def optimize():
        return $eq.optimize(
            objective=$objective,
            constraints=$constraints,
            optimizer=$optimizer,
            ftol=$ftol,
            xtol=$xtol,
            gtol=$gtol,
            maxiter=$maxiter,
            x_scale=$x_scale,
            options=$options,
            verbose=$verbose,
            copy=$copy
        )
    """
    result = py"optimize"()
end

function jl_optimize_fmintr(
    fun, 
    x0, 
    grad;
    hess = "bfs", 
    args = (), 
    method = "dogleg", 
    x_scale = 1, 
    ftol = 1e-6, 
    xtol = 1e-6, 
    gtol = 1e-6, 
    verbose = 1, 
    maxiter = nothing, 
    callback = nothing, 
    options = Dict()
)

py"""
import numpy as np
import desc

def fmintr():
    return desc.optimize.fmintr(
        fun = $fun, 
        x0 = $x0, 
        grad = $grad, 
        hess = $hess, 
        args = $args, 
        method = $method, 
        x_scale = $x_scale, 
        ftol = $ftol, 
        xtol = $xtol, 
        gtol = $gtol, 
        verbose = $verbose, 
        maxiter = $maxiter, 
        callback = $callback, 
        options = $options 
    )
"""
result = py"fmintr"()

end 

function jl_optimize_lsqtr(
    fun, 
    x0, 
    jac; 
    args = (), 
    x_scale = 1, 
    ftol = 1e-6, 
    xtol = 1e-6, 
    gtol = 1e-6, 
    verbose = 1, 
    maxiter = nothing, 
    tr_method = "svd", 
    callback = nothing, 
    options = Dict()
)

py"""
import numpy as np
import desc

def lsqtr():
    return desc.optimize.lsqtr(
        fun = $fun, 
        x0 = $x0, 
        jac = $jac; 
        args = $args, 
        x_scale = $x_scale, 
        ftol = $ftol, 
        xtol = $xtol, 
        gtol = $gtol, 
        verbose = $verbose, 
        maxiter = $maxiter, 
        tr_method = $tr_method, 
        callback = $callback, 
        options = $options
    )
"""
result = py"lsqtr"()

end 
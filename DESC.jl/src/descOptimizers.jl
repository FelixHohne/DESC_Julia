using PyCall

function jl_objective_function(
    objectives;
    eq=nothing, 
    use_jit=true, 
    deriv_mode="batched", 
    verbose=1)
    py"""
    def create_obj_function():
        return desc.objectives.ObjectiveFunction(
            objectives=$objectives,
            eq=$eq, 
            use_jit=$use_jit, 
            deriv_mode=$deriv_mode, 
            verbose=$verbose
    )
    """
    output = py"create_obj_function"()

end 

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

function jl_optimize(optimizer, equilibrium, objective, constraints = (), ftol = nothing, xtol = nothing, gtol = nothing, x_scale = "auto", verbose = 1, maxiter = nothing, options = Dict())

    py"""
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

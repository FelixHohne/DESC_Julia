using PyCall
function jl_equilibrium(descData)
    py"""
    import numpy as np
    import desc
    
    def create_eq_obj():
        eq = desc.equilibrium.Equilibrium(
            Psi = $descData.ψ, 
            NFP = $descData.nfp, 
            L_grid = $descData.L_grid, 
            M_grid = $descData.M_grid, 
            N_grid = $descData.N_grid, 
            # pressure = $descData.pres, 
            # current = $descData.current, 
        )
        return eq
    """
    output = py"create_eq_obj"()
end 



function jl_objective_aspect_ratio(equilibrium::PyObject, target=2, weight=1, normalize=true, normalize_target=true, name="aspect_ratio")
    py"""
    import numpy as np
    import desc
    def create_objective_aspect_ratio():
        # Not supporting grid
        return desc.objectives.AspectRatio(
            eq=$equilibrium, target=$target, weight=$weight, 
            normalize=$normalize, normalize_target=$normalize_target, 
            name=$name)
    """
    output = py"create_objective_aspect_ratio"()
end 


function jl_objective_current_density(equilibrium, target = 0, bounds = nothing, weight = 1, normalize = true, normalize_target = true, grid = nothing, name = "current density")
    py"""
    import numpy as np
    import desc
    def create_objective_current_density():
        # Not supporting grid
        return desc.objectives.CurrentDensity(
            eq=$equilibrium, target=$target, weight=$weight, 
            normalize=$normalize, normalize_target=$normalize_target, 
            name=$name)
    """
    output = py"create_objective_aspect_ratio"()
end 




function jl_objective_function(objective, equilibrium, use_jit = true, deriv_mode = "batched", verbose = 1)
    py"""
    def create_obj_function():
        # Currently only support one objective
        return desc.objectives.ObjectiveFunction(
            objectives=($objective),
            eq=$equilibrium, 
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
            # $objectiveFunction
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

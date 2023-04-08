using PyCall

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

function jl_objective_fix_boundary_r(;
    eq = nothing, 
    target = nothing,  
    bounds = nothing, 
    weight = 1, 
    normalize = true, 
    normalize_target = true, 
    fixed_boundary = false, 
    modes = true, 
    surface_label = nothing, 
    name = "lcfs R"
)

    py"""
    import numpy as np
    import desc
    from desc.objectives.linear_objectives import *
    def create_objective_fix_boundary_r():
        return FixBoundaryR(
            eq=$eq, target=$target, weight=$weight, 
            normalize=$normalize, normalize_target=$normalize_target, 
            name=$name)
    """
    output = py"create_objective_fix_boundary_r"()
end 




function jl_objective_bootstrapRedlConsistency(
    equilibrium::PyObject, 
    target = 0, 
    bounds = nothing, 
    weight = 1, 
    normalize = true, 
    normalize_target = true, 
    grid = nothing, 
    helicity = (1, 0), 
    name = "Bootstrap current self-consistency (Redl)"
)
    py"""
    import numpy as np 
    import desc 
    def createBootstrapRedlConsistency():
        return desc.objectives.BootstrapRedlConsistency(
            eq=$equilibrium, target=$target, bounds=$bounds, weight=$weight, 
            normalize=$normalize, normalize_target=$normalize_target, 
            grid=$grid, helicity=$helicity, name=$name
        )
    output = py"createBootstrapRedlConsistency"()
    """
end 


function jl_objective_current_density(equilibrium, target = 0, bounds = nothing, weight = 1, normalize = true, normalize_target = true, grid = nothing, name = "current density")
    py"""
    import numpy as np
    import desc
    # TODO: Grid
    def create_objective_current_density():
        # Not supporting grid
        return desc.objectives.CurrentDensity(
            eq=$equilibrium, target=$target, weight=$weight, 
            normalize=$normalize, normalize_target=$normalize_target, 
            name=$name)
    """
    output = py"create_objective_aspect_ratio"()
end 


function jl_objective_elongation(
    equilibrium, 
    target = 1, 
    bounds = nothing, 
    weight = 1, 
    normalize = true, 
    normalize_target = true, 
    grid = nothing, 
    name = "elongation"
)
    py"""
    import numpy as np
    import desc
    # TODO: Grid
    def create_objective_elongation():
        # Not supporting grid
        return desc.objectives.Elongation(
            eq=$equilibrium, target=$target, bounds = $bounds, weight=$weight, 
            normalize=$normalize, normalize_target=$normalize_target, 
            grid=$grid, 
            name=$name)
    """
    output = py"create_objective_elongation"()

end 




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

function jl_objective_quasisymmetry_two_term(;
    eq, 
    target=0,  
    bounds = nothing, 
    weight=1, 
    normalize=true, 
    normalize_target=true, 
    grid = nothing, 
    helicity = (1, 0), 
    name="QS two-term"
)

    py"""
    import numpy as np
    import desc
    def create_objective_quasisymmetry_two_term():
        # Not supporting grid
        return desc.objectives.QuasisymmetryTwoTerm(
            eq=$equilibrium, target=$target, weight=$weight, 
            normalize=$normalize, normalize_target=$normalize_target, 
            name=$name)
    """
    output = py"create_objective_quasisymmetry_two_term"()
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


function jl_objective_fix_boundary_z(;
    eq = nothing, 
    target = nothing,  
    bounds = nothing, 
    weight = 1, 
    normalize = true, 
    normalize_target = true, 
    fixed_boundary = false, 
    modes = true, 
    surface_label = nothing, 
    name = "lcfs Z"
)

    py"""
    import numpy as np
    import desc
    from desc.objectives.linear_objectives import *
    def create_objective_fix_boundary_z():
        return FixBoundaryR(
            eq=$eq, target=$target, weight=$weight, 
            normalize=$normalize, normalize_target=$normalize_target, 
            name=$name)
    """
    output = py"create_objective_fix_boundary_z"()
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


function jl_objective_current_density(
    equilibrium = nothing, 
    target = 0, 
    bounds = nothing, 
    weight = 1, 
    normalize = true, 
    normalize_target = true, 
    grid = nothing, 
    name = "current density")
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
    output = py"create_objective_current_density"()
end 


function jl_objective_elongation(
    equilibrium = nothing, 
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


function jl_objective_energy(;
    eq = nothing, 
    target = 0, 
    bounds = nothing, 
    weight = 1, 
    normalize = true, 
    normalize_target = true, 
    grid = nothing, 
    gamma = 0, 
    name = "energy"
)
    py"""
    import numpy as np
    import desc
    # TODO: Grid
    def create_objective_energy():
        # Not supporting grid
        return desc.objectives.Energy(
            eq=$eq, target=$target, bounds = $bounds, weight=$weight, 
            normalize=$normalize, normalize_target=$normalize_target, 
            grid=$grid, gamma=$gamma, 
            name=$name)
    """
    output = py"create_objective_energy"()

end 


function jl_objective_force_balance(;
    eq = nothing,  
    target = 0, 
    bounds = nothing, 
    weight = 1, 
    normalize = true, 
    normalize_target = true, 
    grid = nothing, 
    name = "force"
)
    py"""
    import numpy as np
    import desc
    # TODO: Grid
    def create_objective_force_balance():
        # Not supporting grid
        return desc.objectives.ForceBalance(
            eq=$eq, target=$target, bounds = $bounds, weight=$weight, 
            normalize=$normalize, normalize_target=$normalize_target, 
            grid=$grid, 
            name=$name)
    """
    output = py"create_objective_force_balance"()
end 


function jl_objective_generic_objective(;
    f; 
    eq = nothing,  
    target = 0, 
    bounds = nothing, 
    weight = 1, 
    normalize = true, 
    normalize_target = true, 
    grid = nothing, 
    name = "generic"
)
    py"""
    import numpy as np
    import desc
    # TODO: Grid
    def create_generic_objective():
        # Not supporting grid
        return desc.objectives.GenericObjective(
            f=$f, eq=$eq, target=$target, bounds = $bounds, weight=$weight, 
            normalize=$normalize, normalize_target=$normalize_target, 
            grid=$grid, name=$name)
    """
    output = py"create_generic_objective"()
end 

function jl_objective_helical_force_balance(;
    eq = nothing,  
    target = 0, 
    bounds = nothing, 
    weight = 1, 
    normalize = true, 
    normalize_target = true, 
    grid = nothing, 
    name = "helical force"
)
    py"""
    import numpy as np
    import desc
    # TODO: Grid
    def create_helical_force_objective():
        # Not supporting grid
        return desc.objectives.HelicalForceBalance(
            eq=$eq, target=$target, bounds = $bounds, weight=$weight, 
            normalize=$normalize, normalize_target=$normalize_target, 
            grid=$grid, name=$name)
    """
    output = py"create_helical_force_objective"()
end 


function jl_objective_magnetic_well(;
    eq = nothing,  
    target = nothing, 
    bounds = (0, Inf), 
    weight = 1, 
    normalize = true, 
    normalize_target = true, 
    grid = nothing, 
    name = "Magnetic Well"
)
    py"""
    import numpy as np
    import desc
    # TODO: Grid
    def create_magnetic_well_objective():
        # Not supporting grid
        return desc.objectives.MagneticWelle(
            eq=$eq, target=$target, bounds = $bounds, weight=$weight, 
            normalize=$normalize, normalize_target=$normalize_target, 
            grid=$grid, name=$name)
    """
    output = py"create_magnetic_well_objective"()
end 

function jl_objective_mean_curvature(;
    eq = nothing,  
    target = nothing, 
    bounds = (-Inf, 0), 
    weight = 1, 
    normalize = true, 
    normalize_target = true, 
    grid = nothing, 
    name = "mean-curvature"
)
    py"""
    import numpy as np
    import desc
    # TODO: Grid
    def create_mean_curvature_objective():
        # Not supporting grid
        return desc.objectives.MeanCurvature(
            eq=$eq, target=$target, bounds = $bounds, weight=$weight, 
            normalize=$normalize, normalize_target=$normalize_target, 
            grid=$grid, name=$name)
    """
    output = py"create_mean_curvature_objective"()
end 


function jl_objective_mercier_stability(;
    eq = nothing,  
    target = nothing, 
    bounds = (0, Inf), 
    weight = 1, 
    normalize = true, 
    normalize_target = true, 
    grid = nothing, 
    name = "Mercier Stability"
)
    py"""
    import numpy as np
    import desc
    # TODO: Grid
    def create_mercier_stability_objective():
        # Not supporting grid
        return desc.objectives.MercierStability(
            eq=$eq, target=$target, bounds = $bounds, weight=$weight, 
            normalize=$normalize, normalize_target=$normalize_target, 
            grid=$grid, name=$name)
    """
    output = py"create_mercier_stability_objective"()
end 











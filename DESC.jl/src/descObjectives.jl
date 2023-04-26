using PyCall

function AspectRatio(;
    eq=nothing,
    target=2,
    bounds=nothing,
    weight=1,
    normalize=true,
    normalize_target=true,
    grid=nothing,
    name="aspect ratio",
)

    py"""
    import numpy as np
    import desc
    import desc.objectives
    def create_objective_aspect_ratio():
        return desc.objectives.AspectRatio(
            eq=$eq,
            target=$target,
            bounds=$bounds,
            weight=$weight,
            normalize=$normalize,
            normalize_target=$normalize_target,
            grid = $grid, 
            name=$name)
    """
    output = py"create_objective_aspect_ratio"()
end 

function BootstrapRedlConsistency(
    eq = nothing, 
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


function CurrentDensity(
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
    def create_objective_current_density():
        return desc.objectives.CurrentDensity(
            eq=$equilibrium, target=$target, weight=$weight, 
            normalize=$normalize, normalize_target=$normalize_target, 
            name=$name)
    """
    output = py"create_objective_current_density"()
end 


function Elongation(
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
    def create_objective_elongation():
        return desc.objectives.Elongation(
            eq=$equilibrium, target=$target, bounds = $bounds, weight=$weight, 
            normalize=$normalize, normalize_target=$normalize_target, 
            grid=$grid, 
            name=$name)
    """
    output = py"create_objective_elongation"()

end 


function Energy(;
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
    def create_objective_energy():
        return desc.objectives.Energy(
            eq=$eq, target=$target, bounds = $bounds, weight=$weight, 
            normalize=$normalize, normalize_target=$normalize_target, 
            grid=$grid, gamma=$gamma, 
            name=$name)
    """
    output = py"create_objective_energy"()

end 


function ForceBalance(;
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
    def create_objective_force_balance():
        return desc.objectives.ForceBalance(
            eq=$eq, target=$target, bounds = $bounds, weight=$weight, 
            normalize=$normalize, normalize_target=$normalize_target, 
            grid=$grid, 
            name=$name)
    """
    output = py"create_objective_force_balance"()
end 


function GenericObjective(
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
    def create_generic_objective():
        return desc.objectives.GenericObjective(
            f=$f, eq=$eq, target=$target, bounds = $bounds, weight=$weight, 
            normalize=$normalize, normalize_target=$normalize_target, 
            grid=$grid, name=$name)
    """
    output = py"create_generic_objective"()
end 

function HelicalForceBalance(;
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
    def create_helical_force_objective():
        return desc.objectives.HelicalForceBalance(
            eq=$eq, target=$target, bounds = $bounds, weight=$weight, 
            normalize=$normalize, normalize_target=$normalize_target, 
            grid=$grid, name=$name)
    """
    output = py"create_helical_force_objective"()
end 

function Isodynamicity(;
    eq = nothing,  
    target = 0, 
    bounds = nothing, 
    weight = 1, 
    normalize = false, 
    normalize_target = false, 
    grid = nothing, 
    name = "Isodynamicity"
)
    py"""
    import numpy as np
    import desc
    def create_isodynamicity_objective():
        return desc.objectives.Isodynamicity(
            eq=$eq,
            target=$target,
            bounds=$bounds,
            weight=$weight,
            normalize=$normalize,
            normalize_target=$normalize_target,
            name=$name,
        )
    """
    output = py"create_isodynamicity_objective"()
end 

function MagneticWell(;
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
    def create_magnetic_well_objective():
        return desc.objectives.MagneticWell(
            eq=$eq, target=$target, bounds = $bounds, weight=$weight, 
            normalize=$normalize, normalize_target=$normalize_target, 
            grid=$grid, name=$name)
    """
    output = py"create_magnetic_well_objective"()
end 

function MeanCurvature((;
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
    def create_mean_curvature_objective():
        return desc.objectives.MeanCurvature(
            eq=$eq, target=$target, bounds = $bounds, weight=$weight, 
            normalize=$normalize, normalize_target=$normalize_target, 
            grid=$grid, name=$name)
    """
    output = py"create_mean_curvature_objective"()
end 


function MercierStability(;
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
    def create_mercier_stability_objective():
        return desc.objectives.MercierStability(
            eq=$eq, target=$target, bounds = $bounds, weight=$weight, 
            normalize=$normalize, normalize_target=$normalize_target, 
            grid=$grid, name=$name)
    """
    output = py"create_mercier_stability_objective"()
end 

function ObjectiveFunction(;
    objectives = nothing, 
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


# The objectives field is required, so do not pass in as keyword argument. 
function ObjectiveFunction(
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

function PlasmaVesselDistance(;
    eq = nothing,  
    target = nothing, 
    bounds = (1, Inf), 
    weight = 1, 
    normalize = true, 
    normalize_target = true, 
    surface_grid = nothing, 
    plasma_grid = nothing, 
    name = "plasma vessel distance"
)
    py"""
    import numpy as np
    import desc
    def create_plasma_vessel_distance_objective():
        return desc.objectives.PlasmaVesselDistance(
            eq=$eq, target=$target, bounds = $bounds, weight=$weight, 
            normalize=$normalize, normalize_target=$normalize_target, 
            surface_grid=$surface_grid, plasma_grid=$plasma_grid, 
            name=$name)
    """
    output = py"create_plasma_vessel_distance_objective"()
end 


function PrincipalCurvature(;
    eq = nothing,  
    target = 1, 
    bounds = nothing, 
    weight = 1, 
    normalize = true, 
    normalize_target = true, 
    grid = nothing, 
    name = "principal-curvature"
)
    py"""
    import numpy as np
    import desc
    def create_principal_curvature_objective():
        return desc.objectives.PrincipalCurvature(
            eq=$eq, target=$target, bounds = $bounds, weight=$weight, 
            normalize=$normalize, normalize_target=$normalize_target, 
            grid=$grid, name=$name)
    """
    output = py"create_principal_curvature_objective"()
end 


function QuasisymmetryBoozer(;
    eq = nothing,  
    target = 0, 
    bounds = nothing, 
    weight = 1, 
    normalize = true, 
    normalize_target = true, 
    grid = nothing, 
    helicity=(1,0), 
    M_booz=nothing, 
    N_booz = nothing, 
    name = "QS Boozer"
)
    py"""
    import numpy as np
    import desc
    def create_quasisymmetry_boozer_objective():
        return desc.objectives.QuasisymmetryBoozer(
            eq=$eq, target=$target, bounds = $bounds, weight=$weight, 
            normalize=$normalize, normalize_target=$normalize_target, 
            grid=$grid, helicity=$helicity, M_booz=$M_booz, N_booz=$N_booz, 
            name=$name)
    """
    output = py"create_quasisymmetry_boozer_objective"()
end 

function QuasisymmetryTwoTerm(;
    eq = nothing, 
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
    import desc.objectives
    def create_objective_quasisymmetry_two_term():
        return desc.objectives.QuasisymmetryTwoTerm(
            eq=$eq, target=$target, bounds=$bounds, weight=$weight, 
            normalize=$normalize, normalize_target=$normalize_target, 
            grid=$grid, helicity=$helicity, name=$name)
    """
    output = py"create_objective_quasisymmetry_two_term"()
end 


function QuasisymmetryTripleProduct(;
    eq = nothing,  
    target = 0, 
    bounds = nothing, 
    weight = 1, 
    normalize = true, 
    normalize_target = true, 
    grid = nothing, 
    name = "QS triple product"
)
    py"""
    import numpy as np
    import desc
    def create_quasisymmetry_triple_product_objective():
        return desc.objectives.QuasisymmetryTripleProduct(
            eq=$eq, target=$target, bounds = $bounds, weight=$weight, 
            normalize=$normalize, normalize_target=$normalize_target, 
            grid=$grid, name=$name)
    """
    output = py"create_quasisymmetry_triple_product_objective"()
end 


function RadialForceBalance(;
    eq = nothing,  
    target = 0, 
    bounds = nothing, 
    weight = 1, 
    normalize = true, 
    normalize_target = true, 
    grid = nothing, 
    name = "radial force"
)
    py"""
    import numpy as np
    import desc
    # TODO: Grid
    def create_radial_force_balance_objective():
        return desc.objectives.RadialForceBalance(
            eq=$eq, target=$target, bounds = $bounds, weight=$weight, 
            normalize=$normalize, normalize_target=$normalize_target, 
            grid=$grid, name=$name)
    """
    output = py"create_radial_force_balance_objective"()
end 




function RotationalTransform(;
    eq = nothing,  
    target = 0, 
    bounds = nothing, 
    weight = 1, 
    normalize = true, 
    normalize_target = true, 
    grid = nothing, 
    name = "rotational transform"
)
    py"""
    import numpy as np
    import desc
    def create_rotational_transform_objective():
        return desc.objectives.RotationalTransform(
            eq=$eq, target=$target, bounds = $bounds, weight=$weight, 
            normalize=$normalize, normalize_target=$normalize_target, 
            grid=$grid, name=$name)
    """
    output = py"create_rotational_transform_objective"()
end 


function ToroidalCurrent(;
    eq = nothing,  
    target = 0, 
    bounds = nothing, 
    weight = 1, 
    normalize = true, 
    normalize_target = true, 
    grid = nothing, 
    name = "toroidal current"
)
    py"""
    import numpy as np
    import desc
    def create_toroidal_current_objective():
        return desc.objectives.ToroidalCurrent(
            eq=$eq, target=$target, bounds = $bounds, weight=$weight, 
            normalize=$normalize, normalize_target=$normalize_target, 
            grid=$grid, name=$name)
    """
    output = py"create_toroidal_current_objective"()
end 


function Volume(;
    eq = nothing,  
    target = 1, 
    bounds = nothing, 
    weight = 1, 
    normalize = true, 
    normalize_target = true, 
    grid = nothing, 
    name = "volume"
)
    py"""
    import numpy as np
    import desc
    def create_volume_objective():
        return desc.objectives.Volume(
            eq=$eq, target=$target, bounds = $bounds, weight=$weight, 
            normalize=$normalize, normalize_target=$normalize_target, 
            grid=$grid, name=$name)
    """
    output = py"create_volume_objective"()
end 

function objectives_load(
    load_from;
    file_format = nothing
)    
    py"""
    import desc
    def load():
        return desc.objectives.GenericObjective.load(
           $load_from, 
           file_format = $file_format
        )
    """
    output = py"load"()

end 

# function objectives_compute(
#     args...; 
#     kwargs...
# )    
#     kwargs_dict = Dict(pairs(kwargs))   

#     py"""
#     import numpy as np
#     import desc
#     def compute():
#         return objective.compute(
#             *args, 
#             **kwargs_dict
#         )
#     """
#     output = py"compute"()

# end 

# function objectives_compute_scalar(
#     args...; 
#     kwargs...
# )    
#     kwargs_dict = Dict(pairs(kwargs))   
    
#     py"""
#     import numpy as np
#     import desc
#     def compute():
#         return objective.compute_scalar(
#             *args, 
#             **kwargs_dict
#         )
#     """
#     output = py"compute"()

# end 

# function objectives_compute_scaled(
#     args...; 
#     kwargs...
# )    
#     kwargs_dict = Dict(pairs(kwargs))   
    
#     py"""
#     import numpy as np
#     import desc
#     def compute():
#         return objective.compute_scaled(
#             *args, 
#             **kwargs_dict
#         )
#     """
#     output = py"compute"()

# end 
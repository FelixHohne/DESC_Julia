using PyCall

#= 
Functions implement objectives that are in file 
desc.objectives.linear_objectives. 

Seperated from main objectives file as original API 
did not properly Linear Objectives in API Spec. 
=#

function FixBoundaryR(;
    eq=nothing,
    target=nothing,
    bounds=nothing,
    weight=1,
    normalize=true,
    normalize_target=true,
    modes=true,
    surface_label=nothing,
    name="lcfs R"
)

    py"""
    import numpy as np
    import desc
    from desc.objectives.linear_objectives import *
    
    if not isinstance($modes, bool):
        
        python_modes = $modes 
        python_modes = np.ascontiguousarray(python_modes)
    else:
        python_modes = $modes
    def create_objective_fix_boundary_r():
        return FixBoundaryR(
            eq=$eq,
            target=$target,
            bounds=$bounds,
            weight=$weight,
            normalize=$normalize,
            normalize_target=$normalize_target,
            modes=python_modes,
            surface_label=$surface_label,
            name=$name
        )
    """
    output = py"create_objective_fix_boundary_r"()
end 


function FixBoundaryZ(;
    eq = nothing, 
    target = nothing,  
    bounds = nothing, 
    weight = 1, 
    normalize = true, 
    normalize_target = true, 
    modes = true, 
    surface_label = nothing, 
    name = "lcfs Z"
)

    py"""
    import numpy as np
    import desc
    from desc.objectives.linear_objectives import *

    if not isinstance($modes, bool):
        
        python_modes = $modes 
        python_modes = np.ascontiguousarray(python_modes)
    else:
        python_modes = $modes

    def create_objective_fix_boundary_z():
        return FixBoundaryZ(
            eq=$eq, target=$target, bounds = $bounds, weight=$weight, 
            normalize=$normalize, normalize_target=$normalize_target, 
            modes = python_modes, 
            surface_label = $surface_label, 
            name=$name)
    """
    output = py"create_objective_fix_boundary_z"()
end 

function FixLambdaGauge(;
    eq = nothing, 
    target = 0,  
    bounds = nothing, 
    weight = 1, 
    normalize = false, 
    normalize_target = false, 
    name = "lambda gauge"
)

    py"""
    import numpy as np
    import desc
    from desc.objectives.linear_objectives import *
    def create_objective_fix_lambda_gauge():
        return FixLambdaGauge(
            eq=$eq, target=$target, bounds = $bounds, 
            weight=$weight, normalize=$normalize, 
            normalize_target=$normalize_target, 
            name=$name)
    """
    output = py"create_objective_fix_lambda_gauge"()
end 


function FixThetaSFL(;
    eq=nothing, 
    target=0, 
    weight=1, 
    name="Theta SFL"
)
    py"""
    import numpy as np
    import desc
    from desc.objectives.linear_objectives import *
    def create_objective_ix_theta_sfl():
        return FixThetaSFL(
            eq=$eq, target=$target, weight=$weight, name=$name
        )
    """
    output = py"create_objective_ix_theta_sfl"()
end 

function FixAxisR(;
    eq=nothing,
    target=nothing,
    weight=1,
    modes=true,
    normalize=true,
    normalize_target=true,
    name="axis R",
)
    py"""
    import numpy as np
    import desc
    from desc.objectives.linear_objectives import *
    def create_objective_fix_axis_r():
        return FixAxisR(
            eq=$eq,
            target=$target,
            weight=$weight,
            name=$name,
            normalize=$normalize,
            normalize_target=$normalize_target,
        )
    """
    output = py"create_objective_fix_axis_r"()
end 

#=
For unknown reason, Fix Axiz Z, in both Python and Julia DESC 
cannot take in bounds argument. 
=#
function FixAxisZ(;
    eq=nothing,
    target=nothing,
    bounds=nothing, 
    weight=1,
    normalize=true,
    normalize_target=true,
    modes = true, 
    name="axis Z",
)
    py"""
    import numpy as np
    import desc
    from desc.objectives.linear_objectives import *

    if isinstance($target, np.ndarray):
        new_target = np.ascontiguousarray($target)
        assert new_target.flags['C_CONTIGUOUS']
    else:
        new_target = $target

    if isinstance($weight, np.ndarray):
        new_weight = np.ascontiguousarray($weight)
        assert new_weight.flags['C_CONTIGUOUS']
    else:
        new_weight = $weight
    
    if isinstance($modes, np.ndarray):
        new_modes = np.ascontiguousarray($modes)
        assert new_modes.flags['C_CONTIGUOUS']
    else:
        new_modes = $modes

    def objective():
        return FixAxisZ(
            eq = $eq, 
            target = new_target, 
            # bounds = $bounds, 
            weight = new_weight, 
            normalize = $normalize, 
            normalize_target = $normalize_target, 
            modes = new_modes, 
            name = $name
        )
    """
    output = py"objective"()
end 
function FixPressure(;
    eq=nothing,
    target=nothing,
    bounds=nothing,
    weight=1,
    normalize=true,
    normalize_target=true,
    profile=nothing,
    indices=true,
    name="fixed-pressure",
)
    py"""
    import numpy as np
    import desc
    from desc.objectives.linear_objectives import *
    def create_objective_fix_pressure():
        return FixPressure(
            eq=$eq,
            target=$target,
            bounds=$bounds,
            weight=$weight,
            normalize=$normalize,
            normalize_target=$normalize_target,
            profile=$profile,
            indices=$indices,
            name=$name,
        )
    """
    output = py"create_objective_fix_pressure"()
end 


function FixCurrent(;
    eq=nothing,
    target=nothing,
    bounds=nothing,
    weight=1,
    normalize=true,
    normalize_target=true,
    profile=nothing,
    indices=true,
    name="fixed-current",
)
    py"""
    import numpy as np
    import desc
    from desc.objectives.linear_objectives import *
    def create_objective_fix_current():
        return FixCurrent(
            eq=$eq,
            target=$target,
            bounds=$bounds,
            weight=$weight,
            normalize=$normalize,
            normalize_target=$normalize_target,
            profile=$profile,
            indices=$indices,
            name=$name,
        )
    """
    output = py"create_objective_fix_current"()
end 


function FixPsi(;
    eq=nothing,
    target=nothing,
    bounds=nothing,
    weight=1,
    normalize=true,
    normalize_target=true,
    name="fixed-Psi",
)
    py"""
    import numpy as np
    import desc
    from desc.objectives.linear_objectives import *
    def create_objective_fix_psi():
        return FixPsi(
            eq=$eq,
            target=$target,
            bounds=$bounds,
            weight=$weight,
            normalize=$normalize,
            normalize_target=$normalize_target,
            name=$name,
        )
    """
    output = py"create_objective_fix_psi"()
end 

function FixIota(;
    eq=nothing,
    target=nothing,
    bounds=nothing,
    weight=1,
    normalize=false,
    normalize_target=false,
    profile=nothing, 
    indices=true,
    name="fixed-iota",
)
    py"""
    import numpy as np
    import desc
    from desc.objectives.linear_objectives import *
    def create_objective_fix_iota():
        return FixIota(
            eq=$eq,
            target=$target,
            bounds=$bounds,
            weight=$weight,
            normalize=$normalize,
            normalize_target=$normalize_target,
            profile=$profile, 
            indices=$indices,
            name=$name,
        )
    """
    output = py"create_objective_fix_iota"()
end 


function FixAtomicNumber(; 
    eq = nothing, 
    target = nothing, 
    bounds = nothing, 
    weight = 1, 
    normalize = false, 
    normalize_target = false, 
    profile = nothing, 
    indices = true, 
    name = "fixed-atomic-number"
)
    py"""
    import desc 
    import numpy as np  

    if isinstance($target, np.ndarray):
        new_target = np.ascontiguousarray($target)
        assert new_target.flags['C_CONTIGUOUS']
    else:
        new_target = $target

    if isinstance($weight, np.ndarray):
        new_weight = np.ascontiguousarray($weight)
        assert new_weight.flags['C_CONTIGUOUS']
    else:
        new_weight = $weight

    if isinstance($indices, np.ndarray):
        new_indices = np.ascontiguousarray($indices)
        assert new_indices.flags['C_CONTIGUOUS']
    else:
        new_indices = $indices
    
    def objective():
        return desc.objectives.FixAtomicNumber(
        eq = $eq, 
        target = new_target, 
        bounds = $bounds, 
        weight = new_weight, 
        normalize = $normalize, 
        normalize_target = $normalize_target, 
        profile = $profile, 
        indices = new_indices, 
        name = $name
        )
    """
    output = py"objective"()
end

function FixElectronDensity(; 
    eq = nothing, 
    target = nothing, 
    bounds = nothing, 
    weight = 1, 
    normalize = true, 
    normalize_target = true, 
    profile = nothing, 
    indices = true, 
    name = "fixed-electron-density"
)
    py"""
    import desc 
    import numpy as np  

    if isinstance($target, np.ndarray):
        new_target = np.ascontiguousarray($target)
        assert new_target.flags['C_CONTIGUOUS']
    else:
        new_target = $target

    if isinstance($weight, np.ndarray):
        new_weight = np.ascontiguousarray($weight)
        assert new_weight.flags['C_CONTIGUOUS']
    else:
        new_weight = $weight

    if isinstance($indices, np.ndarray):
        new_indices = np.ascontiguousarray($indices)
        assert new_indices.flags['C_CONTIGUOUS']
    else:
        new_indices = $indices
    
    def objective():
        return desc.objectives.FixElectronDensity(
        eq = $eq, 
        target = new_target, 
        bounds = $bounds, 
        weight = new_weight, 
        normalize = $normalize, 
        normalize_target = $normalize_target, 
        profile = $profile, 
        indices = new_indices, 
        name = $name
        )
    """
    output = py"objective"()

end 

function FixElectronTemperature(; 
    eq = nothing, 
    target = nothing, 
    bounds = nothing, 
    weight = 1, 
    normalize = true, 
    normalize_target = true, 
    profile = nothing, 
    indices = true, 
    name = "fixed-electron-temperature"
)
    py"""
    import desc 
    import numpy as np  

    if isinstance($target, np.ndarray):
        new_target = np.ascontiguousarray($target)
        assert new_target.flags['C_CONTIGUOUS']
    else:
        new_target = $target

    if isinstance($weight, np.ndarray):
        new_weight = np.ascontiguousarray($weight)
        assert new_weight.flags['C_CONTIGUOUS']
    else:
        new_weight = $weight

    if isinstance($indices, np.ndarray):
        new_indices = np.ascontiguousarray($indices)
        assert new_indices.flags['C_CONTIGUOUS']
    else:
        new_indices = $indices
    
    def objective():
        return desc.objectives.FixElectronTemperature(
        eq = $eq, 
        target = new_target, 
        bounds = $bounds, 
        weight = new_weight, 
        normalize = $normalize, 
        normalize_target = $normalize_target, 
        profile = $profile, 
        indices = new_indices, 
        name = $name
        )
    """
    output = py"objective"()

end 


function FixIonTemperature(; 
    eq = nothing, 
    target = nothing, 
    bounds = nothing, 
    weight = 1, 
    normalize = true, 
    normalize_target = true, 
    profile = nothing, 
    indices = true, 
    name = "fixed-ion-temperature"
)
    py"""
    import desc 
    import numpy as np  

    if isinstance($target, np.ndarray):
        new_target = np.ascontiguousarray($target)
        assert new_target.flags['C_CONTIGUOUS']
    else:
        new_target = $target

    if isinstance($weight, np.ndarray):
        new_weight = np.ascontiguousarray($weight)
        assert new_weight.flags['C_CONTIGUOUS']
    else:
        new_weight = $weight

    if isinstance($indices, np.ndarray):
        new_indices = np.ascontiguousarray($indices)
        assert new_indices.flags['C_CONTIGUOUS']
    else:
        new_indices = $indices
    
    def objective():
        return desc.objectives.FixIonTemperature(
        eq = $eq, 
        target = new_target, 
        bounds = $bounds, 
        weight = new_weight, 
        normalize = $normalize, 
        normalize_target = $normalize_target, 
        profile = $profile, 
        indices = new_indices, 
        name = $name
        )
    """
    output = py"objective"()

end 


#=
Bounds argument does not appear to work in Python DESC.
=#
function FixModeR(; 
    eq = nothing, 
    target = nothing, 
    bounds = nothing, 
    weight = 1, 
    normalize = true, 
    normalize_target = true, 
    modes = true, 
    name = "'Fix Mode R"
)
    py"""
    import desc 
    import numpy as np  

    if isinstance($target, np.ndarray):
        new_target = np.ascontiguousarray($target)
        assert new_target.flags['C_CONTIGUOUS']
    else:
        new_target = $target

    if isinstance($weight, np.ndarray):
        new_weight = np.ascontiguousarray($weight)
        assert new_weight.flags['C_CONTIGUOUS']
    else:
        new_weight = $weight

    if isinstance($modes, np.ndarray):
        new_modes = np.ascontiguousarray($modes)
        assert new_modes.flags['C_CONTIGUOUS']
    else:
        new_modes = $modes
    
    def objective():
        return desc.objectives.FixModeR(
        eq = $eq, 
        target = new_target, 
        # bounds = $bounds, 
        weight = new_weight, 
        normalize = $normalize, 
        normalize_target = $normalize_target, 
        modes = new_modes,
        name = $name
        )
    """
    output = py"objective"()

end 

function FixModeZ(; 
    eq = nothing, 
    target = nothing, 
    bounds = nothing, 
    weight = 1, 
    normalize = true, 
    normalize_target = true, 
    modes = true, 
    name = "'Fix Mode Z"
)
    py"""
    import desc 
    import numpy as np  

    if isinstance($target, np.ndarray):
        new_target = np.ascontiguousarray($target)
        assert new_target.flags['C_CONTIGUOUS']
    else:
        new_target = $target

    if isinstance($weight, np.ndarray):
        new_weight = np.ascontiguousarray($weight)
        assert new_weight.flags['C_CONTIGUOUS']
    else:
        new_weight = $weight

    if isinstance($modes, np.ndarray):
        new_modes = np.ascontiguousarray($modes)
        assert new_modes.flags['C_CONTIGUOUS']
    else:
        new_modes = $modes
    
    def objective():
        return desc.objectives.FixModeZ(
        eq = $eq, 
        target = new_target, 
        # bounds = $bounds, 
        weight = new_weight, 
        normalize = $normalize, 
        normalize_target = $normalize_target, 
        modes = new_modes,
        name = $name
        )
    """
    output = py"objective"()

end 


function FixSumModesR(; 
    eq = nothing, 
    target = nothing, 
    bounds = nothing, 
    weight = 1, 
    normalize = true, 
    normalize_target = true, 
    sum_weights = nothing, 
    modes = true, 
    name = "Fix Sum Modes R"
)
    py"""
    import desc 
    import numpy as np  

    if isinstance($target, np.ndarray):
        new_target = np.ascontiguousarray($target)
        assert new_target.flags['C_CONTIGUOUS']
    else:
        new_target = $target

    if isinstance($weight, np.ndarray):
        new_weight = np.ascontiguousarray($weight)
        assert new_weight.flags['C_CONTIGUOUS']
    else:
        new_weight = $weight

    if isinstance($sum_weights, np.ndarray):
        new_sum_weights = np.ascontiguousarray($sum_weights)
        assert new_sum_weights.flags['C_CONTIGUOUS']
    else:
        new_sum_weights = $sum_weights

    if isinstance($modes, np.ndarray):
        new_modes = np.ascontiguousarray($modes)
        assert new_modes.flags['C_CONTIGUOUS']
    else:
        new_modes = $modes
    
    def objective():
        return desc.objectives.FixSumModesR(
        eq = $eq, 
        target = new_target, 
        # bounds = $bounds, 
        weight = new_weight, 
        normalize = $normalize, 
        normalize_target = $normalize_target, 
        sum_weights = new_sum_weights, 
        modes = new_modes,
        name = $name
        )
    """
    output = py"objective"()

end 


function FixSumModesZ(; 
    eq = nothing, 
    target = nothing, 
    bounds = nothing, 
    weight = 1, 
    normalize = true, 
    normalize_target = true, 
    sum_weights = nothing, 
    modes = true, 
    name = "'Fix Sum Modes Z"
)
    py"""
    import desc 
    import numpy as np  

    if isinstance($target, np.ndarray):
        new_target = np.ascontiguousarray($target)
        assert new_target.flags['C_CONTIGUOUS']
    else:
        new_target = $target

    if isinstance($weight, np.ndarray):
        new_weight = np.ascontiguousarray($weight)
        assert new_weight.flags['C_CONTIGUOUS']
    else:
        new_weight = $weight

    if isinstance($sum_weights, np.ndarray):
        new_sum_weights = np.ascontiguousarray($sum_weights)
        assert new_sum_weights.flags['C_CONTIGUOUS']
    else:
        new_sum_weights = $sum_weights

    if isinstance($modes, np.ndarray):
        new_modes = np.ascontiguousarray($modes)
        assert new_modes.flags['C_CONTIGUOUS']
    else:
        new_modes = $modes
    
    def objective():
        return desc.objectives.FixSumModesZ(
        eq = $eq, 
        target = new_target, 
        # bounds = $bounds, 
        weight = new_weight, 
        normalize = $normalize, 
        normalize_target = $normalize_target, 
        sum_weights = new_sum_weights, 
        modes = new_modes,
        name = $name
        )
    """
    output = py"objective"()

end 
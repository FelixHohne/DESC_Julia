using PyCall

function jl_objective_fix_boundary_r(;
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

    # if typeof(modes) != Bool
    #     println(typeof(modes))
    #     modes = PyCall.PyReverseDims(modes)
    # end 

    py"""
    import numpy as np
    import desc
    from desc.objectives.linear_objectives import *
    
    if not isinstance($modes, bool):
        
        python_modes = $modes 
        print("Python modes originally Contiguous:", python_modes.flags['C_CONTIGUOUS'])
        python_modes = np.ascontiguousarray(python_modes)
        print("Python modes now Contiguous:", python_modes.flags['C_CONTIGUOUS'])

        # print("Modes code")
        # print(type($modes))
        # print($modes.dtype)
        # print($modes.shape)
        # print("Originally Contiguous:", $modes.flags['C_CONTIGUOUS'])
        # print("Now: ", np.ascontiguousarray($modes).flags['C_CONTIGUOUS'])
        # new_modes = np.ascontiguousarray($modes)
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


function jl_objective_fix_boundary_z(;
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
        print("Python modes originally Contiguous:", python_modes.flags['C_CONTIGUOUS'])
        python_modes = np.ascontiguousarray(python_modes)
        print("Python modes now Contiguous:", python_modes.flags['C_CONTIGUOUS'])

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

function jl_objective_fix_lambda_gauge(;
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


function jl_objective_fix_theta_sfl(;
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

function jl_objective_fix_axis_r(;
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

function jl_objective_fix_pressure(;
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


function jl_objective_fix_current(;
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


function jl_objective_fix_psi(;
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





function jl_objective_fix_iota(;
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






using PyCall

#=
Functions implement desc.equilibrium API
=#

function Equilibrium(;
    Psi=1.0,
    NFP=nothing,
    L=nothing,
    M=nothing,
    N=nothing,
    L_grid=nothing,
    M_grid=nothing,
    N_grid=nothing,
    node_pattern=nothing,
    pressure=nothing,
    iota=nothing,
    current=nothing,
    electron_temperature=nothing,
    electron_density=nothing,
    ion_temperature=nothing,
    atomic_number=nothing,
    surface=nothing,
    axis=nothing,
    sym=nothing,
    spectral_indexing=nothing, 
    kwargs...
)

    kwargs_dict = Dict(pairs(kwargs))

    py"""
    import numpy as np
    import desc
    import desc.equilibrium

    if isinstance($pressure, np.ndarray):
        new_pressure = np.ascontiguousarray($pressure)
        assert new_pressure.flags['C_CONTIGUOUS']
    else:
        new_pressure = $pressure

    if isinstance($iota, np.ndarray):
        new_iota = np.ascontiguousarray($iota)
        assert new_iota.flags['C_CONTIGUOUS']
    else:
        new_iota = $iota

    if isinstance($current, np.ndarray):
        new_current = np.ascontiguousarray($current)
        assert new_current.flags['C_CONTIGUOUS']
    else:
        new_current = $current

    if isinstance($electron_temperature, np.ndarray):
        new_electron_temperature = np.ascontiguousarray($electron_temperature)
        assert new_electron_temperature.flags['C_CONTIGUOUS']
    else:
        new_electron_temperature = $electron_temperature

    if isinstance($electron_density, np.ndarray):
        new_electron_density = np.ascontiguousarray($electron_density)
        assert new_electron_density.flags['C_CONTIGUOUS']
    else:
        new_electron_density = $electron_density
    
    if isinstance($ion_temperature, np.ndarray):
        new_ion_temperature = np.ascontiguousarray($ion_temperature)
        assert new_ion_temperature.flags['C_CONTIGUOUS']
    else:
        new_ion_temperature = $ion_temperature

    if isinstance($atomic_number, np.ndarray):
        new_atomic_number = np.ascontiguousarray($atomic_number)
        assert new_atomic_number.flags['C_CONTIGUOUS']
    else:
        new_atomic_number = $atomic_number
    
    if isinstance($surface, np.ndarray):
        new_surface = np.ascontiguousarray($surface)
        assert new_surface.flags['C_CONTIGUOUS']
    else:
        new_surface = $surface
    
    if isinstance($axis, np.ndarray):
        new_axis = np.ascontiguousarray($axis)
        assert new_axis.flags['C_CONTIGUOUS']
    else:
        new_axis = $axis

    def create_eq_obj():
        eq = desc.equilibrium.Equilibrium(
            Psi=$Psi,
            NFP=$NFP,
            L=$L,
            M=$M,
            N=$N,
            L_grid=$L_grid,
            M_grid=$M_grid,
            N_grid=$N_grid,
            node_pattern=$node_pattern,
            pressure=new_pressure,
            iota=new_iota,
            current=new_current,
            electron_temperature=new_electron_temperature,
            electron_density=new_electron_density,
            ion_temperature=new_ion_temperature,
            atomic_number=new_atomic_number,
            surface=new_surface,
            axis=new_axis,
            sym=$sym,
            spectral_indexing=$spectral_indexing, 
            **$kwargs_dict
        )
        return eq
    """
    output = py"create_eq_obj"()

end 


function EquilibriaFamily(
    args
)
    py"""
    import numpy as np
    import desc
    import desc.equilibrium

    def create_eq_fam_obj():
        eq = desc.equilibrium.EquilibriaFamily(
            $args
        )
        return eq
    """
    output = py"create_eq_fam_obj"()

end


function equilibrium_family_append(eq_fam, eq) 
    py"""
    import numpy as np
    import desc
    def append():
        $eq_fam.append($eq) 
    """
    py"append"()
end 


function qsc_from_paper(name)
    py"""
    import numpy as np
    import desc
    from qsc import Qsc
    def get_qsc():
        return Qsc.from_paper($name)
    """
    output = py"get_qsc"()
end 

function from_near_axis(
    na_eq; 
    r = 0.1, 
    L = nothing, 
    M = 8, 
    N = nothing, 
    ntheta = nothing, 
    spectral_indexing = "ansi"
)
    py"""
    import numpy as np
    import desc
    from qsc import Qsc
    def from_naxis():
        return desc.equilibrium.Equilibrium.from_near_axis(
            $na_eq, 
            r = $r, 
            L = $L, 
            M = $M, 
            N = $N, 
            ntheta = $ntheta, 
            spectral_indexing = $spectral_indexing
        )
    """
    output = py"from_naxis"()
end 

function equilibrium_compute(
    eq, 
    names;
    grid = nothing, 
    params = nothing, 
    transformers = nothing, 
    profiles = nothing, 
    data = nothing, 
    kwargs...
)

    kwargs_dict = Dict(pairs(kwargs))

    py"""
    import numpy as np
    import desc

    new_params = {}
    for item in $params:
        key = item 
        value = $params[item]

        if isinstance(value, np.ndarray):
            new_value = np.ascontiguousarray(value)
            new_params[key] = new_value
    
        else:
            new_params[key] = value 
    
    new_profiles = {}
    for item in $profiles:
        key = item 
        value = $profiles[item]

        if isinstance(value, np.ndarray):
            new_value = np.ascontiguousarray(value)
            new_profiles[key] = new_value
    
        else:
            new_profiles[key] = value 


    def compute():
        return $eq.compute(
            $names, 
            grid = $grid, 
            params = new_params, 
            transforms = $transforms, 
            profiles = new_profiles, 
            data = $data, 
            **$kwargs_dict

        )
    """
    output = py"compute"()
end

function equilibrium_compute_flux_coords(
    eq, 
    real_coords;
    R_lmn = nothing, 
    Z_lmn = nothing, 
    tol = 1e-6, 
    maxiter = 20, 
    rhomin = 1e-6
)
    py"""
    import numpy as np
    import desc

    if isinstance($real_coords, np.ndarray):
        new_real_coords = np.ascontiguousarray($real_coords)
        assert real_coords.flags['C_CONTIGUOUS']
    else:
        new_real_coords = $real_coords
    
    if isinstance($R_lmn, np.ndarray):
        new_R_lmn = np.ascontiguousarray($R_lmn)
        assert R_lmn.flags['C_CONTIGUOUS']
    else:
        new_R_lmn = $R_lmn

    if isinstance($Z_lmn, np.ndarray):
        new_Z_lmn = np.ascontiguousarray($Z_lmn)
        assert Z_lmn.flags['C_CONTIGUOUS']
    else:
        new_Z_lmn = $Z_lmn
    

    def compute_flux_coords():
        return $eq.compute_flux_coords(
            new_real_coords, 
            R_lmn = new_R_lmn, 
            Z_lmn = new_Z_lmn, 
            tol = $tol, 
            maxiter = $maxiter, 
            rhomin = $rhomin
        )
    """
    output = py"compute_flux_coords"()
end 


function equilibrium_compute_theta_coords(
    eq, 
    flux_coords;
    L_lmn = nothing, 
    tol = 1e-6, 
    maxiter = 20
)
    py"""
    import numpy as np
    import desc

    if isinstance($flux_coords, np.ndarray):
        new_flux_coords = np.ascontiguousarray($flux_coords)
        assert flux_coords.flags['C_CONTIGUOUS']
    else:
        new_flux_coords = $flux_coords
    
    if isinstance($L_lmn, np.ndarray):
        new_L_lmn = np.ascontiguousarray($L_lmn)
        assert L_lmn.flags['C_CONTIGUOUS']
    else:
        new_L_lmn = $L_lmn
    

    def compute_theta_coords():
        return $eq.compute_theta_coords(
            new_flux_coords, 
            L_lmn = nenew_L_lmn, 
            tol = $tol, 
            maxiter = $maxiter
        )
    """
    output = py"compute_theta_coords"()
end

function equilibrium_is_nested(
    eq;
    grid = nothing, 
    R_lmn = nothing, 
    Z_lmn = nothing, 
    msg = nothing
)
    py"""
    import numpy as np
    import desc

    if isinstance($R_lmn, np.ndarray):
        new_R_lmn = np.ascontiguousarray($R_lmn)
        assert R_lmn.flags['C_CONTIGUOUS']
    else:
        new_R_lmn = $R_lmn

    if isinstance($Z_lmn, np.ndarray):
        new_Z_lmn = np.ascontiguousarray($Z_lmn)
        assert Z_lmn.flags['C_CONTIGUOUS']
    else:
        new_Z_lmn = $Z_lmn
    
    def is_nested():
        return $eq.is_nested(
            grid = $grid, 
            R_lmn = new_R_lmn,
            Z_lmn = new_Z_lmn, 
            msg = $msg
        )
    """
    output = py"is_nested"()
end

function equilibrium_load(load_from, file_format) 
    py"""
    import numpy as np
    import desc
    import desc.equilibrium
    def load_eq():
        return desc.equilibrium.Equilibrium.load($load_from, $file_format)
    """
    output = py"load_eq"()
end 

function equilibrium_optimize(
    eq;
    objective = nothing,
    constraints = nothing, 
    optimizer = "proximal-lsq-exact", 
    ftol = nothing, 
    xtol = nothing,
    gtol = nothing, 
    maxiter = 50, 
    x_scale = "auto", 
    options = nothing, 
    verbose = 1, 
    copy = false
)
    py"""
    import numpy as np
    import desc

    if isinstance($x_scale, np.ndarray):
        new_x_scale = np.ascontiguousarray($x_scale)
        assert x_scale.flags['C_CONTIGUOUS']
    else:
        new_x_scale = $x_scale

    def optimize():
        return $eq.optimize(
            objective = $objective,
            constraints = $constraints, 
            optimizer = $optimizer, 
            ftol = $ftol, 
            xtol = $xtol,
            gtol = $gtol, 
            maxiter = $maxiter, 
            x_scale = new_x_scale, 
            options = $options, 
            verbose = $verbose, 
            copy = $copy
        )
    """
    output = py"optimize"()
end


function equilibrium_perturb(
    eq, 
    deltas;
    objective = nothing, 
    constraints = nothing, 
    order = 2, 
    tr_ratio = 0.1, 
    weight = "auto", 
    include_f = true, 
    verbose = 1, 
    copy = false
)

    py"""
    import numpy as np
    import desc

    new_deltas = {}
    for item in $deltas:
        key = item 
        value = $deltas[item]

        if isinstance(value, np.ndarray):
            new_value = np.ascontiguousarray(value)
            new_deltas[key] = new_value
    
        else:
            new_deltas[key] = value 
    
    if isinstance($weight, np.ndarray):
        new_weight = np.ascontiguousarray($weight)
        assert new_weight.flags['C_CONTIGUOUS']
    else:
        new_weight = $weight

    def compute():
        return $eq.perturb(
            new_deltas,
            objective = $objective, 
            constraints = $constraints, 
            order = $order, 
            tr_ratio = $tr_ratio, 
            weight = new_weight,
            include_f = $include_f, 
            verbose = $verbose, 
            copy = $copy
        )
    """
    output = py"compute"()
end


function equilibrium_set_initial_guess(
    eq;
    args...
)

    py"""
    import numpy as np
    import desc

    new_args = []
    for item in $args:
        arg = item

        if isinstance(item, np.ndarray):
            new_item = np.ascontiguousarray(item)    
        else:
            new_item = item 
        new_args.append(new_item)

    def compute():
        return $eq.set_initial_guess(
            *new_args
        )
    """
    output = py"compute"()
end
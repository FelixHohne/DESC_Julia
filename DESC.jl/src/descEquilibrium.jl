using PyCall

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



function jl_save_equilibrium(eq, file_name; file_format = "hdf5")
    py"""
    import numpy as np
    import desc

    import desc.equilibrium
    def save_eq():
        $eq.save($file_name, $file_format) 
    """
    py"save_eq"()
end 


function eq_load(load_from, file_format) 
    py"""
    import numpy as np
    import desc
    import desc.equilibrium
    def load_eq():
        return desc.equilibrium.Equilibrium.load($load_from, $file_format)
    """
    output = py"load_eq"()
end 


function jl_equilibrium_family_append(eq_fam, eq) 
    py"""
    import numpy as np
    import desc
    def append():
        $eq_fam.append($eq) 
    """
    py"append"()
end 

function jl_save_equilibrium_family(eq_fam, file_name; file_format = "hdf5")
    py"""
    import numpy as np
    import desc
    def save_eq_fam():
        $eq_fam.save($file_name, $file_format) 
    """
    py"save_eq_fam"()
end 

function jl_solve_equilibrium(
    eq;
    objective = "force", 
    constraints = nothing, 
    optimizer = "lsq-exact", 
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
    def eq_solve():
        return $eq.solve(
            objective = $objective, 
            constraints = $constraints, 
            optimizer = $optimizer, 
            ftol = $ftol, 
            xtol = $xtol,
            gtol = $gtol, 
            maxiter = $maxiter, 
            x_scale = $x_scale, 
            options = $options, 
            verbose = $verbose, 
            copy = $copy 
        ) 
    """
    result = py"eq_solve"()
end


function jl_change_resolution(
    eq;
    L = nothing, 
    M = nothing, 
    N = nothing, 
    L_grid = nothing, 
    M_grid = nothing, 
    N_grid = nothing, 
    NFP = nothing
)
    py"""
    import numpy as np
    import desc
    def func_call():
        return $eq.change_resolution(
            L = $L, 
            M = $M, 
            N = $N, 
            L_grid = $L_grid, 
            M_grid = $M_grid, 
            N_grid = $N_grid, 
            NFP = $NFP
        ) 
    """
    result = py"func_call"()
end

function jl_perturb(
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
    def func_call():
        return $eq.perturb(
            $deltas,
            objective = $objective, 
            constraints = $constraints, 
            order = $order, 
            tr_ratio = $tr_ratio, 
            weight = $weight, 
            include_f = $include_f, 
            verbose = $verbose, 
            copy = $copy
        ) 
    """
    result = py"func_call"()
end


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



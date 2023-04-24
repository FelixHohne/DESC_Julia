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

function jl_equilibrium(;
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
    spectral_indexing=nothing
)

    py"""
    import numpy as np
    import desc

    import desc.equilibrium

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
            pressure=$pressure,
            iota=$iota,
            current=$current,
            electron_temperature=$electron_temperature,
            electron_density=$electron_density,
            ion_temperature=$ion_temperature,
            atomic_number=$atomic_number,
            surface=$surface,
            axis=$axis,
            sym=$sym,
            spectral_indexing=$spectral_indexing
        )
        return eq
    """
    output = py"create_eq_obj"()

end 


function jl_equilibria_family(
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


function jl_load_equilibrium(load_from, file_format) 
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

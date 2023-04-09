using PyCall

function jl_equilibrium(descData)
    py"""
    import numpy as np
    import desc
    from desc import set_device
    set_device('gpu')
    
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
    from desc import set_device
    set_device('gpu')
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
    from desc import set_device
    set_device('gpu')
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
    from desc import set_device
    set_device('gpu')
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
    from desc import set_device
    set_device('gpu')
    import desc.equilibrium
    def load_eq():
        return desc.equilibrium.Equilibrium.load($load_from, $file_format)
    """
    output = py"load_eq"()
end 



using PyCall 


function jl_fourierRZCurve(
    R_n = 10, 
    Z_n = 0, 
    modes_R = nothing, 
    modes_Z = nothing, 
    NFP = 1, 
    sym = "auto", 
    grid = nothing, 
    name = ""
    ) 
    py"""
    import desc 
    import desc.geometry
    import numpy as np  

    def create_fourier_rz_curve():
        eq = desc.geometry.FourierRZCurve(
            R_n = $R_n, 
            Z_n = $Z_n, 
            modes_R = $modes_R, 
            modes_Z = $modes_Z, 
            NFP = $NFP, 
            sym = $sym, 
            grid = $Grid, 
            name = $name
        )

        return eq 
    """
    output = py"create_fourier_rz_curve"()
end 

function jl_fourierRZToroidalSurface(
    R_lmn = nothing, 
    Z_lmn = nothing, 
    modes_R = nothing, 
    modes_Z = nothing, 
    NFP = 1, 
    sym = "auto", 
    rho = 1, 
    grid = nothing, 
    name ="", 
    check_orientation = true) 
    py"""
    import desc 
    import desc.geometry
    import numpy as np  

    def create_fourier_rz_toroidal_surface():
        eq = desc.geometry.FourierRZToroidalSurface(
            R_lmn = $R_lmn, 
            Z_lmn = $Z_lmn, 
            modes_R = $modes_R, 
            modes_Z = $modes_Z, 
            NFP = $NFP, 
            sym = $sym, 
            rho = $rho, 
            grid = $grid, 
            name = $name, 
            check_orientation = $check_orientation
        )
        return eq 
    """
    output = py"create_fourier_rz_toroidal_surface"()
end 


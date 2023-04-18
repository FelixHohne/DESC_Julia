using PyCall 


function jl_fourierRZCurve(;
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

function jl_geometry_fourier_xyz_curve(
    X_n = [0, 10, 2], 
    Y_n = [0, 0, 0], 
    Z_n = [-2, 0, 0], 
    modes = nothing, 
    grid = nothing, 
    name = ""
)

py"""
import desc 
import desc.geometry
import numpy as np  

    def create_fourier_xyz_curve():
        return desc.geometry.FourierXYZCurve(
            X_n =$X_n,
            Y_n =$Y_n,
            Z_n =$Z_n,
            modes = $modes,
            grid = $grid,
            name = $name
        )

    """
    output = py"create_fourier_xyz_curve"()
end 


function jl_geometry_fourier_planar_curve(
    center = [10, 0, 0], 
    normal = [0, 1, 0], 
    r_n = 2, 
    modes = nothing, 
    grid = nothing, 
    name = ""
)

py"""
import desc 
import desc.geometry
import numpy as np  

    def create_fourier_planar_curve():
        return desc.geometry.FourierPlanarCurve(
            center = $center, 
            normal = $normal, 
            r_n = $r_n, 
            modes = $modes, 
            grid = $grid, 
            name = $name
        )
    """
    output = py"create_fourier_planar_curve"()
end 



function jl_fourierRZToroidalSurface(;
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
        return desc.geometry.FourierRZToroidalSurface(
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
    """
    output = py"create_fourier_rz_toroidal_surface"()
end 

function jl_geometry_zernike_RZ_toroidal_section(
    R_lmn = nothing, 
    Z_lmn = nothing, 
    modes_R = nothing, 
    modes_Z = nothing, 
    spectral_indexing = "fringe", 
    sym = "auto", 
    zeta = 0.0, 
    grid = nothing, 
    name = "", 
    check_orientation = true
)

py"""
import desc 
import desc.geometry
import numpy as np  

    def create_zernike_RZ_toroidal_section():
        return desc.geometry.ZernikeRZToroidalSection(
            R_lmn = $R_lmn, 
            Z_lmn = $Z_lmn, 
            modes_R = $modes_R, 
            modes_Z = $modes_Z, 
            spectral_indexing = $spectral_indexing, 
            sym = $sym, 
            zeta = $zeta, 
            grid = $grid, 
            name = $name, 
            check_orientation = $check_orientation
        )
    """
    output = py"create_zernike_RZ_toroidal_section"()
end 



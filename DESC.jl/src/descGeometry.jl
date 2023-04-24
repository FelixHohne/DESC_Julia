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



function FourierRZToroidalSurface(;
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

    if isinstance($R_lmn, np.ndarray):
        n_R_lmn = np.ascontiguousarray($R_lmn)
        assert n_R_lmn.flags['C_CONTIGUOUS']
    else:
        n_R_lmn = $R_lmn

    if isinstance($Z_lmn, np.ndarray):
        n_Z_lmn = np.ascontiguousarray($Z_lmn)
        assert n_Z_lmn.flags['C_CONTIGUOUS']
    else:
        n_Z_lmn = $Z_lmn
    
    if isinstance($modes_R, np.ndarray):
        n_modes_R = np.ascontiguousarray($modes_R)
        assert n_modes_R.flags['C_CONTIGUOUS']
    else:
        n_modes_R = $modes_R

    if isinstance($modes_Z, np.ndarray):
        n_modes_Z = np.ascontiguousarray($modes_Z)
        assert n_modes_Z.flags['C_CONTIGUOUS']
    else:
        n_modes_Z = $modes_Z


    def create_fourier_rz_toroidal_surface():
        return desc.geometry.FourierRZToroidalSurface(
            R_lmn = n_R_lmn, 
            Z_lmn = n_Z_lmn, 
            modes_R = n_modes_R, 
            modes_Z = n_modes_Z, 
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



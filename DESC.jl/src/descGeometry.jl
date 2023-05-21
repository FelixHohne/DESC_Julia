using PyCall 


#=
Functions implement desc.geometry API
=#
function FourierRZCurve(;
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

    if isinstance($R_n, np.ndarray):
            new_R_n = np.ascontiguousarray($R_n)
            assert new_R_n.flags['C_CONTIGUOUS']
    else:
        new_R_n = $R_n
    
    if isinstance($Z_n, np.ndarray):
        new_Z_n = np.ascontiguousarray($Z_n)
        assert new_Z_n.flags['C_CONTIGUOUS']
    else:
        new_Z_n = $Z_n

    if isinstance($modes_R, np.ndarray):
        new_modes_R = np.ascontiguousarray($modes_R)
        assert new_modes_R.flags['C_CONTIGUOUS']
    else:
        new_modes_R = $modes_R

    if isinstance($modes_Z, np.ndarray):
        new_modes_Z = np.ascontiguousarray($modes_Z)
        assert new_modes_Z.flags['C_CONTIGUOUS']
    else:
        new_modes_Z = $modes_Z

    def create_fourier_rz_curve():
        eq = desc.geometry.FourierRZCurve(
            NFP = $NFP, 
            sym = $sym, 
            grid = $grid, 
            name = $name
        )

        return eq 
    """
    output = py"create_fourier_rz_curve"()
end 

function FourierXYZCurve(;
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

    if isinstance($X_n, np.ndarray):
                new_X_n = np.ascontiguousarray($X_n)
                assert new_X_n.flags['C_CONTIGUOUS']
    else:
        new_X_n = $X_n

    if isinstance($Y_n, np.ndarray):
        new_Y_n = np.ascontiguousarray($Y_n)
        assert new_Y_n.flags['C_CONTIGUOUS']
    else:
        new_Y_n = $Y_n

    if isinstance($Z_n, np.ndarray):
        new_Z_n = np.ascontiguousarray($Z_n)
        assert new_Z_n.flags['C_CONTIGUOUS']
    else:
        new_Z_n = $Z_n

    if isinstance($modes, np.ndarray):
        new_modes = np.ascontiguousarray($modes)
        assert new_modes.flags['C_CONTIGUOUS']
    else:
        new_modes = $modes

    def create_fourier_xyz_curve():
        return desc.geometry.FourierXYZCurve(
            X_n = new_X_n,
            Y_n = new_Y_n,
            Z_n = new_Z_n,
            modes = new_modes,
            grid = $grid,
            name = $name
        )

    """
    output = py"create_fourier_xyz_curve"()
end 


function FourierPlanarCurve(
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

    if isinstance($center, np.ndarray):
        new_center = np.ascontiguousarray($center)
        assert new_center.flags['C_CONTIGUOUS']
    else:
        new_center = $center

    if isinstance($normal, np.ndarray):
        new_normal = np.ascontiguousarray($normal)
        assert new_normal.flags['C_CONTIGUOUS']
    else:
        new_normal = $normal

    if isinstance($r_n, np.ndarray):
        new_r_n = np.ascontiguousarray($r_n)
        assert new_r_n.flags['C_CONTIGUOUS']
    else:
        new_r_n = $r_n

    if isinstance($modes, np.ndarray):
        new_modes = np.ascontiguousarray($modes)
        assert new_modes.flags['C_CONTIGUOUS']
    else:
        new_modes = $modes

    def create_fourier_planar_curve():
        return desc.geometry.FourierPlanarCurve(
            center = new_center, 
            normal = new_normal, 
            r_n = new_r_n, 
            modes = new_modes, 
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

function ZernikeRZToroidalSection(;
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

    def create_zernike_RZ_toroidal_section():
        return desc.geometry.ZernikeRZToroidalSection(
            R_lmn = n_R_lmn, 
            Z_lmn = n_Z_lmn, 
            modes_R = n_modes_R, 
            modes_Z = n_modes_Z, 
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


# Implemented for FourierRZCurve
function geometry_compute_coordinates(
    curve;
    R_n = nothing, 
    Z_n = nothing, 
    grid = nothing, 
    dt = 0, 
    basis = "rpz"
)

py"""
import desc 
import desc.geometry
import numpy as np  

    if isinstance($R_n, np.ndarray):
        n_R_n = np.ascontiguousarray($R_n)
        assert n_R_n.flags['C_CONTIGUOUS']
    else:
        n_R_n = $R_n

    if isinstance($Z_n, np.ndarray):
        n_Z_n = np.ascontiguousarray($Z_n)
        assert n_Z_n.flags['C_CONTIGUOUS']
    else:
        n_Z_n = $Z_n


    def compute_coordinates():
        return $curve.compute_coordinates(
            R_n = n_R_n, 
            Z_n = n_Z_n, 
            grid = $grid, 
            dt = $dt, 
            basis = $basis
        )
    """
    output = py"compute_coordinates"()
end 

# Implemented for FourierXYZCurve
function geometry_compute_coordinates(
    curve;
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

    if isinstance($X_n, np.ndarray):
        n_X_n = np.ascontiguousarray($X_n)
        assert n_X_n.flags['C_CONTIGUOUS']
    else:
        n_X_n = $X_n

    if isinstance($Y_n, np.ndarray):
        n_Y_n = np.ascontiguousarray($Y_n)
        assert n_Y_n.flags['C_CONTIGUOUS']
    else:
        n_Y_n = $Y_n
    
    if isinstance($Z_n, np.ndarray):
        n_Z_n = np.ascontiguousarray($Z_n)
        assert n_Z_n.flags['C_CONTIGUOUS']
    else:
        n_Z_n = $Z_n

    if isinstance($modes, np.ndarray):
        n_modes = np.ascontiguousarray($modes)
        assert n_modes.flags['C_CONTIGUOUS']
    else:
        n_modes = $modes


    def compute_coordinates():
        return $curve.compute_coordinates(
            X_n = n_X_n, 
            Y_n = n_Y_n, 
            Z_n = n_Z_n, 
            modes = n_modes, 
            grid = $grid, 
            name = $name
        )
    """
    output = py"compute_coordinates"()
end 

# Implemented for FourierPlanarCurve
function geometry_compute_coordinates(
    curve;
    center = nothing, 
    normal = nothing, 
    r_n = nothing, 
    grid = nothing, 
    dt = 0, 
    basis = "xyz"
)

py"""
import desc 
import desc.geometry
import numpy as np  

    if isinstance($center, np.ndarray):
        n_center = np.ascontiguousarray($center)
        assert n_center.flags['C_CONTIGUOUS']
    else:
        n_center = $center

    if isinstance($normal, np.ndarray):
        n_normal = np.ascontiguousarray($normal)
        assert n_normal.flags['C_CONTIGUOUS']
    else:
        n_normal = $normal
    
    if isinstance($r_n, np.ndarray):
        n_r_n = np.ascontiguousarray($r_n)
        assert n_r_n.flags['C_CONTIGUOUS']
    else:
        n_r_n = $r_n


    def compute():
        return $curve.compute_coordinates(
            center = n_center, 
            normal = n_normal, 
            r_n = n_r_n, 
            grid = $grid, 
            dt = $dt, 
            basis = $basis
        )
    """
    output = py"compute"()
end 

# Implemented for FourierRZToroidalSurface
function geometry_compute_coordinates(
    curve;
    R_lmn = nothing, 
    Z_lmn = nothing, 
    grid = nothing, 
    dt = 0, 
    dz = 0, 
    basis = "rpz"
)

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


    def compute():
        return $curve.compute_coordinates(
            R_lmn = n_R_lmn, 
            Z_lmn = n_Z_lmn, 
            grid = $grid, 
            dt = $dt, 
            dz = $dz, 
            basis = $basis
        )
    """
    output = py"compute"()
end 

# Implemented for ZernikeRZToroidalSection
function geometry_compute_coordinates(
    curve;
    R_lmn = nothing, 
    Z_lmn = nothing, 
    grid = nothing, 
    dr = 0, 
    dt = 0, 
    basis = "rpz"
)

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

    def compute():
        return $curve.compute_coordinates(
            R_lmn = n_R_lmn, 
            Z_lmn = n_Z_lmn, 
            grid = $grid, 
            dr = $dr, 
            dt = $dt, 
            basis = $basis
        )
    """
    output = py"compute"()
end 


# Implemented for FourierRZCurve
function geometry_compute_curvature(
    curve;
    R_n = nothing, 
    Z_n = nothing, 
    grid = nothing
)

py"""
import desc 
import desc.geometry
import numpy as np  

    if isinstance($R_n, np.ndarray):
        n_R_n = np.ascontiguousarray($R_n)
        assert n_R_n.flags['C_CONTIGUOUS']
    else:
        n_R_n = $R_n

    if isinstance($Z_n, np.ndarray):
        n_Z_n = np.ascontiguousarray($Z_n)
        assert n_Z_n.flags['C_CONTIGUOUS']
    else:
        n_Z_n = $Z_n


    def compute_curvature():
        return $curve.compute_curvature(
            R_n = n_R_n, 
            Z_n = n_Z_n, 
            grid = $grid
        )
    """
    output = py"compute_curvature"()
end 

# Implemented for FourierXYZCurve
function geometry_compute_curvature(
    curve;
    X_n = nothing, 
    Y_n = nothing, 
    Z_n = nothing, 
    grid = nothing
)

py"""
import desc 
import desc.geometry
import numpy as np  

    if isinstance($X_n, np.ndarray):
        n_X_n = np.ascontiguousarray($X_n)
        assert n_X_n.flags['C_CONTIGUOUS']
    else:
        n_X_n = $X_n

    if isinstance($Y_n, np.ndarray):
        n_Y_n = np.ascontiguousarray($Y_n)
        assert n_Y_n.flags['C_CONTIGUOUS']
    else:
        n_Y_n = $Y_n
    
    if isinstance($Z_n, np.ndarray):
        n_Z_n = np.ascontiguousarray($Z_n)
        assert n_Z_n.flags['C_CONTIGUOUS']
    else:
        n_Z_n = $Z_n

    def compute():
        return $curve.compute_curvature(
            X_n = n_X_n, 
            Y_n = n_Y_n, 
            Z_n = n_Z_n, 
            grid = $grid, 
        )
    """
    output = py"compute"()
end 


# Implemented for FourierPlanarCurve
function geometry_compute_curvature(
    curve;
    center = nothing, 
    normal = nothing, 
    r_n = nothing, 
    grid = nothing
)

py"""
import desc 
import desc.geometry
import numpy as np  

    if isinstance($center, np.ndarray):
        n_center = np.ascontiguousarray($center)
        assert n_center.flags['C_CONTIGUOUS']
    else:
        n_center = $center

    if isinstance($normal, np.ndarray):
        n_normal = np.ascontiguousarray($normal)
        assert n_normal.flags['C_CONTIGUOUS']
    else:
        n_normal = $normal
    
    if isinstance($r_n, np.ndarray):
        n_r_n = np.ascontiguousarray($r_n)
        assert n_r_n.flags['C_CONTIGUOUS']
    else:
        n_r_n = $r_n


    def compute():
        return $curve.compute_curvature(
            center = n_center, 
            normal = n_normal, 
            r_n = n_r_n, 
            grid = $grid
        )
    """
    output = py"compute"()
end 


# Implemented for FourierRZToroidalSurface
function geometry_compute_curvature(
    curve;
    R_lmn = nothing, 
    Z_lmn = nothing, 
    grid = nothing, 
    dt = 0, 
    dz = 0, 
    basis = "rpz"
)

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


    def compute():
        return $curve.compute_curvature(
            R_lmn = n_R_lmn, 
            Z_lmn = n_Z_lmn, 
            grid = $grid
        )
    """
    output = py"compute"()
end 

# Implemented for ZernikeRZToroidalSection
function geometry_compute_coordinates(
    curve;
    R_lmn = nothing, 
    Z_lmn = nothing, 
    grid = nothing
)

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

    def compute():
        return $curve.compute_coordinates(
            R_lmn = n_R_lmn, 
            Z_lmn = n_Z_lmn, 
            grid = $grid
            )
    """
    output = py"compute"()
end 


# Implemented for FourierRZCurve
function geometry_compute_frenet_frame(
    curve;
    R_n = nothing, 
    Z_n = nothing, 
    grid = nothing, 
    basis = "rpz"
)

py"""
import desc 
import desc.geometry
import numpy as np  

    if isinstance($R_n, np.ndarray):
        n_R_n = np.ascontiguousarray($R_n)
        assert n_R_n.flags['C_CONTIGUOUS']
    else:
        n_R_n = $R_n

    if isinstance($Z_n, np.ndarray):
        n_Z_n = np.ascontiguousarray($Z_n)
        assert n_Z_n.flags['C_CONTIGUOUS']
    else:
        n_Z_n = $Z_n


    def compute_frenet_frame():
        return $curve.compute_frenet_frame(
            R_n = n_R_n, 
            Z_n = n_Z_n, 
            grid = $grid, 
            basis = $basis
        )
    """
    output = py"compute_frenet_frame"()
end 

# Implemented for FourierXYZCurve
function geometry_compute_frenet_frame(
    curve;
    X_n = nothing, 
    Y_n = nothing, 
    Z_n = nothing, 
    grid = nothing,
    basis = "xyz"
)

py"""
import desc 
import desc.geometry
import numpy as np  

    if isinstance($X_n, np.ndarray):
        n_X_n = np.ascontiguousarray($X_n)
        assert n_X_n.flags['C_CONTIGUOUS']
    else:
        n_X_n = $X_n

    if isinstance($Y_n, np.ndarray):
        n_Y_n = np.ascontiguousarray($Y_n)
        assert n_Y_n.flags['C_CONTIGUOUS']
    else:
        n_Y_n = $Y_n
    
    if isinstance($Z_n, np.ndarray):
        n_Z_n = np.ascontiguousarray($Z_n)
        assert n_Z_n.flags['C_CONTIGUOUS']
    else:
        n_Z_n = $Z_n

    def compute():
        return $curve.compute_frenet_frame(
            X_n = n_X_n, 
            Y_n = n_Y_n, 
            Z_n = n_Z_n, 
            grid = $grid, 
            basis = $basis, 
        )
    """
    output = py"compute"()
end 


# Implemented for FourierPlanarCurve
function geometry_compute_frenet_frame(
    curve;
    center = nothing, 
    normal = nothing, 
    r_n = nothing, 
    grid = nothing, 
    basis = "xyz"
)

py"""
import desc 
import desc.geometry
import numpy as np  

    if isinstance($center, np.ndarray):
        n_center = np.ascontiguousarray($center)
        assert n_center.flags['C_CONTIGUOUS']
    else:
        n_center = $center

    if isinstance($normal, np.ndarray):
        n_normal = np.ascontiguousarray($normal)
        assert n_normal.flags['C_CONTIGUOUS']
    else:
        n_normal = $normal
    
    if isinstance($r_n, np.ndarray):
        n_r_n = np.ascontiguousarray($r_n)
        assert n_r_n.flags['C_CONTIGUOUS']
    else:
        n_r_n = $r_n


    def compute():
        return $curve.compute_frenet_frame(
            center = n_center, 
            normal = n_normal, 
            r_n = n_r_n, 
            grid = $grid, 
            basis = $basis
        )
    """
    output = py"compute"()
end 

# Implemented for FourierRZToroidalSurface
function geometry_compute_normal(
    curve;
    R_lmn = nothing, 
    Z_lmn = nothing, 
    grid = nothing, 
    basis = "rpz"
)

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


    def compute():
        return $curve.compute_normal(
            R_lmn = n_R_lmn, 
            Z_lmn = n_Z_lmn, 
            grid = $grid, 
            basis = $basis
        )
    """
    output = py"compute"()
end 

# Implemented for ZernikeRZToroidalSection
function geometry_compute_normal(
    curve;
    R_lmn = nothing, 
    Z_lmn = nothing, 
    grid = nothing, 
    basis = "rpz"
)

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

    def compute():
        return $curve.compute_normal(
            R_lmn = n_R_lmn, 
            Z_lmn = n_Z_lmn, 
            grid = $grid,
            basis = $basis
        )
    """
    output = py"compute"()
end 

# Implemented for FourierRZCurve
function geometry_compute_length(
    curve;
    R_n = nothing, 
    Z_n = nothing, 
    grid = nothing
)

py"""
import desc 
import desc.geometry
import numpy as np  

    if isinstance($R_n, np.ndarray):
        n_R_n = np.ascontiguousarray($R_n)
        assert n_R_n.flags['C_CONTIGUOUS']
    else:
        n_R_n = $R_n

    if isinstance($Z_n, np.ndarray):
        n_Z_n = np.ascontiguousarray($Z_n)
        assert n_Z_n.flags['C_CONTIGUOUS']
    else:
        n_Z_n = $Z_n


    def compute_length():
        return $curve.compute_length(
            R_n = n_R_n, 
            Z_n = n_Z_n, 
            grid = $grid
        )
    """
    output = py"compute_length"()
end 


# Implemented for FourierXYZCurve
function geometry_compute_length(
    curve;
    X_n = nothing, 
    Y_n = nothing, 
    Z_n = nothing, 
    grid = nothing
)

py"""
import desc 
import desc.geometry
import numpy as np  

    if isinstance($X_n, np.ndarray):
        n_X_n = np.ascontiguousarray($X_n)
        assert n_X_n.flags['C_CONTIGUOUS']
    else:
        n_X_n = $X_n

    if isinstance($Y_n, np.ndarray):
        n_Y_n = np.ascontiguousarray($Y_n)
        assert n_Y_n.flags['C_CONTIGUOUS']
    else:
        n_Y_n = $Y_n
    
    if isinstance($Z_n, np.ndarray):
        n_Z_n = np.ascontiguousarray($Z_n)
        assert n_Z_n.flags['C_CONTIGUOUS']
    else:
        n_Z_n = $Z_n

    def compute():
        return $curve.compute_length(
            X_n = n_X_n, 
            Y_n = n_Y_n, 
            Z_n = n_Z_n, 
            grid = $grid, 
            basis = $basis, 
        )
    """
    output = py"compute"()
end 

# Implemented for FourierPlanarCurve
function geometry_compute_length(
    curve;
    center = nothing, 
    normal = nothing, 
    r_n = nothing, 
    grid = nothing
)

py"""
import desc 
import desc.geometry
import numpy as np  

    if isinstance($center, np.ndarray):
        n_center = np.ascontiguousarray($center)
        assert n_center.flags['C_CONTIGUOUS']
    else:
        n_center = $center

    if isinstance($normal, np.ndarray):
        n_normal = np.ascontiguousarray($normal)
        assert n_normal.flags['C_CONTIGUOUS']
    else:
        n_normal = $normal
    
    if isinstance($r_n, np.ndarray):
        n_r_n = np.ascontiguousarray($r_n)
        assert n_r_n.flags['C_CONTIGUOUS']
    else:
        n_r_n = $r_n


    def compute():
        return $curve.compute_length(
            center = n_center, 
            normal = n_normal, 
            r_n = n_r_n, 
            grid = $grid
        )
    """
    output = py"compute"()
end 

# Implemented for FourierRZToroidalSurface
function geometry_compute_surface_area(
    curve;
    R_lmn = nothing, 
    Z_lmn = nothing, 
    grid = nothing
)

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


    def compute():
        return $curve.compute_surface_area(
            R_lmn = n_R_lmn, 
            Z_lmn = n_Z_lmn, 
            grid = $grid 
        )
    """
    output = py"compute"()
end 

# Implemented for ZernikeRZToroidalSection
function geometry_compute_surface_area(
    curve;
    R_lmn = nothing, 
    Z_lmn = nothing, 
    grid = nothing
)

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

    def compute():
        return $curve.surface_area(
            R_lmn = n_R_lmn, 
            Z_lmn = n_Z_lmn, 
            grid = $grid
        )
    """
    output = py"compute"()
end 

# Implemented for FourierRZCurve
function geometry_compute_torsion(
    curve;
    R_n = nothing, 
    Z_n = nothing, 
    grid = nothing
)

py"""
import desc 
import desc.geometry
import numpy as np  

    if isinstance($R_n, np.ndarray):
        n_R_n = np.ascontiguousarray($R_n)
        assert n_R_n.flags['C_CONTIGUOUS']
    else:
        n_R_n = $R_n

    if isinstance($Z_n, np.ndarray):
        n_Z_n = np.ascontiguousarray($Z_n)
        assert n_Z_n.flags['C_CONTIGUOUS']
    else:
        n_Z_n = $Z_n


    def compute_torsion():
        return $curve.compute_torsion(
            R_n = n_R_n, 
            Z_n = n_Z_n, 
            grid = $grid
        )
    """
    output = py"compute_torsion"()
end 

# Implemented for FourierXYZCurve
function geometry_compute_torsion(
    curve;
    X_n = nothing, 
    Y_n = nothing, 
    Z_n = nothing, 
    grid = nothing
)

py"""
import desc 
import desc.geometry
import numpy as np  

    if isinstance($X_n, np.ndarray):
        n_X_n = np.ascontiguousarray($X_n)
        assert n_X_n.flags['C_CONTIGUOUS']
    else:
        n_X_n = $X_n

    if isinstance($Y_n, np.ndarray):
        n_Y_n = np.ascontiguousarray($Y_n)
        assert n_Y_n.flags['C_CONTIGUOUS']
    else:
        n_Y_n = $Y_n
    
    if isinstance($Z_n, np.ndarray):
        n_Z_n = np.ascontiguousarray($Z_n)
        assert n_Z_n.flags['C_CONTIGUOUS']
    else:
        n_Z_n = $Z_n

    def compute():
        return $curve.compute_torsion(
            X_n = n_X_n, 
            Y_n = n_Y_n, 
            Z_n = n_Z_n, 
            grid = $grid, 
            basis = $basis, 
        )
    """
    output = py"compute"()
end 

# Implemented for FourierPlanarCurve
function geometry_compute_torsion(
    curve;
    center = nothing, 
    normal = nothing, 
    r_n = nothing, 
    grid = nothing
)

py"""
import desc 
import desc.geometry
import numpy as np  

    if isinstance($center, np.ndarray):
        n_center = np.ascontiguousarray($center)
        assert n_center.flags['C_CONTIGUOUS']
    else:
        n_center = $center

    if isinstance($normal, np.ndarray):
        n_normal = np.ascontiguousarray($normal)
        assert n_normal.flags['C_CONTIGUOUS']
    else:
        n_normal = $normal
    
    if isinstance($r_n, np.ndarray):
        n_r_n = np.ascontiguousarray($r_n)
        assert n_r_n.flags['C_CONTIGUOUS']
    else:
        n_r_n = $r_n


    def compute():
        return $curve.compute_torsion(
            center = n_center, 
            normal = n_normal, 
            r_n = n_r_n, 
            grid = $grid
        )
    """
    output = py"compute"()
end 


# Implemented for FourierXYZCurve
function geometry_set_coeffs(
    curve;
    X = nothing, 
    Y = nothing, 
    Z = nothing
)

py"""
import desc 
import desc.geometry
import numpy as np  

    if isinstance($X, np.ndarray):
        n_X = np.ascontiguousarray($X)
        assert n_X.flags['C_CONTIGUOUS']
    else:
        n_X = $X

    if isinstance($Y, np.ndarray):
        n_Y = np.ascontiguousarray($Y)
        assert n_Y.flags['C_CONTIGUOUS']
    else:
        n_Y = $Y
    
    if isinstance($Z, np.ndarray):
        n_Z = np.ascontiguousarray($Z)
        assert n_Z.flags['C_CONTIGUOUS']
    else:
        n_Z = $Z

    def compute():
        return $curve.compute_torsion(
            X = n_X, 
            Y = n_Y, 
            Z = n_Z
        )
    """
    output = py"compute"()
end 


function geometry_load(
    load_from; 
    file_format = nothing
)

py"""
import desc 
import desc.geometry
import numpy as np  

    def load():
        return desc.geometry.FourierRZCurve.load(
            $load_from, 
            file_format = $file_format
        )
    """
    output = py"load"()
end 
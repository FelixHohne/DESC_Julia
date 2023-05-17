using PyCall 

#=
Functions implement desc.profiles API
=#

function PowerSeriesProfile(;
    params = [0], 
    modes = nothing, 
    grid = nothing, 
    sym = "auto", 
    name = ""
) 
    py"""
    import desc 
    import desc.profiles
    import numpy as np  

    if isinstance($params, np.ndarray):
        new_params = np.ascontiguousarray($params)
        assert new_params.flags['C_CONTIGUOUS']
    else:
        new_params = $params

    if isinstance($modes, np.ndarray):
        new_modes = np.ascontiguousarray($modes)
        assert new_modes.flags['C_CONTIGUOUS']
    else:
        new_modes = $modes

    def create_power_series_profile():
        return desc.profiles.PowerSeriesProfile(
           params=$params,
           modes=$modes, 
           grid=$grid, 
           sym = $sym, 
           name = $name
        )
    """
    output = py"create_power_series_profile"()
end 


function SplineProfile(;    
    values = [0, 0, 0], 
    knots = nothing, 
    grid = nothing, 
    method = "cubic2", 
    name = ""
) 


    py"""
    import desc 
    import desc.profiles
    import numpy as np  

    if isinstance($values, np.ndarray):
        new_values = np.ascontiguousarray($values)
        assert new_values.flags['C_CONTIGUOUS']
    else:
        new_values = $values

    if isinstance($knots, np.ndarray):
        new_knots = np.ascontiguousarray($knots)
        assert new_knots.flags['C_CONTIGUOUS']
    else:
        new_knots = $knots

    def create_spline_profile():
        return desc.profiles.SplineProfile(
           values = new_values, 
           knots = new_knots, 
           grid= $grid, 
           method = $method, 
           name = $name
        )
    """
    output = py"create_spline_profile"()
end 


function MTanhProfile(;    
    params = [0, 0, 1, 1, 0], 
    grid = nothing, 
    name = ""
) 
    py"""
    import desc 
    import desc.profiles
    import numpy as np  

    if isinstance($params, np.ndarray):
        new_params = np.ascontiguousarray($params)
        assert new_params.flags['C_CONTIGUOUS']
    else:
        new_params = $params
    

    def create_MTanhProfile_profile():
        return desc.profiles.MTanhProfile(
            params = new_params, 
            grid = $grid, 
            name = $name
        )
    """
    output = py"create_MTanhProfile_profile"()
end 

function FourierZernikeProfile(;    
    params = [0],
    modes = nothing, 
    grid = nothing, 
    sym = "auto", 
    NFP = 1,  
    name = ""
) 
    py"""
    import desc 
    import desc.profiles
    import numpy as np  

    if isinstance($params, np.ndarray):
        new_params = np.ascontiguousarray($params)
        assert new_params.flags['C_CONTIGUOUS']
    else:
        new_params = $params

    if isinstance($modes, np.ndarray):
        new_modes = np.ascontiguousarray($modes)
        assert new_modes.flags['C_CONTIGUOUS']
    else:
        new_modes = $modes

    def create_FourierZernike_profile():
        return desc.profiles.FourierZernikeProfile(
            params = new_params, 
            modes = new_modes, 
            grid = $grid, 
            sym = $sym, 
            NFP = $NFP, 
            name = $name
        )
    """
    output = py"create_FourierZernike_profile"()
end 

function ScaledProfile(;    
    profile, 
    scale
) 
    py"""
    import desc 
    import desc.profiles
    import numpy as np  

    def create_scaled_profile():
        return desc.profiles.ScaledProfile(
            profile = $profile, 
            scale = $scale
        )
    """
    output = py"create_scaled_profile"()
end 

function SumProfile(;    
    profiles
) 
    py"""
    import desc 
    import desc.profiles
    import numpy as np  

    def create_sum_profile():
        return desc.profiles.SumProfile(
            profiles = $profiles
        )
    """
    output = py"create_sum_profile"()
end 

function ProductProfile(;    
    profiles
) 
    py"""
    import desc 
    import desc.profiles
    import numpy as np  

    def create_product_profile():
        return desc.profiles.ProductProfile(
            profiles = $profiles
        )
    """
    output = py"create_product_profile"()
end 


function profiles_compute(
    profile;
    params = nothing, 
    grid = nothing, 
    dr = 0, 
    dt = 0, 
    dz = 0
) 
    py"""
    import desc 
    import desc.profiles
    import numpy as np  

    if isinstance($params, np.ndarray):
        new_params = np.ascontiguousarray($params)
        assert new_params.flags['C_CONTIGUOUS']
    else:
        new_params = $params

    def compute():
        return $profile.compute(
           params=new_params,
           grid = $grid, 
           dr = $dr, 
           dt = $dt, 
           dz = $dz
        )
    """
    output = py"compute"()
end 

function profiles_from_values(
    x, 
    y; 
    order = 6, 
    rcond = nothing, 
    w = nothing, 
    grid = nothing, 
    sym = "auto", 
    name = ""
) 
    py"""
    import desc 
    import desc.profiles
    import numpy as np  

    if isinstance($x, np.ndarray):
        new_x = np.ascontiguousarray($x)
        assert new_x.flags['C_CONTIGUOUS']
    else:
        new_x = $x

    if isinstance($y, np.ndarray):
        new_y = np.ascontiguousarray($y)
        assert new_y.flags['C_CONTIGUOUS']
    else:
        new_y = $y

    if isinstance($w, np.ndarray):
        new_w = np.ascontiguousarray($w)
        assert new_w.flags['C_CONTIGUOUS']
    else:
        new_w = $w

    def compute():
        return desc.profiles.PowerSeriesProfile.from_values(
           new_x, 
           new_y,
           order = $order, 
           rcond = $rcond, 
           w = new_w, 
           grid = $grid, 
           sym = $sym, 
           name = $name
        )
    """
    output = py"compute"()
end 


function profiles_load(
    load_from; 
    file_format = nothing
)

py"""
import desc 
import desc.geometry
import numpy as np  

    def load():
        return desc.profiles.PowerSeriesProfile.load(
            $load_from, 
            file_format = $file_format
        )
    """
    output = py"load"()
end 

function profiles_to_fourierzernike(
    profile; 
    L = 6, 
    M = 0, 
    N = 0, 
    NFP = 1, 
    xs = 100, 
    w = nothing
) 
    py"""
    import desc 
    import desc.profiles
    import numpy as np  

    if isinstance($xs, np.ndarray):
        new_xs = np.ascontiguousarray($xs)
        assert new_xs.flags['C_CONTIGUOUS']
    else:
        new_x = $xs

    if isinstance($w, np.ndarray):
        new_w = np.ascontiguousarray($w)
        assert new_w.flags['C_CONTIGUOUS']
    else:
        new_w = $w

    def compute():
        return $profile.to_fourierzernike(
            L = $L, 
            M = $M, 
            N = $N, 
            NFP = $NFP, 
            xs = new_xs, 
            w = new_w
        )
    """
    output = py"compute"()
end 


function profiles_to_mtanh(
    profile; 
    order = 4, 
    xs = 100, 
    w = nothing, 
    p0 = nothing, 
    pmax = nothing, 
    pmin = nothing, 
    kwargs...
) 

    kwargs_dict = Dict(pairs(kwargs))

    py"""
    import desc 
    import desc.profiles
    import numpy as np  

    if isinstance($xs, np.ndarray):
        new_xs = np.ascontiguousarray($xs)
        assert new_xs.flags['C_CONTIGUOUS']
    else:
        new_xs = $xs

    if isinstance($w, np.ndarray):
        new_w = np.ascontiguousarray($w)
        assert new_w.flags['C_CONTIGUOUS']
    else:
        new_w = $w

    if isinstance($p0, np.ndarray):
        new_p0 = np.ascontiguousarray($p0)
        assert new_p0.flags['C_CONTIGUOUS']
    else:
        new_p0 = $p0

    if isinstance($pmin, np.ndarray):
        new_pmin = np.ascontiguousarray($pmin)
        assert new_pmin.flags['C_CONTIGUOUS']
    else:
        new_pmin = $pmin
    
        if isinstance($pmax, np.ndarray):
        new_pmax = np.ascontiguousarray($pmax)
        assert new_pmax.flags['C_CONTIGUOUS']
    else:
        new_pmax = $pmax

    def compute():
        return $profile.to_mtanh(
            order = $order,
            xs = new_xs, 
            w = new_w, 
            p0 = new_p0, 
            pmin = new_pmin, 
            pmax = new_pmax, 
            **kwargs_dict
        )
    """
    output = py"compute"()
end 



function profiles_to_powerseries(
    profile; 
    order = 6, 
    xs = 100, 
    sym = "auto", 
    rcond = nothing, 
    w = nothing
) 
    py"""
    import desc 
    import desc.profiles
    import numpy as np  

    if isinstance($xs, np.ndarray):
        new_xs = np.ascontiguousarray($xs)
        assert new_xs.flags['C_CONTIGUOUS']
    else:
        new_xs = $xs

    if isinstance($w, np.ndarray):
        new_w = np.ascontiguousarray($w)
        assert new_w.flags['C_CONTIGUOUS']
    else:
        new_w = $w
    
    def compute():
        return $profile.to_powerseries(
            order = $order,
            xs = new_xs, 
            sym = $sym, 
            rcond = $rcond, 
            w = new_w 
        )
    """
    output = py"compute"()
end 



function profiles_to_spline(
    profile; 
    knots = 20, 
    method = "cubic2"
) 
    py"""
    import desc 
    import desc.profiles
    import numpy as np  

    if isinstance($knots, np.ndarray):
        new_knots = np.ascontiguousarray($knots)
        assert new_knots.flags['C_CONTIGUOUS']
    else:
        new_knots = $knots
    
    def compute():
        return $profile.to_spline(
            knots = $knots, 
            method = $method
        )
    """
    output = py"to_spline"()
end 

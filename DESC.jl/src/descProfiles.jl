using PyCall 

function jl_profiles_power_series_profile(;
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


function jl_profiles_spline_profile(;    
    # TODO: Where did params go
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

    def create_spline_profile():
        return desc.profiles.SplineProfile(
           values = $values, 
           knots = $knots, 
           grid=$grid, 
           method = $cubic2, 
           name = $name
        )
    """
    output = py"create_spline_profile"()
end 


function jl_profiles_MTanhProfile(;    
    params = [0, 0, 1, 1, 0], 
    grid = nothing, 
    name = ""
) 
    py"""
    import desc 
    import desc.profiles
    import numpy as np  

    def create_MTanhProfile_profile():
        return desc.profiles.MTanhProfile(
            params = $params, 
            grid = $grid, 
            name = $name
        )
    """
    output = py"create_MTanhProfile_profile"()
end 

function jl_profiles_fourier_zernike_profile(;    
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

    def create_FourierZernike_profile():
        return desc.profiles.FourierZernikeProfile(
            params = $params, 
            modes = $modes, 
            grid = $grid, 
            sym = $sym, 
            NFP = $NFP, 
            name = $name
        )
    """
    output = py"create_FourierZernike_profile"()
end 

function jl_profiles_scaled_profile(;    
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

function jl_profiles_sum_profile(;    
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

function jl_profiles_product_profile(;    
    profiles
) 
    py"""
    import desc 
    import desc.profiles
    import numpy as np  

    def create_product_profile():
        return desc.profiles.ProductFile(
            profiles = $profiles
        )
    """
    output = py"create_product_profile"()
end 
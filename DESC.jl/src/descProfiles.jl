using PyCall 

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

    def create_MTanhProfile_profile():
        return desc.profiles.MTanhProfile(
            params = $params, 
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
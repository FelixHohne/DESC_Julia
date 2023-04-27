using PyCall 

function Transform(
    grid, 
    basis;
    derivs = 0, 
    rcond = "auto", 
    build = true, 
    build_pinv = false, 
    method = "auto"
)
    py"""
    import desc 
    import numpy as np  

    if isinstance($derivs, np.ndarray):
        new_derivs = np.ascontiguousarray($derivs)
        assert new_derivs.flags['C_CONTIGUOUS']
    else:
        new_derivs = $derivs

    def transform():
        return desc.transform.Transform(
        grid = $grid, 
        basis = $basis, 
        derivs = new_derivs, 
        rcond = $rcond, 
        build = $build, 
        build_pinv = $build_pinv, 
        method = $method
        )
    """
    output = py"transform"()

end 

function transform_change_derivatives(
    transform, 
    derivs; 
    build = true
)
    py"""
    import desc 
    import desc.profiles
    import numpy as np  

    if isinstance($derivs, np.ndarray):
        new_derivs = np.ascontiguousarray($derivs)
        assert new_derivs.flags['C_CONTIGUOUS']
    else:
        new_derivs = $derivs

    def transform():
        return $transform.change_derivatives(
            new_derivs, 
            build = $build
        )
    """
    output = py"transform"()

end 

function transform_fit(
    transform,
    x
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

    def transform():
        return $transform.fit(
            new_x
        )
    """
    output = py"transform"()

end 

function transform_load(
    load_from; 
    file_format
)
    py"""
    import desc 
    import desc.profiles
    import numpy as np  

    def load():
        return Transform.load(
            $load_from, 
            file_format = $file_format
        )
    """
    output = py"load"()

end 


function transform_project(
    transform, 
    y
)
    py"""
    import desc 
    import numpy as np  

    if isinstance($y, np.ndarray):
        new_y = np.ascontiguousarray($y)
        assert new_y.flags['C_CONTIGUOUS']
    else:
        new_y = $y
    
    def load():
        return $transform.project(
            new_y
        )
    """
    output = py"load"()

end 



function transform_transform(
    transform, 
    c; 
    dr = 0, 
    dt = 0, 
    dz = 0
)
    py"""
    import desc 
    import numpy as np  

    if isinstance($c, np.ndarray):
        new_c = np.ascontiguousarray($c)
        assert new_c.flags['C_CONTIGUOUS']
    else:
        new_c = $c
    
    def load():
        return $transform.transform(
            new_c, 
            dr = $dr, 
            dt = $dt, 
            dz = $dz
        )
    """
    output = py"load"()

end 

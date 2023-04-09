function jl_solve_continuation(
    eqfam;
    objective = "force", 
    optimizer = "lsq-exact", 
    pert_order = 2, 
    ftol = nothing, 
    xtol = nothing, 
    nfev = 100, 
    verbose = 1, 
    checkpoint_path = nothing
)
    # TODO: Natively support kwargs
    py"""
    import numpy as np
    import desc
    import desc.continuation 

    from desc import set_device
    set_device('gpu')
    
    def create_solve_continuation():
        return desc.continuation.solve_continuation(
            eq = $eq, 
            objective = $objective, 
            optimizer = $optimizer, 
            pert_order = $pert_order, 
            ftol = $ftol, 
            xtol = $xtol, 
            nfev = $nfev, 
            verbose = $verbose, 
            checkpoint_path = $checkpoint_path
        )
    """
    output = py"create_solve_continuation"()
    return output 

end 



function jl_solve_continuation_automatic(
    eq;
    objective = "force", 
    optimizer = "lsq-exact", 
    pert_order = 2, 
    ftol = nothing, 
    xtol = nothing, 
    nfev = 100, 
    verbose = 1, 
    checkpoint_path = nothing,
    mres_step = 6, # implemented kwargs 
    pres_step = 0.5, 
    bdry_step = 0.25
)
    # TODO: Natively support kwargs
    py"""
    import numpy as np
    import desc
    import desc.continuation 

    from desc import set_device
    set_device('gpu')
    

    def create_solve_continuation_automatic():
        return desc.continuation.solve_continuation_automatic(
            eq = $eq, 
            objective = $objective, 
            optimizer = $optimizer, 
            pert_order = $pert_order, 
            ftol = $ftol, 
            xtol = $xtol, 
            nfev = $nfev, 
            verbose = $verbose, 
            checkpoint_path = $checkpoint_path, 
            mres_step = $mres_step, 
            pres_step = $pres_step, 
            bdry_step = $bdry_step
        )
    """
    output = py"create_solve_continuation_automatic"()
    return output 

end 
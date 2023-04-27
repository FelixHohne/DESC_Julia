function solve_continuation(
    eqfam;
    objective = "force", 
    optimizer = "lsq-exact", 
    pert_order = 2, 
    ftol = nothing, 
    xtol = nothing,
    gtol = nothing,  
    nfev = 100, 
    verbose = 1, 
    checkpoint_path = nothing
)
    py"""
    import numpy as np
    import desc
    import desc.continuation 
    
    def create_solve_continuation():
        return desc.continuation.solve_continuation(
            $eqfam, 
            objective = $objective, 
            optimizer = $optimizer, 
            pert_order = $pert_order, 
            ftol = $ftol, 
            xtol = $xtol, 
            gtol = $gtol, 
            nfev = $nfev, 
            verbose = $verbose, 
            checkpoint_path = $checkpoint_path
        )
    """
    output = py"create_solve_continuation"()
    return output 

end 

function solve_continuation_automatic(
    eq;
    objective = "force", 
    optimizer = "lsq-exact", 
    pert_order = 2, 
    ftol = nothing, 
    xtol = nothing, 
    gtol = nothing, 
    nfev = 100, 
    verbose = 1, 
    checkpoint_path = nothing,
    mres_step = 6, 
    pres_step = 0.5, 
    bdry_step = 0.25
)
    py"""
    import numpy as np
    import desc
    import desc.continuation 

    def create_solve_continuation_automatic():
        return desc.continuation.solve_continuation_automatic(
            eq = $eq, objective = $objective, optimizer = $optimizer, 
            pert_order = $pert_order, ftol = $ftol, xtol = $xtol, gtol=$gtol, 
            nfev = $nfev, verbose = $verbose, checkpoint_path = $checkpoint_path, 
            mres_step = $mres_step, pres_step = $pres_step, bdry_step = $bdry_step
        )
    """
    output = py"create_solve_continuation_automatic"()
    return output 

end 
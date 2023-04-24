function AutoDiffDerivative(
    fun;
    argnum = 0, 
    mode = "fwd", 
    kwargs...
)

    kwargs_dict = Dict(pairs(kwargs))

    py"""
    import numpy as np
    import desc
    import desc.continuation 
    
    def create_solve_continuation():
        return desc.derivatives.AutoDiffDerivative(
            $fun, 
            argnum = $argnum, 
            mode = $mode, 
            **$kwargs_dict
        )
    """
    output = py"create_solve_continuation"()
    return output 

end 


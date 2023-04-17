using PyCall
function desc_jl_use_gpu_if_available()
    py"""
    def use_gpu():
        import jax 
        if jax.default_backend() == "gpu":
            from desc import set_device
            set_device('gpu')
        from desc.equilibrium import Equilibrium
    """
    py"use_gpu"()
end


function jl_get_fixed_boundary_constraints(; 
    profiles=true, iota = true, kinetic=false, normalize=true)

    py"""
    import desc 
    from desc.objectives.utils import *

    def call_get_fixed_boundary_constraints():
        return get_fixed_boundary_constraints()
    """
    result = py"call_get_fixed_boundary_constraints"()
    
end 

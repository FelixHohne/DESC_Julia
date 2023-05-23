using PyCall

#=
Various useful utility functions. 
=#
function use_gpu_if_available()
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

function get_fixed_boundary_constraints(; 
    profiles=true, iota = true, kinetic=false, normalize=true)

    py"""
    import desc 
    from desc.objectives.utils import *

    def call_get_fixed_boundary_constraints():
        return get_fixed_boundary_constraints(
            profiles = $profiles, 
            iota = $iota, 
            kinetic = $kinetic, 
            normalize = $normalize
        )
    """
    result = py"call_get_fixed_boundary_constraints"()
    
end 

function get_equilibrium_objective(; 
    mode = "force", normalize = true)

    py"""
    import desc 

    def get_equilibrium_objective():
        return desc.objectives.get_equilibrium_objective(
            mode = $mode, 
            normalize = $normalize
        )
    """
    result = py"get_equilibrium_objective"()
end 

function get_NAE_constraints(
    desc_eq, 
    qsc_eq; 
    profiles = true, 
    iota = false, 
    order = 1
)
    py"""
    import desc 

    def get_NAE_constraints():
        return desc.objectives.get_NAE_constraints(
            desc_eq = $desc_eq, 
            qsc_eq = $qsc_eq, 
            profiles = $profiles, 
            iota = $iota, 
            order = $order
        )
    """
    result = py"get_NAE_constraints"()
end 


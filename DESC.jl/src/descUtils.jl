using PyCall
function desc_jl_use_gpu_if_available()
    py"""
    def use_gpu():
        import jax 
        print(jax.default_backend())
        if jax.default_backend() == "gpu":
            from desc import set_device
            set_device('gpu')
        from desc.equilibrium import Equilibrium
    """
    py"use_gpu"()
end


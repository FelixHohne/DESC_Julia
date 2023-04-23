using PyCall 

function jl_plotting_plot_boundaries(
    eqs; 
    labels = nothing, 
    zeta = nothing, 
    ax = nothing, 
    return_data = false, 
    kwargs...
)
    kwargs_dict = Dict(pairs(kwargs))

    py"""
    import desc 
    import desc.plotting
    import numpy as np  

    def plot():
        return desc.plotting.plot_boundaries(
            eqs; 
            labels = $labels, 
            zeta = $zeta, 
            ax = $ax, 
            return_data = $return_data, 
            **$kwargs_dict
            )
    """
    output = py"plot"()
    plt = pyimport("matplotlib.pyplot")
    plt.show()
end 


function jl_plotting_plot_boundary(
    eq; 
    zeta = nothing, 
    plot_axis = false, 
    ax = nothing, 
    return_data = false, 
    kwargs...
)
    kwargs_dict = Dict(pairs(kwargs))

    py"""
    import desc 
    import desc.plotting
    import numpy as np  

    def plot():
        return desc.plotting.plot_boundary(
            eq; 
            zeta = $zeta, 
            plot_axis = $plot_axis, 
            ax = $ax, 
            return_data = $return_data, 
            **$kwargs_dict
            )
    """
    output = py"plot"()
    plt = pyimport("matplotlib.pyplot")
    plt.show()
end 

function jl_plotting_plot_coefficients(
    eq; 
    L = true, 
    M = true, 
    N = true, 
    ax = nothing, 
    kwargs...
)
    kwargs_dict = Dict(pairs(kwargs))

    py"""
    import desc 
    import desc.plotting
    import numpy as np  

    def plot():
        return desc.plotting.plot_coefficients(
            eq = $eq; 
            L = $L, 
            M = $M, 
            N = $N, 
            ax = $ax, 
            **$kwargs_dict
            )
    """
    output = py"plot"()
    plt = pyimport("matplotlib.pyplot")
    plt.show()
end 




function jl_plotting_plot_coils(
    coils, 
    grid = nothing, 
    ax = nothing, 
    return_data = false, 
    kwargs...
)
    kwargs_dict = Dict(pairs(kwargs))

    py"""
    import desc 
    import desc.plotting
    import numpy as np  

    def plot():
        return desc.plotting.plot_coils(
            coils = $coils, 
            grid = $grid, 
            ax = $ax, 
            return_data = $return_data, 
            **$kwargs_dict
            )
    """
    output = py"plot"()
    plt = pyimport("matplotlib.pyplot")
    plt.show()
end 



function jl_plotting_plot_comparison(
    eqs; 
    rho = 8,
    theta = 8, 
    zeta = nothing, 
    cmap = "rainbow", 
    colors = nothing, 
    lws = nothing, 
    linestyles = nothing, 
    labels = nothing, 
    return_data = false, 
    kwargs...
)
    kwargs_dict = Dict(pairs(kwargs))

    py"""
    import desc 
    import desc.plotting
    import numpy as np  

    def plot():
        return desc.plotting.plot_comparison(
            eqs; 
            rho = $rho;
            theta = $theta, 
            zeta = $zeta, 
            cmap = $cmap, 
            colors = $colors, 
            lws = $lws, 
            linestyles = $linestyles, 
            labels = $labels, 
            return_data = $return_data, 
            **$kwargs_dict
            )
    """
    output = py"plot"()
    plt = pyimport("matplotlib.pyplot")
    plt.show()
end 

function jl_plotting_plot_field_lines_sfl(
    eq, 
    rho;
    seed_thetas = 0, 
    phi_start = 0, 
    phi_end = 6.283185307179586, 
    dphi = 0.01, 
    ax = nothing, 
    return_data = nothing, 
    kwargs...
)
    kwargs_dict = Dict(pairs(kwargs))

    py"""
    import desc 
    import desc.plotting
    import numpy as np  

    def plot():
        return desc.plotting.plot_field_lines_sfl(
            eq = $eq, 
            rho = $rho;
            seed_thetas = $seed_thetas, 
            phi_start = $phi_start, 
            phi_end = $phi_end, 
            dphi = $dphi, 
            ax = $ax, 
            return_data = $return_data, 
            **$kwargs_dict
            )
    """
    output = py"plot"()
    plt = pyimport("matplotlib.pyplot")
    plt.show()
end 


function jl_plotting_plot_fsa(
    eq, 
    name, 
    wirth_sqrt_g = true, 
    log = false, 
    rho = 20, 
    M = nothing, 
    N = nothing, 
    norm_f = false, 
    ax = nothing, 
    return_data = false, 
    kwargs...
)
    kwargs_dict = Dict(pairs(kwargs))

    py"""
    import desc 
    import desc.plotting
    import numpy as np  

    def plot():
        return desc.plotting.plot_fsa(
            eq = $eq, 
            name = $name, 
            wirth_sqrt_g = $wirth_sqrt_g, 
            log = $log, 
            rho = $rho, 
            M = $M, 
            N = $N, 
            norm_f = $norm_f, 
            ax = $ax, 
            return_data = $return_data, 
            **$kwargs_dict
            )
    """
    output = py"plot"()
    plt = pyimport("matplotlib.pyplot")
    plt.show()
end 



function jl_plotting_plot_grid(
    grid;
    return_data = false, 
    kwargs...
)
    kwargs_dict = Dict(pairs(kwargs))

    py"""
    import desc 
    import desc.plotting
    import numpy as np  

    def plot():
        return desc.plotting.plot_grid(
            grid$grid, 
            return_data = $return_data, 
            **$kwargs_dict
        )
    """
    output = py"plot"()
    plt = pyimport("matplotlib.pyplot")
    plt.show()
end 



function jl_plotting_plot_qs_error(
    eq; 
    log=true, 
    fB = true, 
    fC = true, 
    fT = true, 
    helicity = (1, 0), 
    rho = nothing, 
    ax = nothing, 
    return_data = false, 
    kwargs...
)

kwargs_dict = Dict(pairs(kwargs))

py"""
import desc 
import desc.plotting
import numpy as np  

def plot():
    return desc.plotting.plot_qs_error(
        eq; 
        log=$log, 
        fB = $fB, 
        fC = $fC, 
        fT = $fT, 
        helicity = $helicity, 
        rho = $rho, 
        ax = $ax, 
        return_data = $return_data, 
        **$kwargs_dict
    )
"""

output = py"plot"()
plt = pyimport("matplotlib.pyplot")
plt.show()

end 



function jl_plotting_plot_logo(;
    savepath = nothing, 
    kwargs...
)
    kwargs_dict = Dict(pairs(kwargs))
    
    py"""
    import desc 
    import desc.plotting
    import numpy as np  

    def plot():
        return desc.plotting.plot_logo(
            savepath, 
            **$kwargs_dict
        )
    """
    output = py"plot"()
    plt = pyimport("matplotlib.pyplot")
    plt.show()
end 

function jl_plotting_plot_section(
    eq, 
    name, 
    grid = nothing, 
    log = false, 
    norm_f = false, 
    ax = nothing, 
    return_data = false, 
    kwargs...
)
    kwargs_dict = Dict(pairs(kwargs))

    py"""
    import desc 
    import desc.plotting
    import numpy as np  

    def plot():
        return desc.plotting.plot_section(
            eq=$eq, 
            name=$name, 
            grid = $grid, 
            log = $log, 
            norm_f = $norm_f, 
            ax = $ax, 
            return_data = $return_data, 
            **$kwargs_dict
        )
    """
    output = py"plot"()
    plt = pyimport("matplotlib.pyplot")
    plt.show()

end 


function jl_plotting_plot_surfaces(
    eq; 
    rho = 8, 
    theta = 8, 
    zeta = nothing, 
    ax = nothing,
    return_data = false, 
    kwargs...
)

kwargs_dict = Dict(pairs(kwargs))

py"""
import desc 
import desc.plotting
import numpy as np  

def plot_surface():
    return desc.plotting.plot_surfaces(
        eq = $eq, 
        rho = $rho, 
        theta = $theta, 
        zeta = $zeta, 
        ax = $ax,
        return_data = $return_data, 
        **$kwargs_dict
    )
"""
output = py"plot_surface"()
plt = pyimport("matplotlib.pyplot")
plt.show()
end 
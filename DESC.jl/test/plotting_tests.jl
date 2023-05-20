#=
These tests are designed to ensure that the plotting functions can be correctly initialized 
=#
@testset "DESC Plotting Function Initialization Tests" begin 

    surface = DESC.FourierRZToroidalSurface(
        R_lmn=[10, 1],
        modes_R=[[0, 0], [1, 0]],  # modes given as [m,n] for each coefficient
        Z_lmn=[0, -1],
        modes_Z=[[0, 0], [-1, 0]],
    )

    pressure = DESC.PowerSeriesProfile(params=[0, 0], modes=[0, 2])
    iota = DESC.PowerSeriesProfile(params=[1, 1.5], modes=[0, 2])

    eq = DESC.Equilibrium(
        surface=surface,
        pressure=pressure,
        iota=iota,
        Psi=1.0,  
        NFP=1,  
        L=7,  
        M=7, 
        N=0, 
        L_grid=12, 
        M_grid=12,  
        N_grid=0, 
        sym=true,  
        output_plots = false
    )

    #=
    Make sure kwargs work
    =#
    DESC.plot_surfaces(
        eq, 
        figsize=(16, 16), 
        title_font_size = 32,
        xlabel_fontsize = 20, 
        ylabel_fontsize = 20, 
        output_plots = false
    )


    obj = DESC.plot_1d(
        eq, 
        "p"; 
        output_plots = false, 
        figsize = (6,6), 

    )
    obj_rep = string(obj)
    @assert occursin("Figure size 862x862", obj_rep)

    obj = DESC.plot_2d(
        eq, 
        "sqrt(g)"; 
        output_plots = false, 
        figsize = (6,6)
    )

    obj_rep = string(obj)
    @assert occursin("Figure size 862x862", obj_rep)


    obj = DESC.plot_3d(
        eq, 
        "|F|"; 
        output_plots = false, 
        figsize = (6,6)
    )

    obj_rep = string(obj)
    @assert occursin("Figure size 862x862", obj_rep)

    obj = DESC.plot_boozer_modes(
        eq;
        output_plots = false, 
        figsize = (6,6)
    )

    obj_rep = string(obj)
    @assert occursin("Figure size 862x862", obj_rep)

    obj = DESC.plot_boundaries(
        [eq];
        output_plots = false, 
        figsize = (6,6)
    )

    obj_rep = string(obj)
    @assert occursin("Figure size 862x862", obj_rep)

    obj = DESC.plot_boundary(
        eq;
        output_plots = false, 
        figsize = (6,6)
    )

    obj_rep = string(obj)
    @assert occursin("Figure size 862x862", obj_rep)

    obj = DESC.plot_coefficients(
        eq, 
        output_plots = false, 
        figsize = (6,6)
    )

    obj_rep = string(obj)
    @assert occursin("Figure size 862x862", obj_rep)


    # Construction of coils not in the API Documentation, 
    # not currently supported. 
    #=  
     obj = DESC.plot_coils(
        eq, 
        output_plots = false, 
        figsize = (6,6)
    )

    obj_rep = string(obj)
    @assert occursin("Figure size 862x862", obj_rep)
    =#

    obj = DESC.plot_comparison(
        [eq, eq], 
        output_plots = false, 
        figsize = (6,6)
    )

    obj_rep = string(obj)
    @assert occursin("Figure size 862x862", obj_rep)


    obj = DESC.plot_field_lines_sfl(
        eq, 
        1, 
        output_plots = false, 
        figsize = (6,6)
    )

    obj_rep = string(obj)
    @assert occursin("Figure size 862x862", obj_rep)


    obj = DESC.plot_fsa(
        eq, 
        "B_theta", 
        with_sqrt_g = false, 
        output_plots = false, 
        figsize = (6,6)
    )

    obj_rep = string(obj)
    @assert occursin("Figure size 862x862", obj_rep)


    grid = DESC.ConcentricGrid(20, 10, 1, node_pattern = "jacobi")
    
    obj = DESC.plot_grid(
        grid,
        output_plots = false, 
        figsize=(6,6))
    
    obj_rep = string(obj)
    
    @assert occursin("Figure size 862x862", obj_rep)

    obj = DESC.plot_qs_error(
        eq; 
        helicity = (1, eq.NFP), 
        log=true, 
        output_plots = false, 
        figsize=(6,6))
    
    obj_rep = string(obj)

    @assert occursin("Figure size 862x862", obj_rep)

    obj = DESC.plot_section(
        eq,
        "J^rho"; 
        output_plots = false, 
        figsize=(6,6))
    
    obj_rep = string(obj)

    @assert occursin("Figure size 862x862", obj_rep)

    obj = DESC.plot_surfaces(
        eq;
        output_plots = false, 
        figsize=(6,6))
    
    obj_rep = string(obj)

    @assert occursin("Figure size 862x862", obj_rep)
   
end 


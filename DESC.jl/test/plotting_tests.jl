@testset "DESC Interactive example" begin 

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
    )
    DESC.plot_3d(
        eq, 
        "|F|"
    )

    DESC.plot_surfaces(
        eq, 
        figsize=(16, 16), 
        title_font_size = 32,
        xlabel_fontsize = 20, 
        ylabel_fontsize = 20, 
        return_data = true
    )

    DESC.plot_boozer_surface(
        eq
    )
end 

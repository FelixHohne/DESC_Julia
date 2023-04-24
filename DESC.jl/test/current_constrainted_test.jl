@testset "current constraint test" begin 

    surface_2D = DESC.FourierRZToroidalSurface(
        R_lmn=[10, -1], 
        Z_lmn=[1],
        modes_R=[0 0; 1 0],
        modes_Z=[-1 0],
        NFP=5
    )
    current = DESC.PowerSeriesProfile(modes=[0], params=[0])
    pressure = DESC.PowerSeriesProfile(modes=[0], params=[0])
    
    eq = DESC.Equilibrium(current=current, pressure=pressure, surface=surface_2D, sym=true)
    eq.change_resolution(L=6, M=6, L_grid=12, M_grid=12)

    # objective = DESC.ObjectiveFunction(DESC.ForceBalance())





end 

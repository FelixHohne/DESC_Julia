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

    objective = DESC.ObjectiveFunction(objectives=(DESC.ForceBalance()))

    constraints = DESC.get_fixed_boundary_constraints(profiles=true, iota=false)

    optimizer = DESC.Optimizer("lsq-exact")

    eq, solver_outputs = eq.solve(
    objective=objective, constraints=constraints, optimizer=optimizer, verbose=3
    )

    DESC.plot_1d(eq, "current")
    DESC.plot_section(eq, "|F|", norm_F=true, log=true);
    DESC.plot_1d(eq, "p")
    DESC.plot_1d(eq, "iota")

    surface_3D = DESC.FourierRZToroidalSurface(
        R_lmn=[10, -1, -0.3, 0.3],  
        Z_lmn=[1, -0.3, -0.3],
        modes_R=[0 0; 1 0; 1 1; -1 -1], 
        modes_Z=[-1 0; -1 1; 1 -1],
        NFP=5
    )

    eq.change_resolution(L=10, M=10, N=6, L_grid=20, M_grid=20, N_grid=12)
    surface_3D.change_resolution(eq.L, eq.M, eq.N)

    objective = DESC.ObjectiveFunction(objectives=(DESC.ForceBalance()))
    constraints = get_fixed_boundary_constraints(profiles=true, iota=false)





end 

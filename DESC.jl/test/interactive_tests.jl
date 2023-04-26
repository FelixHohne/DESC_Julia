@testset "DESC Interactive example" begin 
    DESC.desc_jl_use_gpu_if_available()
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
        Psi=1.0,  # flux (in Webers) within the last closed flux surface
        NFP=1,  # number of field periods
        L=7,  # radial spectral resolution
        M=7,  # poloidal spectral resolution
        N=0,  # toroidal spectral resolution (axisymmetric case, so we don't need any toroidal modes)
        L_grid=12,  # real space radial resolution, slightly oversampled
        M_grid=12,  # real space poloidal resolution, slightly oversampled
        N_grid=0,  # real space toroidal resolution (axisymmetric, so we don't need any grid points toroidally)
        sym=true,  # explicitly enforce stellarator symmetry
    )

    optimizer = DESC.Optimizer("lsq-exact")

    constraints = (
        DESC.FixBoundaryR(), 
        DESC.FixBoundaryZ(), 
        DESC.FixPressure(), 
        DESC.FixIota(), 
        DESC.FixPsi()
    )

    objectives = DESC.ForceBalance()
    obj = DESC.ObjectiveFunction(objectives)
    constraints2 = DESC.get_fixed_boundary_constraints()

    eq.solve(
        verbose=3, 
        ftol=1e-8, 
        maxiter=50, 
        constraints=constraints, 
        optimizer=optimizer, 
        objective=obj
    )


end 
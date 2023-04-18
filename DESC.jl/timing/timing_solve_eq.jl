using Pkg
Pkg.activate("../") 
using DESC 

DESC.desc_jl_use_gpu_if_available()
    surface = DESC.jl_fourierRZToroidalSurface(
        R_lmn=[10, 1],
        modes_R=[[0, 0], [1, 0]],  # modes given as [m,n] for each coefficient
        Z_lmn=[0, -1],
        modes_Z=[[0, 0], [-1, 0]],
    )

    pressure = DESC.jl_profiles_power_series_profile(params=[0, 0], modes=[0, 2])
    iota = DESC.jl_profiles_power_series_profile(params=[1, 1.5], modes=[0, 2])

    eq = DESC.jl_equilibrium(
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

    optimizer = DESC.jl_create_optimizer("lsq-exact")

    constraints = (
        DESC.jl_objective_fix_boundary_r(), 
        DESC.jl_objective_fix_boundary_z(), 
        DESC.jl_objective_fix_pressure(), 
        DESC.jl_objective_fix_iota(), 
        DESC.jl_objective_fix_psi()
    )

    objectives = DESC.jl_objective_force_balance()
    obj = DESC.jl_objective_function(objectives)
    constraints2 = DESC.jl_get_fixed_boundary_constraints()

    DESC.jl_solve_equilibrium(
        eq,
        verbose=3, 
        ftol=1e-8, 
        maxiter=50, 
        constraints=constraints, 
        optimizer=optimizer, 
        objective=obj
    )
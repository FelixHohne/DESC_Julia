using Pkg
Pkg.activate("../") 
using DESC 

function solve_eq_timing()

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
end 

function solve_continuation_automatic_timing()
    DESC.desc_jl_use_gpu_if_available()
    surf = DESC.jl_fourierRZToroidalSurface(
        R_lmn = [1, 0.125, 0.1],
        Z_lmn = [-0.125, -0.1],
        modes_R = [[0, 0], [1, 0], [0, 1]],
        modes_Z = [[-1, 0], [0, -1]],
        NFP = 4
      )
  
    eq = DESC.jl_equilibrium(
        M = 8, 
        N = 8, 
        Psi=0.04, 
        surface=surf
    )

    DESC.jl_solve_continuation_automatic(eq, objective = "force", verbose=3, bdry_step=0.5)
end 

function solve_full_set()



@testset "set eq for QAS_output.h5" begin 
    DESC.desc_jl_use_gpu_if_available()
    surf = DESC.jl_fourierRZToroidalSurface(
    R_lmn = [1, 0.125, 0.1],
    Z_lmn = [-0.125, -0.1],
    modes_R = [[0, 0], [1, 0], [0, 1]],
    modes_Z = [[-1, 0], [0, -1]],
    NFP = 4
    )
  
    eq = DESC.jl_equilibrium(
        M = 8, 
        N = 8, 
        Psi=0.04, 
        surface=surf
    )

    eq = last(DESC.jl_solve_continuation_automatic(eq, objective = "force", verbose=3, bdry_step=0.5))
    DESC.jl_save_equilibrium(eq, "QAS_output_continuation_automatic.hdf5")
    end 
  end 
  
  
  @testset "QAS_output.h5" begin 
      DESC.desc_jl_use_gpu_if_available()
      eq = DESC.jl_load_equilibrium("QAS_output_continuation_automatic.hdf5", "hdf5")
      println(eq)
      eq_fam = DESC.jl_equilibria_family(eq)
      grid = DESC.jl_linear_grid(M=eq.M, N=eq.N, NFP=eq.NFP, rho=[0.6, 0.8, 1.0], sym=true)
  
        for n in 1:eq.M
          objective = DESC.jl_objective_function(
            (DESC.jl_objective_quasisymmetry_two_term(helicity = (1, eq.NFP), grid=grid, normalize=false),
            DESC.jl_objective_aspect_ratio(target=8, weight=1e1, normalize=false)), 
            verbose = 0
          )
  
          R_abs = abs.(eq.surface.R_basis.modes)
          max_Rabs = maximum(R_abs, dims=2)
          # equivalent to np.squeeze. See https://stackoverflow.com/questions/52505760/dropping-singleton-dimensions-in-julia
          # ensures indexing with findall broadcasts correctly
          max_Rabs = dropdims(max_Rabs, dims = tuple(findall(size(max_Rabs) .== 1)...))
          R_elem_mode = eq.surface.R_basis.modes[findall(>(n), max_Rabs), :]
          R_modes = vcat([0 0 0], R_elem_mode)
        
          Z_abs = abs.(eq.surface.Z_basis.modes)
          max_Zabs = maximum(Z_abs, dims=2)
          max_Zabs = dropdims(max_Zabs, dims = tuple(findall(size(max_Zabs) .== 1)...))
          Z_modes = eq.surface.Z_basis.modes[findall(>(n), max_Zabs), :]
  
          constraints = (
            DESC.jl_objective_force_balance(), 
            DESC.jl_objective_fix_boundary_r(modes=R_modes), 
            DESC.jl_objective_fix_boundary_z(modes=Z_modes), 
            DESC.jl_objective_fix_pressure(), 
            DESC.jl_objective_fix_current(), 
            DESC.jl_objective_fix_psi()
          )
  
          println("Constructed DESC constraints")
  
          optimizer = DESC.jl_create_optimizer("lsq-exact")
  
          println("Defined optimizer")
  
          println("Beginning equilibrium optimization")
          eq_new, out = DESC.jl_optimize_equilibrium(
            last(eq_fam); 
            objective=objective, 
            constraints=constraints, 
            optimizer=optimizer, 
            maxiter=1, 
            verbose=3, 
            copy=true, 
            options = Dict(
              "initial_trust_radius" => 0.5,
              "perturb_options" => Dict("verbose" => 0),
              "solve_options" => Dict("verbose" => 0)
            )
          )
          DESC.jl_equilibrium_family_append(eq_fam, eq_new)
          
      end 
  
      DESC.jl_save_equilibrium_family(eq_fam, "qas_julia_test_results.hdf5")
    end 
  end 
  

end 


println("Interactive example timing")
@time solve_eq_timing()

println("Solve continuation automatic timing")
@time solve_continuation_automatic_timing()

println("Solve full set")
@time solve_full_optimization()
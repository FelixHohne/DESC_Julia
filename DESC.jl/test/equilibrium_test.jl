@testset "Eq Construction Test" begin 
eq = DESC.Equilibrium()
new_eq1 = eq.change_resolution(L = 10, M = 3, L_grid = 9, M_grid = 1, NFP = 12)
eq = DESC.Equilibrium(Psi=1.0)
eq = DESC.Equilibrium(Psi=1.0, NFP = 100)
eq = DESC.Equilibrium(Psi=1.0, NFP = 100, L_grid =1, M_grid = 2, N_grid = 3)
eq = DESC.Equilibrium(Psi=1.0, NFP = 100, L_grid =1, M_grid = 2, axis = [0 10 3], pressure = [0 10; 2 5])
eq = DESC.Equilibrium(Psi=1.0, NFP = 100, L_grid =1, M_grid = 2, axis = [0 10 3], pressure = [0 10; 2 5], iota = [0 1; 2 3])
eq = DESC.Equilibrium(Psi=1.0, NFP = 100, L_grid =1, M_grid = 2, axis = [0 10 3], pressure = [0 10; 2 5], iota = [0 1; 2 3], surface = [0 0 0 10 0; 1 1 0 1 1])
end 

plotting = false 
saving = false
@testset "Near Axis Eq Test" begin 
    qsc_eq = DESC.qsc_from_paper("precise QA")
    ntheta = 75
    r = 0.35
    desc_eq = DESC.from_near_axis(
        qsc_eq,  
        r=r,  
        L=8,  
        M=8, 
        N=8, 
        ntheta=ntheta
    )

    eq_fit = desc_eq.copy()  
    constraints = DESC.get_fixed_boundary_constraints(iota=false)
    println(constraints)

    desc_eq.solve(
        verbose=3,
        ftol=1e-2,
        objective="force",
        maxiter=100,
        xtol=1e-6,
        constraints=constraints
    )   

    if saving 
        desc_eq.save("DESC_from_NAE_precise_QA_output.h5")
    end 

    if plotting 
        DESC.plot_comparison(
            [desc_eq, eq_fit], 
            labels=["DESC-solved equilibrium", "NAE surfaces"],
            figsize = (12, 12), 
            theta = 0, 
            colors = ["k", "r"], 
            linestyles = ["-", ":"], 
            lws = [1, 2],
        )
    end 
end 

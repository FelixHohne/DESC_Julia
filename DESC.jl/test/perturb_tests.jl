@testset "Perturb tests" begin 

    eq = DESC.Equilibrium()
    objective = DESC.get_equilibrium_objective()
    constraints = DESC.get_fixed_boundary_constraints(iota = false)

    deltas = Dict("Psi" => float(eq.Psi))
    eq = DESC.perturb(
        eq, 
        objective, 
        constraints, 
        deltas, 
        order = 0, 
        verbose = 2, 
        copy = true
    )

    println("Perturb done")
end 

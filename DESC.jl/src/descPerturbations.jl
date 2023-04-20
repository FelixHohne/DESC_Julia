using PyCall 

function jl_perturbations_get_deltas(
    things1, 
    things2
)
    py"""
    import desc 
    import desc.perturbations
    import numpy as np  

    def get_deltas():
        return desc.perturbations.get_deltas(
            things1 = $things1, 
            things2 = $things2
        )
    """
    output = py"get_deltas"()
end 

function jl_perturbations_perturb(
    eq, 
    objective, 
    constraints, 
    deltas;
    order = 2, 
    tr_ratio = 0.1, 
    weight = "auto", 
    include_f = true, 
    verbose = 1, 
    copy = true
)

py"""
import desc 
import desc.perturbations
import numpy as np  

def perturb():
    return desc.perturbations.perturb(
       eq = $eq, 
       objective = $objective, 
       constraints = $constraints, 
       deltas = $deltas, 
       order = $order, 
       tr_ratio = $tr_ratio, 
       weight = $weight, 
       include_f = $include_f, 
       verbose = $verbose, 
       copy = $copy
    )
"""
output = py"perturb"()
end 


function jl_perturbations_optimal_perturb(
    eq, 
    objective_f, 
    objective_g; 
    dR = false, 
    dZ = false, 
    dL = false, 
    dp = false, 
    di = false, 
    dPsi = false, 
    dRb = false, 
    dZb = false, 
    subspace = nothing, 
    order = 2, 
    tr_ratio = [0.1, 0.25], 
    cutoff = nothing, 
    verbose = 1, 
    copy = true
)

py"""
import desc 
import desc.perturbations
import numpy as np  

def perturb():
    return desc.perturbations.optimal_perturb(
        eq = $eq, 
        objective_f = $objective_f, 
        objective_g = $objective_g, 
        dR = $dR, 
        dZ = $dZ, 
        dL = $dL, 
        dp = $dp, 
        di = $di, 
        dPsi = $dPsi, 
        dRb = $dRb, 
        dZb = $dZb, 
        subspace = $subspace, 
        order = $order, 
        tr_ratio = $tr_ratio, 
        cutoff = $cutoff, 
        verbose = $verbose, 
        copy = $copy
    )
"""
output = py"perturb"()
end 


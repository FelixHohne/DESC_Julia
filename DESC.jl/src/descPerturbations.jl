using PyCall 

function get_deltas(
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

function perturb(
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

    new_deltas = {}
    for item in $deltas:
        key = item 
        value = $deltas[item]

        if isinstance(value, np.ndarray):
            new_value = np.ascontiguousarray(value)
            new_deltas[key] = new_value

        else:
            new_deltas[key] = value 

    if isinstance($weight, np.ndarray):
        new_weight = np.ascontiguousarray($weight)
        assert new_weight.flags['C_CONTIGUOUS']
    else:
        new_weight = $weight

    if isinstance($tr_ratio, np.ndarray):
        new_tr_ratio = np.ascontiguousarray($tr_ratio)
        assert new_tr_ratio.flags['C_CONTIGUOUS']
    else:
        new_tr_ratio = $tr_ratio


    def perturb():
        return desc.perturbations.perturb(
        eq = $eq, 
        objective = $objective, 
        constraints = $constraints, 
        deltas = new_deltas, 
        order = $order, 
        tr_ratio = $new_tr_ratio, 
        weight = new_weight, 
        include_f = $include_f, 
        verbose = $verbose, 
        copy = $copy
        )
    """
    output = py"perturb"()
end 


function optimal_perturb(
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

if isinstance($dR , np.ndarray):
    new_dR = np.ascontiguousarray($dR )
    assert new_dR.flags['C_CONTIGUOUS']
else:
    new_dR = $dR 

if isinstance($dZ, np.ndarray):
    new_dZ = np.ascontiguousarray($dZ)
    assert new_dZ.flags['C_CONTIGUOUS']
else:
    new_dZ = $dZ 

if isinstance($dL , np.ndarray):
    new_dL = np.ascontiguousarray($dL)
    assert new_dL.flags['C_CONTIGUOUS']
else:
    new_dL = $dL 

if isinstance($dp, np.ndarray):
    new_dp = np.ascontiguousarray($dp)
    assert new_dp.flags['C_CONTIGUOUS']
else:
    new_dp = $dp 

if isinstance($di, np.ndarray):
    new_di = np.ascontiguousarray($di)
    assert new_di.flags['C_CONTIGUOUS']
else:
    new_di = $di 

if isinstance($dPsi, np.ndarray):
    new_dPsi = np.ascontiguousarray($dPsi)
    assert new_dPsi.flags['C_CONTIGUOUS']
else:
    new_dPsi = $dPsi

if isinstance($dRb, np.ndarray):
    new_dRb = np.ascontiguousarray($dRb)
    assert new_dRb.flags['C_CONTIGUOUS']
else:
    new_dRb = $dRb 

if isinstance($dZb, np.ndarray):
    new_dZb = np.ascontiguousarray($dZb)
    assert new_dZb.flags['C_CONTIGUOUS']
else:
    new_dZb = $dZb 

if isinstance($subspace, np.ndarray):
    new_subspace = np.ascontiguousarray($subspace)
    assert new_subspace.flags['C_CONTIGUOUS']
else:
    new_subspace = $subspace

def perturb():
    return desc.perturbations.optimal_perturb(
        eq = $eq, 
        objective_f = $objective_f, 
        objective_g = $objective_g, 
        dR = new_dR, 
        dZ = new_dZ, 
        dL = new_dL, 
        dp = new_dp, 
        di = new_di, 
        dPsi = new_dPsi, 
        dRb = new_dRb, 
        dZb = new_dZb, 
        subspace = new_subspace, 
        order = $order, 
        tr_ratio = $tr_ratio, 
        cutoff = $cutoff, 
        verbose = $verbose, 
        copy = $copy
    )
"""
output = py"perturb"()
end 


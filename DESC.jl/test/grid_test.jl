@testset "Grid Construction Test" begin 

    obj = DESC.LinearGrid(M=2, N=3, NFP=4, rho=[0.6, 0.8, 1.0], sym=true)
    obj_rep = string(obj)
    @assert occursin("LinearGrid", obj_rep)

    obj = DESC.QuadratureGrid(2, 3, 3, NFP=4)
    obj_rep = string(obj)
    @assert occursin("QuadratureGrid", obj_rep)

    obj = DESC.ConcentricGrid(2, 3, 3)
    obj_rep = string(obj)
    @assert occursin("ConcentricGrid", obj_rep)

end 

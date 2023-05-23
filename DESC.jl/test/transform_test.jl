@testset "Basis Transform Test" begin 

    grid_1 = DESC.LinearGrid(L = 10, N = 1)
    grid_2 = DESC.LinearGrid(M = 2, N = 2)
    grid_3 = DESC.ConcentricGrid(4, 2, 2)

    basis_1 = DESC.DoubleFourierSeries(1, 1)
    basis_2 = DESC.FourierZernikeBasis(-1, 1, 1)

    transf_11 = DESC.Transform(grid_1, basis_1)
    transf_21 = DESC.Transform(grid_2, basis_1)
    transf_31 = DESC.Transform(grid_3, basis_1)
    transf_32 = DESC.Transform(grid_3, basis_2)
    transf_32b = DESC.Transform(grid_3, basis_2)

    @assert transf_11.eq(transf_21) == false
    @assert transf_31.eq(transf_32) == false 
    @assert transf_32.eq(transf_32b)

end 
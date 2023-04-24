@testset "Basis Construction Test" begin 
    DESC.PowerSeries(1, sym="even")
    DESC.FourierSeries(10, NFP=12, sym="cos")
    DESC.DoubleFourierSeries(10, 9)
    DESC.DoubleFourierSeries(1, 1, NFP=123, sym=false)
    DESC.ZernikePolynomial(10, 1, sym="cos", spectral_indexing="ansi")
    DESC.ZernikePolynomial(-1, 1,  sym="cos", spectral_indexing="ansi")
    DESC.FourierZernikeBasis(1, 2, 3, sym=false, spectral_indexing="ansi")
end 

@testset "Basis Operations Test" begin 
    basis = DESC.PowerSeries(1, sym="even")
    basis.change_resolution(10)
    @assert basis.L == 10
    new_basis = basis.copy(deepcopy = true)
    @assert new_basis.L == 10
    @assert new_basis.eq(basis)
    DESC.basis_evaluate(new_basis, randn(100, 3))
    DESC.basis_evaluate(new_basis, randn(100, 3); modes = [0 0 0])
    new_basis.get_idx(L= 1, M = 2, N = 3, error = false)
    new_basis.save("test_basis.hdf5")
    basis.load("test_basis.hdf5")
end 

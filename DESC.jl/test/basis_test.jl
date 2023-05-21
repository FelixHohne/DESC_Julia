@testset "Basis Construction Test" begin 

    obj = DESC.PowerSeries(1, sym="even")
    obj_rep = string(obj)
    @assert occursin("PowerSeries", obj_rep)

    obj = DESC.FourierSeries(10, NFP=12, sym="cos")
    obj_rep = string(obj)
    @assert occursin("FourierSeries", obj_rep)


    obj = DESC.DoubleFourierSeries(10, 9)
    obj_rep = string(obj)
    @assert occursin("DoubleFourierSeries", obj_rep)

    obj =  DESC.DoubleFourierSeries(1, 1, NFP=123, sym=false)
    obj_rep = string(obj)
    @assert occursin("DoubleFourierSeries", obj_rep)

    obj = DESC.ZernikePolynomial(10, 1, sym="cos", spectral_indexing="ansi")
    obj_rep = string(obj)
    @assert occursin("ZernikePolynomial", obj_rep)

   
    obj = DESC.ZernikePolynomial(-1, 1,  sym="cos", spectral_indexing="ansi")
    obj_rep = string(obj)
    @assert occursin("ZernikePolynomial", obj_rep)

    obj = DESC.FourierZernikeBasis(1, 2, 3, sym=false, spectral_indexing="ansi")
    obj_rep = string(obj)
    @assert occursin("FourierZernikeBasis", obj_rep)

    
    
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

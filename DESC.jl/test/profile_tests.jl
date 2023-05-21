@testset "Grid Construction Test" begin 

    obj = DESC.PowerSeriesProfile(modes=[0, 1, 2, 4], params=[1, 0, -2, 1], sym="auto")
    obj_rep = string(obj)
    @assert occursin("PowerSeriesProfile", obj_rep)

    obj = DESC.SplineProfile(values=[1,2,3])
    obj_rep = string(obj)
    @assert occursin("SplineProfile", obj_rep)

    obj1 = DESC.MTanhProfile()
    obj_rep = string(obj1)
    @assert occursin("MTanhProfile", obj_rep)

    obj = DESC.ScaledProfile(1, obj)
    obj_rep = string(obj)
    @assert occursin("ScaledProfile", obj_rep)

    obj = DESC.FourierZernikeProfile()
    obj_rep = string(obj)
    @assert occursin("FourierZernikeProfile", obj_rep)

    new_obj = DESC.SumProfile(obj1, obj)
    obj_rep = string(new_obj)
    @assert occursin("SumProfile", obj_rep)

    new_obj = DESC.ProductProfile(obj1, obj)
    obj_rep = string(new_obj)
    @assert occursin("ProductProfile", obj_rep)
end 

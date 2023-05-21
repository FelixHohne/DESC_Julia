@testset "Geometry Construction Test" begin 
    obj = DESC.FourierRZCurve(R_n = [-0.2, 10, 0.3], Z_n = [0.3, 0.0, -0.2], NFP=2, sym=false)
    obj_rep = string(obj)
    @assert occursin("FourierRZCurve", obj_rep)

    obj = DESC.FourierXYZCurve(modes=[-1, 0, 1])
    obj_rep = string(obj)
    @assert occursin("FourierXYZCurve", obj_rep)

    obj = DESC.FourierPlanarCurve()
    obj_rep = string(obj)
    @assert occursin("FourierPlanarCurve", obj_rep)

    R0 = 10.0
    a_s = 2.0
    obj = DESC.FourierRZToroidalSurface(
        R_lmn=[R0, a_s], Z_lmn=[-a_s], modes_R=[[0, 0], [1, 0]], modes_Z=[[-1, 0]]
    )

    obj_rep = string(obj)
    @assert occursin("FourierRZToroidalSurface", obj_rep)

    obj = DESC.ZernikeRZToroidalSection(spectral_indexing="ansi")
    obj_rep = string(obj)
    @assert occursin("ZernikeRZToroidalSection", obj_rep)


end 

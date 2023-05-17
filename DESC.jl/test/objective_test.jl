@testset "Construct objectives" begin 

    obj = DESC.AspectRatio()
    obj_rep = string(obj)
    @assert occursin("AspectRatio", obj_rep)

    obj = DESC.BootstrapRedlConsistency()
    obj_rep = string(obj)
    @assert occursin("BootstrapRedlConsistency", obj_rep)

    obj = DESC.CurrentDensity()
    obj_rep = string(obj)
    @assert occursin("CurrentDensity", obj_rep)

    obj = DESC.Elongation()
    obj_rep = string(obj)
    @assert occursin("Elongation", obj_rep)

    obj = DESC.Energy()
    obj_rep = string(obj)
    @assert occursin("Energy", obj_rep)

    obj = DESC.FixAtomicNumber()
    obj_rep = string(obj)
    @assert occursin("FixAtomicNumber", obj_rep)

    obj = DESC.FixAxisR()
    obj_rep = string(obj)
    @assert occursin("FixAxisR", obj_rep)

    obj = DESC.FixAxisZ()
    obj_rep = string(obj)
    @assert occursin("FixAxisZ", obj_rep)

    obj = DESC.FixBoundaryR()
    obj_rep = string(obj)
    @assert occursin("FixBoundaryR", obj_rep)

    obj = DESC.FixBoundaryZ()
    obj_rep = string(obj)
    @assert occursin("FixBoundaryZ", obj_rep)

    obj = DESC.FixCurrent()
    obj_rep = string(obj)
    @assert occursin("FixCurrent", obj_rep)

    obj = DESC.FixElectronDensity()
    obj_rep = string(obj)
    @assert occursin("FixElectronDensity", obj_rep)

    obj = DESC.FixElectronTemperature()
    obj_rep = string(obj)
    @assert occursin("FixElectronTemperature", obj_rep)

    obj = DESC.FixIonTemperature()
    obj_rep = string(obj)
    @assert occursin("FixIonTemperature", obj_rep)

    obj = DESC.FixIota()
    obj_rep = string(obj)
    @assert occursin("FixIota", obj_rep)

    obj = DESC.FixModeR()
    obj_rep = string(obj)
    @assert occursin("FixModeR", obj_rep)

    obj = DESC.FixModeR()
    obj_rep = string(obj)
    @assert occursin("FixModeR", obj_rep)

    obj = DESC.FixModeZ()
    obj_rep = string(obj)
    @assert occursin("FixModeZ", obj_rep)

    obj = DESC.FixPressure()
    obj_rep = string(obj)
    @assert occursin("FixPressure", obj_rep)

    obj = DESC.FixPsi()
    obj_rep = string(obj)
    @assert occursin("FixPsi", obj_rep)

    obj = DESC.FixSumModesR()
    obj_rep = string(obj)
    @assert occursin("FixSumModesR", obj_rep)

    obj = DESC.FixSumModesZ()
    obj_rep = string(obj)
    @assert occursin("FixSumModesZ", obj_rep)

    obj = DESC.FixThetaSFL()
    obj_rep = string(obj)
    @assert occursin("FixThetaSFL", obj_rep)

    obj = DESC.ForceBalance()
    obj_rep = string(obj)
    @assert occursin("ForceBalance", obj_rep)

    obj = DESC.HelicalForceBalance()
    obj_rep = string(obj)
    @assert occursin("HelicalForceBalance", obj_rep)

    obj = DESC.Isodynamicity()
    obj_rep = string(obj)
    @assert occursin("Isodynamicity", obj_rep)

    obj = DESC.MagneticWell()
    obj_rep = string(obj)
    @assert occursin("MagneticWell", obj_rep)

    obj = DESC.MeanCurvature()
    obj_rep = string(obj)
    @assert occursin("MeanCurvature", obj_rep)

    obj = DESC.MercierStability()
    obj_rep = string(obj)
    @assert occursin("MercierStability", obj_rep)

    obj = DESC.PrincipalCurvature()
    obj_rep = string(obj)
    @assert occursin("PrincipalCurvature", obj_rep)

    obj = DESC.QuasisymmetryBoozer()
    obj_rep = string(obj)
    @assert occursin("QuasisymmetryBoozer", obj_rep)

    obj = DESC.QuasisymmetryTwoTerm()
    obj_rep = string(obj)
    @assert occursin("QuasisymmetryTwoTerm", obj_rep)

    obj = DESC.QuasisymmetryTripleProduct()
    obj_rep = string(obj)
    @assert occursin("QuasisymmetryTripleProduct", obj_rep)

    obj = DESC.RadialForceBalance()
    obj_rep = string(obj)
    @assert occursin("RadialForceBalance", obj_rep)

    obj = DESC.ToroidalCurrent()
    obj_rep = string(obj)
    @assert occursin("ToroidalCurrent", obj_rep)

    obj = DESC.Volume()
    obj_rep = string(obj)
    @assert occursin("Volume", obj_rep)

end 
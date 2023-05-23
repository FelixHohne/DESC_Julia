@testset "VMec tests" begin 

    eq = DESC.vmecio_load("wout_SOLOVEV.nc", profile = "current")
    @assert isnothing(eq.iota)
    eq = DESC.vmecio_load("wout_SOLOVEV.nc", profile = "iota")
    @assert isnothing(eq.current)
    DESC.vmecio_save(eq, "desc_wout_solovev")

end 

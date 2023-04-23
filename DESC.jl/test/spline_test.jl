@testset "extract R modes" begin 

  DESC.desc_jl_use_gpu_if_available()

  surf = DESC.jl_fourierRZToroidalSurface(
    R_lmn = [10, 1],
    modes_R = [[0, 0], [1, 0]],
    Z_lmn = [0, -1],
    modes_Z = [[0, 0], [-1, 0]],
  )

  knots = LinRange(0, 1, 20)
  pressure_values = zeros(20)
  iota_values = 1 .+ 1.5 * (knots.^2)

  pressure = DESC.jl_profiles_spline_profile(
    values = pressure_values, knots = knots
  )

  iota = DESC.jl_profiles_spline_profile(
    values = iota_values, 
    knots = knots
  )


end 

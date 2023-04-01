# values are benchmarked against DESC for the same eq

@testset "basic python call from julia" begin 
  py"""
  import numpy as np
  from desc.equilibrium import Equilibrium
  """
end 

@testset "check julia - python data conversion" begin 
  d = DESC.DescData(joinpath(@__DIR__, "QAS_output.h5"));
  out = DESC.jl_equilibrium(d)
end 

# @testset "blah blah lbah" begin 
#   d = DESC.DescData(joinpath(@__DIR__, "QAS_output.h5"));
#   eq = DESC.jl_equilibrium(d)
#   objective = DESC.test_jl_objective_aspect_ratio(eq)
# end 


@testset "construct Desc Optimizer" begin 
  d = DESC.DescData(joinpath(@__DIR__, "QAS_output.h5"));
  eq = DESC.jl_equilibrium(d)
  objective_aspect_ratio = DESC.jl_objective_aspect_ratio(eq)
  objective_current_density = DESC.jl_obective_current_density(eq)
  objective_function = DESC.jl_objective_function((objective_aspect_ratio, objective_current_density), eq, true, "batched", 2)
  optimizer = DESC.jl_create_optimizer("lsq-exact")
  constraints = DESC.jl_create_example_constraints()
  optimize_result = DESC.jl_optimize(optimizer, eq, objective_function, constraints)
end 


@testset "check desc utils" begin
  rtol = 1.0E-6
  d = DESC.DescData(joinpath(@__DIR__, "QAS_output.h5"));

  funcdict = Dict()
  funcdict["R"] = desc_R
  funcdict["Z"] = desc_Z
  funcdict["lambda"] = desc_λ
  funcdict["R_r"] = desc_dRdρ
  funcdict["R_t"] = desc_dRdθ
  funcdict["R_z"] = desc_dRdζ
  funcdict["Z_r"] = desc_dZdρ
  funcdict["Z_t"] = desc_dZdθ
  funcdict["Z_z"] = desc_dZdζ
  funcdict["lambda_r"] = desc_dλdρ
  funcdict["lambda_t"] = desc_dλdθ
  funcdict["lambda_z"] = desc_dλdζ
  funcdict["sqrt(g)"] = desc_jacobian
  funcdict["g_rr"] = desc_g_ρρ
  funcdict["g_tt"] = desc_g_θθ
  funcdict["g_zz"] = desc_g_ζζ
  funcdict["g_rt"] = desc_g_ρθ
  funcdict["g_rz"] = desc_g_ρζ
  funcdict["g_tz"] = desc_g_θζ
  funcdict["R_rr"] = DESC.desc_d2Rdρ2
  funcdict["R_tt"] = DESC.desc_d2Rdθ2
  funcdict["R_zz"] = DESC.desc_d2Rdζ2
  funcdict["R_rt"] = DESC.desc_d2Rdρdθ
  funcdict["R_rz"] = DESC.desc_d2Rdρdζ
  funcdict["R_tz"] = DESC.desc_d2Rdθdζ
  funcdict["Z_rr"] = DESC.desc_d2Zdρ2
  funcdict["Z_tt"] = DESC.desc_d2Zdθ2
  funcdict["Z_zz"] = DESC.desc_d2Zdζ2
  funcdict["Z_rt"] = DESC.desc_d2Zdρdθ
  funcdict["Z_rz"] = DESC.desc_d2Zdρdζ
  funcdict["Z_tz"] = DESC.desc_d2Zdθdζ
  funcdict["lambda_rr"] = DESC.desc_d2λdρ2
  funcdict["lambda_tt"] = DESC.desc_d2λdθ2
  funcdict["lambda_zz"] = DESC.desc_d2λdζ2
  funcdict["lambda_rt"] = DESC.desc_d2λdρdθ
  funcdict["lambda_rz"] = DESC.desc_d2λdρdζ
  funcdict["lambda_tz"] = DESC.desc_d2λdθdζ
  funcdict["sqrt(g)_r"] = DESC.desc_djacobiandρ
  


  @testset "check desc scalar parameters" begin
    rf = readlines(joinpath(@__DIR__, "desc_scalar_test.txt"))
    for line in rf
      dum = split(line)
      func = funcdict[dum[1]]
      ρ = parse(Float64, dum[2])
      θ = parse(Float64, dum[3])
      ζ = parse(Float64, dum[4])
      zc = ZernikeCoordinates(ρ,θ,ζ)
      @test isapprox(func(zc, d), parse(Float64, dum[5]), atol=rtol)
    end
  end

  @testset "check second derivatives" begin
    rf = readlines(joinpath(@__DIR__, "desc_second_deriv.txt"))
    for line in rf
      dum = split(line)
      func = funcdict[dum[1]]
      ρ = parse(Float64, dum[2])
      θ = parse(Float64, dum[3])
      ζ = parse(Float64, dum[4])
      zc = ZernikeCoordinates(ρ,θ,ζ)
  #    println(dum[1])
      @test isapprox(func(zc, d), parse(Float64, dum[5]), atol=rtol)
    end
  end
    
  @testset "test desc vector parameters" begin
    zc = ZernikeCoordinates(0.3, 0.0, 0.251327412)
    @test isapprox(desc_eρ(zc, d),[0.18425619,  0.,         -0.13898179] , rtol=rtol)
    @test isapprox(desc_eθ(zc, d),[ 0.03116248,  0.,         -0.20558328] , rtol=rtol)
    @test isapprox(desc_eζ(zc, d),[-0.18010022,  1.63049995,  0.01049823] , rtol=rtol)
  end

  @testset "test desc flux surface parameters" begin
    zc = ZernikeCoordinates(0.3, 0.0, 0.251327412)
    @test isapprox(desc_ψ(zc, d), 0.00712001, atol=rtol);
    @test isapprox(desc_dψdρ(zc, d), 0.04746671, atol=rtol);
  end
end

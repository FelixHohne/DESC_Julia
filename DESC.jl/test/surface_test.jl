@testset "fourier test" begin
  funcdict = Dict()
  funcdict["|B|"] = (:b, :none)
  funcdict["|B|_r"] = (:b, :ds)
  funcdict["B^theta"] = (:bsupu, :none)
  funcdict["B^zeta"] = (:bsupv, :none)
  funcdict["B_rho"] = (:bsubs, :none)
  funcdict["B_theta"] = (:bsubu, :none)
  funcdict["B_zeta"] = (:bsubv, :none)
  funcdict["B^theta_r"] = (:bsupu, :ds)
  funcdict["B^zeta_r"] = (:bsupv, :ds)
  funcdict["B_rho_r"] = (:bsubs, :ds)
  funcdict["B_theta_r"] = (:bsubu, :ds)
  funcdict["B_zeta_r"] = (:bsubv, :ds)
  funcdict["J^theta"] = (:jsupu, :none)
  funcdict["J^zeta"] = (:jsupv, :none)
  funcdict["J^rho"] = (:jsups, :none)

  rtol = 1.0E-6
  rtol_lo = 1.0E-3
  atol = 1.0E-6
  d = DESC.DescData(joinpath(@__DIR__, "QAS_output.h5"));
  ζmax = 2π/d.nfp
  for ρ in range(0.01,1.0,3)
    Rfour = DESC.descFourierDataArray(ρ, d, :R_lmn, :R_modes)
    Zfour = DESC.descFourierDataArray(ρ, d, :Z_lmn, :Z_modes)
    λfour = DESC.descFourierDataArray(ρ, d, :λ_lmn, :λ_modes)
    for θ in range(0.01,2π,6), ζ in range(0.0, ζmax, 5)
      zc = ZernikeCoordinates(ρ, θ, ζ)
      @test isapprox(desc_R(zc, d), inverseTransform(zc, Rfour), atol=rtol)
      @test isapprox(desc_Z(zc, d), inverseTransform(zc, Zfour), atol=rtol)
      @test isapprox(desc_λ(zc, d), inverseTransform(zc, λfour), atol=rtol)
      @test isapprox(desc_dRdρ(zc, d), inverseTransform(zc, Rfour, deriv=:ds)*(2*ρ), atol=rtol)
      @test isapprox(desc_dZdρ(zc, d), inverseTransform(zc, Zfour, deriv=:ds)*(2*ρ), atol=rtol)
      @test isapprox(desc_dλdρ(zc, d), inverseTransform(zc, λfour, deriv=:ds)*(2*ρ), atol=rtol)
      @test isapprox(desc_dRdθ(zc, d), inverseTransform(zc, Rfour, deriv=:dθ), atol=rtol)
      @test isapprox(desc_dZdθ(zc, d), inverseTransform(zc, Zfour, deriv=:dθ), atol=rtol)
      @test isapprox(desc_dλdθ(zc, d), inverseTransform(zc, λfour, deriv=:dθ), atol=rtol)
      @test isapprox(desc_dRdζ(zc, d), inverseTransform(zc, Rfour, deriv=:dζ), atol=rtol)
      @test isapprox(desc_dZdζ(zc, d), inverseTransform(zc, Zfour, deriv=:dζ), atol=rtol)
      @test isapprox(desc_dλdζ(zc, d), inverseTransform(zc, λfour, deriv=:dζ), atol=rtol)
      #2nd derivative tests
      @test isapprox(DESC.desc_d2Rdρdθ(zc, d), 
                     inverseTransform(zc, Rfour, deriv=:dsdθ)*(2*ρ), atol=rtol)
      @test isapprox(DESC.desc_d2Zdρdθ(zc, d), 
                     inverseTransform(zc, Zfour, deriv=:dsdθ)*(2*ρ), atol=rtol)
      @test isapprox(DESC.desc_d2λdρdθ(zc, d), 
                     inverseTransform(zc, λfour, deriv=:dsdθ)*(2*ρ), atol=rtol)
      @test isapprox(DESC.desc_d2Rdρdζ(zc, d), 
                     inverseTransform(zc, Rfour, deriv=:dsdζ)*(2*ρ), atol=rtol)
      @test isapprox(DESC.desc_d2Zdρdζ(zc, d), 
                     inverseTransform(zc, Zfour, deriv=:dsdζ)*(2*ρ), atol=rtol)
      @test isapprox(DESC.desc_d2λdρdζ(zc, d), 
                     inverseTransform(zc, λfour, deriv=:dsdζ)*(2*ρ), atol=rtol)
      @test isapprox(DESC.desc_d2Rdθ2(zc, d), 
                     inverseTransform(zc, Rfour, deriv=:dθdθ), atol=rtol)
      @test isapprox(DESC.desc_d2Zdθ2(zc, d), 
                     inverseTransform(zc, Zfour, deriv=:dθdθ), atol=rtol)
      @test isapprox(DESC.desc_d2λdθ2(zc, d), 
                     inverseTransform(zc, λfour, deriv=:dθdθ), atol=rtol)
      @test isapprox(DESC.desc_d2Rdθdζ(zc, d), 
                     inverseTransform(zc, Rfour, deriv=:dθdζ), atol=rtol)
      @test isapprox(DESC.desc_d2Zdθdζ(zc, d), 
                     inverseTransform(zc, Zfour, deriv=:dθdζ), atol=rtol)
      @test isapprox(DESC.desc_d2λdθdζ(zc, d), 
                     inverseTransform(zc, λfour, deriv=:dθdζ), atol=rtol)
      @test isapprox(DESC.desc_d2Rdζ2(zc, d), 
                     inverseTransform(zc, Rfour, deriv=:dζdζ), atol=rtol)
      @test isapprox(DESC.desc_d2Zdζ2(zc, d), 
                     inverseTransform(zc, Zfour, deriv=:dζdζ), atol=rtol)
      @test isapprox(DESC.desc_d2λdζ2(zc, d), 
                     inverseTransform(zc, λfour, deriv=:dζdζ), atol=rtol)
    end
  end
  @testset "surface tests" begin
    for ρ in range(0.01, 1.0, 3)
      descsurf = DescSurface(ρ^2, d, simple=true)#note surface uses s
      for θ in range(0.01,π,3), ζ in range(0.0, ζmax, 3)
        zc = ZernikeCoordinates(ρ, θ, ζ)
        @test isapprox(desc_eρ(zc, d), desc_eρ(zc, descsurf), atol=rtol)
        @test isapprox(desc_eθ(zc, d), desc_eθ(zc, descsurf), atol=rtol)
        @test isapprox(desc_eζ(zc, d), desc_eζ(zc, descsurf), atol=rtol)
        @test isapprox(desc_g_ρρ(zc, d), desc_g_ρρ(zc, descsurf), atol=rtol)
        @test isapprox(desc_g_θθ(zc, d), desc_g_θθ(zc, descsurf), atol=rtol)
        @test isapprox(desc_g_ζζ(zc, d), desc_g_ζζ(zc, descsurf), atol=rtol)
        @test isapprox(desc_g_ρθ(zc, d), desc_g_ρθ(zc, descsurf), atol=rtol)
        @test isapprox(desc_g_ρζ(zc, d), desc_g_ρζ(zc, descsurf), atol=rtol)
        @test isapprox(desc_g_θζ(zc, d), desc_g_θζ(zc, descsurf), atol=rtol)
        @test isapprox(desc_jacobian(zc, d), desc_jacobian(zc, descsurf), atol=rtol)
        @test isapprox(DESC.desc_d2Rdρ2(zc, d), 
                       DESC.desc_dρdρ(zc, :rmn, descsurf), atol=rtol)
        @test isapprox(DESC.desc_d2Zdρ2(zc, d), 
                       DESC.desc_dρdρ(zc, :zmn, descsurf), atol=rtol)
        @test isapprox(DESC.desc_d2λdρ2(zc, d), 
                       DESC.desc_dρdρ(zc, :λmn, descsurf), atol=rtol)
        @test isapprox(DESC.desc_djacobiandρ(zc,d),
                       DESC.desc_djacobiandρ(zc, descsurf), atol=rtol)
      end
    end
  end
  #these values benchmarked against DESC
  @testset "full surface test" begin
    ρ = 0.3
    descsurf = DescSurface(ρ^2, d, simple=false)

    #iota = desc_iota(descsurf)
    @test isapprox(descsurf.iota[1], -0.41639632, rtol=rtol_lo)
    @test isapprox(descsurf.iota[2]*2*ρ, -0.20832255,rtol=rtol_lo)
    @test isapprox(descsurf.vol[1], 0.2798715051806796,rtol=rtol_lo)
    @test isapprox(descsurf.vol[2]*2*ρ, 1.859755881358821,rtol=rtol_lo)
    for θ in range(0.01,π,15), ζ in range(0.0, ζmax, 15)
      fc = FluxCoordinates(ρ^2, θ, ζ)
      @test isapprox(DESC.desc_jacobian(fc, descsurf), surface_get(fc, descsurf, :g), atol=rtol_lo)
    end

    rf = readlines(joinpath(@__DIR__, "desc_b_and_j.txt"))
    for line in rf
      dum = split(line)
      func = funcdict[dum[1]]
      θ = parse(Float64, dum[3])
      ζ = parse(Float64, dum[4])
      k = parse(Float64, dum[5])
      fc = FluxCoordinates(ρ^2,θ,ζ)
      zc = ZernikeCoordinates(ρ,θ,ζ)
      q = func[1]
      deriv = func[2]
      if abs(k) < 1.0E-6
        tolname = :atol
        tolval = atol
      else
        tolname = :rtol
        tolval = rtol_lo
      end
      #the jsups values are not accurate enough, and it's not expected for them to be
      if q == :jsups 
        continue
      end

      k2 = 0
      if q == :bsubs && deriv == :ds
        k2 = DESC.desc_dbsubρdρ(zc, descsurf)
      elseif q == :jsups  
        k2 = surface_get(fc, descsurf, q, deriv=deriv)/2/ρ
      elseif deriv == :ds || q == :bsubs
        k2 = surface_get(fc, descsurf, q, deriv=deriv)*2*ρ
      else 
        k2 = surface_get(fc, descsurf, q, deriv=deriv)
      end
#      println(func," ",θ," ",ζ," ",k," ",k2)
      @test eval(:(isapprox($k, $k2, $tolname = $tolval)))
    end

  end
end

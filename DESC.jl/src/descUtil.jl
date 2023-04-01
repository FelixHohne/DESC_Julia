#may want to preevaluate some of these
function desc_R(zc::ZernikeCoordinates, desc::DescData)
  return zernike_eval(zc, desc.R_lmn, desc.R_modes, desc.nfp)
end

function desc_Z(zc::ZernikeCoordinates, desc::DescData)
  return zernike_eval(zc, desc.Z_lmn, desc.Z_modes, desc.nfp)
end

function desc_λ(zc::ZernikeCoordinates, desc::DescData)
  return zernike_eval(zc, desc.λ_lmn, desc.λ_modes, desc.nfp)
end

function desc_dRdρ(zc::ZernikeCoordinates, desc::DescData)
  return zernike_eval(zc, desc.R_lmn, desc.R_modes, desc.nfp, deriv="dρ")
end

function desc_dRdθ(zc::ZernikeCoordinates, desc::DescData)
  return zernike_eval(zc, desc.R_lmn, desc.R_modes, desc.nfp, deriv="dθ")
end

function desc_dRdζ(zc::ZernikeCoordinates, desc::DescData)
  return zernike_eval(zc, desc.R_lmn, desc.R_modes, desc.nfp, deriv="dζ")
end

function desc_dZdρ(zc::ZernikeCoordinates, desc::DescData)
  return zernike_eval(zc, desc.Z_lmn, desc.Z_modes, desc.nfp, deriv="dρ")
end

function desc_dZdθ(zc::ZernikeCoordinates, desc::DescData)
  return zernike_eval(zc, desc.Z_lmn, desc.Z_modes, desc.nfp, deriv="dθ")
end

function desc_dZdζ(zc::ZernikeCoordinates, desc::DescData)
  return zernike_eval(zc, desc.Z_lmn, desc.Z_modes, desc.nfp, deriv="dζ")
end

function desc_dλdρ(zc::ZernikeCoordinates, desc::DescData)
  return zernike_eval(zc, desc.λ_lmn, desc.λ_modes, desc.nfp, deriv="dρ")
end

function desc_dλdθ(zc::ZernikeCoordinates, desc::DescData)
  return zernike_eval(zc, desc.λ_lmn, desc.λ_modes, desc.nfp, deriv="dθ")
end

function desc_dλdζ(zc::ZernikeCoordinates, desc::DescData)
  return zernike_eval(zc, desc.λ_lmn, desc.λ_modes, desc.nfp, deriv="dζ")
end

### second R derivatives

function desc_d2Rdρ2(zc::ZernikeCoordinates, desc::DescData)
  return zernike_eval(zc, desc.R_lmn, desc.R_modes, desc.nfp, deriv="dρdρ")
end

function desc_d2Rdθ2(zc::ZernikeCoordinates, desc::DescData)
  return zernike_eval(zc, desc.R_lmn, desc.R_modes, desc.nfp, deriv="dθdθ")
end

function desc_d2Rdζ2(zc::ZernikeCoordinates, desc::DescData)
  return zernike_eval(zc, desc.R_lmn, desc.R_modes, desc.nfp, deriv="dζdζ")
end

function desc_d2Rdρdθ(zc::ZernikeCoordinates, desc::DescData)
  return zernike_eval(zc, desc.R_lmn, desc.R_modes, desc.nfp, deriv="dρdθ")
end

function desc_d2Rdρdζ(zc::ZernikeCoordinates, desc::DescData)
  return zernike_eval(zc, desc.R_lmn, desc.R_modes, desc.nfp, deriv="dρdζ")
end

function desc_d2Rdθdζ(zc::ZernikeCoordinates, desc::DescData)
  return zernike_eval(zc, desc.R_lmn, desc.R_modes, desc.nfp, deriv="dθdζ")
end

### second Z derivatives

function desc_d2Zdρ2(zc::ZernikeCoordinates, desc::DescData)
  return zernike_eval(zc, desc.Z_lmn, desc.Z_modes, desc.nfp, deriv="dρdρ")
end

function desc_d2Zdθ2(zc::ZernikeCoordinates, desc::DescData)
  return zernike_eval(zc, desc.Z_lmn, desc.Z_modes, desc.nfp, deriv="dθdθ")
end

function desc_d2Zdζ2(zc::ZernikeCoordinates, desc::DescData)
  return zernike_eval(zc, desc.Z_lmn, desc.Z_modes, desc.nfp, deriv="dζdζ")
end

function desc_d2Zdρdθ(zc::ZernikeCoordinates, desc::DescData)
  return zernike_eval(zc, desc.Z_lmn, desc.Z_modes, desc.nfp, deriv="dρdθ")
end

function desc_d2Zdρdζ(zc::ZernikeCoordinates, desc::DescData)
  return zernike_eval(zc, desc.Z_lmn, desc.Z_modes, desc.nfp, deriv="dρdζ")
end

function desc_d2Zdθdζ(zc::ZernikeCoordinates, desc::DescData)
  return zernike_eval(zc, desc.Z_lmn, desc.Z_modes, desc.nfp, deriv="dθdζ")
end

### second λ derivatives

function desc_d2λdρ2(zc::ZernikeCoordinates, desc::DescData)
  return zernike_eval(zc, desc.λ_lmn, desc.λ_modes, desc.nfp, deriv="dρdρ")
end

function desc_d2λdθ2(zc::ZernikeCoordinates, desc::DescData)
  return zernike_eval(zc, desc.λ_lmn, desc.λ_modes, desc.nfp, deriv="dθdθ")
end

function desc_d2λdζ2(zc::ZernikeCoordinates, desc::DescData)
  return zernike_eval(zc, desc.λ_lmn, desc.λ_modes, desc.nfp, deriv="dζdζ")
end

function desc_d2λdρdθ(zc::ZernikeCoordinates, desc::DescData)
  return zernike_eval(zc, desc.λ_lmn, desc.λ_modes, desc.nfp, deriv="dρdθ")
end

function desc_d2λdρdζ(zc::ZernikeCoordinates, desc::DescData)
  return zernike_eval(zc, desc.λ_lmn, desc.λ_modes, desc.nfp, deriv="dρdζ")
end

function desc_d2λdθdζ(zc::ZernikeCoordinates, desc::DescData)
  return zernike_eval(zc, desc.λ_lmn, desc.λ_modes, desc.nfp, deriv="dθdζ")
end

### Jacobian and metric tensor elements
function desc_eρ(zc::ZernikeCoordinates, desc::DescData)
  return SVector(desc_dRdρ(zc, desc),
                  0,
                  desc_dZdρ(zc, desc))
end

function desc_eθ(zc::ZernikeCoordinates, desc::DescData)
  return SVector(desc_dRdθ(zc, desc),
                  0,
                  desc_dZdθ(zc, desc))
end

function desc_eζ(zc::ZernikeCoordinates, desc::DescData)
  return SVector(desc_dRdζ(zc, desc),
                  desc_R(zc, desc),
                  desc_dZdζ(zc, desc))
end

#rho derivatives of basis vectors
function desc_deρdρ(zc::ZernikeCoordinates, desc::DescData)
  return SVector(desc_d2Rdρ2(zc, desc), 0.0, desc_d2Zdρ2(zc, desc))
end

function desc_deθdρ(zc::ZernikeCoordinates, desc::DescData)
  return SVector(desc_d2Rdρdθ(zc, desc), 0.0, desc_d2Zdρdθ(zc, desc))
end

function desc_deζdρ(zc::ZernikeCoordinates, desc::DescData)
  return SVector(desc_d2Rdρdζ(zc, desc), 
                 desc_dRdρ(zc, desc),  
                 desc_d2Zdρdζ(zc, desc))
end


function desc_jacobian(zc::ZernikeCoordinates, desc::DescData)
  eρ = desc_eρ(zc, desc)
  eθ = desc_eθ(zc, desc)
  eζ = desc_eζ(zc, desc)
  return dot(eρ, cross(eθ, eζ))
end

#jacobian derivative
#note stuff like this can be refactored away eventually, since
#it's identical to the function that uses the surface
function desc_djacobiandρ(zc::ZernikeCoordinates, desc::DescData)
  deρdρ = desc_deρdρ(zc, desc)
  deθdρ = desc_deθdρ(zc, desc)
  deζdρ = desc_deζdρ(zc, desc)
  eρ = desc_eρ(zc, desc)
  eθ = desc_eθ(zc, desc)
  eζ = desc_eζ(zc, desc)
  v = dot(deρdρ, cross(eθ, eζ))
  v += dot(eρ, cross(deθdρ, eζ))
  v += dot(eρ, cross(eθ, deζdρ))
  return v
end
  

function desc_g_ρρ(zc::ZernikeCoordinates, desc::DescData)
  eρ = desc_eρ(zc, desc)
  return dot(eρ, eρ)
end

function desc_g_θθ(zc::ZernikeCoordinates, desc::DescData)
  eθ = desc_eθ(zc, desc)
  return dot(eθ, eθ)
end

function desc_g_ζζ(zc::ZernikeCoordinates, desc::DescData)
  eζ = desc_eζ(zc, desc)
  return dot(eζ, eζ)
end

function desc_g_ρθ(zc::ZernikeCoordinates, desc::DescData)
  return dot(desc_eρ(zc, desc), desc_eθ(zc, desc))
end

function desc_g_ρζ(zc::ZernikeCoordinates, desc::DescData)
  return dot(desc_eρ(zc, desc), desc_eζ(zc, desc))
end

function desc_g_θζ(zc::ZernikeCoordinates, desc::DescData)
  return dot(desc_eθ(zc, desc), desc_eζ(zc, desc))
end


### Flux functions
function desc_ψ(ρ::Real, desc::DescData)
  return ρ^2 * desc.ψ /(2π)
end
function desc_ψ(zc::ZernikeCoordinates, desc::DescData)
  return desc_ψ(zc.ρ, desc)
end

function desc_dψdρ(ρ::Real, desc::DescData)
  return ρ * desc.ψ / π
end

function desc_dψdρ(zc::ZernikeCoordinates, desc::DescData)
  return desc_dψdρ(zc.ρ, desc)
end

function desc_current(zc::ZernikeCoordinates, desc::DescData)
  ρ = zc.ρ
  return desc_current(ρ, desc)
end

function desc_current(ρ::Real, desc::DescData)
  v = sum( desc.current .* ρ.^desc.current_modes)
  return v
end

function desc_dcurrentdρ(zc::ZernikeCoordinates, desc::DescData)
  ρ = zc.ρ
  return desc_dcurrentdρ(ρ, desc)
end

function desc_dcurrentdρ(ρ::Real, desc::DescData)
  v = sum( desc.current .* desc.current_modes.*ρ.^(desc.current_modes .- 1))
  return v
end
  

function desc_pressure(zc::ZernikeCoordinates, desc::DescData)
  ρ = zc.ρ
  return desc_pressure(ρ, desc)
end

function desc_pressure(ρ::Real, desc::DescData)
  v = sum( desc.pres .* ρ.^desc.pres_modes)
  return v
end

function desc_dpressuredρ(zc::ZernikeCoordinates, desc::DescData)
  ρ = zc.ρ
  return desc_dpressuredρ(ρ, desc)
end

function desc_dpressuredρ(ρ::Real, desc::DescData)
  v = sum( desc.pres .* desc.pres_modes.*ρ.^(desc.pres_modes .-1 ))
  return v
end


#iota is a mess, we need to get various quantities on a grid
function desc_iota(zc::ZernikeCoordinates, desc::DescData; θres = 128, ζres = 128)
  ρ = zc.ρ
  return desc_iota(ρ, desc; θres = θres, ζres = ζres)
end

function desc_iota(ρ::Real, desc::DescData; θres = 128, ζres = 128)
  println("Warning: it is preferred to create a surface first which will calculate iota")
  dθ = 2π/(θres+1)
  dζ = 2π/desc.nfp/(ζres+1)
  dA = dθ*dζ
  nmat = zeros(Float64, θres, ζres)
  dmat = zeros(Float64, θres, ζres)
  #note when surface averaging, the division by the jacobian drops out
  for (θi,θ) in enumerate(range(0,2π-dθ,θres))
    for (ζi, ζ) in enumerate(range(0, (2π/desc.nfp-dζ), ζres))
      zc = ZernikeCoordinates(ρ,θ,ζ)
      g = desc_jacobian(zc, desc)
      λ_ζ = desc_dλdζ(zc, desc)
      g_θθ = desc_g_θθ(zc, desc)
      λ_θ = desc_dλdθ(zc, desc)
      g_θζ = desc_g_θζ(zc, desc)
      #the formula divides by the jacobian, but then multiplies by it again to do the surface average
      nmat[θi,ζi] = ((λ_ζ * g_θθ) - (1 + λ_θ)*g_θζ)/g
      dmat[θi,ζi] = g_θθ/g
    end
  end
  #return nmat, dmat
  return sum(nmat) / (sum(dmat))
end


mutable struct DescData{FT,IT} <: PlasmaEquilibriumToolkit.AbstractMagneticEquilibrium
  #solved::IT
  L_grid::IT
  M_grid::IT
  N_grid::IT
  nfp::IT
  R_lmn::Array{FT,1}
  Z_lmn::Array{FT,1}
  λ_lmn::Array{FT,1}
  R_modes::Array{IT,2}
  Z_modes::Array{IT,2}
  λ_modes::Array{IT,2}
  current::Array{FT,1} #in vmec this is jcuru and jcurv, not sure how to best do this
  current_modes::Array{FT,1} #in vmec this is jcuru and jcurv, not sure how to best do this
  pres::Array{FT,1} #name consistent with vmec
  pres_modes::Array{FT,1} #name consistent with vmec
  ψ::FT #total toroidal flux
end


function DescData(filename::AbstractString; eq = -1)

  desc_h5 = h5open(filename)
  if "_equilibria" in keys(desc_h5)
    ml = length(desc_h5["_equilibria"]) - 2
    if eq < 0 || eq >= ml
      eql = string(ml)
    else
      eql = string(eq)
    end
    desc_h5 = desc_h5["_equilibria"][eql]
  end

  #solved = read(desc_h5["_solved"])
  L_grid = read(desc_h5["_L_grid"])
  M_grid = read(desc_h5["_M_grid"])
  N_grid = read(desc_h5["_N_grid"])
  IT = typeof(L_grid)
  nfp = read(desc_h5["_NFP"])
  nfp = IT(nfp)
  R_lmn = read(desc_h5["_R_lmn"])
  Z_lmn = read(desc_h5["_Z_lmn"])
  λ_lmn = read(desc_h5["_L_lmn"])
  R_modes = read(desc_h5["_R_basis"]["_modes"])
  Z_modes = read(desc_h5["_Z_basis"]["_modes"])
  λ_modes = read(desc_h5["_L_basis"]["_modes"])
  current = read(desc_h5["_current"]["_params"])
  current_modes = read(desc_h5["_current"]["_basis"]["_modes"])[1,:]
  pres = read(desc_h5["_pressure"]["_params"])
  pres_modes = read(desc_h5["_pressure"]["_basis"]["_modes"])[1,:]
  ψ = read(desc_h5["_Psi"])
  
  
  return DescData(L_grid, M_grid, N_grid, nfp,
                  R_lmn, Z_lmn, λ_lmn, 
                  R_modes, Z_modes, λ_modes,
                  current, current_modes,
                  pres, pres_modes, ψ)
end

struct ZernikeCoordinates{T <: Real, A <: Real} <: PlasmaEquilibriumToolkit.AbstractMagneticCoordinates
  ρ::T
  θ::A
  ζ::A
  ZernikeCoordinates{T,A}(ρ::T, θ::A, ζ::A) where {T,A} = new(ρ, θ, ζ)
end

function ZernikeCoordinates(ρ,θ,ζ)
  ρ2, θ2, ζ2 = promote(ρ, θ, ζ)
  return ZernikeCoordinates{typeof(ρ2), typeof(θ2)}(ρ2, θ2, ζ2)
end

mutable struct DescSurface{T} <: PlasmaEquilibriumToolkit.AbstractMagneticSurface
  rmn::SurfaceFourierArray{T}
  zmn::SurfaceFourierArray{T}
  λmn::SurfaceFourierArray{T}


  r::Union{Nothing, Interpolations.Extrapolation}  
  z::Union{Nothing, Interpolations.Extrapolation}  
  λ::Union{Nothing, Interpolations.Extrapolation} 
  g::Union{Nothing, Interpolations.Extrapolation}
  b::Union{Nothing, Interpolations.Extrapolation}
  bsubs::Union{Nothing, Interpolations.Extrapolation}
  bsubu::Union{Nothing, Interpolations.Extrapolation}
  bsubv::Union{Nothing, Interpolations.Extrapolation}
  bsupu::Union{Nothing, Interpolations.Extrapolation}
  bsupv::Union{Nothing, Interpolations.Extrapolation}
  curru::Union{Nothing, Interpolations.Extrapolation}
  currv::Union{Nothing, Interpolations.Extrapolation}
  jsupu::Union{Nothing, Interpolations.Extrapolation}
  jsupv::Union{Nothing, Interpolations.Extrapolation}
  jsups::Union{Nothing, Interpolations.Extrapolation}

  drds::Union{Nothing, Interpolations.Extrapolation}
  dzds::Union{Nothing, Interpolations.Extrapolation}
  dλds::Union{Nothing, Interpolations.Extrapolation}
  dbds::Union{Nothing, Interpolations.Extrapolation}
  dgds::Union{Nothing, Interpolations.Extrapolation}
  dbsubsds::Union{Nothing, Interpolations.Extrapolation}
  dbsubuds::Union{Nothing, Interpolations.Extrapolation}
  dbsubvds::Union{Nothing, Interpolations.Extrapolation}
  dbsupuds::Union{Nothing, Interpolations.Extrapolation}
  dbsupvds::Union{Nothing, Interpolations.Extrapolation}
  dcurruds::Union{Nothing, Interpolations.Extrapolation}
  dcurrvds::Union{Nothing, Interpolations.Extrapolation}
  djsupuds::Union{Nothing, Interpolations.Extrapolation}
  djsupvds::Union{Nothing, Interpolations.Extrapolation}
  djsupsds::Union{Nothing, Interpolations.Extrapolation}



  phi::SVector{3,T}
  pres::SVector{2,T}
  current::SVector{2,T}
  iota::SVector{2,T}
  vol::SVector{2,T}

  s::T
  ρ::T

  nfp::Int
  mpol::Int
  ntor::Int
  mnmax::Int
  mnmax_nyq::Int #note there is no nyquist sampling difference in desc

  xm::Vector{T}
  xn::Vector{T}
  xm_nyq::Vector{T}
  xn_nyq::Vector{T}
end

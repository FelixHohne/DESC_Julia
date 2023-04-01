#Note, right now there is no accurate transformation between flux and zernike coordinates
#in the ideal case we would convert everything between the two that will take
#the angle into account. Some functions may need to be rewritten in this case
#right now the assumption is that the Flux θ is equal to the Zernike θ which may
#not be the case for something like VMEC coordinates


function descFourierDataArray(ρ::T, desc::DescData{T}, field::Symbol,
                             desc_modes::Symbol) where {T <: AbstractFloat}
  data = getfield(desc, field)
  modes = getfield(desc, desc_modes)

  #construct a dictionary to relate m,n number to array.  Multiple
  #zernike modes will interact with each array, and each
  #zernike mode will contribute to the + and - n modes
  modeDict = Dict()
  modeDict[(0,0)] = 1#force this always so everything aligns
  larr = modes[1,:]
  marr = modes[2,:]
  narr = modes[3,:]
  count = 1
  for i in 1:length(larr)
    m = abs(marr[i])
    n = abs(narr[i])*desc.nfp
    if !((m,n) in keys(modeDict))
      count += 1
      modeDict[(m,n)] = count
    end
    if !((m,-n) in keys(modeDict))
      count += 1
      modeDict[(m,-n)] = count
    end
  end
  n_modes = count
      
  m_vec = Vector{T}(undef, n_modes)
  n_vec = Vector{T}(undef, n_modes)
  cos_vec = zeros(T, n_modes)
  sin_vec = zeros(T, n_modes)
  dcos_vec = zeros(T, n_modes)
  dsin_vec = zeros(T, n_modes)
  d2cos_vec = zeros(T, n_modes)
  d2sin_vec = zeros(T, n_modes)

  #create the m and n arrays
  for k in modeDict
    m_vec[k[2]] = k[1][1]
    n_vec[k[2]] = k[1][2]
  end
  
  #four options based on signs of m and n
  for i in 1:length(marr)
    l = larr[i]
    msigned = marr[i]
    nsigned = narr[i]
    m = abs(msigned)
    n = abs(nsigned)*desc.nfp
    q = 0.5 * DESC.zernike_r(ρ,l,m)*data[i]
    dqdρ = 0.5 * DESC.zernike_r(ρ,l,m,deriv=1)*data[i]
    dq2 = 0.5 * DESC.zernike_r(ρ,l,m,deriv=2)*data[i]
    # ds = 2ρdρ
    if ρ == 0
      ρ = 0.001
    end
    dq = dqdρ / (2 * ρ) #turn into dq/ds
    #todo at some point verify this against a vmec eq or something
    dq2 = dq2/(4*ρ^2) - dqdρ/4/ρ^3

    j1 = modeDict[(m,n)]
    j2 = modeDict[(m,-n)]
    #case: m pos and n pos
    #corresponds to cos(θ)*cos(ζ)
    #1/2 cos(θ + ζ) + cos(θ - ζ)
    if msigned >= 0 && nsigned >= 0
      cos_vec[j1] += q
      cos_vec[j2] += q
      dcos_vec[j1] += dq
      dcos_vec[j2] += dq
      d2cos_vec[j1] += dq2
      d2cos_vec[j2] += dq2
    #case: m neg and n neg
    #corresonds to sin(θ)*sin(ζ)
    #1/2 cos(θ + ζ) - cos(θ - ζ)
    #note default is -ζ so it's the positive m,n that is negative here
    #somehow the signs are opposite of what I thought, but this
    #causes the tests to pass.
    elseif msigned < 0 && nsigned < 0
      cos_vec[j1] += q
      cos_vec[j2] -= q
      dcos_vec[j1] += dq
      dcos_vec[j2] -= dq
      d2cos_vec[j1] += dq2
      d2cos_vec[j2] -= dq2
    #case: m pos and n neg
    #corresonds to cos(θ)*sin(ζ)
    #1/2 sin(θ + ζ) - sin(θ - ζ)
    elseif msigned >= 0 && nsigned < 0
      sin_vec[j1] -= q
      sin_vec[j2] += q
      dsin_vec[j1] -= dq
      dsin_vec[j2] += dq
      d2sin_vec[j1] -= dq2
      d2sin_vec[j2] += dq2
    #case: m neg and n pos
    #corresonds to sin(θ)*cos(ζ)
    #1/2 sin(θ + ζ) + sin(θ - ζ)
    else
      sin_vec[j1] += q
      sin_vec[j2] += q
      dsin_vec[j1] += dq
      dsin_vec[j2] += dq
      d2sin_vec[j1] += dq2
      d2sin_vec[j2] += dq2
    end
  end
  return StructArray{SurfaceFourierData{T}}((m_vec, n_vec, cos_vec, sin_vec,
                                            dcos_vec, dsin_vec, 
                                            d2cos_vec, d2sin_vec))
end

#can really be a multiple dispatch call to "make_surface_interpolation" in PET
function make_desc_surface_interpolation(f::Function,
                                         dfds::Function,
                                         descsurf::DescSurface,
                                         θres::Int,
                                         ζres::Int)

  nfp = descsurf.nfp
  θs = range(0, length = θres, step = 2π/θres)
  ζs = range(0, length = ζres, step = 2π/nfp/ζres)
  knots = (θs, ζs)
  itp_types = (BSpline(Cubic(Periodic(OnCell()))),
                 BSpline(Cubic(Periodic(OnCell()))))
  itp = (f) -> scale(interpolate(f, itp_types), knots...)
  extp = (f) -> extrapolate(itp(f), (Periodic(), Periodic()))
  T = typeof(f(FluxCoordinates(0.1, 0.1, 0.1), descsurf))
  field = Array{T, 2}(undef, (length(θs), length(ζs)))
  deriv = similar(field)
  for (j, ζ) in enumerate(ζs)
    for (i, θ) in enumerate(θs)
      field[i, j] = f(FluxCoordinates(descsurf.s, θ, ζ), descsurf) 
      deriv[i, j] = dfds(FluxCoordinates(descsurf.s, θ, ζ), descsurf) 
    end
  end
  return (extp(field), extp(deriv))
end

# Helper function for surface_get, will add to it as needed
function get_function_and_deriv(q::Symbol)
  if q == :r
    return(DESC.desc_R, DESC.desc_dRds)
  elseif q == :z
    return(DESC.desc_Z, DESC.desc_dZds)
  elseif q == :λ
    return(DESC.desc_λ, DESC.desc_dλds)
  elseif q == :g
    return(DESC.desc_jacobian, DESC.desc_djacobiands)
  elseif q == :b
    return(DESC.desc_bmag, DESC.desc_dbmagds)
  elseif q == :bsupu
    return(DESC.desc_bsupu, DESC.desc_dbsupuds)
  elseif q == :bsupv
    return(DESC.desc_bsupv, DESC.desc_dbsupvds)
  elseif q == :bsubu
    return(DESC.desc_bsubu, DESC.desc_dbsubuds)
  elseif q == :bsubv
    return(DESC.desc_bsubv, DESC.desc_dbsubvds)
  elseif q == :bsubs
    return(DESC.desc_bsubs, DESC.desc_dbsubsds)
  elseif q == :jsups
    return(DESC.desc_jsups, DESC.unknown_deriv)
  elseif q == :jsupu
    return(DESC.desc_jsupu, DESC.unknown_deriv)
  elseif q == :jsupv
    return(DESC.desc_jsupv, DESC.unknown_deriv)
  end
end

function PlasmaEquilibriumToolkit.surface_get(x::C,
                                              descsurf::DescSurface,
                                              quantity::Symbol;
                                              minθres = 128,
                                              minζres = 128,
                                              deriv::Symbol=:none) where {C <: AbstractMagneticCoordinates}
  q = nothing
  try
    q = getfield(descsurf, quantity)
  catch e
    quantity_string = string(quantity)
    println("vmecsurf has no quantity, "*quantity_string)
    return nothing
  end

  function set_surf_field()
    (f, dfds) = get_function_and_deriv(quantity) 
    (q, dqds) = make_desc_surface_interpolation(f, dfds, descsurf, minθres, minζres)
    quantity_ds = Symbol("d"*string(quantity)*"ds")
    
    setfield!(descsurf, quantity, q)
    setfield!(descsurf, quantity_ds, dqds)
  end
    
  if q == nothing
    set_surf_field()
    q = getfield(descsurf, quantity)
  end

  #check if size is too small
  if size(q.itp)[1] < minθres || size(q.itp)[2] < minζres
    set_surf_field()
    q = getfield(descsurf, quantity)
  end
  
  if deriv == :none
    return q(x.θ, x.ζ)
  elseif deriv == :ds
    quantity_ds = Symbol("d"*string(quantity)*"ds")
    dqds = getfield(descsurf, quantity_ds)
    return dqds(x.θ, x.ζ)
  elseif deriv == :dθ
    return Interpolations.gradient(q, x.θ, x.ζ)[1]
  elseif deriv == :dζ
    return Interpolations.gradient(q, x.θ, x.ζ)[2]
  elseif deriv == :dsdθ
    quantity_ds = Symbol("d"*string(quantity)*"ds")
    dqds = getfield(descsurf, quantity_ds)
    return Interpolations.gradient(dqds, x.θ, x.ζ)[1]
  elseif deriv == :dsdζ
    quantity_ds = Symbol("d"*string(quantity)*"ds")
    dqds = getfield(descsurf, quantity_ds)
    return Interpolations.gradient(dqds, x.θ, x.ζ)[2]
  end
end

function desc_λ(c::Union{ZernikeCoordinates, FluxCoordinates}, descsurf::DescSurface)
  return inverseTransform(c, descsurf.λmn)
end

function desc_dλds(c::Union{ZernikeCoordinates, FluxCoordinates}, descsurf::DescSurface)
  return inverseTransform(c, descsurf.λmn, deriv=:ds)
end

function desc_R(c::Union{ZernikeCoordinates, FluxCoordinates}, descsurf::DescSurface)
  return inverseTransform(c, descsurf.rmn)
end

function desc_dRds(c::Union{ZernikeCoordinates, FluxCoordinates}, descsurf::DescSurface)
  return inverseTransform(c, descsurf.rmn, deriv=:ds)
end

function desc_Z(c::Union{ZernikeCoordinates, FluxCoordinates}, descsurf::DescSurface)
  return inverseTransform(c, descsurf.zmn)
end

function desc_dZds(c::Union{ZernikeCoordinates, FluxCoordinates}, descsurf::DescSurface)
  return inverseTransform(c, descsurf.zmn, deriv=:ds)
end



#calculate 2nd derivative functions 
function desc_dρdρ(zc::ZernikeCoordinates, value::Symbol, descsurf::DescSurface)
  ρ = descsurf.ρ
  qmn = getfield(descsurf, value)
  dqds = inverseTransform(zc, qmn, deriv=:ds)
  d2qds2 = inverseTransform(zc, qmn, deriv=:dsds)
  d2qdρ2 = d2qds2*4*ρ^2 + dqds * 2
  return d2qdρ2
end

function desc_eρ(zc::ZernikeCoordinates, descsurf::DescSurface)
  dRds = inverseTransform(zc, descsurf.rmn, deriv=:ds)
  dRdρ = dRds * 2 * zc.ρ
  dZds = inverseTransform(zc, descsurf.zmn, deriv=:ds)
  dZdρ = dZds * 2 * zc.ρ
  return SVector(dRdρ, 0, dZdρ)
end

function desc_es(fc::FluxCoordinates, descsurf::DescSurface)
  dRds = inverseTransform(fc, descsurf.rmn, deriv=:ds)
  dZds = inverseTransform(fc, descsurf.zmn, deriv=:ds)
  return SVector(dRds, 0, dZds)
end

#note for θ and ζ derivs the radial coordinate is not used
function desc_eθ(c::Union{ZernikeCoordinates, FluxCoordinates}, 
                 descsurf::DescSurface)
  dRdθ = inverseTransform(c, descsurf.rmn, deriv=:dθ)
  dZdθ = inverseTransform(c, descsurf.zmn, deriv=:dθ)
  return SVector(dRdθ, 0, dZdθ)
end

function desc_eζ(c::Union{ZernikeCoordinates, FluxCoordinates}, 
                 descsurf::DescSurface)
  dRdζ = inverseTransform(c, descsurf.rmn, deriv=:dζ)
  dZdζ = inverseTransform(c, descsurf.zmn, deriv=:dζ)
  R = inverseTransform(c, descsurf.rmn)
  return SVector(dRdζ, R, dZdζ)
end

#ρ derivatives of basis vectors
function desc_deθdρ(zc::ZernikeCoordinates, descsurf::DescSurface)

  dR2dsdθ = inverseTransform(zc, descsurf.rmn, deriv=:dsdθ)
  dR2dρdθ = dR2dsdθ * 2 * descsurf.ρ
  dZ2dsdθ = inverseTransform(zc, descsurf.zmn, deriv=:dsdθ)
  dZ2dρdθ = dZ2dsdθ * 2 * descsurf.ρ
  return SVector(dR2dρdθ, 0, dZ2dρdθ)
end

function desc_deθds(fc::FluxCoordinates, descsurf::DescSurface)

  dR2dsdθ = inverseTransform(fc, descsurf.rmn, deriv=:dsdθ)
  dZ2dsdθ = inverseTransform(fc, descsurf.zmn, deriv=:dsdθ)
  return SVector(dR2dsdθ, 0, dZ2dsdθ)
end

function desc_deζdρ(zc::ZernikeCoordinates, descsurf::DescSurface)

  dR2dsdζ = inverseTransform(zc, descsurf.rmn, deriv=:dsdζ)
  dR2dρdζ = dR2dsdζ * 2 * descsurf.ρ
  dRds = inverseTransform(zc, descsurf.rmn, deriv=:ds)
  dRdρ = dRds * 2 * descsurf.ρ
  dZ2dsdζ = inverseTransform(zc, descsurf.zmn, deriv=:dsdζ)
  dZ2dρdζ = dZ2dsdζ * 2 * descsurf.ρ
  return SVector(dR2dρdζ, dRdρ, dZ2dρdζ)
end

function desc_deζds(fc::FluxCoordinates, descsurf::DescSurface)

  dR2dsdζ = inverseTransform(fc, descsurf.rmn, deriv=:dsdζ)
  dRds = inverseTransform(fc, descsurf.rmn, deriv=:ds)
  dZ2dsdζ = inverseTransform(fc, descsurf.zmn, deriv=:dsdζ)
  return SVector(dR2dsdζ, dRds, dZ2dsdζ)
end

function desc_deρdρ(zc::ZernikeCoordinates, descsurf::DescSurface)
  d2Rdρ2 = desc_dρdρ(zc, :rmn, descsurf)
  d2Zdρ2 = desc_dρdρ(zc, :zmn, descsurf)
  return SVector(d2Rdρ2, 0.0, d2Zdρ2)
end

function desc_desds(fc::FluxCoordinates, descsurf::DescSurface)
  d2Rds2 = inverseTransform(fc, descsurf.rmn, deriv=:dsds)
  d2Zds2 = inverseTransform(fc, descsurf.zmn, deriv=:dsds)
  return SVector(d2Rds2, 0.0, d2Zds2)
end

# Metric tensor elements
function desc_g_ρρ(c::Union{ZernikeCoordinates, FluxCoordinates},
                  descsurf::DescSurface)
  return dot(desc_eρ(c,descsurf), desc_eρ(c,descsurf))
end

function desc_g_θθ(c::Union{ZernikeCoordinates, FluxCoordinates},
                  descsurf::DescSurface)
  return dot(desc_eθ(c,descsurf), desc_eθ(c,descsurf))
end

function desc_g_ζζ(c::Union{ZernikeCoordinates, FluxCoordinates},
                  descsurf::DescSurface)
  return dot(desc_eζ(c,descsurf), desc_eζ(c,descsurf))
end

function desc_g_ρθ(c::Union{ZernikeCoordinates, FluxCoordinates},
                  descsurf::DescSurface)
  return dot(desc_eρ(c,descsurf), desc_eθ(c,descsurf))
end

function desc_g_ρζ(c::Union{ZernikeCoordinates, FluxCoordinates},
                  descsurf::DescSurface)
  return dot(desc_eρ(c,descsurf), desc_eζ(c,descsurf))
end

function desc_g_θζ(c::Union{ZernikeCoordinates, FluxCoordinates},
                  descsurf::DescSurface)
  return dot(desc_eθ(c,descsurf), desc_eζ(c,descsurf))
end

#two derivatives we need for calcs
function desc_dg_θθdρ(c::ZernikeCoordinates,
                      descsurf::DescSurface)
  return 2 * dot(desc_eθ(c, descsurf), desc_deθdρ(c, descsurf))
end

function desc_dg_θθds(c::FluxCoordinates,
                      descsurf::DescSurface)
  return 2 * dot(desc_eθ(c, descsurf), desc_deθds(c, descsurf))
end

function desc_dg_θζdρ(c::Union{ZernikeCoordinates, FluxCoordinates},
                      descsurf::DescSurface)
  return (dot(desc_deθdρ(c, descsurf), desc_eζ(c, descsurf)) + 
          dot(desc_eθ(c, descsurf), desc_deζdρ(c, descsurf)))
end

function desc_jacobian(c::ZernikeCoordinates, descsurf::DescSurface)
  eρ = desc_eρ(c, descsurf)
  eθ = desc_eθ(c, descsurf)
  eζ = desc_eζ(c, descsurf)
  return dot(eρ, cross(eθ, eζ))
end

function desc_jacobian(c::FluxCoordinates, descsurf::DescSurface)
  es = desc_es(c, descsurf)
  eθ = desc_eθ(c, descsurf)
  eζ = desc_eζ(c, descsurf)
  return dot(es, cross(eθ, eζ))
end

function desc_djacobiandρ(zc::ZernikeCoordinates, descsurf::DescSurface)
  deρdρ = desc_deρdρ(zc, descsurf)
  deθdρ = desc_deθdρ(zc, descsurf)
  deζdρ = desc_deζdρ(zc, descsurf)
  eρ = desc_eρ(zc, descsurf)
  eθ = desc_eθ(zc, descsurf)
  eζ = desc_eζ(zc, descsurf)
  v = dot(deρdρ, cross(eθ, eζ))
  v += dot(eρ, cross(deθdρ, eζ))
  v += dot(eρ, cross(eθ, deζdρ))
  return v
end

function desc_djacobiands(fc::FluxCoordinates, descsurf::DescSurface)
  #this is wrong
  desds = desc_desds(fc, descsurf)
  deθds = desc_deθds(fc, descsurf)
  deζds = desc_deζds(fc, descsurf)
  es = desc_es(fc, descsurf)
  eθ = desc_eθ(fc, descsurf)
  eζ = desc_eζ(fc, descsurf)
  v = dot(desds, cross(eθ, eζ))
  v += dot(es, cross(deθds, eζ))
  v += dot(es, cross(eθ, deζds))
  return v
#  zc = ZernikeCoordinates(descsurf.ρ, fc.θ, fc.ζ)
#  return desc_djacobiandρ(zc, descsurf) / (2 * zc.ρ)
end

function desc_iota( descsurf::DescSurface;)
  #current term needs to be benchmarked
  μ0 = 1.25663706E-6
  ρ = descsurf.ρ
  current_term = μ0/(2π) * descsurf.current[1] / (descsurf.phi[2]*(2*ρ))
  function numerator(v::SVector{T}) where {T}

    θ = v[1]
    ζ = v[2]
    zc = ZernikeCoordinates(descsurf.ρ,θ,ζ)
    λ_ζ = surface_get(zc, descsurf, :λ, deriv=:dζ)
    λ_θ = surface_get(zc, descsurf, :λ, deriv=:dθ)
    g_θθ = desc_g_θθ(zc, descsurf)
    g_θζ = desc_g_θζ(zc, descsurf)
    sqrtg = surface_get(zc, descsurf, :g)*(2*ρ)
    return ((λ_ζ * g_θθ) - (1 + λ_θ)*g_θζ)/sqrtg
  end

  function denominator(v::SVector{T}) where {T}
    θ = v[1]
    ζ = v[2]
    zc = ZernikeCoordinates(ρ,θ,ζ)
    sqrtg = surface_get(zc, descsurf, :g)*(2*ρ)
    return desc_g_θθ(zc, descsurf)/sqrtg
  end

  ζmax = 2π/descsurf.nfp

  num = hcubature(numerator,[0.0,0.0],[2π,ζmax],rtol=1.0E-6)[1]*descsurf.nfp/4/π^2 
  den = hcubature(denominator,[0.0,0.0],[2π,ζmax],rtol=1.0E-6)[1]*descsurf.nfp/4/π^2
  return (current_term + num)/den
end 

function desc_diotadρ(descsurf::DescSurface, iota::T) where {T}
  μ0 = 1.25663706E-6
  ρ = descsurf.ρ
  current_term = μ0/2π * descsurf.current[1] / (descsurf.phi[2] * 2*ρ)
  ψ_rr = descsurf.phi[3]/π

  #note the 2ρ factors cancel
  current_term_r = (μ0/2π) * descsurf.current[2]/(descsurf.phi[2])
  current_term_r -= current_term * ψ_rr / (descsurf.phi[2]*2*ρ)

  function numerator_r(v::SVector{T}) where {T}

    θ = v[1]
    ζ = v[2]
    ρ = descsurf.ρ
    zc = ZernikeCoordinates(ρ,θ,ζ)
    λ_ζ = surface_get(zc, descsurf, :λ, deriv=:dζ)
    λ_θ = surface_get(zc, descsurf, :λ, deriv=:dθ)
    g_θθ = desc_g_θθ(zc, descsurf)
    g_θζ = desc_g_θζ(zc, descsurf)
    sqrtg = surface_get(zc, descsurf, :g) #function of s
    gρ = sqrtg * 2 * ρ
    numerator = ((λ_ζ * g_θθ) - (1 + λ_θ)*g_θζ)/gρ
    λ_ρζ = surface_get(zc, descsurf, :λ, deriv=:dsdζ)*2*ρ
    λ_ρθ = surface_get(zc, descsurf, :λ, deriv=:dsdθ)*2*ρ
    dg_θθdρ = desc_dg_θθdρ(zc, descsurf)
    dg_θζdρ = desc_dg_θζdρ(zc, descsurf)
    dgds = surface_get(zc, descsurf, :g, deriv=:ds)
    dgdρ = 2*(sqrtg + 2*ρ^2* dgds)
    t1 = (λ_ρζ * g_θθ) + (λ_ζ * dg_θθdρ) - (λ_ρθ * g_θζ)
    t1 -= (1 + λ_θ) * dg_θζdρ
    t1 /= gρ
    t2 = numerator * dgdρ/gρ
    return t1 - t2
  end

  function denominator(v::SVector{T}) where {T}
    θ = v[1]
    ζ = v[2]
    zc = ZernikeCoordinates(descsurf.ρ,θ,ζ)
    sqrtg = surface_get(zc, descsurf, :g)
    gρ = sqrtg * 2 * ρ
    return desc_g_θθ(zc, descsurf)/gρ
  end
  
  function denominator_r(v::SVector{T}) where {T}
    θ = v[1]
    ζ = v[2]
    zc = ZernikeCoordinates(descsurf.ρ,θ,ζ)
    g_θθ = desc_g_θθ(zc, descsurf)
    dg_θθdρ = desc_dg_θθdρ(zc, descsurf)
    sqrtg = surface_get(zc, descsurf, :g)
    gρ = sqrtg * 2 * ρ
    dgds = surface_get(zc, descsurf, :g, deriv=:ds)
    dgdρ = 2*(sqrtg + 2*ρ^2* dgds)
#    dgdρ = desc_djacobiandρ(zc, descsurf)
    return (dg_θθdρ - (g_θθ)*dgdρ/gρ)/gρ
  end
  

  ζmax = 2π/descsurf.nfp
  num_r = hcubature(numerator_r,[0.0,0.0],[2π,ζmax],rtol=1.0E-6)[1]*descsurf.nfp/4/π^2 
  den = hcubature(denominator,[0.0,0.0],[2π,ζmax],rtol=1.0E-6)[1]*descsurf.nfp/4/π^2 
  den_r = hcubature(denominator_r,[0.0,0.0],[2π,ζmax],rtol=1.0E-6)[1]*descsurf.nfp/4/π^2 
  diotadρ = (current_term_r + num_r - (iota * den_r))/den
  return diotadρ

end    

#Benchmark for testing trapezoidal sum
function desc_iota_rough(descsurf::DescSurface; θres = 128, ζres = 128)
  dθ = 2π/(θres + 1)
  dζ = 2π/descsurf.nfp/(ζres + 1)
  nmat = zeros(Float64, θres, ζres)
  dmat = zeros(Float64, θres, ζres)
  for (θi, θ) in enumerate(range(0, 2π-dθ, θres))
    for (ζi, ζ) in enumerate(range(0, 2π/descsurf.nfp - dζ, ζres))
      zc = ZernikeCoordinates(descsurf.ρ,θ,ζ)
      g = desc_jacobian(zc, descsurf)
      λ_ζ = inverseTransform(zc,descsurf.λmn,deriv=:dζ)
      λ_θ = inverseTransform(zc,descsurf.λmn,deriv=:dθ)
      g_θθ = desc_g_θθ(zc, descsurf)
      g_θζ = desc_g_θζ(zc, descsurf)
      nmat[θi,ζi] = ((λ_ζ * g_θθ) - (1 + λ_θ)*g_θζ)/g
      dmat[θi,ζi] = g_θθ/g
    end
  end
  return sum(nmat) / sum(dmat)
end

#B^θ - using vmec convention
function desc_bsupu(c::Union{FluxCoordinates, ZernikeCoordinates}, descsurf::DescSurface)
  g = surface_get(c, descsurf, :g)
  dλdζ = surface_get(c, descsurf, :λ, deriv = :dζ)
  return (descsurf.phi[2] / g) * (descsurf.iota[1] - dλdζ)
end

function desc_dbsupuds(c::Union{FluxCoordinates, ZernikeCoordinates}, descsurf::DescSurface)
  dψds = descsurf.phi[2]
  g = surface_get(c, descsurf, :g)
  dgds = surface_get(c, descsurf, :g, deriv=:ds)
  ι = descsurf.iota[1]
  dιds = descsurf.iota[2]
  dλdζ = surface_get(c, descsurf, :λ, deriv = :dζ)
  dλ2dsdζ = surface_get(c, descsurf, :λ, deriv = :dsdζ)
  B0 = dψds/g

  return B0*((dιds - dλ2dsdζ)-(dgds*(ι - dλdζ)/g))
end

#B^ζ - using vmec convention
function desc_bsupv(c::Union{FluxCoordinates, ZernikeCoordinates}, descsurf::DescSurface)
  g = surface_get(c, descsurf, :g)
  dλdθ = surface_get(c, descsurf, :λ, deriv = :dθ)
  return (descsurf.phi[2] / g) * (1+ dλdθ)
end

function desc_dbsupvds(c::Union{FluxCoordinates, ZernikeCoordinates}, descsurf::DescSurface)
  dψds = descsurf.phi[2]
  g = surface_get(c, descsurf, :g)
  dgds = surface_get(c, descsurf, :dgds)
  dλdθ = surface_get(c, descsurf, :λ, deriv = :dθ)
  dλ2dsdθ = surface_get(c, descsurf, :λ, deriv = :dsdθ)
  B0 = dψds/g
  return B0*(dλ2dsdθ - (dgds*(1 + dλdθ)/g))
end

function desc_bfield(c::Union{FluxCoordinates, ZernikeCoordinates}, descsurf::DescSurface)
  eθ = desc_eθ(c, descsurf)
  eζ = desc_eζ(c, descsurf)
  Bsupθ = surface_get(c, descsurf, :bsupu)
  Bsupζ = surface_get(c, descsurf, :bsupv)
  return Bsupθ * eθ .+ Bsupζ * eζ
end

function desc_dBds(fc::FluxCoordinates, descsurf::DescSurface)
  t1 = surface_get(fc, descsurf, :bsupu, deriv=:ds) * desc_eθ(fc, descsurf)
  t2 = surface_get(fc, descsurf, :bsupu) * desc_deθds(fc, descsurf)
  t3 = surface_get(fc, descsurf, :bsupv, deriv=:ds) * desc_eζ(fc, descsurf)
  t4 = surface_get(fc, descsurf, :bsupv) * desc_deζds(fc, descsurf)
  return t1 .+ t2 .+ t3 .+ t4
end

function desc_dBdρ(zc::ZernikeCoordinates, descsurf::DescSurface)
  fc = FluxCoordinates(descsurf.s, zc.θ, zc.ζ)
  return desc_dBds(fc, descsurf) * 2 * zc.ρ
end

function desc_bsubu(c::Union{FluxCoordinates, ZernikeCoordinates}, descsurf::DescSurface)
  B = desc_bfield(c, descsurf)
  eθ = desc_eθ(c, descsurf)
  return dot(B, eθ)
end

function desc_bsubv(c::Union{FluxCoordinates, ZernikeCoordinates}, descsurf::DescSurface)
  B = desc_bfield(c, descsurf)
  eζ = desc_eζ(c, descsurf)
  return dot(B, eζ)
end

function desc_bsubs(c::FluxCoordinates, descsurf::DescSurface)
  B = desc_bfield(c, descsurf)
  es = desc_es(c, descsurf)
  return dot(B, es)
end

#probably will never use this function
function desc_bsubρ(c::ZernikeCoordinates, descsurf::DescSurface)
  B = desc_bfield(c, descsurf)
  eρ = desc_es(c, descsurf)
  return dot(B, eρ)
end

function desc_dbsubuds(c::FluxCoordinates, descsurf::DescSurface)
  t1 = dot(desc_dBds(c, descsurf), desc_eθ(c, descsurf)) 
  t2 = dot(desc_bfield(c, descsurf), desc_deθds(c, descsurf)) 
  return t1 + t2
end

function desc_dbsubvds(c::FluxCoordinates, descsurf::DescSurface)
  t1 = dot(desc_dBds(c, descsurf), desc_eζ(c, descsurf)) 
  t2 = dot(desc_bfield(c, descsurf), desc_deζds(c, descsurf)) 
  return t1 + t2
end

function desc_dbsubsds(c::FluxCoordinates, descsurf::DescSurface)
  t1 = dot(desc_dBds(c, descsurf), desc_es(c, descsurf)) 
  t2 = dot(desc_bfield(c, descsurf), desc_desds(c, descsurf)) 
  return t1 + t2
end

#function used for benchmarking against desc
function desc_dbsubρdρ(zc::ZernikeCoordinates, descsurf::DescSurface)
  fc = FluxCoordinates(descsurf.s, zc.θ, zc.ζ)
  dbsubsds = surface_get(fc, descsurf, :bsubs, deriv=:ds)
  bsubs = surface_get(fc, descsurf, :bsubs)
  ρ = descsurf.ρ
  return (dbsubsds * 4 * ρ^2 + 2*bsubs)
end

function desc_bmag(c::Union{ZernikeCoordinates, FluxCoordinates}, descsurf::DescSurface)
  return norm(desc_bfield(c, descsurf))
end

function desc_dbmagds(fc::FluxCoordinates, descsurf::DescSurface)
  dbsupuds = surface_get(fc, descsurf, :bsupu, deriv=:ds)
  bsupu = surface_get(fc, descsurf, :bsupu)
  dbsupvds = surface_get(fc, descsurf, :bsupv, deriv=:ds)
  bsupv = surface_get(fc, descsurf, :bsupv)
  dbsubuds = surface_get(fc, descsurf, :bsubu, deriv=:ds)
  bsubu = surface_get(fc, descsurf, :bsubu)
  dbsubvds = surface_get(fc, descsurf, :bsubv, deriv=:ds)
  bsubv = surface_get(fc, descsurf, :bsubv)
  v = dbsupuds * bsubu + bsupu * dbsubuds
  v += dbsupvds * bsubv + bsupv * dbsubvds
  v /= 2*desc_bmag(fc, descsurf)
  return v
end

function desc_volume(descsurf::DescSurface)
  function integrand(v::SVector{T}) where {T}
    θ = v[1]
    ζ = v[2]
    fc = FluxCoordinates(descsurf.s,θ,ζ)
    eθ = desc_eθ(fc, descsurf)
    eζ = desc_eζ(fc, descsurf)
    t1 = cross(eθ, eζ)[3]
    z = surface_get(fc, descsurf, :z)
    return z*t1
  end
  ζmax = 2π/descsurf.nfp
  return hcubature(integrand,[0.0,0.0],[2π,ζmax],rtol=1.0E-6)[1] * descsurf.nfp
end


function desc_dVds(descsurf::DescSurface)
  function integrand(v::SVector{T}) where {T}
    θ = v[1]
    ζ = v[2]
    fc = FluxCoordinates(descsurf.s,θ,ζ)
    return abs(surface_get(fc, descsurf, :g))
  end
  ζmax = 2π/descsurf.nfp
  return hcubature(integrand,[0.0,0.0],[2π,ζmax],rtol=1.0E-6)[1]*descsurf.nfp
end

#Not sure if these are the correct quantities.  I think the vmec ones
#may be the covariant quantities.
function desc_jsups(c::Union{ZernikeCoordinates, FluxCoordinates}, descsurf::DescSurface)
  dbsubζdθ = surface_get(c, descsurf, :bsubv, deriv=:dθ)
  dbsubθdζ = surface_get(c, descsurf, :bsubu, deriv=:dζ)
  g = surface_get(c, descsurf, :g)
  μ0 = 1.25663706E-6
  return (dbsubζdθ - dbsubθdζ) / (μ0 * g)
end

function desc_jsupu(c::Union{ZernikeCoordinates, FluxCoordinates}, descsurf::DescSurface)
  dbsubζds = surface_get(c, descsurf, :bsubv, deriv=:ds)
  dbsubsdζ = surface_get(c, descsurf, :bsubs, deriv=:dζ)
  g = surface_get(c, descsurf, :g)
  μ0 = 1.25663706E-6
  return (dbsubsdζ - dbsubζds) / (μ0 * g)
end

function desc_jsupv(c::Union{ZernikeCoordinates, FluxCoordinates}, descsurf::DescSurface)
  dbsubuds = surface_get(c, descsurf, :bsubu, deriv=:ds)
  dbsubsdθ = surface_get(c, descsurf, :bsubs, deriv=:dθ)
  g = surface_get(c, descsurf, :g)
  μ0 = 1.25663706E-6
  return (dbsubuds - dbsubsdθ) / (μ0 * g)
end

function unknown_deriv(c::Union{ZernikeCoordinates, FluxCoordinates}, descsurf::DescSurface)
  return 0.0
end

#calculation for adding all the missing quantities for the desc surface
function desc_finish_surface(descsurf::DescSurface)
  #calculate sqrt(g)
  (g, dgds) = make_desc_surface_interpolation(desc_jacobian,
                                              desc_djacobiands,
                                              descsurf, 128, 128)
  setfield!(descsurf,:g,g)
  setfield!(descsurf,:dgds,dgds)


  #calculate iota, time consuming...
  T = typeof(descsurf.ρ)
  iota = desc_iota(descsurf)
  diotads = desc_diotadρ(descsurf, iota)/(2*descsurf.ρ)
  descsurf.iota = SVector{2,T}(iota, diotads)
  v = desc_volume(descsurf)
  dVds = desc_dVds(descsurf)
  descsurf.vol = SVector{2,T}(v, dVds)
  #all other quantities should be calculated and added
  #as soon as the first surface_get is called
  return descsurf
end



#constructor for DescSurface
function DescSurface(s::T, desc::DescData; 
                     simple::Bool=false) where {T <: AbstractFloat}
  ρ = sqrt(s)
  if ρ == 0
    ρ = 0.0001
  end
  rmn = descFourierDataArray(ρ,desc,:R_lmn, :R_modes)
  zmn = descFourierDataArray(ρ,desc,:Z_lmn, :Z_modes)
  λmn = descFourierDataArray(ρ,desc,:λ_lmn, :λ_modes)

  #sort to have a consistent m and n
  sortfunc(sfa::SurfaceFourierData) = (sfa.m, sfa.n)
  sort!(rmn, by=sortfunc)
  sort!(zmn, by=sortfunc)
  sort!(λmn, by=sortfunc)


  phi1 = desc_ψ(ρ, desc)
  dphids = desc_dψdρ(ρ, desc)/(2*ρ)
  phi = SVector{3,T}(phi1, dphids, desc.ψ)
  pres = SVector{2,T}(desc_pressure(ρ, desc), desc_dpressuredρ(ρ, desc)/(2*ρ))
  current = SVector{2,T}(desc_current(ρ, desc), desc_dcurrentdρ(ρ, desc)/(2*ρ))

  #construct the m and n vectors
  xm = [sfa.m for sfa in rmn]
  xn = [sfa.n for sfa in rmn]

  mpol = Int64(maximum(xm))
  ntor = Int64(maximum(xn))
  mnmax = length(xm)
  dummy = SVector{2,T}(0.0, 0.0)

  descsurf = DescSurface(rmn, zmn, λmn,
                     nothing, nothing, nothing, nothing, nothing, nothing,
                     nothing, nothing, nothing, nothing, nothing, nothing,
                     nothing, nothing, nothing,
                     nothing, nothing, nothing, nothing, nothing, nothing,
                     nothing, nothing, nothing, nothing, nothing, nothing,
                     nothing, nothing, nothing,
                     phi, pres, current, dummy, dummy, s, ρ,
                     desc.nfp, mpol, ntor, mnmax, mnmax,
                     xm,xn,xm,xn) #nyquist the same as non-nyquist


  
  if simple == true
    return descsurf
  else
    return desc_finish_surface(descsurf)
  end
end

#utility functions for zernike polynomials
function zernike_prefactor(l::Integer,m::Integer)
  m = abs(m)
  s_max = (l - m) ÷ 2
  prefactor_array = zeros(Float64, l+1)#note zero index
  for s in 0:s_max
    numerator = (-1)^s * factorial(l-s)
    denom1 = factorial(((l + m)÷2) - s)
    denom2 = factorial(((l - m)÷2) - s)
    denom = factorial(s)*denom1*denom2
    prefactor_array[l-(2*s)+1] = numerator / denom
  end
  return prefactor_array
end

function zernike_r(ρ::Float64, l::Integer, m::Integer;
                   deriv=0)
  z_prefactor = zernike_prefactor(l, m)
  v = 0
  for (i, c) in enumerate(z_prefactor)
    if c == 0
      continue
    end
    exponent = i-1
    for d in 1:deriv
      c *= exponent
      exponent -= 1
    end
    if exponent < 0 
      continue
    end
    v += c * ρ^exponent
  end
  return v
end

#generalized fourier used in Zernike polynomials
#If arg is ζ then n should be n * NFP
function zernike_fourier(arg::Float64, n::Integer; deriv=0)
  if n < 0
    feval = 0
    n = abs(n)
  else
    feval = 1
  end
  mult = 1
  for d in 1:deriv
    mult *= n
    feval += 1
  end
  feval %= 4
  if feval == 0
    return mult*sin(n*arg)
  elseif feval == 1
    return mult*cos(n*arg)
  elseif feval == 2
    return -mult*sin(n*arg)
  elseif feval == 3
    return -mult*cos(n*arg)
  end
end

#derivoptions
#dρdρ
#dρdθ
#dρdζ
#dθdθ
#dθdζ
#dζdζ
#dρ
#dθ
#dζ

function zernike_eval(zc::ZernikeCoordinates, 
                      coeffs::Vector{Float64},
                      modes::Matrix{IT}, 
                      nfp::Integer; 
                      deriv::String="") where {IT}
  v = 0
  larr = modes[1,:]
  marr = modes[2,:]
  narr = modes[3,:]
  ρ = zc.ρ
  θ = zc.θ
  ζ = zc.ζ
  #probably should check to make sure all the lengths are correct
  for i in 1:length(coeffs)
    l = larr[i]
    m = marr[i]
    n = narr[i]
    if deriv == ""
      v += coeffs[i]*zernike_r(ρ,l,m)*zernike_fourier(θ, m)*zernike_fourier(ζ, nfp*n)
    elseif deriv == "dρ"
      v += coeffs[i]*zernike_r(ρ,l,m,deriv=1)*zernike_fourier(θ, m)*zernike_fourier(ζ, nfp*n)
    elseif deriv == "dρdρ"
      v += coeffs[i]*zernike_r(ρ,l,m,deriv=2)*zernike_fourier(θ, m)*zernike_fourier(ζ, nfp*n)
    elseif deriv == "dρdθ" || deriv == "dθdρ"
      v += coeffs[i]*zernike_r(ρ,l,m,deriv=1)*zernike_fourier(θ, m,deriv=1)*zernike_fourier(ζ, nfp*n)
    elseif deriv == "dρdζ" || deriv == "dζdρ"
      v += coeffs[i]*zernike_r(ρ,l,m,deriv=1)*zernike_fourier(θ, m)*zernike_fourier(ζ, nfp*n, deriv=1)
    elseif deriv == "dθ"
      v += coeffs[i]*zernike_r(ρ,l,m)*zernike_fourier(θ, m,deriv=1)*zernike_fourier(ζ, nfp*n)
    elseif deriv == "dθdθ"
      v += coeffs[i]*zernike_r(ρ,l,m)*zernike_fourier(θ, m,deriv=2)*zernike_fourier(ζ, nfp*n)
    elseif deriv == "dθdζ" || deriv == "dζdθ"
      v += coeffs[i]*zernike_r(ρ,l,m)*zernike_fourier(θ, m,deriv=1)*zernike_fourier(ζ, nfp*n, deriv=1)
    elseif deriv == "dζ"
      v += coeffs[i]*zernike_r(ρ,l,m)*zernike_fourier(θ, m)*zernike_fourier(ζ, nfp*n, deriv = 1)
    elseif deriv == "dζdζ"
      v += coeffs[i]*zernike_r(ρ,l,m)*zernike_fourier(θ, m)*zernike_fourier(ζ, nfp*n, deriv = 2)
    end
  end
  return v
end

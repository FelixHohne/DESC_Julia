#get values of R,Z for a given ζ, ρ value
function desc_RZ(ρ::Float64, ζ::Float64,desc::DescData; θres = 100)
  θs = range(0,2π,θres)
  Rs = [desc_R(ZernikeCoordinates(ρ, θ, ζ), desc) for θ in θs]
  Zs = [desc_Z(ZernikeCoordinates(ρ, θ, ζ), desc) for θ in θs]
  return Rs, Zs
end
  

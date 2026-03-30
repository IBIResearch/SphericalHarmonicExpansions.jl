import Base.setindex!, Base.getindex, Base.isapprox, Base.==, Base.!=,
        Base.+, Base.-, Base.*, Base./, Base.write

mutable struct SphericalHarmonicCoefficients{T<:Real, vecT<:AbstractVector{T}}
  c::vecT
  L::Int64
  R::Float64
  solid::Bool

  function SphericalHarmonicCoefficients(c::AbstractVector{T}, L::Real, R::Real,solid::Bool) where T<:Real
    if L<0 || !isinteger(L)
      throw(DomainError(L,"L ∈ ℕ₀ required, but is $L"))
    end
    if (!isinteger(sqrt(length(c))) || length(c)==0)
      throw(DomainError(sqrt(length(c)),"Input vector needs to be of size (L+1)², where L ∈ ℕ₀."))
    end 
    return new{T, typeof(c)}(c,L,R,solid)
  end

  SphericalHarmonicCoefficients(c::AbstractVector{<:Real}) = SphericalHarmonicCoefficients(c, sqrt(length(c))-1, 1.0, false)
  
  SphericalHarmonicCoefficients(L::Integer) = SphericalHarmonicCoefficients(zeros((L+1)^2), L, 1.0, false)

  SphericalHarmonicCoefficients(L::Int,R::Real,solid::Bool) = SphericalHarmonicCoefficients(zeros((L+1)^2), L, R, solid)

  SphericalHarmonicCoefficients(c::AbstractVector{<:Real},R::Float64,solid::Bool) = SphericalHarmonicCoefficients(c,sqrt(length(c))-1,R,solid)

  function SphericalHarmonicCoefficients(c::AbstractArray{<:AbstractVector{T}},R::AbstractArray{<:Real},solid::AbstractArray{<:Bool}) where T<:Real
    if size(c) != size(R) || size(c) != size(solid)
      throw(DomainError(size(c),"Arrays do not have the same size."))
    end

    L = [sqrt(length(co))-1 for co in c]
    return reshape([SphericalHarmonicCoefficients(c[n],L[n],R[n],solid[n]) for n in eachindex(c)],size(c)...)
  end

  SphericalHarmonicCoefficients(c::AbstractArray{<:AbstractVector{<:Real}},R::T,solid::Bool)  where T<:Real = SphericalHarmonicCoefficients(c,fill(R, size(c)),fill(solid, size(c)) )

  # write and read coefficients to/from an HDF5-file
  function SphericalHarmonicCoefficients(path::String)
    coeffs, R, solid = h5open(path,"r") do file
      coeffsArray = read(file, "/coeffs")
      R = read(file, "/normalization")
      solid = (read(file, "/solid") .== 1)

      if size(coeffsArray)[1:end-1] == ()
        coeffs = Array{Vector{Float64}}(undef,1)
      else
        coeffs = Array{Vector{Float64}}(undef,size(coeffsArray)[1:end-1])
      end
      coeffsArray = reshape(coeffsArray,(Int(length(coeffsArray)/size(coeffsArray)[end]),size(coeffsArray)[end]))
      for n in axes(coeffsArray,1)
        coeffs[n] = coeffsArray[n,:]
      end

      return coeffs, R, solid
    end

    return SphericalHarmonicCoefficients(coeffs, R, solid)
  end
end

function write(path::String, coeffs::AbstractArray{<:SphericalHarmonicCoefficients})

  if size(coeffs) != (1,)
      coeffsArray = coeffs[1].c'
      for (n,co) in enumerate(coeffs[2:end])
          coeffsArray = vcat(coeffsArray,co.c')
      end
      coeffsArray = reshape(coeffsArray,(size(coeffs)...,(coeffs[1].L+1)^2))
  else
      coeffsArray = coeffs[1].c
  end

  R = Array{Float64}(undef,size(coeffs))
  solid = Array{Int}(undef,size(coeffs))
  for (n,co) in enumerate(coeffs)
      R[n] = co.R
      solid[n] = co.solid ? 1 : 0
  end

  h5open(path,"w") do file
      write(file, "/coeffs", coeffsArray)
      write(file, "/normalization", R)
      write(file, "/solid", solid)
  end
end

isapprox(shc1::SphericalHarmonicCoefficients, shc2::SphericalHarmonicCoefficients; kargs...) = isapprox(shc1.c,shc2.c;kargs...) && isapprox(shc1.R,shc2.R;kargs...) && isapprox(shc1.solid,shc2.solid;kargs...)
getindex(shc::SphericalHarmonicCoefficients,I) = getindex(shc.c,I)
getindex(shc::SphericalHarmonicCoefficients,l,m) = getindex(shc.c, l*(l+1)+m+1)
setindex!(shc::SphericalHarmonicCoefficients,value,I) = setindex!(shc.c,value,I)
setindex!(shc::SphericalHarmonicCoefficients,value,l,m) = setindex!(shc.c,value,l*(l+1)+m+1)
==(shca::SphericalHarmonicCoefficients, shcb::SphericalHarmonicCoefficients) = shca.c == shcb.c
!=(shca::SphericalHarmonicCoefficients, shcb::SphericalHarmonicCoefficients) = shca.c != shcb.c
+(shc::SphericalHarmonicCoefficients, value) = SphericalHarmonicCoefficients(shc.c .+ value,shc.R,shc.solid)
-(shc::SphericalHarmonicCoefficients, value) = SphericalHarmonicCoefficients(shc.c .- value,shc.R,shc.solid)
*(shc::SphericalHarmonicCoefficients, value) = SphericalHarmonicCoefficients(shc.c .* value,shc.R,shc.solid)
/(shc::SphericalHarmonicCoefficients, value) = SphericalHarmonicCoefficients(shc.c ./ value,shc.R,shc.solid)
+(value, shc::SphericalHarmonicCoefficients) = +(shc::SphericalHarmonicCoefficients, value)
-(value, shc::SphericalHarmonicCoefficients) = SphericalHarmonicCoefficients(value .- shc.c,shc.R,shc.solid)
*(value, shc::SphericalHarmonicCoefficients) = *(shc::SphericalHarmonicCoefficients, value)
LinearAlgebra.normalize(shc::SphericalHarmonicCoefficients,R::AbstractFloat) = SphericalHarmonicCoefficients([shc[l,m] * 1/(R^l) for l = 0:shc.L for m = -l:l],shc.R/R,shc.solid)
LinearAlgebra.normalize!(shc::SphericalHarmonicCoefficients,R::AbstractFloat) = SphericalHarmonicCoefficients([shc[l,m] *= 1/(R^l) for l = 0:shc.L for m = -l:l],shc.R /= R,shc.solid)

function +(shca::SphericalHarmonicCoefficients, shcb::SphericalHarmonicCoefficients)
  if shca.R != shcb.R
	  throw(DomainError(shca.R,"Coefficients do not have the same normalization factor."))
  end
  if shca.solid != shcb.solid
    throw(DomainError(shca.solid,"Coefficients do not have the same type."))
  end
  return SphericalHarmonicCoefficients(shca.c + shcb.c,shca.R,shca.solid)
end
function -(shca::SphericalHarmonicCoefficients, shcb::SphericalHarmonicCoefficients)
  if shca.R != shcb.R
    throw(DomainError(shca.R,"Coefficients do not have the same normalization factor."))
  end
  if shca.solid != shcb.solid
    throw(DomainError(shca.solid,"Coefficients do not have the same type."))
  end
  return SphericalHarmonicCoefficients(shca.c - shcb.c,shca.R,shca.solid)
end

## Convertions ##
# convert spherical coefficients to solid coefficients
function solid!(shc::SphericalHarmonicCoefficients)
  if !(shc.solid)
    for l = 0:shc.L
      for m = -l:l
        shc[l,m] *= sqrt((2*l+1)/(4*pi))
      end
    end
    shc.solid = true
  end
  return shc
end
# convert solid coefficients to spherical coefficients
function spherical!(shc::SphericalHarmonicCoefficients)
  if shc.solid
    for l = 0:shc.L
      for m = -l:l
        shc[l,m] *= sqrt(4*pi/(2*l+1))
      end
    end
    shc.solid = false
  end
  return shc
end

"""
    sphericalHarmonicsExpansion(Clm::SphericalHarmonicCoefficients, x::Variable, y::Variable, z::Variable)
*Description:*  Calculation of the spherical or solid harmonics expansion in Cartesian coordinates
                for given coefficients which define the maximum degree of the spherical harmonics\\
*Input:*  `Clm`       - Array with coefficients (length = (l+1)², l = max. deg. of the spherical harmonics)\\
          `x, y, z`   - Cartesian coordinates\\
*Output:*  Spherical/Solid harmonics expansion
"""
function sphericalHarmonicsExpansion(Clm::SphericalHarmonicCoefficients, x::T, y::T, z::T) where {T<:AbstractVariable}

  sum = 0
  for l in 0:Clm.L, m in -l:l
    if Clm[l,m] == 0
      continue
    end
    if Clm.solid
      # solid expansion
      sum += Clm[l,m] * zlm(l,m,x,y,z)
    else
      # spherical expansion
      sum += Clm[l,m] * rlylm(l,m,x,y,z)
    end
  end
  return sum
end

function sphericalHarmonicsExpansion(Clm::Vector{PolyVar{true}}, x::T, y::T, z::T, solid::Bool) where {T<:AbstractVariable}
  L = sqrt(length(Clm)) - 1
  if !isinteger(L)
    throw(DomainError(L,"Input vector needs to be of size (L+1)², where L ∈ ℕ₀."))
  else
    L = convert(Int,L)
  end

  sum = 0

  for l in 0:L
    for m in -l:l
      if Clm[l * (l+1) + m + 1] != 0
          if solid
              # solid expansion
              sum += Clm[l * (l+1) + m + 1] * zlm(l,m,x,y,z)
          else
              # spherical expansion
              sum += Clm[l * (l+1) + m + 1] * rlylm(l,m,x,y,z)
          end
      end
    end
  end
  return sum
end

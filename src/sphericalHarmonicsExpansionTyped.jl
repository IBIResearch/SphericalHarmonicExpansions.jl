# Type-generic expansion construction for the research branch.
# The legacy `sphericalHarmonicsExpansion` API remains unchanged in
# sphericalHarmonicsExpansion.jl and acts as the Float64 baseline.

"""
    sphericalHarmonicsExpansion_typed(T, Clm, x, y, z; scaling=:distributed)

Construct the complete spherical- or solid-harmonic expansion using floating
coefficient type `T`. The basis is selected from `Clm.solid`:
- `false`: radial spherical harmonics `r^l Y_l^m`
- `true`: solid harmonics `Z_l^m`

This method is intended for controlled Float64/BigFloat comparisons. It does
not silently change precision; for BigFloat use `setprecision` or the wrapper
`sphericalHarmonicsExpansion_bigfloat`.
"""
function sphericalHarmonicsExpansion_typed(::Type{T}, Clm::SphericalHarmonicCoefficients,
                                            x::Variable, y::Variable, z::Variable;
                                            scaling::Symbol=:distributed) where {T<:AbstractFloat}
    out = zero(T) * (x + y + z)
    for l in 0:Int(Clm.L), m in -l:l
        coefficient = T(Clm[l,m])
        iszero(coefficient) && continue
        basis = if Clm.solid
            zlm_typed(T, l, m, x, y, z; scaling=scaling)
        else
            rlylm_typed(T, l, m, x, y, z; scaling=scaling)
        end
        out += coefficient * basis
    end
    return out
end

"""
    sphericalHarmonicsExpansion_bigfloat(Clm, x, y, z; precision=256, scaling=:adaptive)

Construct the complete expansion with an explicit BigFloat precision.
"""
function sphericalHarmonicsExpansion_bigfloat(Clm::SphericalHarmonicCoefficients,
                                               x::Variable, y::Variable, z::Variable;
                                               precision::Integer=256,
                                               scaling::Symbol=:adaptive)
    return setprecision(BigFloat, Int(precision)) do
        sphericalHarmonicsExpansion_typed(BigFloat, Clm, x, y, z; scaling=scaling)
    end
end

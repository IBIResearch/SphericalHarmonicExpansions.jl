module SphericalHarmonicExpansions

using LinearAlgebra, Reexport, HDF5, Combinatorics
import StaticPolynomials
@reexport using DynamicPolynomials

include("sphericalHarmonic.jl")
export ylm, rlylm

include("sphericalHarmonicsExpansion.jl")
export SphericalHarmonicCoefficients, sphericalHarmonicsExpansion, solidHarmonicsExpansion,
       solid!, spherical!, normalize, normalize!

include("fastfunc.jl")
export @fastfunc, fastfunc

include("sphericalQuadrature.jl")

include("translation.jl")
export translation

using WignerD # Wigner d matrix
include("rotation.jl")
export rotation, pointReflection

end # module

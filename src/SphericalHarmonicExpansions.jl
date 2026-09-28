module SphericalHarmonicExpansions

using LinearAlgebra, Reexport, HDF5, Combinatorics
import MultivariatePolynomials: terms, degree, polynomial, coefficients
import StaticPolynomials
@reexport using TypedPolynomials

include("sphericalHarmonic.jl")
export ylm, rlylm, zlm, ylm_typed, ylm_bigfloat, rlylm_typed, zlm_typed, ylmCoefficientLog, ylmCoefficient, binomial_as

include("sphericalHarmonicsExpansion.jl")
include("sphericalHarmonicsExpansionTyped.jl")
export SphericalHarmonicCoefficients, sphericalHarmonicsExpansion,
       sphericalHarmonicsExpansion_typed, sphericalHarmonicsExpansion_bigfloat,
       solid!, spherical!, normalize, normalize!

include("PrecisionPolicy.jl")
export recommend_precision, evaluate_with_precision

include("fastfunc.jl")
export @fastfunc, fastfunc

include("sphericalQuadrature.jl")

include("translation.jl")
export translation

using WignerD # Wigner d matrix
include("rotation.jl")
export rotation, pointReflection

end # module

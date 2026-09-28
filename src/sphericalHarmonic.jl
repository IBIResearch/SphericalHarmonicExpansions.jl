# Spherical harmonics in Cartesian polynomial form.
#
# This working copy preserves the original public Float64 API while adding
# type-generic and adaptively scaled variants for the research project.

"""Convert an exact binomial coefficient to the selected floating type."""
binomial_as(::Type{T}, n::Integer, k::Integer) where {T<:AbstractFloat} = T(binomial(big(n), big(k)))
binomial_float(n::Integer, k::Integer) = binomial_as(Float64, n, k)

"""Logarithm of the normalisation factor for Y_l^m."""
function ylmCoefficientLog(::Type{T}, l::Integer, m::Integer) where {T<:AbstractFloat}
    l >= 0 || throw(DomainError(l, "l must be non-negative"))
    0 <= m <= l || throw(DomainError(m, "expected 0 <= m <= l"))
    logkpi = log(T(4) * T(pi))
    for i in (l-m+1):(l+m)
        logkpi += log(T(i))
    end
    return T(0.5) * (log(T(2*l + 1)) - logkpi)
end

ylmCoefficientLog(l::Integer, m::Integer) = ylmCoefficientLog(Float64, l, m)

function ylmCoefficient(::Type{T}, l::Integer, m::Integer) where {T<:AbstractFloat}
    out = exp(ylmCoefficientLog(T, l, m))
    if out != zero(T) && isfinite(out)
        return out
    end
    error("ylmCoefficient could not be represented in $(T)")
end

ylmCoefficient(l::Integer, m::Integer) = ylmCoefficient(Float64, l, m)

function ylmCosSinPolynomial(::Type{T}, m::Integer, x, y) where {T<:AbstractFloat}
    terms = [
        (isodd(j) ? -one(T) : one(T)) * binomial_as(T, m, 2*j) * y^(2*j) * x^(m-2*j)
        for j in 0:fld(m, 2)
    ]
    return polynomial(terms)
end

ylmCosSinPolynomial(m::Integer, x, y) = ylmCosSinPolynomial(Float64, m, x, y)

function ylmSinSinPolynomial(::Type{T}, m::Integer, x, y) where {T<:AbstractFloat}
    terms = [
        (isodd(j) ? -one(T) : one(T)) * binomial_as(T, m, 2*j + 1) * y^(2*j + 1) * x^(m-2*j-1)
        for j in 0:fld(m - 1, 2)
    ]
    return polynomial(terms)
end

ylmSinSinPolynomial(m::Integer, x, y) = ylmSinSinPolynomial(Float64, m, x, y)

function _normalize_poly(p, ::Type{T}) where {T<:AbstractFloat}
    coeffs = collect(coefficients(p))
    isempty(coeffs) && error("zero polynomial detected during adaptive scaling")
    all(isfinite, coeffs) || error("non-finite coefficient detected during adaptive scaling")
    maxc = maximum(abs, coeffs)
    (maxc == zero(T) || !isfinite(maxc)) && error("invalid maximum coefficient during adaptive scaling")
    return (one(T) / maxc) * p, log(maxc)
end

function _initial_z_polynomial(::Type{T}, l::Integer, m::Integer, x, y, z) where {T<:AbstractFloat}
    terms = [
        (isodd(k) ? -one(T) : one(T)) * binomial_as(T, l, k) * z^(2*(l-k))
        for k in 0:cld(l - abs(m), 2)
    ]
    return polynomial(terms) + zero(T) * (x + y)
end

function _xy_factor(::Type{T}, m::Integer, x, y) where {T<:AbstractFloat}
    if m > 0
        return sqrt(T(2)) * ylmCosSinPolynomial(T, m, x, y)
    elseif m < 0
        return sqrt(T(2)) * ylmSinSinPolynomial(T, abs(m), x, y)
    end
    return nothing
end

"""
    ylm_typed(T, l, m, x, y, z; scaling=:distributed)

Type-generic polynomial construction. Supported scaling modes:
- `:none`: logarithmic normalisation applied once at the end.
- `:distributed`: negative log-normalisation distributed over derivative steps.
- `:adaptive`: polynomial normalised after every derivative step.
"""
function ylm_typed(::Type{T}, l::Integer, m::Integer, x, y, z; scaling::Symbol=:distributed) where {T<:AbstractFloat}
    abs(m) <= l || throw(DomainError(m, "-l <= m <= l expected, but m = $m and l = $l"))
    l >= 0 || throw(DomainError(l, "l must be non-negative"))
    scaling in (:none, :distributed, :adaptive) || throw(ArgumentError("unsupported scaling mode: $scaling"))

    p = _initial_z_polynomial(T, l, m, x, y, z)
    nsteps = l + abs(m)

    if scaling == :adaptive
        logscale = zero(T)
        for i in 1:nsteps
            c = i <= l ? inv(T(2*i)) : one(T)
            p = c * differentiate(p, z, Val{1}())
            p, increment = _normalize_poly(p, T)
            logscale += increment
        end
        logscale += ylmCoefficientLog(T, l, abs(m))
        xy = _xy_factor(T, m, x, y)
        if xy !== nothing
            p = xy * p
        end
        p, increment = _normalize_poly(p, T)
        logscale += increment
        return exp(logscale) * p
    end

    logout = ylmCoefficientLog(T, l, abs(m))
    if scaling == :distributed && nsteps > 0 && logout < zero(T)
        step_scale = exp(logout / T(nsteps))
        final_scale = one(T)
    else
        step_scale = one(T)
        final_scale = exp(logout)
    end

    for i in 1:nsteps
        c = i <= l ? inv(T(2*i)) : one(T)
        p = step_scale * c * differentiate(p, z, Val{1}())
    end
    p = final_scale * p
    xy = _xy_factor(T, m, x, y)
    return xy === nothing ? p : xy * p
end

"""Backwards-compatible Float64 API using distributed scaling."""
ylm(l, m, x, y, z) = ylm_typed(Float64, l, m, x, y, z; scaling=:distributed)

"""Construct a BigFloat polynomial with an explicit precision."""
function ylm_bigfloat(l, m, x, y, z; precision::Integer=256, scaling::Symbol=:adaptive)
    return setprecision(BigFloat, Int(precision)) do
        ylm_typed(BigFloat, l, m, x, y, z; scaling=scaling)
    end
end

function _multinomial_big(i::Integer, j::Integer, k::Integer)
    n = big(i + j + k)
    return factorial(n) ÷ (factorial(big(i)) * factorial(big(j)) * factorial(big(k)))
end

function trinomialExpansion(::Type{T}, n::Integer, x, y, z) where {T<:AbstractFloat}
    multiindices = [(i, j, n-i-j) for i in 0:n for j in 0:n-i]
    terms = [T(_multinomial_big(i, j, k)) * x^(2*i) * y^(2*j) * z^(2*k) for (i,j,k) in multiindices]
    return polynomial(terms)
end

trinomialExpansion(n, x, y, z) = trinomialExpansion(Float64, n, x, y, z)

function rlylm_typed(::Type{T}, l, m, x, y, z; scaling::Symbol=:distributed) where {T<:AbstractFloat}
    p = ylm_typed(T, l, m, x, y, z; scaling=scaling)
    pout = zero(T) * (x + y + z)
    for term in terms(p)
        deg = degree(term)
        degR = l - deg
        degR >= 0 || error("negative radial degree")
        iseven(degR) || error("radial degree must be even")
        pout += trinomialExpansion(T, div(degR, 2), x, y, z) * term
    end
    return pout
end

rlylm(l, m, x, y, z) = rlylm_typed(Float64, l, m, x, y, z; scaling=:distributed)

function zlm_typed(::Type{T}, l, m, x, y, z; scaling::Symbol=:distributed) where {T<:AbstractFloat}
    return sqrt(T(4) * T(pi) / T(2*l + 1)) * rlylm_typed(T, l, m, x, y, z; scaling=scaling)
end

zlm(l, m, x, y, z) = zlm_typed(Float64, l, m, x, y, z; scaling=:distributed)

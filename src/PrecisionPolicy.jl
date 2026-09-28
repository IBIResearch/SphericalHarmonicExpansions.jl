const _PRECISION_POLICY_VERSION = "l-worst-case-v1-2026-06-15"
const _PRECISION_POLICY_PATH = normpath(joinpath(@__DIR__, "..", "data", "precision_policy_v1.csv"))
const _PRECISION_POLICY_HEADER = (
    "basis,variant,l,tolerance,recommended_numeric_type,precision_bits," *
    "measured_bound,scope,source_reference,policy_version"
)

struct _PrecisionPolicyRow
    basis::Symbol
    variant::Symbol
    l::Int
    tolerance_exponent::Int
    numeric_type::Symbol
    precision_bits::Int
    measured_bound::String
    scope::String
    source_reference::String
    policy_version::String
end

const _PRECISION_POLICY_CACHE = Ref{Union{Nothing,Vector{_PrecisionPolicyRow}}}(nothing)
const _PRECISION_POLICY_LOCK = ReentrantLock()

_csv_unquote(value::AbstractString) = replace(strip(value, ['"']), "\"\"" => "\"")

function _policy_tolerance_exponent(value::AbstractString)
    matched = match(r"^1e-(\d+)$", value)
    matched === nothing && error("unsupported policy tolerance format: $value")
    return parse(Int, only(matched.captures))
end

function _load_precision_policy()
    cached = _PRECISION_POLICY_CACHE[]
    cached === nothing || return cached

    lock(_PRECISION_POLICY_LOCK) do
        cached = _PRECISION_POLICY_CACHE[]
        cached === nothing || return cached

        lines = readlines(_PRECISION_POLICY_PATH)
        isempty(lines) && error("precision policy file is empty: $(_PRECISION_POLICY_PATH)")
        header = replace(lines[1], '\ufeff' => "", '"' => "")
        header == _PRECISION_POLICY_HEADER || error(
            "unexpected precision policy header in $(_PRECISION_POLICY_PATH)",
        )
        rows = _PrecisionPolicyRow[]
        sizehint!(rows, length(lines) - 1)

        for (offset, line) in enumerate(Iterators.drop(lines, 1))
            line_number = offset + 1
            isempty(strip(line)) && continue
            fields = _csv_unquote.(split(line, ','))
            length(fields) == 10 || error(
                "expected 10 policy columns at line $line_number, found $(length(fields))",
            )
            push!(rows, _PrecisionPolicyRow(
                Symbol(fields[1]),
                Symbol(fields[2]),
                parse(Int, fields[3]),
                _policy_tolerance_exponent(fields[4]),
                Symbol(fields[5]),
                parse(Int, fields[6]),
                fields[7],
                fields[8],
                fields[9],
                fields[10],
            ))
        end

        _PRECISION_POLICY_CACHE[] = rows
        return rows
    end
end

function _requested_tolerance_exponent(tolerance::Real)
    isfinite(tolerance) && tolerance > zero(tolerance) || return nothing
    value = BigFloat(tolerance)
    exponent = round(Int, -log10(value))
    reference = BigFloat(10)^(-exponent)
    relative_distance = abs(value - reference) / reference
    return relative_distance <= BigFloat(128) * eps(Float64) ? exponent : nothing
end

function _policy_result(status::Symbol, tolerance, basis, variant, mode;
                        numeric_type=nothing, precision_bits=nothing,
                        measured_bound=nothing, scope=nothing,
                        source_reference=nothing,
                        policy_version=_PRECISION_POLICY_VERSION)
    return (
        numeric_type=numeric_type,
        precision_bits=precision_bits,
        measured_bound=measured_bound,
        requested_tolerance=tolerance,
        basis=basis,
        variant=variant,
        mode=mode,
        status=status,
        scope=scope,
        source_reference=source_reference,
        policy_version=policy_version,
    )
end

"""
    recommend_precision(; l, tolerance, basis=:ylm,
                          variant=:distributed, mode=:worst_case)

Return the smallest measured arithmetic that satisfies `tolerance` for the
exact requested order `l`. The bundled policy covers only measured
`ylm`/`distributed` all-`m` worst cases and decimal tolerances `1e-1` through
`1e-30`. It does not interpolate orders or tolerances.

The result is a named tuple. A recommendation is usable only when
`status == :supported`; otherwise `numeric_type` and `precision_bits` are
`nothing`. The measured bound applies to the tested points and is not a global
mathematical guarantee over the sphere.
"""
function recommend_precision(; l::Integer, tolerance::Real,
                             basis::Symbol=:ylm,
                             variant::Symbol=:distributed,
                             mode::Symbol=:worst_case)
    mode == :worst_case || return _policy_result(
        :unsupported_mode, tolerance, basis, variant, mode,
    )

    rows = _load_precision_policy()
    configuration_rows = filter(row -> row.basis == basis && row.variant == variant, rows)
    isempty(configuration_rows) && return _policy_result(
        :configuration_not_measured, tolerance, basis, variant, mode,
    )

    order_rows = filter(row -> row.l == l, configuration_rows)
    isempty(order_rows) && return _policy_result(
        :order_not_measured, tolerance, basis, variant, mode,
    )

    isfinite(tolerance) && tolerance > zero(tolerance) || return _policy_result(
        :invalid_tolerance, tolerance, basis, variant, mode,
    )

    requested_exponent = _requested_tolerance_exponent(tolerance)
    max_exponent = maximum(row.tolerance_exponent for row in order_rows)
    if requested_exponent === nothing
        minimum_tolerance = BigFloat(10)^(-max_exponent)
        status = BigFloat(tolerance) < minimum_tolerance ?
                 :not_supported_by_measured_grid : :tolerance_not_measured
        return _policy_result(status, tolerance, basis, variant, mode)
    end

    requested_exponent > max_exponent && return _policy_result(
        :not_supported_by_measured_grid, tolerance, basis, variant, mode,
    )

    index = findfirst(row -> row.tolerance_exponent == requested_exponent, order_rows)
    index === nothing && return _policy_result(
        :tolerance_not_measured, tolerance, basis, variant, mode,
    )
    row = order_rows[index]

    numeric_type = row.numeric_type == :Float64 ? Float64 : BigFloat
    measured_bound = setprecision(BigFloat, 384) do
        parse(BigFloat, row.measured_bound)
    end
    return _policy_result(
        :supported, tolerance, basis, variant, mode;
        numeric_type=numeric_type,
        precision_bits=row.precision_bits,
        measured_bound=measured_bound,
        scope=row.scope,
        source_reference=row.source_reference,
        policy_version=row.policy_version,
    )
end

@polyvar _precision_policy_x _precision_policy_y _precision_policy_z

function _evaluate_ylm_at(::Type{T}, l, m, x, y, z, variant) where {T<:AbstractFloat}
    polynomial_value = ylm_typed(
        T, l, m,
        _precision_policy_x, _precision_policy_y, _precision_policy_z;
        scaling=variant,
    )
    point = (T(x), T(y), T(z))
    return polynomial_value(
        (_precision_policy_x, _precision_policy_y, _precision_policy_z) => point,
    )
end

"""
    evaluate_with_precision(basis, l, m, x, y, z;
                            tolerance, variant=:distributed, mode=:worst_case)

Select a measured precision with [`recommend_precision`](@ref) and evaluate
`Y_l^m` at a numeric Cartesian point. Unsupported requests return the policy
result with `value=nothing`; supported requests return the same provenance
fields and the computed `value`.
"""
function evaluate_with_precision(basis::Symbol, l::Integer, m::Integer,
                                 x::Real, y::Real, z::Real;
                                 tolerance::Real,
                                 variant::Symbol=:distributed,
                                 mode::Symbol=:worst_case)
    recommendation = recommend_precision(
        l=l,
        tolerance=tolerance,
        basis=basis,
        variant=variant,
        mode=mode,
    )
    recommendation.status == :supported || return merge(
        recommendation, (m=Int(m), value=nothing),
    )

    basis == :ylm || return merge(
        recommendation, (status=:configuration_not_measured, m=Int(m), value=nothing),
    )

    value = if recommendation.numeric_type === Float64
        _evaluate_ylm_at(Float64, l, m, x, y, z, variant)
    else
        setprecision(BigFloat, recommendation.precision_bits) do
            _evaluate_ylm_at(BigFloat, l, m, x, y, z, variant)
        end
    end
    return merge(recommendation, (m=Int(m), value=value))
end


#!/usr/bin/env julia
# Benchmark script comparing TypedPolynomials and DynamicPolynomials
# Usage: from the repo root or this folder run:
#   julia --project=benchmark benchmark.jl

using Pkg

# Activate a local project in this folder so the benchmark has its own environment
Pkg.activate(@__DIR__)

# Ensure required packages are available. This adds packages if they are missing.

function ensure_packages(pkgs::Vector{String})
	for pkg in pkgs
		try
			# try to load the package module without using `import` syntax that must be top-level
			Base.require(Main, Symbol(pkg))
		catch
			println("Adding package: $pkg")
			Pkg.add(pkg)
			Base.require(Main, Symbol(pkg))
		end
	end
end

ensure_packages(["BenchmarkTools", "Random", "DynamicPolynomials", "TypedPolynomials", "MultivariatePolynomials"])

# Import commonly-used modules at top-level so `import` is not executed in local scope
import DynamicPolynomials
import TypedPolynomials
import MultivariatePolynomials

using BenchmarkTools, Random

const RNG = MersenneTwister(1234)

# Helper: create nvars symbolic variables for a given package using the @polyvar macro
function create_polyvars(pkg::String, nvars::Int)
	names = ["x$(i)" for i in 1:nvars]
	mac = "@polyvar " * join(names, " ")
	# Load the package module (without bringing exported names into Main)
	Base.require(Main, Symbol(pkg))
	# Invoke the module-qualified macro to avoid ambiguity when multiple packages export @polyvar
	modmac = "$(pkg).@polyvar " * join(names, " ")
	Core.eval(Main, Meta.parse(modmac))
	return Symbol.(names)
end

# Build a random polynomial (as an expression evaluated in Main) using the package's variables.
# degree controls max exponent per variable, nterms is number of random monomials.
function build_random_poly(varsyms::Vector{Symbol}, degree::Int, nterms::Int)
	nvars = length(varsyms)
	coeffs = randn(RNG, nterms)
	function random_monomial()
		parts = String[]
		for i in 1:nvars
			p = rand(RNG, 0:degree)
			if p > 0
				push!(parts, string(varsyms[i], "^", p))
			end
		end
		isempty(parts) ? "1" : join(parts, "*")
	end
	terms = ["($(coeffs[i]))*$(random_monomial())" for i in 1:nterms]
	expr = join(terms, " + ")
	return Core.eval(Main, Meta.parse(expr))
end

# Try to evaluate polynomial p at numeric values vals (vector matching varsyms order).
# First try MultivariatePolynomials.evaluate with a Dict, then fall back to string substitution.
function try_evaluate(p, varsyms::Vector{Symbol}, vals::Vector)
	try
		# attempt dict-based evaluate
		mapping = Dict{Any,Any}()
		for (s,v) in zip(varsyms, vals)
			mapping[s] = v
		end
		return MultivariatePolynomials.evaluate(p, mapping)
	catch err1
		try
			# try positional evaluate (some implementations accept a tuple/vector)
			return MultivariatePolynomials.evaluate(p, vals)
		catch err2
			# fallback: stringify and replace variable names with numeric literals then eval
			s = string(p)
			for (sym, val) in zip(varsyms, vals)
				s = replace(s, string(sym) => "($(val))")
			end
			return Core.eval(Main, Meta.parse(s))
		end
	end
end

## Substitution tests have been removed per user request.

function bench_package(pkg::String; nvars=3, degree=6, nterms=80)
	println("\n=== Benchmarking $pkg ===")
	# create variables
	varsyms = create_polyvars(pkg, nvars)

	println("Building polynomials (definition)...")
	# Benchmark the construction (definition) of a polynomial, but assign the results
	# separately so the created polynomials are available for later benchmarks.
	println("(timing polynomial construction — result will be created afterwards)")
	@btime build_random_poly($varsyms, $degree, $nterms)
	p1 = build_random_poly(varsyms, degree, nterms)
	@btime build_random_poly($varsyms, $degree, $nterms)
	p2 = build_random_poly(varsyms, degree, nterms)

	println("Running basic operation benchmarks:")
	println("Addition:")
	@btime $p1 + $p2
	println("Multiplication:")
	@btime $p1 * $p2

	# Prepare numeric values for evaluation
	vals = randn(RNG, nvars)
	println("Evaluation at random numeric point:")
	@btime try_evaluate($p1, $varsyms, $vals)

	# (Substitution tests removed)
end

function main()
	# Run benchmarks across several size categories (small -> xlarge).
	# Each tuple is (nvars, degree, nterms).
	sizes = [
		(3, 2, 20),        # tiny
		(4, 5, 60),        # small
		(5, 10, 300),      # medium
		(6, 20, 1200),     # large
		(8, 40, 5000)      # xlarge (may be slow)
	]

	pkgs = ["DynamicPolynomials", "TypedPolynomials"]

	for (nvars, degree, nterms) in sizes
		println("\n===== SIZE: nvars=$nvars degree=$degree nterms=$nterms =====")
		for pkg in pkgs
			bench_package(pkg, nvars=nvars, degree=degree, nterms=nterms)
		end
	end
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end

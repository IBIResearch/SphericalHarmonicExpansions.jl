#!/usr/bin/env julia
# Parse benchmark `results.txt` (produced by `benchmark.jl`) and plot timings

using Pkg

# Use the benchmark folder as a project so packages are local
Pkg.activate(@__DIR__)

function ensure(pkgs::Vector{String})
    for p in pkgs
        try
            Base.require(Main, Symbol(p))
        catch
            println("Adding package: $p")
            Pkg.add(p)
            Base.require(Main, Symbol(p))
        end
    end
end

ensure(["DataFrames", "Plots", "Measures"])  # CSV not needed; file is plain text

import DataFrames
import Plots
import Measures

const RESULTS_PATH = joinpath(@__DIR__, "results.txt")

function parse_time(line::AbstractString)
    # match number and unit (ns, μs/us, ms, s)
    m = match(r"([0-9]*\.?[0-9]+)\s*(ns|μs|us|ms|s)", line)
    if m === nothing
        return missing
    end
    val = parse(Float64, m.captures[1])
    unit = m.captures[2]
    if unit == "ns"
        return val * 1e-9
    elseif unit == "μs" || unit == "us"
        return val * 1e-6
    elseif unit == "ms"
        return val * 1e-3
    elseif unit == "s"
        return val
    else
        return missing
    end
end

function parse_results(path::AbstractString)
    lines = readlines(path)
    n = length(lines)
    i = 1
    records = Vector{Dict{String,Any}}()

    size_re = r"^===== SIZE: nvars=(\d+) degree=(\d+) nterms=(\d+) ====="
    pkg_re = r"^=== Benchmarking\s+(\S+)\s+==="

    # initialize size variables so they are always defined in this scope
    nvars = 0
    degree = 0
    nterms = 0
    size_label = ""

    while i <= n
        line = lines[i]
        if occursin("===== SIZE:", line)
            sm = match(size_re, line)
            if sm === nothing
                i += 1
                continue
            end
            nvars = parse(Int, sm.captures[1])
            degree = parse(Int, sm.captures[2])
            nterms = parse(Int, sm.captures[3])
            size_label = "nvars=$(nvars) deg=$(degree) nterms=$(nterms)"
            i += 1
            continue
        elseif occursin("=== Benchmarking", line)
            pm = match(pkg_re, line)
            if pm === nothing
                i += 1
                continue
            end
            pkg = pm.captures[1]
            # initialize record with size info
            rec = Dict("pkg"=>pkg, "nvars"=>nvars, "degree"=>degree, "nterms"=>nterms, "size_label"=>size_label,
                "build"=>missing, "addition"=>missing, "multiplication"=>missing, "evaluation"=>missing)

            # look ahead until next package or size header or EOF
            j = i+1
            while j <= n && !occursin("===== SIZE:", lines[j]) && !occursin("=== Benchmarking", lines[j])
                l = strip(lines[j])
                if startswith(l, "(timing") || startswith(l, "Building polynomials")
                    # find next line that contains a timing (parse_time != missing)
                    k = j+1
                    while k <= n && parse_time(lines[k]) === missing
                        k += 1
                    end
                    if k <= n
                        t = parse_time(lines[k])
                        if t !== missing && rec["build"] === missing
                            rec["build"] = t
                        end
                    end
                    j = k
                elseif startswith(l, "Addition:")
                    # next non-empty line is timing
                    k = j+1
                    while k <= n && isempty(strip(lines[k]))
                        k += 1
                    end
                    if k <= n
                        rec["addition"] = parse_time(lines[k])
                    end
                    j = k
                elseif startswith(l, "Multiplication:")
                    k = j+1
                    while k <= n && isempty(strip(lines[k]))
                        k += 1
                    end
                    if k <= n
                        rec["multiplication"] = parse_time(lines[k])
                    end
                    j = k
                elseif occursin("Evaluation", l)
                    # next non-empty line is timing
                    k = j+1
                    while k <= n && isempty(strip(lines[k]))
                        k += 1
                    end
                    if k <= n
                        rec["evaluation"] = parse_time(lines[k])
                    end
                    j = k
                end
                j += 1
            end
            push!(records, rec)
            i = j
            continue
        end
        i += 1
    end

    return records
end

function records_to_df(records)
    cols = Dict(:pkg=>String[], :nvars=>Int[], :degree=>Int[], :nterms=>Int[], :size_label=>String[],
        :build_s=>Float64[], :addition_s=>Float64[], :multiplication_s=>Float64[], :evaluation_s=>Float64[])

    for r in records
        push!(cols[:pkg], r["pkg"])
        push!(cols[:nvars], r["nvars"])
        push!(cols[:degree], r["degree"])
        push!(cols[:nterms], r["nterms"])
        push!(cols[:size_label], r["size_label"])
        push!(cols[:build_s], isnothing(r["build"]) || r["build"]===missing ? NaN : r["build"]) 
        push!(cols[:addition_s], isnothing(r["addition"]) || r["addition"]===missing ? NaN : r["addition"]) 
        push!(cols[:multiplication_s], isnothing(r["multiplication"]) || r["multiplication"]===missing ? NaN : r["multiplication"]) 
        push!(cols[:evaluation_s], isnothing(r["evaluation"]) || r["evaluation"]===missing ? NaN : r["evaluation"]) 
    end

    return DataFrames.DataFrame(cols)
end

function plot_from_df(df)
    Plots.gr()
    size_labels = unique(df.size_label)
    size_labels_ml = [replace(s, " " => "\n") for s in size_labels]
    pkgs = unique(df.pkg)
    x = 1:length(size_labels)

    # ✅ Use Measures.Length units (e.g. 10mm) to avoid type errors
    p = Plots.plot(layout=(2,2), size=(1200, 800),
                   left_margin=10Measures.mm, right_margin=5Measures.mm,
                   top_margin=5Measures.mm, bottom_margin=12Measures.mm)

    metrics = [(:build_s, "Build time / s", true),
               (:addition_s, "Addition time / s", true),
               (:multiplication_s, "Multiplication time / s", true),
               (:evaluation_s, "Evaluation time / s", true)]

    for (idx, (col, ylabel, logy)) in enumerate(metrics)
        # collect all values for this metric to compute tick ranges
        all_vals = Float64[]
        for pkg in pkgs
            row = df[df.pkg .== pkg, :]
            vals = [row[row.size_label .== sl, col][1] for sl in size_labels]
            append!(all_vals, filter(!isnan, vals))
            # plot each package's line; request legend in top-left
            Plots.plot!(p, x, vals, label=pkg, marker=:o, subplot=idx, legend=:topleft)
        end

        Plots.xlabel!(p, "size category", subplot=idx)
        Plots.xticks!(p, x, size_labels_ml, subplot=idx)
        Plots.ylabel!(p, ylabel, subplot=idx)

        if logy && !isempty(all_vals)
            # compute decade ticks with step=1 (every order of magnitude)
            mn = minimum(all_vals)
            mx = maximum(all_vals)
            emin = floor(Int, log10(max(mn, 1e-12)))
            emax = ceil(Int, log10(max(mx, 1e-12)))
            ticks = 10.0 .^ (emin:emax)
            Plots.plot!(p, yscale=:log10, yticks=ticks, subplot=idx)
        end
    end

    out = joinpath(@__DIR__, "results_plot.png")
    println("Saving plot to: ", out)
    Plots.savefig(p, out)
end

function main()
    if !isfile(RESULTS_PATH)
        error("results.txt not found at $(RESULTS_PATH). Run benchmark.jl first and save output to results.txt")
    end
    records = parse_results(RESULTS_PATH)
    df = records_to_df(records)
    println(df)
    plot_from_df(df)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

using Pkg
Pkg.activate(".")

using Dates
using Logging
using Printf
using Random
using Statistics

include("../generate_three_case_curves.jl")
include("run_upper_bound_audit.jl")

safe_num(x::Real) = isfinite(x) ? float(x) : NaN

parse_bool(s::AbstractString) = lowercase(strip(s)) in ("1", "true", "yes", "y", "on")

function eval_with_logging(thunk::Function, suppress_warnings::Bool)
    if suppress_warnings
        return with_logger(ConsoleLogger(stderr, Logging.Error)) do
            thunk()
        end
    end
    return thunk()
end

function parse_args(args::Vector{String})
    opts = Dict{String, String}()
    i = 1
    while i <= length(args)
        arg = args[i]
        if startswith(arg, "--")
            if occursin("=", arg)
                key, value = split(arg[3:end], "=", limit = 2)
                opts[key] = value
                i += 1
            else
                key = arg[3:end]
                if i == length(args)
                    error("Missing value for option --$key")
                end
                opts[key] = args[i + 1]
                i += 2
            end
        else
            opts["output-dir"] = arg
            i += 1
        end
    end
    return opts
end

function maybe_git_commit()
    try
        return strip(read(`git rev-parse --short HEAD`, String))
    catch
        return "unknown"
    end
end

function summarize_runs(runs)
    values = [run.r for run in runs]
    if isempty(values)
        return (
            best_idx = 0,
            best_seed = -1,
            best_r = NaN,
        )
    end

    best_idx = argmax(values)
    return (
        best_idx = best_idx,
        best_seed = runs[best_idx].seed,
        best_r = runs[best_idx].r,
    )
end

function sample_case1_restart(eps::Real, seed::Int; tol::Real)
    Random.seed!(seed)
    u0 = UMIN + rand() * (UMAX - UMIN)
    r, u, x = optimize_case1_at_epsilon(eps, u0; tol = tol)
    return (r = safe_num(r), u = u, x = x, seed = seed, start_u = u0)
end

function dense_scan_case1(
    eps::Real;
    u_min::Real = -3.0,
    u_max::Real = 1.0,
    u_step::Real = 0.01,
    tol::Real = 1e-7,
    suppress_warnings::Bool = true,
)
    best_r = -Inf
    best_u = float(u_min)
    best_x = 10.0^best_u

    for u in u_min:u_step:u_max
        cand = eval_with_logging(suppress_warnings) do
            safe_case1_score(eps, u; tol = tol)
        end
        r, uu, x = cand
        if r > best_r
            best_r = r
            best_u = uu
            best_x = x
        elseif r == best_r && abs(uu) < abs(best_u)
            best_u = uu
            best_x = x
        end
    end

    return (r = safe_num(best_r), u = best_u, x = best_x)
end

"""
    n1_theorem_formula(eps)

Draft analytical expression used as a falsifiable theorem candidate:
R(ϵ) = R0 * max(1 - ϵ / ϵc, 0)^α

The constants are fixed in this sprint and documented in manuscript TODO-PROOF
markers. This function is intentionally explicit and deterministic.
"""
function n1_theorem_formula(eps::Real)
    r0 = 0.893925772022
    eps_cutoff = 0.07
    alpha = 8.0

    e = clamp(float(eps), 0.0, 0.25)
    t = max(1.0 - e / eps_cutoff, 0.0)
    return r0 * t^alpha
end

function write_points_csv(path::AbstractString, points)
    open(path, "w") do io
        println(
            io,
            "epsilon,numeric_R,opt_x_numeric,scan_R,opt_x_scan,formula_R," *
            "abs_diff_numeric_scan,rel_diff_numeric_scan," *
            "abs_diff_numeric_formula,rel_diff_numeric_formula,best_seed",
        )
        for p in points
            row = [
                fmt_csv_value(p.eps),
                fmt_csv_value(p.numeric_r),
                fmt_csv_value(p.opt_x_numeric),
                fmt_csv_value(p.scan_r),
                fmt_csv_value(p.opt_x_scan),
                fmt_csv_value(p.formula_r),
                fmt_csv_value(p.abs_diff_ns),
                fmt_csv_value(p.rel_diff_ns),
                fmt_csv_value(p.abs_diff_nf),
                fmt_csv_value(p.rel_diff_nf),
                string(p.best_seed),
            ]
            println(io, join(row, ","))
        end
    end
end

function write_summary_csv(path::AbstractString, metrics::Vector{Pair{String, String}})
    open(path, "w") do io
        println(io, "metric,value")
        for (k, v) in metrics
            println(io, k, ",", v)
        end
    end
end

function write_meta_json(path::AbstractString; kwargs...)
    kv = Dict{String, Any}(string(k) => v for (k, v) in pairs(kwargs))
    keys_sorted = sort(collect(keys(kv)))
    open(path, "w") do io
        println(io, "{")
        for (idx, k) in enumerate(keys_sorted)
            v = kv[k]
            comma = idx == length(keys_sorted) ? "" : ","
            if v isa AbstractString
                vv = replace(v, "\\" => "\\\\")
                vv = replace(vv, '"' => "\\\"")
                println(io, "  \"", k, "\": \"", vv, "\"", comma)
            elseif v isa Bool
                println(io, "  \"", k, "\": ", lowercase(string(v)), comma)
            elseif v isa Integer || v isa AbstractFloat
                println(io, "  \"", k, "\": ", v, comma)
            else
                vv = replace(string(v), "\\" => "\\\\")
                vv = replace(vv, '"' => "\\\"")
                println(io, "  \"", k, "\": \"", vv, "\"", comma)
            end
        end
        println(io, "}")
    end
end

function run_n1_theorem_validation(;
    output_dir::AbstractString = "results/sprint2_n1",
    eps_min::Real = 0.0,
    eps_max::Real = 0.25,
    eps_step::Real = 0.005,
    restarts::Int = 5,
    tol::Real = 1e-7,
    seed::Int = 20260319,
    suppress_warnings::Bool = true,
    u_min::Real = -3.0,
    u_max::Real = 1.0,
    u_step::Real = 0.01,
    violation_tol::Real = 1e-7,
)
    started = time()
    created_at = Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS.s")
    mkpath(output_dir)

    eps_grid = epsilon_values(eps_min, eps_max, eps_step)
    points = NamedTuple[]

    for (idx, eps) in enumerate(eps_grid)
        seeds = [seed + idx * 10_000 + j for j in 1:restarts]
        runs = eval_with_logging(suppress_warnings) do
            [sample_case1_restart(eps, s; tol = tol) for s in seeds]
        end
        s = summarize_runs(runs)
        best = runs[s.best_idx]

        scan = dense_scan_case1(
            eps;
            u_min = u_min,
            u_max = u_max,
            u_step = u_step,
            tol = tol,
            suppress_warnings = suppress_warnings,
        )

        numeric_r = best.r
        numeric_x = best.x
        numeric_seed = s.best_seed
        if isfinite(scan.r) && (!isfinite(numeric_r) || scan.r > numeric_r)
            numeric_r = scan.r
            numeric_x = scan.x
            numeric_seed = -1
        end

        formula_r = n1_theorem_formula(eps)

        abs_diff_ns = abs(numeric_r - scan.r)
        abs_diff_nf = abs(numeric_r - formula_r)
        denom = max(abs(numeric_r), 1e-12)
        rel_diff_ns = abs_diff_ns / denom
        rel_diff_nf = abs_diff_nf / denom

        push!(
            points,
            (
                eps = float(eps),
                numeric_r = numeric_r,
                opt_x_numeric = numeric_x,
                scan_r = scan.r,
                opt_x_scan = scan.x,
                formula_r = formula_r,
                abs_diff_ns = abs_diff_ns,
                rel_diff_ns = rel_diff_ns,
                abs_diff_nf = abs_diff_nf,
                rel_diff_nf = rel_diff_nf,
                best_seed = numeric_seed,
            ),
        )

        @printf(
            "[%d/%d] eps=%.6f numeric=%.6g scan=%.6g formula=%.6g |ns|=%.3g |nf|=%.3g\n",
            idx,
            length(eps_grid),
            eps,
            numeric_r,
            scan.r,
            formula_r,
            abs_diff_ns,
            abs_diff_nf,
        )
    end

    points_csv = joinpath(output_dir, "n1_points.csv")
    write_points_csv(points_csv, points)

    audit = run_upper_bound_audit(
        output_dir = output_dir,
        eps_min = eps_min,
        eps_max = eps_max,
        eps_step = eps_step,
        tol = tol,
        suppress_warnings = suppress_warnings,
        violation_tol = violation_tol,
        verbose = true,
    )

    abs_ns = [p.abs_diff_ns for p in points]
    rel_ns = [p.rel_diff_ns for p in points]
    abs_nf = [p.abs_diff_nf for p in points]
    rel_nf = [p.rel_diff_nf for p in points]
    max_abs_ns = maximum(abs_ns)
    mean_abs_ns = mean(abs_ns)
    max_rel_ns = maximum(rel_ns)
    mean_rel_ns = mean(rel_ns)
    max_abs_nf = maximum(abs_nf)
    mean_abs_nf = mean(abs_nf)
    max_rel_nf = maximum(rel_nf)
    mean_rel_nf = mean(rel_nf)
    threshold_ok = max_abs_ns <= 5e-4

    runtime = time() - started
    summary_csv = joinpath(output_dir, "n1_summary.csv")
    metrics = Pair{String, String}[
        "eps_count" => string(length(eps_grid)),
        "eps_min" => fmt_csv_value(eps_min),
        "eps_max" => fmt_csv_value(eps_max),
        "eps_step" => fmt_csv_value(eps_step),
        "restarts" => string(restarts),
        "tol" => fmt_csv_value(tol),
        "seed" => string(seed),
        "u_min" => fmt_csv_value(u_min),
        "u_max" => fmt_csv_value(u_max),
        "u_step" => fmt_csv_value(u_step),
        "max_abs_numeric_scan" => fmt_csv_value(max_abs_ns),
        "mean_abs_numeric_scan" => fmt_csv_value(mean_abs_ns),
        "max_rel_numeric_scan" => fmt_csv_value(max_rel_ns),
        "mean_rel_numeric_scan" => fmt_csv_value(mean_rel_ns),
        "max_abs_numeric_formula" => fmt_csv_value(max_abs_nf),
        "mean_abs_numeric_formula" => fmt_csv_value(mean_abs_nf),
        "max_rel_numeric_formula" => fmt_csv_value(max_rel_nf),
        "mean_rel_numeric_formula" => fmt_csv_value(mean_rel_nf),
        "numeric_scan_within_5e-4" => string(threshold_ok),
        "bound_violation_tol" => fmt_csv_value(violation_tol),
        "bound_violation_count" => string(audit.violation_count),
        "bound_max_violation" => fmt_csv_value(audit.max_violation),
        "bound_audit_csv" => "upper_bound_audit.csv",
        "runtime_seconds" => fmt_csv_value(runtime),
    ]
    write_summary_csv(summary_csv, metrics)

    meta_json = joinpath(output_dir, "n1_meta.json")
    write_meta_json(
        meta_json;
        created_at = created_at,
        script = "src/analysis/run_n1_theorem_validation.jl",
        output_dir = output_dir,
        points_csv = basename(points_csv),
        summary_csv = basename(summary_csv),
        meta_file = basename(meta_json),
        upper_bound_audit_csv = basename(audit.audit_csv),
        eps_min = float(eps_min),
        eps_max = float(eps_max),
        eps_step = float(eps_step),
        eps_count = length(eps_grid),
        restarts = restarts,
        tol = float(tol),
        seed = seed,
        suppress_warnings = suppress_warnings,
        u_min = float(u_min),
        u_max = float(u_max),
        u_step = float(u_step),
        theorem_formula = "R0*max(1-eps/eps_cutoff,0)^alpha",
        theorem_formula_R0 = 0.893925772022,
        theorem_formula_eps_cutoff = 0.07,
        theorem_formula_alpha = 8.0,
        bound_violation_tol = float(violation_tol),
        bound_violation_count = audit.violation_count,
        bound_max_violation = audit.max_violation,
        numeric_scan_max_abs_threshold = 5e-4,
        numeric_scan_within_threshold = threshold_ok,
        git_commit = maybe_git_commit(),
        runtime_seconds = runtime,
    )

    return (
        points_csv = points_csv,
        summary_csv = summary_csv,
        meta_json = meta_json,
        upper_bound_audit_csv = audit.audit_csv,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    opts = parse_args(ARGS)
    out_dir = get(opts, "output-dir", "results/sprint2_n1")
    eps_min = parse(Float64, get(opts, "eps-min", "0.0"))
    eps_max = parse(Float64, get(opts, "eps-max", "0.25"))
    eps_step = parse(Float64, get(opts, "eps-step", "0.005"))
    restarts = parse(Int, get(opts, "restarts", "5"))
    tol = parse(Float64, get(opts, "tol", "1e-7"))
    seed = parse(Int, get(opts, "seed", "20260319"))
    suppress_warnings = parse_bool(get(opts, "suppress-warnings", "true"))
    u_min = parse(Float64, get(opts, "u-min", "-3.0"))
    u_max = parse(Float64, get(opts, "u-max", "1.0"))
    u_step = parse(Float64, get(opts, "u-step", "0.01"))
    violation_tol = parse(Float64, get(opts, "violation-tol", "1e-7"))

    outputs = run_n1_theorem_validation(
        output_dir = out_dir,
        eps_min = eps_min,
        eps_max = eps_max,
        eps_step = eps_step,
        restarts = restarts,
        tol = tol,
        seed = seed,
        suppress_warnings = suppress_warnings,
        u_min = u_min,
        u_max = u_max,
        u_step = u_step,
        violation_tol = violation_tol,
    )

    println("N=1 theorem validation complete:")
    println("  ", outputs.points_csv)
    println("  ", outputs.summary_csv)
    println("  ", outputs.meta_json)
    println("  ", outputs.upper_bound_audit_csv)
end

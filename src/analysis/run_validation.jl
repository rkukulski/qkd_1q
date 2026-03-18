using Pkg
Pkg.activate(".")

using Dates
using Logging
using Printf
using Random
using Statistics

include("../generate_three_case_curves.jl")
include("../bounds.jl")

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
            # Single positional argument is treated as output-dir.
            opts["output-dir"] = arg
            i += 1
        end
    end
    return opts
end

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

function choose_winner(r1::Real, r2::Real, r3::Real)
    v1 = isfinite(r1) ? r1 : -Inf
    v2 = isfinite(r2) ? r2 : -Inf
    v3 = isfinite(r3) ? r3 : -Inf
    if v1 >= v2 && v1 >= v3
        return "case1"
    elseif v2 >= v1 && v2 >= v3
        return "case2"
    else
        return "case3"
    end
end

function sample_case1_restart(eps::Real, seed::Int; tol::Real)
    Random.seed!(seed)
    u0 = UMIN + rand() * (UMAX - UMIN)
    r, u, x = optimize_case1_at_epsilon(eps, u0; tol = tol)
    return (r = safe_num(r), u = u, x = x, seed = seed, start_u = u0)
end

function sample_case2_restart(eps::Real, seed::Int; tol::Real)
    Random.seed!(seed)
    u0 = UMIN + rand() * (UMAX - UMIN)
    v0 = VMIN + rand() * (VMAX - VMIN)
    r, u, v, x, y = optimize_case2_at_epsilon(eps, (u0, v0); tol = tol)
    return (r = safe_num(r), u = u, v = v, x = x, y = y, seed = seed, start_u = u0, start_v = v0)
end

function summarize_runs(runs)
    values = [run.r for run in runs]
    m = isempty(values) ? NaN : mean(values)
    s = length(values) <= 1 ? 0.0 : std(values)
    min_v = isempty(values) ? NaN : minimum(values)
    max_v = isempty(values) ? NaN : maximum(values)

    if isempty(values)
        return (
            mean = NaN,
            std = NaN,
            min = NaN,
            max = NaN,
            best_seed = -1,
            best_idx = 0,
            best_r = NaN,
        )
    end

    best_idx = argmax(values)
    best_seed = runs[best_idx].seed
    best_r = runs[best_idx].r

    return (
        mean = m,
        std = s,
        min = min_v,
        max = max_v,
        best_seed = best_seed,
        best_idx = best_idx,
        best_r = best_r,
    )
end

function case_value(point, case_name::String)
    if case_name == "case1"
        return point.r1
    elseif case_name == "case2"
        return point.r2
    else
        return point.r3
    end
end

function estimate_crossovers(points)
    crossovers = NamedTuple[]
    if length(points) < 2
        return crossovers
    end

    for i in 2:length(points)
        prev = points[i - 1]
        curr = points[i]
        if prev.winner == curr.winner
            continue
        end

        from_case = prev.winner
        to_case = curr.winner
        prev_diff = case_value(prev, from_case) - case_value(prev, to_case)
        curr_diff = case_value(curr, from_case) - case_value(curr, to_case)

        est_eps = (prev.eps + curr.eps) / 2
        denom = curr_diff - prev_diff
        if isfinite(prev_diff) && isfinite(curr_diff) && abs(denom) > 1e-12
            est_eps = prev.eps + (-prev_diff) * (curr.eps - prev.eps) / denom
            est_eps = clamp(est_eps, min(prev.eps, curr.eps), max(prev.eps, curr.eps))
        end

        push!(
            crossovers,
            (
                from_case = from_case,
                to_case = to_case,
                est_eps = est_eps,
                left_eps = prev.eps,
                right_eps = curr.eps,
            ),
        )
    end

    return crossovers
end

function json_escape(s::AbstractString)
    t = replace(s, "\\" => "\\\\")
    t = replace(t, '"' => "\\\"")
    t = replace(t, "\n" => "\\n")
    return t
end

function maybe_git_commit()
    try
        return strip(read(`git rev-parse --short HEAD`, String))
    catch
        return "unknown"
    end
end

function write_points_csv(path::AbstractString, points)
    open(path, "w") do io
        println(
            io,
            "epsilon,case1_R,case2_R,case3_R,upper_bound,gap_case1,gap_case2,gap_case3,best_case,case1_x_opt,case2_x_opt,case2_y_opt,case1_best_seed,case2_best_seed",
        )
        for p in points
            row = [
                fmt_csv_value(p.eps),
                fmt_csv_value(p.r1),
                fmt_csv_value(p.r2),
                fmt_csv_value(p.r3),
                fmt_csv_value(p.upper),
                fmt_csv_value(p.g1),
                fmt_csv_value(p.g2),
                fmt_csv_value(p.g3),
                p.winner,
                fmt_csv_value(p.case1_x),
                fmt_csv_value(p.case2_x),
                fmt_csv_value(p.case2_y),
                string(p.case1_best_seed),
                string(p.case2_best_seed),
            ]
            println(io, join(row, ","))
        end
    end
end

function write_summary_csv(path::AbstractString, stats_rows, crossovers)
    open(path, "w") do io
        println(
            io,
            "row_type,case,epsilon,mean,std,min,max,best_seed,best_R,crossover_from,crossover_to,crossover_est_eps,bracket_eps_left,bracket_eps_right",
        )

        for row in stats_rows
            out = [
                "stats",
                row.case,
                fmt_csv_value(row.eps),
                fmt_csv_value(row.mean),
                fmt_csv_value(row.std),
                fmt_csv_value(row.min),
                fmt_csv_value(row.max),
                string(row.best_seed),
                fmt_csv_value(row.best_r),
                "",
                "",
                "",
                "",
                "",
            ]
            println(io, join(out, ","))
        end

        for c in crossovers
            out = [
                "crossover",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                c.from_case,
                c.to_case,
                fmt_csv_value(c.est_eps),
                fmt_csv_value(c.left_eps),
                fmt_csv_value(c.right_eps),
            ]
            println(io, join(out, ","))
        end
    end
end

function write_meta_json(path::AbstractString; kwargs...)
    kv = Dict{String, Any}(string(k) => v for (k, v) in pairs(kwargs))
    open(path, "w") do io
        println(io, "{")
        keys_sorted = sort(collect(keys(kv)))
        for (idx, key) in enumerate(keys_sorted)
            value = kv[key]
            comma = idx < length(keys_sorted) ? "," : ""
            if value isa AbstractString
                println(io, "  \"$(json_escape(key))\": \"$(json_escape(value))\"$comma")
            elseif value isa Bool
                println(io, "  \"$(json_escape(key))\": $(value ? "true" : "false")$comma")
            elseif value isa Real
                println(io, "  \"$(json_escape(key))\": $(fmt_csv_value(value))$comma")
            else
                println(io, "  \"$(json_escape(key))\": \"$(json_escape(string(value)))\"$comma")
            end
        end
        println(io, "}")
    end
end

function run_validation(
    ;
    output_dir::AbstractString = ".",
    eps_min::Real = 0.0,
    eps_max::Real = 0.25,
    eps_step::Real = 0.01,
    restarts::Int = 3,
    tol::Real = 1e-7,
    seed::Int = 20260318,
    suppress_warnings::Bool = true,
)
    if restarts <= 0
        error("restarts must be positive.")
    end

    mkpath(output_dir)

    eps_grid = epsilon_values(eps_min, eps_max, eps_step)

    points = NamedTuple[]
    stats_rows = NamedTuple[]

    t_start = time()
    for (idx, eps) in enumerate(eps_grid)
        seeds = [seed + 1000 * (idx - 1) + i for i in 1:restarts]

        case1_runs = eval_with_logging(suppress_warnings) do
            [sample_case1_restart(eps, s; tol = tol) for s in seeds]
        end
        case2_runs = eval_with_logging(suppress_warnings) do
            [sample_case2_restart(eps, s; tol = tol) for s in seeds]
        end
        r3 = eval_with_logging(suppress_warnings) do
            safe_num(score_six_state(eps; tol = tol))
        end

        s1 = summarize_runs(case1_runs)
        s2 = summarize_runs(case2_runs)

        best_case1 = case1_runs[s1.best_idx]
        best_case2 = case2_runs[s2.best_idx]

        # case3 has no explicit restart optimization; use same value for summary rows.
        s3 = (
            mean = r3,
            std = 0.0,
            min = r3,
            max = r3,
            best_seed = seeds[1],
            best_r = r3,
        )

        upper = safe_num(upper_bound(eps))
        g1 = isfinite(upper) && isfinite(best_case1.r) ? upper - best_case1.r : NaN
        g2 = isfinite(upper) && isfinite(best_case2.r) ? upper - best_case2.r : NaN
        g3 = isfinite(upper) && isfinite(r3) ? upper - r3 : NaN

        winner = choose_winner(best_case1.r, best_case2.r, r3)

        push!(
            points,
            (
                eps = eps,
                r1 = best_case1.r,
                r2 = best_case2.r,
                r3 = r3,
                upper = upper,
                g1 = g1,
                g2 = g2,
                g3 = g3,
                winner = winner,
                case1_x = best_case1.x,
                case2_x = best_case2.x,
                case2_y = best_case2.y,
                case1_best_seed = s1.best_seed,
                case2_best_seed = s2.best_seed,
            ),
        )

        push!(stats_rows, (case = "case1", eps = eps, mean = s1.mean, std = s1.std, min = s1.min, max = s1.max, best_seed = s1.best_seed, best_r = s1.best_r))
        push!(stats_rows, (case = "case2", eps = eps, mean = s2.mean, std = s2.std, min = s2.min, max = s2.max, best_seed = s2.best_seed, best_r = s2.best_r))
        push!(stats_rows, (case = "case3", eps = eps, mean = s3.mean, std = s3.std, min = s3.min, max = s3.max, best_seed = s3.best_seed, best_r = s3.best_r))

        println(
            @sprintf(
                "[%d/%d] eps=%.6f | case1=%.6g case2=%.6g case3=%.6g upper=%.6g winner=%s",
                idx,
                length(eps_grid),
                eps,
                best_case1.r,
                best_case2.r,
                r3,
                upper,
                winner,
            ),
        )
    end

    crossovers = estimate_crossovers(points)

    points_csv = joinpath(output_dir, "validation_points.csv")
    summary_csv = joinpath(output_dir, "validation_summary.csv")
    meta_json = joinpath(output_dir, "validation_meta.json")

    write_points_csv(points_csv, points)
    write_summary_csv(summary_csv, stats_rows, crossovers)

    write_meta_json(
        meta_json;
        created_at = string(now()),
        script = "src/analysis/run_validation.jl",
        eps_min = float(eps_min),
        eps_max = float(eps_max),
        eps_step = float(eps_step),
        eps_count = length(eps_grid),
        restarts = restarts,
        tol = float(tol),
        seed = seed,
        suppress_warnings = suppress_warnings,
        git_commit = maybe_git_commit(),
        runtime_seconds = time() - t_start,
        output_dir = output_dir,
        points_csv = basename(points_csv),
        summary_csv = basename(summary_csv),
        meta_file = basename(meta_json),
    )

    println("Saved: $points_csv")
    println("Saved: $summary_csv")
    println("Saved: $meta_json")

    return (points_csv = points_csv, summary_csv = summary_csv, meta_json = meta_json)
end

if abspath(PROGRAM_FILE) == @__FILE__
    opts = parse_args(ARGS)

    output_dir = get(opts, "output-dir", ".")
    eps_min = parse(Float64, get(opts, "eps-min", "0.0"))
    eps_max = parse(Float64, get(opts, "eps-max", "0.25"))
    eps_step = parse(Float64, get(opts, "eps-step", "0.01"))
    restarts = parse(Int, get(opts, "restarts", "3"))
    tol = parse(Float64, get(opts, "tol", "1e-7"))
    seed = parse(Int, get(opts, "seed", "20260318"))
    suppress_warnings = parse_bool(get(opts, "suppress-warnings", "true"))

    run_validation(
        output_dir = output_dir,
        eps_min = eps_min,
        eps_max = eps_max,
        eps_step = eps_step,
        restarts = restarts,
        tol = tol,
        seed = seed,
        suppress_warnings = suppress_warnings,
    )
end

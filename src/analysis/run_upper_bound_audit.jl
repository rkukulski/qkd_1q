using Pkg
Pkg.activate(".")

using Logging
using Printf

if !isdefined(Main, :QKDProtocol)
    include("../generate_three_case_curves.jl")
end
if !isdefined(Main, :upper_bound)
    include("../bounds.jl")
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

function run_upper_bound_audit(;
    output_dir::AbstractString = "results/sprint2_n1",
    eps_min::Real = 0.0,
    eps_max::Real = 0.25,
    eps_step::Real = 0.005,
    tol::Real = 1e-7,
    suppress_warnings::Bool = true,
    violation_tol::Real = 1e-7,
    verbose::Bool = true,
)
    mkpath(output_dir)
    eps_grid = epsilon_values(eps_min, eps_max, eps_step)
    audit_csv = joinpath(output_dir, "upper_bound_audit.csv")

    prev_u_case1 = nothing
    prev_uv_case2 = nothing

    violation_count = 0
    max_violation = 0.0

    open(audit_csv, "w") do io
        println(
            io,
            "epsilon,case1_R,case2_R,case3_R,max_case_R,upper_bound,violation,violates",
        )

        for (idx, eps) in enumerate(eps_grid)
            r1, u1, _ = eval_with_logging(suppress_warnings) do
                optimize_case1_at_epsilon(eps, prev_u_case1; tol = tol)
            end
            prev_u_case1 = u1

            r2, u2, v2, _, _ = eval_with_logging(suppress_warnings) do
                optimize_case2_at_epsilon(eps, prev_uv_case2; tol = tol)
            end
            prev_uv_case2 = (u2, v2)

            r3 = eval_with_logging(suppress_warnings) do
                score_six_state(eps; tol = tol)
            end

            v1 = isfinite(r1) ? float(r1) : -Inf
            v2 = isfinite(r2) ? float(r2) : -Inf
            v3 = isfinite(r3) ? float(r3) : -Inf
            max_case = max(v1, v2, v3)
            ub = safe_num(upper_bound(eps))

            violation = (isfinite(max_case) && isfinite(ub)) ? (max_case - ub) : NaN
            violates = isfinite(violation) && violation > violation_tol

            if violates
                violation_count += 1
                max_violation = max(max_violation, violation)
            end

            row = [
                fmt_csv_value(eps),
                fmt_csv_value(r1),
                fmt_csv_value(r2),
                fmt_csv_value(r3),
                fmt_csv_value(max_case),
                fmt_csv_value(ub),
                fmt_csv_value(violation),
                violates ? "true" : "false",
            ]
            println(io, join(row, ","))

            if verbose
                @printf(
                    "[%d/%d] eps=%.6f max_case=%.6g upper=%.6g violation=%.6g\n",
                    idx,
                    length(eps_grid),
                    eps,
                    max_case,
                    ub,
                    violation,
                )
            end
        end
    end

    return (
        audit_csv = audit_csv,
        eps_count = length(eps_grid),
        violation_count = violation_count,
        max_violation = max_violation,
        violation_tol = float(violation_tol),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    opts = parse_args(ARGS)
    out_dir = get(opts, "output-dir", "results/sprint2_n1")
    eps_min = parse(Float64, get(opts, "eps-min", "0.0"))
    eps_max = parse(Float64, get(opts, "eps-max", "0.25"))
    eps_step = parse(Float64, get(opts, "eps-step", "0.005"))
    tol = parse(Float64, get(opts, "tol", "1e-7"))
    suppress_warnings = parse_bool(get(opts, "suppress-warnings", "true"))
    violation_tol = parse(Float64, get(opts, "violation-tol", "1e-7"))

    stats = run_upper_bound_audit(
        output_dir = out_dir,
        eps_min = eps_min,
        eps_max = eps_max,
        eps_step = eps_step,
        tol = tol,
        suppress_warnings = suppress_warnings,
        violation_tol = violation_tol,
        verbose = true,
    )

    @printf(
        "Upper-bound audit complete: rows=%d, violations=%d, max=%.6g\n",
        stats.eps_count,
        stats.violation_count,
        stats.max_violation,
    )
    println("Audit CSV: ", stats.audit_csv)
end

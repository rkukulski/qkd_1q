using Pkg
Pkg.activate(".")

using Logging
using Printf
using Statistics

include("../generate_three_case_curves.jl")

parse_bool(s::AbstractString) = lowercase(strip(s)) in ("1", "true", "yes", "y", "on")

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
            opts["input-points"] = arg
            i += 1
        end
    end
    return opts
end

function eval_with_logging(thunk::Function, suppress_warnings::Bool)
    if suppress_warnings
        return with_logger(ConsoleLogger(stderr, Logging.Error)) do
            thunk()
        end
    end
    return thunk()
end

function h2_local(a::Real)
    x = float(a)
    if x <= 1e-9 || x >= 1 - 1e-9
        return 0.0
    end
    return -x * log2(x) - (1 - x) * log2(1 - x)
end

function read_points(path::AbstractString)
    eps = Float64[]
    numeric_r = Float64[]
    opt_x = Float64[]

    open(path, "r") do io
        header = readline(io)
        expected = "epsilon,numeric_R,opt_x_numeric,scan_R,opt_x_scan,formula_R," *
                   "abs_diff_numeric_scan,rel_diff_numeric_scan," *
                   "abs_diff_numeric_formula,rel_diff_numeric_formula,best_seed"
        if strip(header) != expected
            error("Unexpected header in $path")
        end

        for line in eachline(io)
            parts = split(line, ",")
            length(parts) == 11 || continue
            push!(eps, parse(Float64, parts[1]))
            push!(numeric_r, parse(Float64, parts[2]))
            push!(opt_x, parse(Float64, parts[3]))
        end
    end

    return (eps = eps, numeric_r = numeric_r, opt_x = opt_x)
end

function write_diag_csv(path::AbstractString, rows)
    open(path, "w") do io
        println(
            io,
            "epsilon,numeric_R,opt_x_numeric,local_r_left,local_r_center," *
            "local_r_right,local_max_ok,active_x,p_active,pe_active," *
            "entropy_gap,recomposed_R,recompose_abs_err,p_at_x_minus," *
            "p_at_x_plus,active_sensitivity_ok,active_sensitivity_strict," *
            "lambda_dual,lambda_fd_proxy,lambda_positive,lambda_agreement_err",
        )
        for r in rows
            out = [
                fmt_csv_value(r.eps),
                fmt_csv_value(r.numeric_r),
                fmt_csv_value(r.opt_x),
                fmt_csv_value(r.r_left),
                fmt_csv_value(r.r_center),
                fmt_csv_value(r.r_right),
                r.local_max_ok ? "true" : "false",
                fmt_csv_value(r.active_x),
                fmt_csv_value(r.p_active),
                fmt_csv_value(r.pe_active),
                fmt_csv_value(r.entropy_gap),
                fmt_csv_value(r.recomposed_r),
                fmt_csv_value(r.recompose_abs_err),
                fmt_csv_value(r.p_minus),
                fmt_csv_value(r.p_plus),
                r.active_sensitivity_ok ? "true" : "false",
                r.active_sensitivity_strict ? "true" : "false",
                fmt_csv_value(r.lambda_dual),
                fmt_csv_value(r.lambda_fd_proxy),
                r.lambda_positive ? "true" : "false",
                fmt_csv_value(r.lambda_agreement_err),
            ]
            println(io, join(out, ","))
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

function run_n1_proof_diagnostics(;
    input_points::AbstractString = "results/sprint2_n1/n1_points.csv",
    output_dir::AbstractString = "results/sprint2_n1",
    local_u_step::Real = 0.01,
    x_sensitivity_step::Real = 0.001,
    x_sensitivity_tol::Real = 1e-8,
    x_sensitivity_strict_drop::Real = 1e-5,
    lambda_tol::Real = 1e-6,
    tol::Real = 1e-7,
    suppress_warnings::Bool = true,
)
    mkpath(output_dir)
    d = read_points(input_points)

    rows = NamedTuple[]
    positive_rows = 0
    local_ok_count = 0
    sensitivity_ok_count = 0
    sensitivity_strict_count = 0
    lambda_positive_count = 0
    first_zero_eps = NaN
    pre_zero_x = Float64[]
    recompose_errs = Float64[]
    p_plus_drop = Float64[]
    lambda_errors = Float64[]
    sensitivity_ok_interior = 0
    lambda_positive_interior = 0
    interior_count = 0
    lambda_errors_interior = Float64[]
    recompose_errs_interior = Float64[]

    for idx in eachindex(d.eps)
        eps = d.eps[idx]
        numeric_r = d.numeric_r[idx]
        x_opt = d.opt_x[idx]

        if isnan(first_zero_eps) && numeric_r <= 1e-12
            first_zero_eps = eps
        end
        if numeric_r > 1e-12
            push!(pre_zero_x, x_opt)
        end

        if numeric_r <= 1e-12 || !isfinite(x_opt) || x_opt <= 0
            push!(
                rows,
                (
                    eps = eps,
                    numeric_r = numeric_r,
                    opt_x = x_opt,
                    r_left = NaN,
                    r_center = NaN,
                    r_right = NaN,
                    local_max_ok = true,
                    active_x = NaN,
                    p_active = NaN,
                    pe_active = NaN,
                    entropy_gap = NaN,
                    recomposed_r = NaN,
                    recompose_abs_err = NaN,
                    p_minus = NaN,
                    p_plus = NaN,
                    active_sensitivity_ok = true,
                    active_sensitivity_strict = true,
                    lambda_dual = NaN,
                    lambda_fd_proxy = NaN,
                    lambda_positive = true,
                    lambda_agreement_err = NaN,
                ),
            )
            continue
        end

        u = log10(x_opt)
        u_left = clamp(u - local_u_step, UMIN, UMAX)
        u_right = clamp(u + local_u_step, UMIN, UMAX)

        r_left = eval_with_logging(suppress_warnings) do
            safe_case1_score(eps, u_left; tol = tol)[1]
        end
        r_right = eval_with_logging(suppress_warnings) do
            safe_case1_score(eps, u_right; tol = tol)[1]
        end

        qkd = build_case1_protocol(x_opt)
        r_center = numeric_r
        active_x = eval_with_logging(suppress_warnings) do
            PAB_eps(qkd, eps; tol = tol)
        end
        p_diag = eval_with_logging(suppress_warnings) do
            P_eps_diagnostics(qkd, active_x; tol = tol)
        end
        p_active = p_diag.p
        lambda_dual = p_diag.lambda_dual
        x_minus = clamp(active_x - x_sensitivity_step, 0.5 + 1e-9, 1 - 1e-9)
        x_plus = clamp(active_x + x_sensitivity_step, 0.5 + 1e-9, 1 - 1e-9)
        p_minus = eval_with_logging(suppress_warnings) do
            P_eps(qkd, x_minus; tol = tol)
        end
        p_plus = eval_with_logging(suppress_warnings) do
            P_eps(qkd, x_plus; tol = tol)
        end
        pe_active = eval_with_logging(suppress_warnings) do
            min_PE_eps(qkd, active_x; tol = tol)
        end
        entropy_gap = h2_local(pe_active) - h2_local(active_x)
        recomposed_r = p_active * max(entropy_gap, 0.0)
        recompose_abs_err = abs(recomposed_r - numeric_r)

        local_max_ok = r_center + 1e-9 >= max(r_left, r_right)
        active_sensitivity_ok =
            (p_plus + x_sensitivity_tol >= p_active) &&
            (p_minus <= p_active + x_sensitivity_tol)
        active_sensitivity_strict =
            (p_plus - p_active) > x_sensitivity_strict_drop
        denom = 2 * x_sensitivity_step * max(p_active, 1e-12)
        lambda_fd_proxy = (p_plus - p_minus) / denom
        lambda_positive = isfinite(lambda_dual) && lambda_dual > lambda_tol
        lambda_agreement_err = abs(lambda_dual - lambda_fd_proxy)
        positive_rows += 1
        local_ok_count += local_max_ok ? 1 : 0
        sensitivity_ok_count += active_sensitivity_ok ? 1 : 0
        sensitivity_strict_count += active_sensitivity_strict ? 1 : 0
        lambda_positive_count += lambda_positive ? 1 : 0
        push!(recompose_errs, recompose_abs_err)
        push!(p_plus_drop, p_plus - p_active)
        push!(lambda_errors, lambda_agreement_err)
        if eps > 1e-12
            interior_count += 1
            sensitivity_ok_interior += active_sensitivity_ok ? 1 : 0
            lambda_positive_interior += lambda_positive ? 1 : 0
            push!(lambda_errors_interior, lambda_agreement_err)
            push!(recompose_errs_interior, recompose_abs_err)
        end

        push!(
            rows,
            (
                eps = eps,
                numeric_r = numeric_r,
                opt_x = x_opt,
                r_left = r_left,
                r_center = r_center,
                r_right = r_right,
                local_max_ok = local_max_ok,
                active_x = active_x,
                p_active = p_active,
                pe_active = pe_active,
                entropy_gap = entropy_gap,
                recomposed_r = recomposed_r,
                recompose_abs_err = recompose_abs_err,
                p_minus = p_minus,
                p_plus = p_plus,
                active_sensitivity_ok = active_sensitivity_ok,
                active_sensitivity_strict = active_sensitivity_strict,
                lambda_dual = lambda_dual,
                lambda_fd_proxy = lambda_fd_proxy,
                lambda_positive = lambda_positive,
                lambda_agreement_err = lambda_agreement_err,
            ),
        )

        @printf(
            "[%d/%d] eps=%.6f local_ok=%s sens_ok=%s lambda=%.4g lambda_fd=%.4g p_rise=%.3g recompose_err=%.3g\n",
            idx,
            length(d.eps),
            eps,
            local_max_ok ? "true" : "false",
            active_sensitivity_ok ? "true" : "false",
            lambda_dual,
            lambda_fd_proxy,
            p_plus - p_active,
            recompose_abs_err,
        )
    end

    monotone_pre_zero = true
    for i in 2:length(pre_zero_x)
        if pre_zero_x[i] > pre_zero_x[i - 1] + 1e-12
            monotone_pre_zero = false
            break
        end
    end

    diag_csv = joinpath(output_dir, "n1_proof_diagnostics.csv")
    write_diag_csv(diag_csv, rows)

    max_err = isempty(recompose_errs) ? NaN : maximum(recompose_errs)
    mean_err = isempty(recompose_errs) ? NaN : mean(recompose_errs)
    ratio = positive_rows == 0 ? NaN : local_ok_count / positive_rows
    sensitivity_ratio = positive_rows == 0 ? NaN : sensitivity_ok_count / positive_rows
    strict_ratio =
        positive_rows == 0 ? NaN : sensitivity_strict_count / positive_rows
    min_drop = isempty(p_plus_drop) ? NaN : minimum(p_plus_drop)
    mean_drop = isempty(p_plus_drop) ? NaN : mean(p_plus_drop)
    lambda_ratio = positive_rows == 0 ? NaN : lambda_positive_count / positive_rows
    max_lambda_err = isempty(lambda_errors) ? NaN : maximum(lambda_errors)
    mean_lambda_err = isempty(lambda_errors) ? NaN : mean(lambda_errors)
    sensitivity_ratio_interior =
        interior_count == 0 ? NaN : sensitivity_ok_interior / interior_count
    lambda_ratio_interior =
        interior_count == 0 ? NaN : lambda_positive_interior / interior_count
    max_lambda_err_interior =
        isempty(lambda_errors_interior) ? NaN : maximum(lambda_errors_interior)
    mean_lambda_err_interior =
        isempty(lambda_errors_interior) ? NaN : mean(lambda_errors_interior)
    max_err_interior =
        isempty(recompose_errs_interior) ? NaN : maximum(recompose_errs_interior)
    mean_err_interior =
        isempty(recompose_errs_interior) ? NaN : mean(recompose_errs_interior)

    summary_csv = joinpath(output_dir, "n1_proof_diagnostics_summary.csv")
    metrics = Pair{String, String}[
        "input_points" => basename(input_points),
        "positive_rows" => string(positive_rows),
        "local_max_ok_count" => string(local_ok_count),
        "local_max_ok_ratio" => fmt_csv_value(ratio),
        "active_sensitivity_ok_count" => string(sensitivity_ok_count),
        "active_sensitivity_ok_ratio" => fmt_csv_value(sensitivity_ratio),
        "active_sensitivity_ok_ratio_eps_gt0" => fmt_csv_value(sensitivity_ratio_interior),
        "active_sensitivity_strict_count" => string(sensitivity_strict_count),
        "active_sensitivity_strict_ratio" => fmt_csv_value(strict_ratio),
        "min_p_plus_minus_p_active" => fmt_csv_value(min_drop),
        "mean_p_plus_minus_p_active" => fmt_csv_value(mean_drop),
        "lambda_positive_count" => string(lambda_positive_count),
        "lambda_positive_ratio" => fmt_csv_value(lambda_ratio),
        "lambda_positive_ratio_eps_gt0" => fmt_csv_value(lambda_ratio_interior),
        "max_lambda_agreement_err" => fmt_csv_value(max_lambda_err),
        "mean_lambda_agreement_err" => fmt_csv_value(mean_lambda_err),
        "max_lambda_agreement_err_eps_gt0" => fmt_csv_value(max_lambda_err_interior),
        "mean_lambda_agreement_err_eps_gt0" => fmt_csv_value(mean_lambda_err_interior),
        "first_zero_eps" => fmt_csv_value(first_zero_eps),
        "x_monotone_nonincreasing_pre_zero" => string(monotone_pre_zero),
        "max_recompose_abs_err" => fmt_csv_value(max_err),
        "mean_recompose_abs_err" => fmt_csv_value(mean_err),
        "max_recompose_abs_err_eps_gt0" => fmt_csv_value(max_err_interior),
        "mean_recompose_abs_err_eps_gt0" => fmt_csv_value(mean_err_interior),
        "positive_rows_eps_gt0" => string(interior_count),
        "local_u_step" => fmt_csv_value(local_u_step),
        "x_sensitivity_step" => fmt_csv_value(x_sensitivity_step),
        "x_sensitivity_tol" => fmt_csv_value(x_sensitivity_tol),
        "x_sensitivity_strict_drop" => fmt_csv_value(x_sensitivity_strict_drop),
        "lambda_tol" => fmt_csv_value(lambda_tol),
        "tol" => fmt_csv_value(tol),
    ]
    write_summary_csv(summary_csv, metrics)

    return (
        diagnostics_csv = diag_csv,
        summary_csv = summary_csv,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    opts = parse_args(ARGS)
    input_points = get(opts, "input-points", "results/sprint2_n1/n1_points.csv")
    output_dir = get(opts, "output-dir", "results/sprint2_n1")
    local_u_step = parse(Float64, get(opts, "local-u-step", "0.01"))
    x_sensitivity_step = parse(Float64, get(opts, "x-sensitivity-step", "0.001"))
    x_sensitivity_tol = parse(Float64, get(opts, "x-sensitivity-tol", "1e-8"))
    x_sensitivity_strict_drop = parse(
        Float64,
        get(opts, "x-sensitivity-strict-drop", "1e-5"),
    )
    lambda_tol = parse(Float64, get(opts, "lambda-tol", "1e-6"))
    tol = parse(Float64, get(opts, "tol", "1e-7"))
    suppress_warnings = parse_bool(get(opts, "suppress-warnings", "true"))

    out = run_n1_proof_diagnostics(
        input_points = input_points,
        output_dir = output_dir,
        local_u_step = local_u_step,
        x_sensitivity_step = x_sensitivity_step,
        x_sensitivity_tol = x_sensitivity_tol,
        x_sensitivity_strict_drop = x_sensitivity_strict_drop,
        lambda_tol = lambda_tol,
        tol = tol,
        suppress_warnings = suppress_warnings,
    )

    println("N=1 proof diagnostics complete:")
    println("  ", out.diagnostics_csv)
    println("  ", out.summary_csv)
end

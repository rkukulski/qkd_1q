using Pkg
Pkg.activate(".")

using Printf

include("struct.jl")
include("skr.jl")

const UMIN = -3.0
const UMAX = 1.0
const VMIN = -3.0
const VMAX = 1.0

function build_case1_protocol(x::Real)
    p_array = [1.0 1.0]
    psi_array = ComplexF64[
        1.0 1.0;
        0.0 x;;;
    ]
    q_array = [1.0]
    return QKDProtocol("case1", p_array, psi_array, q_array)
end

function build_case2_protocol(x::Real, y::Real)
    p_array = [
        1.0 1.0;
        x x
    ]
    psi_array = ComplexF64[
        1.0 0.0;
        0.0 1.0;;;
        1.0 1.0;
        1.0 -1.0
    ]
    q_array = [1.0, y]
    return QKDProtocol("case2", p_array, psi_array, q_array)
end

function safe_case1_score(eps::Real, u::Real; tol::Real = 1e-7)
    uu = clamp(float(u), UMIN, UMAX)
    x = 10.0^uu
    try
        qkd = build_case1_protocol(x)
        r = R_eps(qkd, eps; tol = tol)
        if isfinite(r)
            return r, uu, x
        end
    catch
    end
    return -Inf, uu, x
end

function safe_case2_score(eps::Real, u::Real, v::Real; tol::Real = 1e-7)
    uu = clamp(float(u), UMIN, UMAX)
    vv = clamp(float(v), VMIN, VMAX)
    x = 10.0^uu
    y = 10.0^vv
    try
        qkd = build_case2_protocol(x, y)
        r = R_eps(qkd, eps; tol = tol)
        if isfinite(r)
            return r, uu, vv, x, y
        end
    catch
    end
    return -Inf, uu, vv, x, y
end

function score_six_state(eps::Real; tol::Real = 1e-7)
    try
        r = R_eps(six_state, eps; tol = tol)
        return isfinite(r) ? r : NaN
    catch
        return NaN
    end
end

function pick_best_index(scores::Vector{Float64}, preferred_idx::Int)
    best_score = maximum(scores)
    ties = findall(s -> s == best_score, scores)
    return preferred_idx in ties ? preferred_idx : first(ties)
end

function global_scan_case1(eps::Real; tol::Real = 1e-7)
    candidates = collect(UMIN:0.5:UMAX)
    vals = [safe_case1_score(eps, u; tol = tol) for u in candidates]
    scores = [v[1] for v in vals]
    idx = pick_best_index(scores, Int(cld(length(candidates), 2)))
    return vals[idx]
end

function global_scan_case2(eps::Real; tol::Real = 1e-7)
    grid = collect(UMIN:1.0:UMAX)
    best = (-Inf, 0.0, 0.0, 1.0, 1.0)
    for u in grid
        for v in grid
            cand = safe_case2_score(eps, u, v; tol = tol)
            if cand[1] > best[1]
                best = cand
            elseif cand[1] == best[1]
                # Prefer less extreme parameters on ties.
                if abs(cand[2]) + abs(cand[3]) < abs(best[2]) + abs(best[3])
                    best = cand
                end
            end
        end
    end
    return best
end

function optimize_case1_at_epsilon(eps::Real, prev_u::Union{Nothing, Float64}; tol::Real = 1e-7)
    if prev_u === nothing
        return global_scan_case1(eps; tol = tol)
    end

    u = clamp(prev_u, UMIN, UMAX)
    delta = 0.4
    best_r, best_u, best_x = safe_case1_score(eps, u; tol = tol)

    for _ in 1:6
        u_left = clamp(u - delta, UMIN, UMAX)
        u_mid = clamp(u, UMIN, UMAX)
        u_right = clamp(u + delta, UMIN, UMAX)
        u_anchor = 0.0

        cands = unique([u_left, u_mid, u_right, u_anchor])
        vals = [safe_case1_score(eps, cand; tol = tol) for cand in cands]
        scores = [v[1] for v in vals]
        preferred_idx = something(findfirst(==(u_mid), cands), 1)

        idx = pick_best_index(scores, preferred_idx)
        chosen = vals[idx]
        best_r, best_u, best_x = chosen
        u = best_u

        if cands[idx] == u_mid
            delta /= 2
            if delta < 0.02
                break
            end
        end
    end

    # Rescue step for numerical plateaus and boundary lock-in.
    if best_r <= 1e-12
        global_best = global_scan_case1(eps; tol = tol)
        if global_best[1] > best_r
            return global_best
        end
    end

    return best_r, best_u, best_x
end

function optimize_case2_at_epsilon(
    eps::Real,
    prev_uv::Union{Nothing, Tuple{Float64, Float64}};
    tol::Real = 1e-7
)
    if prev_uv === nothing
        return global_scan_case2(eps; tol = tol)
    end

    u = clamp(prev_uv[1], UMIN, UMAX)
    v = clamp(prev_uv[2], VMIN, VMAX)
    delta = 0.5
    best = safe_case2_score(eps, u, v; tol = tol)

    # If the warm-start position is degenerate at this epsilon, use the global
    # scan result as the new starting point and continue with the warm-start loop.
    if best[1] <= 1e-6
        gs = global_scan_case2(eps; tol = tol)
        best = gs
        u = gs[2]
        v = gs[3]
    end

    for _ in 1:8
        center = (clamp(u, UMIN, UMAX), clamp(v, VMIN, VMAX))
        cands = [
            center,
            (clamp(u - delta, UMIN, UMAX), clamp(v, VMIN, VMAX)),
            (clamp(u + delta, UMIN, UMAX), clamp(v, VMIN, VMAX)),
            (clamp(u, UMIN, UMAX), clamp(v - delta, VMIN, VMAX)),
            (clamp(u, UMIN, UMAX), clamp(v + delta, VMIN, VMAX)),
            (clamp(u - delta, UMIN, UMAX), clamp(v - delta, VMIN, VMAX)),
            (clamp(u + delta, UMIN, UMAX), clamp(v + delta, VMIN, VMAX)),
            (0.0, 0.0),
            (0.0, clamp(v, VMIN, VMAX)),
            (clamp(u, UMIN, UMAX), 0.0)
        ]
        cands = unique(cands)

        vals = [safe_case2_score(eps, cand[1], cand[2]; tol = tol) for cand in cands]
        scores = [val[1] for val in vals]
        preferred_idx = something(findfirst(==(center), cands), 1)

        idx = pick_best_index(scores, preferred_idx)
        best = vals[idx]
        u = best[2]
        v = best[3]

        if cands[idx] == center
            delta /= 2
            if delta < 0.03
                break
            end
        end
    end

    if best[1] <= 1e-12
        global_best = global_scan_case2(eps; tol = tol)
        if global_best[1] > best[1]
            return global_best
        end
    end

    return best
end

function fmt_csv_value(x::Real)
    if isnan(x)
        return "NaN"
    elseif isinf(x)
        return x > 0 ? "Inf" : "-Inf"
    else
        return @sprintf("%.12g", float(x))
    end
end

function write_row(
    io::IO,
    case_name::AbstractString,
    eps::Real,
    r::Real,
    x_opt::Real,
    y_opt::Real
)
    println(
        io,
        string(
            case_name, ",",
            fmt_csv_value(eps), ",",
            fmt_csv_value(r), ",",
            fmt_csv_value(x_opt), ",",
            fmt_csv_value(y_opt)
        )
    )
end

function epsilon_values(eps_min::Real, eps_max::Real, eps_step::Real)
    if eps_step <= 0
        error("eps_step must be positive.")
    end
    if eps_min > eps_max
        error("eps_min must be <= eps_max.")
    end
    span = eps_max - eps_min
    n = Int(floor(span / eps_step + 1e-12))
    values = [eps_min + i * eps_step for i in 0:n]
    if values[end] < eps_max - 1e-12
        push!(values, float(eps_max))
    else
        values[end] = float(eps_max)
    end
    return values
end

function generate_three_case_curves(
    output_csv::AbstractString = "three_case_curves.csv";
    eps_step::Real = 0.0005,
    eps_min::Real = 0.0,
    eps_max::Real = 0.25,
    tol::Real = 1e-7
)
    eps_grid = epsilon_values(eps_min, eps_max, eps_step)

    println("Generating data for $(length(eps_grid)) epsilon values.")
    println("Output CSV: $output_csv")

    prev_u_case1 = nothing
    prev_uv_case2 = nothing

    open(output_csv, "w") do io
        println(io, "case,epsilon,R,x_opt,y_opt")
        for (idx, eps) in enumerate(eps_grid)
            r1, u1, x1 = optimize_case1_at_epsilon(eps, prev_u_case1; tol = tol)
            prev_u_case1 = u1
            row_r1 = isfinite(r1) ? r1 : NaN

            r2, u2, v2, x2, y2 = optimize_case2_at_epsilon(eps, prev_uv_case2; tol = tol)
            prev_uv_case2 = (u2, v2)
            row_r2 = isfinite(r2) ? r2 : NaN

            r3 = score_six_state(eps; tol = tol)

            write_row(io, "case1", eps, row_r1, x1, NaN)
            write_row(io, "case2", eps, row_r2, x2, y2)
            write_row(io, "case3", eps, r3, NaN, NaN)

            println(
                @sprintf(
                    "[%d/%d] eps=%.6f | case1 R=%.6g x=%.6g | case2 R=%.6g x=%.6g y=%.6g | case3 R=%.6g",
                    idx,
                    length(eps_grid),
                    eps,
                    row_r1,
                    x1,
                    row_r2,
                    x2,
                    y2,
                    r3
                )
            )
            flush(io)
        end
    end

    println("Done.")
    return output_csv
end

if abspath(PROGRAM_FILE) == @__FILE__
    output_csv = length(ARGS) >= 1 ? ARGS[1] : "three_case_curves.csv"
    eps_step = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 0.0005
    generate_three_case_curves(output_csv; eps_step = eps_step)
end

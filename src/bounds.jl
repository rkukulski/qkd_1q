function binary_entropy(a::Real)
    x = float(a)
    if x <= 1e-9 || x >= 1 - 1e-9
        return 0.0
    end
    return -x * log2(x) - (1 - x) * log2(1 - x)
end

"""
    upper_bound(eps)

Upper reference curve used in legacy figures for eps in [0, 0.25].
Returns `NaN` outside this interval.
"""
function upper_bound(eps::Real)
    e = float(eps)
    if e < 0 || e > 0.25
        return NaN
    end

    p = (4.0 / 3.0) * e
    c1 = sqrt(p)
    c0 = (sqrt(4.0 - 3.0 * p) - sqrt(p)) / 2.0
    term = 0.5 + 0.5 * (c1^2 / 2.0 + c1 * abs(c0 + c1 / 2.0))

    return binary_entropy(term) - binary_entropy(1.0 - 2.0 * e / 3.0)
end

function upper_bound_curve(eps_values)
    return [upper_bound(eps) for eps in eps_values]
end

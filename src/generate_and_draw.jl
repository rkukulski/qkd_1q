include("bounds.jl")

function upper_value(; delta = 0.001)
    interval = 0:delta:0.25
    return upper_bound_curve(interval)
end

function generate(qkd::QKDProtocol; delta = 0.001)
    interval = 0:delta:0.25 
    return map(eps -> R_eps(qkd, eps), interval)
end

function draw_eps_R(reses, filename; delta = 0.001)
    interval = 0:delta:0.25  
    
    Pl=plot(aspect_ratio=1/4)
    
    for i in eachindex(reses)
        plot!(Pl, interval, reses[i][1], color = palette(:tab10)[mod(i-1,10)+1], 
        label = "$(reses[i][2])")
    end

    xlims!(Pl, 0, 0.25)
    ylims!(Pl, 0, 1)
    xlabel!("\$ϵ\$", fontsize=20)
    ylabel!("\$R_ϵ\$", fontsize=20)
    savefig(Pl,"$(filename).pdf")
end

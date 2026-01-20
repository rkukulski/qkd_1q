function upper_value(; delta = 0.001)
    function h2(a)
        if a<=1e-9 || a>=1-1e-9
            return 0
        end
        return -a*log2(a) - (1-a)*log2(1-a)
    end
    interval = 0:delta:0.25
    p(eps) = 4/3*eps
    c1(p) = sqrt(p)
    c0(p) = (sqrt(4-3*p)-sqrt(p))/2
    return map(eps -> h2(1/2+1/2*(c1(p(eps))^2/2+c1(p(eps))*abs(c0(p(eps))+c1(p(eps))/2))) - h2(1-2/3*eps), interval)
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
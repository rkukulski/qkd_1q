function draw_eps_R(qkds::Vector{QKDProtocol}, filename; delta = 0.005)
    interval = 0:delta:1
    qkds_eps_R = [map(eps -> R_eps(qkd, eps), interval) for qkd in qkds]

    Pl=plot(aspect_ratio=1)
    for i in eachindex(qkds)
        plot!(Pl, interval, qkds_eps_R[i], color = palette(:tab10)[mod(i-1,10)+1], 
        label = "$(qkds[i].name)")
    end

    xlims!(Pl, 0, 1)
    ylims!(Pl, 0, 1)
    xlabel!("\$ϵ\$", fontsize=20)
    ylabel!("\$R_ϵ\$", fontsize=20)
    savefig(Pl,"$(filename).pdf")
end
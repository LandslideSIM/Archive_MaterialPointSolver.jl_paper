using CairoMakie
using MaterialPointSolver

let 
    figfont = MaterialPointSolver.fonttnr
    fig = Figure(size=(400, 280), fonts=(; regular=figfont, bold=figfont), fontsize=12)
    ax = Axis(fig[1, 1], xlabel=L"\lambda", ylabel=L"(r_\infty-r_0)/r_0")
    datax = [0, 0.5, 1, 2, 4]
    datay = [0, 0.7, 1.26, 2.1, 3.45]
    expx = [0, 0.28478964401294504, 0.7508090614886731, 0.8025889967637541, 
        1.2427184466019419, 1.631067961165048, 1.6828478964401297, 2.122977346278317, 
        2.2394822006472497, 2.498381877022654, 3.067961165048543, 3.4692556634304212, 
        4.0129449838187705]
    expy = [0, 0.4948453608247423, 0.8453608247422677, 0.9484536082474229, 
        1.3298969072164946, 1.5979381443298966, 1.7216494845360826, 2.1134020618556706, 
        2.2061855670103094, 2.3814432989690726, 2.6597938144329896, 2.9587628865979383, 
        3.247422680412371]
    refx = [0, 0.181229773, 0.530744337, 0.724919094, 0.957928803, 1.25566343, 1.462783172, 
        1.682847896, 1.928802589, 2.122977346, 2.57605178, 2.873786408, 3.1197411, 
        3.585760518, 3.702265372, 3.909385113]
    refy = [0, 0.18556701, 0.525773196, 0.670103093, 0.896907216, 1.226804124, 1.360824742, 
        1.536082474, 1.731958763, 1.896907216, 2.164948454, 2.360824742, 2.484536082, 
        2.75257732, 2.81443299, 2.917525773]
    p1 = scatterlines!(ax, datax, datay, markersize=8, color=:blue)
    p2 = scatter!(ax, expx, expy, markersize=8, marker=:utriangle, color=:red)
    p3 = scatterlines!(ax, refx, refy, markersize=0, color=(:green, 0.3), linewidth=10)

    axislegend(ax, [p1, p2, p3], ["uGIMP", "experiments", "reference"], position=:lt)
    limits!(ax, -0.2, 4.2, -0.2, 3.8)
    display(fig)
    save(joinpath(@__DIR__, "granularcompare.pdf"), fig)
end
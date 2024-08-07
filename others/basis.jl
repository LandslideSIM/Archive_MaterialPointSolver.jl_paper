using CairoMakie
using MaterialPointSolver

function uGIMP(Δx)
    lp = 0.5
    h = 1.0
    T2 = Float64
    T1 = Int64
    if abs(Δx) < T2(0.5)*lp # 0.25
        Ni = T2(1.0) - ((T2(4.0) * Δx * Δx + lp * lp) / (T2(4.0) * h * lp))
        dN = -((T2(8.0) * Δx) / (T2(4.0) * h * lp))
    elseif (T2(0.5) * lp) ≤ abs(Δx) < (h - T2(0.5) * lp) # 0.25 ~ 0.75
        Ni = T2(1.0) - (abs(Δx) / h)
        dN = sign(Δx) * (T2(-1.0) / h)
    elseif (h - T2(0.5) * lp) ≤ abs(Δx) < (h + T2(0.5) * lp) # 0.75 ~ 1.25
        Ni = ((h + T2(0.5) * lp - abs(Δx)) ^ T1(2)) / (T2(2.0) * h * lp)
        dN = -sign(Δx) * ((h + T2(0.5) * lp - abs(Δx)) / (h * lp))
    else
        Ni = T2(0.0)
        dN = T2(0.0)
    end
    return Ni, dN
end

function smpm(Δx)
    h = 1.0
    T2 = Float64
    if abs(Δx) ≤ h
        Ni = T2(1.0) - abs(Δx) / h
        dN = -sign(Δx) / h
    else
        Ni = T2(0.0)
        dN = T2(0.0)
    end
    return Ni, dN
end

let 
    figfont = MaterialPointSolver.fonttnr
    fig = Figure(size=(502, 180), fonts=(; regular=figfont, bold=figfont), fontsize=8, 
        figure_padding=0)

    a01 = fig[1, 1] = GridLayout()
    a02 = fig[1, 2] = GridLayout()

    ax1 = Axis(a01[1, 1], xlabel=L"x", ylabel=L"N_{i}", xgridvisible=false,
        ygridvisible=false, xticks=(-1.25:0.5:1.25), yticks=([-1.5, -0.5, 0, 1.0], 
        [L"0", L"1", L"0", L"1"]), aspect=1.5)
    ax2 = Axis(a02[1, 1], aspect=DataAspect())

    #=======================================================================================
    |  figure a                                                                            |
    =======================================================================================#
    figlinewidth = 1
    figmarkersize = 3
    figstrokewidth = 1

    ds = 1000
    xdots = range(-1.5, 1.5, length=ds) |> collect
    y = smpm.(xdots)
    y1 = [y[i, 1][1] for i in 1:ds]
    y2 = [y[i, 1][2] for i in 1:ds]
    p5 = scatter!(ax1, xdots, y1.-1.5, markersize=figmarkersize, color=:black)

    ds = 1000
    xdots = range(-1.5, 1.5, length=ds) |> collect
    y = uGIMP.(xdots)
    y1 = [y[i, 1][1] for i in 1:ds]
    y2 = [y[i, 1][2] for i in 1:ds]

    sta1 = findall(i->abs(xdots[i])<0.25, 1:ds)
    sta2 = findall(i->0.25≤abs(xdots[i])<0.75, 1:ds)
    sta3 = findall(i->0.75≤abs(xdots[i])<1.25, 1:ds)
    sta4 = findall(i->1.25≤abs(xdots[i]), 1:ds)

    p1 = scatter!(ax1, xdots[sta1], y1[sta1], markersize=figmarkersize, color=:red)
    p2 = scatter!(ax1, xdots[sta2], y1[sta2], markersize=figmarkersize, color=:blue)
    p3 = scatter!(ax1, xdots[sta3], y1[sta3], markersize=figmarkersize, color=:purple)
    p4 = scatter!(ax1, xdots[sta4], y1[sta4], markersize=figmarkersize, color=:green)

    text!(ax1, 0.85, 0.65, text="uGIMP")
    text!(ax1, 0.85,-0.85, text="sMPM")

    vlines!(ax1, [-1.25, -0.75, -0.25, 0.25, 0.75, 1.25], color=:gray, linestyle=:dash, 
        linewidth=0.5)
    hlines!(ax1, [1, -0.5, 0, -1.5], color=:gray, linestyle=:dash, linewidth=0.5)
    limits!(ax1, -1.6, 1.6, -1.7, 1.2)

    #=======================================================================================
    |  figure b                                                                            |
    =======================================================================================#
    figlinewidth = 1
    figmarkersize = 15
    figstrokewidth = 1

    lines!(ax2, [-1, -1], [-1.5, 2.5], color=:black, linewidth=figlinewidth)
    lines!(ax2, [ 0,  0], [-1.5, 2.5], color=:black, linewidth=figlinewidth)
    lines!(ax2, [ 1,  1], [-1.5, 2.5], color=:black, linewidth=figlinewidth)
    lines!(ax2, [ 2,  2], [-1.5, 2.5], color=:black, linewidth=figlinewidth)
    lines!(ax2, [-1.5, 2.5], [-1, -1], color=:black, linewidth=figlinewidth)
    lines!(ax2, [-1.5, 2.5], [ 0,  0], color=:black, linewidth=figlinewidth)
    lines!(ax2, [-1.5, 2.5], [ 1,  1], color=:black, linewidth=figlinewidth)
    lines!(ax2, [-1.5, 2.5], [ 2,  2], color=:black, linewidth=figlinewidth)

    Δ = 0.5
    p1 = poly!(ax2, Point2f[(-1.25+Δ, -1.25+Δ), (1.25+Δ, -1.25+Δ), (1.25+Δ, 1.25+Δ), 
                       (-1.25+Δ,  1.25+Δ)], 
        color=(:purple, 0.2), strokewidth=0, transparency=true, visible=true)
    p2 = scatter!(ax2, [-1 -1; -1 0; -1 1; -1 2; 0 -1; 0 2; 1 -1; 1 2; 2 -1; 2 0; 2 1; 2 2], 
        markersize=figmarkersize, marker=:xcross, color=(:gray, 0.5), strokecolor=:white, 
        strokewidth=figstrokewidth)
    p3 = scatter!(ax2, [0+Δ 0+Δ], markersize=figmarkersize, strokecolor=:white, 
        strokewidth=figstrokewidth, color=:purple)
    p4 = scatter!(ax2, [0 0; 0 1; 1 0; 1 1], 
        markersize=figmarkersize, marker=:xcross, color=(:green, 0.5), strokecolor=:white, 
        strokewidth=figstrokewidth)

    l1 = Legend(a02[1, 2], [p2, p4, p3, p1], [" inactive grid nodes", " active grid nodes", 
        " material particle", " particle effect range"], framevisible=false, rowgap=1, 
        labelsize=8)

    hidespines!(ax2)
    hidedecorations!(ax2)

    #=======================================================================================
    |  layout configuration                                                                |
    =======================================================================================#
    colgap!(a02, 1)
    colgap!(fig.layout, 50)
    colsize!(fig.layout, 1, Auto(0.8))
    Label(fig[1, 1, Bottom()], "(a)", padding=(5, 5, 0, 16), halign=:center)
    Label(fig[1, 2, Bottom()], "(b)", padding=(5, 5, 0, 16), halign=:center)

    display(fig)
    save(joinpath(@__DIR__, "outputs/basis.pdf"), fig, px_per_unit=1)
end

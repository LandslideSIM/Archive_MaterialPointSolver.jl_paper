#= post-processing video
@views let
    figfont = MaterialPointSolver.fontcmu
    figuretheme = Theme(
        Axis = (titlesize=30, titlefont=figfont, xticklabelfont=figfont, 
            yticklabelfont=figfont),
        Colorbar = (spinewidth=0, ticklabelfont=figfont, labelfont=figfont)
    )
    hdf5_path  = joinpath(args.project_path, "$(args.project_name).h5")
    fid        = h5open(hdf5_path, "r")
    video_step = 5
    video_fps  = 30
    itr        = (read(fid["FILE_NUM"])-1) |> Int64
    mp_pos     = Observable(fid["mp_init"] |> read)
    mp_num     = length(mp_pos[][:, 1])
    mp_Vsc     = Observable(zeros(mp_num))
    mp_σm      = Observable(zeros(mp_num))
    e1         = Observable(Point2f[(0.0, 0.0)])
    e2         = Observable(Point2f[(0.0, 0.0)])
    e3         = Observable(Point2f[(0.0, 0.0)])
    anno       = Observable(0.0)

    with_theme(figuretheme) do
        fig = Figure(size=(1080, 820), fontsize=25, font=figfont)
        Label(fig[0, 1:4], "Two elastic bodies collision", font=figfont, fontsize=30)
        ax1 = Axis(fig[1, 1], aspect=DataAspect(), xticks=(-0.5:0.2:0.5), 
            yticks=(-0.5:0.2:0.5))
        ax2 = Axis(fig[1, 3], aspect=DataAspect(), xticks=(-0.5:0.2:0.5), 
            yticks=(-0.5:0.2:0.5))
        ax3 = Axis(fig[2, 1:4], aspect=AxisAspect(4), xticks=(0:0.5:3.0), 
            yticks=(0.5:1:2.5))
        
        p1 = scatter!(ax1, mp_pos, color=mp_σm, colormap=:darktest, markersize=2, 
            colorrange=(-300, 100))
        Colorbar(fig[1, 2], p1, width=20, spinewidth=0, ticklabelfont=figfont, 
            label="Mean stress")
        p2 = scatter!(ax2, mp_pos, color=mp_Vsc, colormap=:darktest, markersize=2 ,
            colorrange=(0, 0.3))
        Colorbar(fig[1, 4], p2, width=20, spinewidth=0, ticklabelfont=figfont, 
            label="Σ Veloclty")

        p3 = scatterlines!(ax3, e1, markersize=0, linewidth=3, color=:blue )
        p4 = scatterlines!(ax3, e2, markersize=0, linewidth=3, color=:green)
        p5 = scatterlines!(ax3, e3, markersize=0, linewidth=3, color=:red  )
        axislegend(ax3, [p3, p4, p5], ["kinetic", "strain", "total"], "Energy", 
            labelsize=20, labelfont=figfont, titlefont=figfont)
        vlines!(ax3, anno  , color=:black, linewidth=1)
        vlines!(ax3, 1.5858, color=:red  , linewidth=2, linestyle=:dash)
        text!(ax3, 0.8, 1.3, text="analytical contact", font=figfont, fontsize=20)

        lines!(ax1, [0, 0], [-0.2, 0.2], color=:red, linewidth=3)
        lines!(ax1, [-0.2, 0.2], [0, 0], color=:red, linewidth=3)
        lines!(ax2, [0, 0], [-0.2, 0.2], color=:red, linewidth=3)
        lines!(ax2, [-0.2, 0.2], [0, 0], color=:red, linewidth=3)
        xlims!(ax1, -0.60, 0.6)
        ylims!(ax1, -0.60, 0.6)
        xlims!(ax2, -0.60, 0.6)
        ylims!(ax2, -0.60, 0.6)
        xlims!(ax3, -0.1, 3.5)
        ylims!(ax3, -0.1, 3.0)
        p = Progress(length(1:video_step:itr)-1; 
            desc="\e[1;36m[ Info:\e[0m $(lpad("video", 7))", color=:white, barlen=12, 
            barglyphs=BarGlyphs(" ◼◼  "))
        CairoMakie.record(fig, joinpath(@__DIR__, "outputs/video.mp4"),
            1:video_step:itr; framerate=video_fps) do i
            mp_σij = fid["group$(i)/sig"]   |> read
            mp_Ms  = fid["group$(i)/mass"]  |> read
            mp_ϵij = fid["group$(i)/eps_s"] |> read
            mp_Vol = fid["group$(i)/vol"]   |> read
            mp_Vs  = fid["group$(i)/v_s"]   |> read
            time   = fid["group$(i)/time"]  |> read
            Vsc  = sqrt.(mp_Vs[:, 1].^2 .+mp_Vs[:, 2].^2)
            σm   = (mp_σij[:, 1].+mp_σij[:, 2].+mp_σij[:, 3])./3
            tmp1 = Point2f(time, sum(sum(0.5.*mp_Vs.^2 .*mp_Ms, dims=2)))
            tmp2 = Point2f(time, sum(0.5.*mp_σij.*mp_ϵij.*mp_Vol))
            tmp3 = Point2f(time, tmp1[2]+tmp2[2])
            mp_pos[] = fid["group$(i)/mp_pos"] |> read
            mp_Vsc[] = Vsc
            mp_σm[]  = σm
            anno[]   = time
            push!(e1[], tmp1)
            push!(e2[], tmp2)
            push!(e3[], tmp3)
            next!(p)
        end
        close(fid)
    end
end
=#

# post-processing plots
@views let 
    hdf5_path = joinpath(args.project_path, "$(args.project_name).h5")
    fid       = h5open(hdf5_path, "r")
    itr       = (read(fid["FILE_NUM"])-1) |> Int64
    mp_pos    = fid["mp_init"] |> read
    data_step = 5
    x_time    = Float64[]
    e1        = Float64[]
    e2        = Float64[]
    e3        = Float64[]
    for i in 1:data_step:itr
        time   = fid["group$(i)/time"] |> read
        mp_Vs  = fid["group$(i)/v_s"]  |> read
        mp_σij = fid["group$(i)/sig"]  |> read
        mp_ϵij = fid["group$(i)/eps_s"]|> read
        mp_vol = fid["group$(i)/vol"]  |> read
        mp_Ms  = fid["group$(i)/mass"] |> read
        tmp1   = sum(sum(0.5.*mp_Vs.^2 .*mp_Ms, dims=2))
        tmp2   = sum(0.5.*mp_σij.*mp_ϵij.*mp_vol)
        tmp3   = tmp1+tmp2
        push!(e1, tmp1)
        push!(e2, tmp2)
        push!(e3, tmp3)
        push!(x_time, time)
    end
    close(fid)

    figfont=MaterialPointSolver.fonttnr
    fig  = Figure(size=(650, 300), fonts=(; regular=figfont, bold=figfont), fontsize=16)
    axis = Axis(fig[1, 1], xlabel="Time (s)", ylabel="Energy (J)", xticks=0:0.5:3.0, 
        yticks=0:0.5:3.0)
    p1 = scatterlines!(axis, x_time, e1, markersize=0, linewidth=3, color=:blue )
    p2 = scatterlines!(axis, x_time, e2, markersize=0, linewidth=3, color=:green)
    p3 = scatterlines!(axis, x_time, e3, markersize=0, linewidth=3, color=:red  )
    axislegend(axis, [p1, p2, p3], ["kinetic", "strain", "total"], "Energy", 
        labelsize=16)
    vlines!(axis, 1.5858, color=:red, linewidth=2, linestyle=:dash)
    text!(axis, 0.5, 1.3, text="analytical contact")
    limits!(axis, -0.1, 3.7, -0.1, 2.7)
    display(fig)
    save(joinpath(@__DIR__, "outputs/energy.pdf"), fig)
end

@views let 
    hdf5_path = joinpath(args.project_path, "$(args.project_name).h5")
    fid       = h5open(hdf5_path, "r")
    itr       = [602, 801, 1001, 1200] # 1.5s 2.0s 2.5s 3.0s
    mp_pos    = fid["mp_init"] |> read
    t1_pos    = fid["group602/mp_pos" ] |> read
    t2_pos    = fid["group801/mp_pos" ] |> read
    t3_pos    = fid["group1001/mp_pos"] |> read
    t4_pos    = fid["group1200/mp_pos"] |> read
    t1_σij    = fid["group602/sig"    ] |> read
    t2_σij    = fid["group801/sig"    ] |> read
    t3_σij    = fid["group1001/sig"   ] |> read
    t4_σij    = fid["group1200/sig"   ] |> read
    t1_Vs     = fid["group602/v_s"    ] |> read
    t2_Vs     = fid["group801/v_s"    ] |> read
    t3_Vs     = fid["group1001/v_s"   ] |> read
    t4_Vs     = fid["group1200/v_s"   ] |> read
    t1_σm     = (t1_σij[:, 1].+t1_σij[:, 2].+t1_σij[:, 3])./3
    t2_σm     = (t2_σij[:, 1].+t2_σij[:, 2].+t2_σij[:, 3])./3
    t3_σm     = (t3_σij[:, 1].+t3_σij[:, 2].+t3_σij[:, 3])./3
    t4_σm     = (t4_σij[:, 1].+t4_σij[:, 2].+t4_σij[:, 3])./3
    t1_Vst    = sqrt.(t1_Vs[:, 1].^2 .+t1_Vs[:, 2].^2)
    t2_Vst    = sqrt.(t2_Vs[:, 1].^2 .+t2_Vs[:, 2].^2)
    t3_Vst    = sqrt.(t3_Vs[:, 1].^2 .+t3_Vs[:, 2].^2)
    t4_Vst    = sqrt.(t4_Vs[:, 1].^2 .+t4_Vs[:, 2].^2)
    close(fid)

    figfont       = MaterialPointSolver.fonttnr
    figmarkersize = 1
    figlinewidth  = 1.5
    figcolormap   = :gist_rainbow # :darktest

    fig = Figure(size=(900, 410), fonts=(; regular=figfont, bold=figfont), fontsize=14)
    ax1 = Axis(fig[1, 1], aspect=DataAspect(), xticks=(-0.5:0.5:0.5), yticks=(-0.5:0.5:0.5))
    ax2 = Axis(fig[1, 2], aspect=DataAspect(), xticks=(-0.5:0.5:0.5), yticks=(-0.5:0.5:0.5))
    ax3 = Axis(fig[1, 3], aspect=DataAspect(), xticks=(-0.5:0.5:0.5), yticks=(-0.5:0.5:0.5))
    ax4 = Axis(fig[1, 4], aspect=DataAspect(), xticks=(-0.5:0.5:0.5), yticks=(-0.5:0.5:0.5))
    ax5 = Axis(fig[2, 1], aspect=DataAspect(), xticks=(-0.5:0.5:0.5), yticks=(-0.5:0.5:0.5))
    ax6 = Axis(fig[2, 2], aspect=DataAspect(), xticks=(-0.5:0.5:0.5), yticks=(-0.5:0.5:0.5))
    ax7 = Axis(fig[2, 3], aspect=DataAspect(), xticks=(-0.5:0.5:0.5), yticks=(-0.5:0.5:0.5))
    ax8 = Axis(fig[2, 4], aspect=DataAspect(), xticks=(-0.5:0.5:0.5), yticks=(-0.5:0.5:0.5))

    Label(fig[1, 1, Top()], "Time: 1.5s", padding=(0, 5, 5, 0), halign=:center)
    Label(fig[1, 2, Top()], "Time: 2.0s", padding=(0, 5, 5, 0), halign=:center)
    Label(fig[1, 3, Top()], "Time: 2.5s", padding=(0, 5, 5, 0), halign=:center)
    Label(fig[1, 4, Top()], "Time: 3.0s", padding=(0, 5, 5, 0), halign=:center)

    p1 = scatter!(ax1, t1_pos, color=t1_σm, colormap=figcolormap, markersize=figmarkersize, 
        colorrange=(-300, 100))
    p2 = scatter!(ax2, t2_pos, color=t2_σm, colormap=figcolormap, markersize=figmarkersize, 
        colorrange=(-300, 100))
    p3 = scatter!(ax3, t3_pos, color=t3_σm, colormap=figcolormap, markersize=figmarkersize,
        colorrange=(-300, 100))
    p4 = scatter!(ax4, t4_pos, color=t4_σm, colormap=figcolormap, markersize=figmarkersize,
        colorrange=(-300, 100))
    Colorbar(fig[1, 5], p1, width=10, spinewidth=0, ticks=-300:200:100, tickformat="{:.1e}",
        label="Mean stress (Pa)")

    p5 = scatter!(ax5, t1_pos, color=t1_Vst, colormap=figcolormap, markersize=figmarkersize,
        colorrange=(0, 0.3))
    p6 = scatter!(ax6, t2_pos, color=t2_Vst, colormap=figcolormap, markersize=figmarkersize,
        colorrange=(0, 0.3))
    p7 = scatter!(ax7, t3_pos, color=t3_Vst, colormap=figcolormap, markersize=figmarkersize,
        colorrange=(0, 0.3))
    p8 = scatter!(ax8, t4_pos, color=t4_Vst, colormap=figcolormap, markersize=figmarkersize,
        colorrange=(0, 0.3))
    Colorbar(fig[2, 5], p5, width=10, spinewidth=0, ticks=0:0.1:0.3, tickformat="{:.1e}",
        label="Veloclty magnitude (m/s)")

    lines!(ax1, [-0.1 0.1; 0 0], color=:red, linewidth=figlinewidth)
    lines!(ax1, [0 0; -0.1 0.1], color=:red, linewidth=figlinewidth)
    lines!(ax2, [-0.1 0.1; 0 0], color=:red, linewidth=figlinewidth)
    lines!(ax2, [0 0; -0.1 0.1], color=:red, linewidth=figlinewidth)
    lines!(ax3, [-0.1 0.1; 0 0], color=:red, linewidth=figlinewidth)
    lines!(ax3, [0 0; -0.1 0.1], color=:red, linewidth=figlinewidth)
    lines!(ax4, [-0.1 0.1; 0 0], color=:red, linewidth=figlinewidth)
    lines!(ax4, [0 0; -0.1 0.1], color=:red, linewidth=figlinewidth)
    lines!(ax5, [-0.1 0.1; 0 0], color=:red, linewidth=figlinewidth)
    lines!(ax5, [0 0; -0.1 0.1], color=:red, linewidth=figlinewidth)
    lines!(ax6, [-0.1 0.1; 0 0], color=:red, linewidth=figlinewidth)
    lines!(ax6, [0 0; -0.1 0.1], color=:red, linewidth=figlinewidth)
    lines!(ax7, [-0.1 0.1; 0 0], color=:red, linewidth=figlinewidth)
    lines!(ax7, [0 0; -0.1 0.1], color=:red, linewidth=figlinewidth)
    lines!(ax8, [-0.1 0.1; 0 0], color=:red, linewidth=figlinewidth)
    lines!(ax8, [0 0; -0.1 0.1], color=:red, linewidth=figlinewidth)
    
    limits!(ax1, -0.50, 0.5, -0.50, 0.5); limits!(ax2, -0.50, 0.5, -0.50, 0.5)
    limits!(ax3, -0.50, 0.5, -0.50, 0.5); limits!(ax4, -0.50, 0.5, -0.50, 0.5)
    limits!(ax5, -0.50, 0.5, -0.50, 0.5); limits!(ax6, -0.50, 0.5, -0.50, 0.5)
    limits!(ax7, -0.50, 0.5, -0.50, 0.5); limits!(ax8, -0.50, 0.5, -0.50, 0.5)
    hidespines!(ax1); hidespines!(ax2); hidespines!(ax3); hidespines!(ax4)
    hidespines!(ax5); hidespines!(ax6); hidespines!(ax7); hidespines!(ax8)

    display(fig)
    save(joinpath(@__DIR__, "outputs/vec_mean.pdf"), fig)
end

let
    figfont = MaterialPointSolver.fonttnr
    fig = Figure(size=(500, 500), fonts=(; regular=figfont, bold=figfont), fontsize=20)
    axis = Axis(fig[1, 1], aspect=DataAspect(), xticks=(-0.4:0.2:0.4), 
        yticks=(-0.4:0.2:0.4), xlabel=L"x\ (m)", ylabel=L"y\ (m)")
    poly!(axis, Circle(Point2f(-0.3, -0.3), 0.2), color="#cbeae0")
    poly!(axis, Circle(Point2f( 0.3,  0.3), 0.2), color="#ddd5fb")
    lines!(axis, [-0.3 -0.3+0.2*cos(-π/4); -0.3 -0.3+0.2*sin(-π/4)], color=:red, 
        linewidth=2)
    lines!(axis, [ 0.3  0.3+0.2*cos(3π/4);  0.3  0.3+0.2*sin(3π/4)], color=:red, 
        linewidth=2)
    lines!(axis, [-0.05 0.05; 0 0], color=:red, linewidth=3)
    lines!(axis, [0 0; -0.05 0.05], color=:red, linewidth=3)

    arrows!(axis, [0.3, -0.3], [0.3, -0.3], [-0.7, 0.7], [-0.7, 0.7], arrowsize=16, 
        lengthscale=0.3, linewidth=2, color=:blue)
    scatter!(axis, [0.3, -0.3], [0.3, -0.3], color=:black, markersize=10)
    text!(axis, -0.15, -0.50, text=L"r")
    text!(axis,  0.10,  0.45, text=L"r")
    text!(axis, -0.08, -0.2, text=L"v")
    text!(axis,  0.01,  0.15, text=L"-v")
    text!(axis,  0.34,  0.23, text=L"O_{2}")
    text!(axis, -0.41, -0.35, text=L"O_{1}")
    xlims!(axis, -0.6, 0.6)
    ylims!(axis, -0.6, 0.6)
    display(fig)
    save(joinpath(@__DIR__, "outputs/model.pdf"), fig)
end

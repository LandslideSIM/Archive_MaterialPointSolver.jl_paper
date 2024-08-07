# visualization
let 
    figfont = MaterialPointSolver.fontcmu
    fig = Figure(size=(500, 180), fonts=(; regular=figfont, bold=figfont), fontsize=16)
    axis = Axis(fig[1, 1], xlabel=L"x\ (mm)", ylabel=L"y\ (mm)", aspect=DataAspect(),
        xticks=(0:0.1:0.5, string.(0:100:500)), yticks=(0:0.05:0.1, string.(0:50:100)))
    p1 = scatter!(axis, mp.pos, color=log10.(mp.epII.+1), colorrange=(0.0, 1.4), 
        colormap=:turbo, markersize=3)
    Colorbar(fig[1, 2], p1, label=L"log_{10}(\epsilon_{II}+1)", size=10,
        ticks=0.1:0.4:1.4, ticklabelsize=16, spinewidth=0)
    rowsize!(fig.layout, 1, 102)
    limits!(axis, -0.012, 0.46, -0.012, 0.14)
    display(fig)
    save(joinpath(@__DIR__, "outputs/runout.pdf"), fig)
end

let 
    figfont = MaterialPointSolver.fontcmu
    fig = Figure(size=(500, 180), fonts=(; regular=figfont, bold=figfont), fontsize=16)
    assets = joinpath(assetsdir, "data")
    failure = readdlm(joinpath(assets, "2d_collapse_experiments/failure.csv"), ',', Float64)
    surface = readdlm(joinpath(assets, "2d_collapse_experiments/surface.csv"), ',', Float64)
    ux = mp.pos[:, 1].-mp.init[:, 1]
    u1 = findall(i->ux[i]>0.001, 1:mp.num)
    u2 = findall(i->ux[i]≤0.001, 1:mp.num)
    ux[u1] .= -1; ux[u2] .= 1
    colors = ["#ef0000", "#00008e"]
    cmap = cgrad(colors, 2; categorical=true)
    axis = Axis(fig[1, 1], xlabel=L"x\ (mm)", ylabel=L"y\ (mm)", aspect=DataAspect(),
        xticks=(0:0.1:0.5, string.(0:100:500)), yticks=(0:0.05:0.1, string.(0:50:100)))
    p1 = scatter!(axis, mp.pos, color=ux, colormap=cmap, marker=:circle, markersize=3)
    scatterlines!(axis, failure[:, 1], failure[:, 2], color=:green, marker=:utriangle,
        markersize=10, strokewidth=0, linewidth=2)
    scatterlines!(axis, surface[:, 1], surface[:, 2], color=:green, marker=:diamond,
        markersize=10, strokewidth=0, linewidth=2)
    limits!(axis, -0.012, 0.46, -0.012, 0.14)
    Colorbar(fig[1, 2], p1, label="\n ", size=10, ticklabelsize=16, spinewidth=0,
        ticks=(-1:1:1, [L">1mm", L"\Delta u", L"≤1mm"]))
    rowsize!(fig.layout, 1, 102)
    display(fig)
    save(joinpath(@__DIR__, "outputs/comparison.pdf"), fig)
end

let 
    figfont = MaterialPointSolver.fontcmu
    fig = Figure(size=(400, 180), fonts=(; regular=figfont, bold=figfont), fontsize=16)

    mp_line = Int64[]
    for i in 1:mp.num
        for j in collect(0:0.02:0.98)
            isapprox(mp.init[i, 2], j, atol=mp.space_y/1.9) ? push!(mp_line, i) : nothing
        end
        for k in collect(0:0.02:0.2)
            isapprox(mp.init[i, 1], k, atol=mp.space_x/1.9) ? push!(mp_line, i) : nothing
        end
        unique!(mp_line)
    end
    colors = ["#F9807d", "#00C1C8"]
    cmap = cgrad(colors, 2; categorical=true)
    axis = Axis(fig[1, 1], xlabel=L"x\ (mm)", ylabel=L"y\ (mm)", aspect=DataAspect(),
        xticks=(0:0.1:0.5, string.(0:100:500)), yticks=(0:0.05:0.1, string.(0:50:100)))
    scatter!(axis, mp.pos, color=:gray, marker=:circle, markersize=3)
    scatter!(axis, mp.pos[mp_line, :], color=:red, marker=:circle, markersize=2)
    limits!(axis, -0.012, 0.46, -0.012, 0.14)
    display(fig)
    save(joinpath(@__DIR__, "outputs/refline.pdf"), fig)
end

let 
    figfont = MaterialPointSolver.fontcmu
    fig = Figure(size=(400, 180), fonts=(; regular=figfont, bold=figfont), fontsize=16)

    axis = Axis(fig[1, 1], xlabel=L"x\ (mm)", ylabel=L"y\ (mm)", aspect=DataAspect(),
        xticks=(0:0.1:0.5, string.(0:100:500)), yticks=(0:0.05:0.1, string.(0:50:100)), 
        xgridvisible=false, ygridvisible=false)

    poly!(axis, Point2f[(0, 0), (0, 0.1), (0.2, 0.1), (0.2, 0)], color=:gray)
    poly!(axis, Point2f[(0, -1e4), (-1e4, -1e4), (-1e4, -1e4), (0, 1e4)], color=:black)
    poly!(axis, Point2f[(0, 0), (0, -1e4), (1e4, -1e4), (1e4, 0)], color=:black)
    
    text!(axis, 0.09, 0.105, text=L"l")
    text!(axis, 0.21, 0.04, text=L"h")

    limits!(axis, -0.012, 0.46, -0.012, 0.14)
    resize_to_layout!(fig)
    display(fig)
    save(joinpath(@__DIR__, "outputs/model.pdf"), fig)
end

@inbounds let
    # prepare data =========================================================================
    area = findall(i->mp.init[i, 1]≥0.15&&mp.init[i, 2]≥0.05, 1:mp.num)
    coord = mp.init[area, :]
    x_num = unique(coord[:, 1]) |> length
    y_num = unique(coord[:, 2]) |> length
    sorted_coord = copy(coord)
    for i in 2:2:x_num
        c_local = (i-1)*y_num
        for j in 1:y_num
            sorted_coord[c_local+j, :] .= coord[c_local+y_num+1-j, :]
        end
    end
    datalength = size(sorted_coord, 1)
    idx = Int[]
    for i in 1:datalength
        for j in 1:mp.num
            if sorted_coord[i, :] == mp.init[j, :]
                push!(idx, j)
            end
        end
    end
    fid = h5open(joinpath(args.project_path, "$(args.project_name).h5"), "r")
    itr = (read(fid["FILE_NUM"])-1) |> Int64
    data_x = zeros(datalength*itr, 3)
    data_y = zeros(datalength*itr, 3)
    Threads.@threads for i in 1:itr
        c_time = fid["group$(i)/time"] |> read
        c_v_s  = fid["group$(i)/v_s" ] |> read
        c_local = (i-1)*datalength+1
        data_x[c_local:c_local+datalength-1, 1] .= collect(1:1:datalength)
        data_x[c_local:c_local+datalength-1, 2] .= c_time
        data_x[c_local:c_local+datalength-1, 3] .= c_v_s[idx, 1]
        data_y[c_local:c_local+datalength-1, 1] .= collect(1:1:datalength)
        data_y[c_local:c_local+datalength-1, 2] .= c_time
        data_y[c_local:c_local+datalength-1, 3] .= c_v_s[idx, 2]
    end
    close(fid)
    # plots ================================================================================
    figfont = MaterialPointSolver.fonttnr
    fig1 = Figure(size=(750, 400), fonts=(; regular=figfont, bold=figfont), fontsize=22)
    fig2 = Figure(size=(750, 400), fonts=(; regular=figfont, bold=figfont), fontsize=22)
    
    ax1 = Axis(fig1[1, 1], xlabel="Particle ID", ylabel="Time (s)", #title=L"x\ vvelocity",
        aspect=1.8)
    ax2 = Axis(fig2[1, 1], xlabel="Particle ID", ylabel="Time (s)", #title=L"y\ vvelocity",
        aspect=1.8)
    pl1 = scatter!(ax1, data_x[:, 1], data_x[:, 2], color=data_x[:, 3], markersize=3, 
        colormap=:rainbow1, colorrange=(0, 1.0))
    pl2 = scatter!(ax2, data_y[:, 1], data_y[:, 2], color=data_y[:, 3], markersize=3, 
        colormap=:rainbow1, colorrange=(-0.5, 0))
    hidespines!(ax1)
    hidespines!(ax2)
    Colorbar(fig1[1, 2], pl1, spinewidth=0, label=L"V_{x}\ \text{(m/s)}", height=Relative(1/1))
    Colorbar(fig2[1, 2], pl2, spinewidth=0, label=L"V_{y}\ \text{(m/s)}", height=Relative(1/1))
    limits!(ax1, 0, 1600, 0, 1)
    limits!(ax2, 0, 1600, 0, 1)
    save(joinpath(@__DIR__, "outputs/velocity_x.pdf"), fig1)
    save(joinpath(@__DIR__, "outputs/velocity_y.pdf"), fig2)
    display(fig1)
    display(fig2)
end
using Base.Threads, Plots, LazyGrids, Interpolations

include("functions_detector.jl")

default(aspect_ratio=:auto)

@time begin

    xStart = -1e-4
    yStart = -1e-4
    xEnd = -xStart
    yEnd = -yStart
    count = Int(1e2)
    tCount = Int(1e2)

    xData = range(xStart, xEnd, length=count)
    yData = range(yStart, yEnd, length=count)

    dtSensor = 1.3e-14
    t = range(start=0, step=dtSensor, length=tCount)
    c = 3e8
    e0 = 8.8541878128e-12
    mu0 = 1.25663706212e-6

    sensor = createGridArray(xStart, xEnd, count, yStart, yEnd, count)

    trajectory = createTrajectory(xStart, xEnd, tCount)

    vecData = initDistanceArray(sensor, trajectory)

    @threads for i in 1:tCount
        vecData[i, :, :, :] = calculateEField(trajectory[:, i], sensor, vecData[i, :, :, :], t[i])
    end

    xUnsynced = vecData[:, :, :, 1]
    yUnsynced = vecData[:, :, :, 2]
    t_data = vecData[:, :, :, 3]

    xSynced = Array{Float64}(undef, size(xUnsynced)...)
    ySynced = Array{Float64}(undef, size(yUnsynced)...)

    @threads for i in 1:count
        for j in 1:count
            itpx = linear_interpolation(t_data[:, i, j], xUnsynced[:, i, j], extrapolation_bc=Line())
            xSynced[:, i, j] .= itpx(t)
            itpy = linear_interpolation(t_data[:, i, j], yUnsynced[:, i, j], extrapolation_bc=Line())
            ySynced[:, i, j] .= itpy(t)
        end
    end

    distSynced = (xSynced .^ 2 + ySynced .^ 2) .^ 0.5

    #= 
        retard = true

        if retard == true
            p = contourf(xData, yData, permutedims(interp_dist_data[1, :, :], [2, 1]), linewidth=0, aspect_ratio=:equal, colormap=:jet, levels=100,) #=ylims=[yStart,yEnd]./5=#

            for i in 1:count
                p.series_list[1].plotattributes[:z] = permutedims(interp_dist_data[i, :, :], [1, 2])
                display(p)
                sleep(1 / 10)
            end

        else
            p = contourf(xData, yData, permutedims(dist_data[1, :, :], [2, 1]), linewidth=0, aspect_ratio=:equal, colormap=:jet, levels=100,) #=ylims=[yStart,yEnd]./5=#

            for i in 1:count
                p.series_list[1].plotattributes[:z] = permutedims(dist_data[i, :, :], [1, 2])
                display(p)
                sleep(1 / 10)
            end

        end
     =#
end
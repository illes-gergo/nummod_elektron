using Base.Threads, LazyGrids, Interpolations, PlotlyJS

include("functions_detector.jl")


@time begin

    xStart = -2e-4
    yStart = -2e-4
    xEnd = -xStart
    yEnd = -yStart
    count = Int(1e2)
    tCount = Int(5e2)

    xData = range(xStart, xEnd, length=count)
    yData = range(yStart, yEnd, length=count)

    dtSensor = 1e-14
    t = range(start=0, step=dtSensor, length=tCount)
    c = 3e8
    e0 = 8.8541878128e-12
    mu0 = 1.25663706212e-6
    qe = -1.60217663e-19

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


    EFieldX, EFieldY = calculateRadiationField(xSynced, ySynced, distSynced, dtSensor)

    trace1 = scatter(x=t, y=EFieldX[:,50,1], name="X")
    trace2 = scatter(x=t, y=EFieldY[:,50,1], name="Y")

    plot([trace1,trace2])

end
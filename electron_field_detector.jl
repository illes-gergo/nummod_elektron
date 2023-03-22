using Base.Threads, Plots, LazyGrids, Interpolations

include("functions_detector.jl")

default(aspect_ratio=:auto)

xStart = -1e-4;
yStart = -1e-4;
xEnd = -xStart;
yEnd = -yStart;
count = Int(1e2);

xData = range(xStart, xEnd, length=count)
yData = range(yStart, yEnd, length=count)

dtSensor = 1.3e-14;
t = range(start=0, step=dtSensor, length=count)
c = 3e8;

sensor = createGridArray(xStart, xEnd, count, yStart, yEnd, count)

trajectory = createTrajectory(xStart, xEnd, count)

data = initDataArray(sensor, trajectory)

@threads for i in 1:count
    data[i, :, :, :] .= calculateEField(trajectory[:, i], sensor, data[i, :, :, :], t[i])
end

t_data = data[:, :, :, 2];
dist_data = data[:, :, :, 1]
interp_dist_data = copy(dist_data)

for i in 1:count
    for j in 1:count
        itp = linear_interpolation(t_data[:, i, j], dist_data[:, i, j], extrapolation_bc=Line())
        interp_dist_data[:, i, j] .= itp(t)
    end
end

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
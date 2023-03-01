using Base.Threads, Plots

include("functions_detector.jl")

xStart = -1e-4;
yStart = -2e-4;
xEnd = -xStart;
yEnd = -yStart;
count = Int(1e3);

dtSensor = 1e-9;
t = range(start=0, step=dtSensor, length=count)
c = 3e8;

sensor = createGridArray(xStart, xEnd, count, yStart, yEnd, count)

trajectory = createTrajectory(xStart, xEnd, count)

data = initDataArray(sensor, trajectory)

@threads for i in 1:count
    calculateEField(trajectory[:, i], sensor, data[i, :, :, :], t[i])
end
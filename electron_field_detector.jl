using Base.Threads

include("functions_detector.jl")

xStart = -1e-6;
yStart = -2e-6;
xEnd = -xStart;
yEnd = -yStart;
count = 1e3;

dtSensor = 1e-9;
c = 3e8;

sensor = createGridArray(xStart,xEnd,count,yStart,yEnd,count)



trajectory = createTrajectory(xStart,xEnd,count)


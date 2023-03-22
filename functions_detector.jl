function createGridArray(xStart, xEnd, xCount, yStart, yEnd, yCount)
    xGrid, yGrid = ndgrid(range(xStart, xEnd, length=xCount), range(yStart, yEnd, length=yCount))
    output = cat(xGrid, yGrid, dims=3)
    return output
end

function createTrajectory(xStart, xEnd, count)
    xVals = range(xStart, xEnd, Int(count))'
    return cat(xVals, zeros(size(xVals)), dims=1)
end

function initDataArray(sensor, trajectory)
    return Array{Float64}(undef, (size(trajectory, 2), size(sensor, 1), size(sensor, 2), 2))
end

function calculateEField(electron_position, sensor, data, t)
    c0 = 3e8
    XDiff = electron_position[1] .- sensor[:, :, 1]
    YDiff = electron_position[2] .- sensor[:, :, 2]
    distances = (XDiff .^ 2 .+ YDiff .^ 2) .^ 0.5
    data[:, :, 1] .= distances .^ 0.5
    data[:, :, 2] .= t .+ distances ./ c0
    return data
end
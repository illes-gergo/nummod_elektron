function createGridArray(xStart, xEnd, xCount, yStart, yEnd, yCount)
    xGrid, yGrid = ndgrid(range(xStart, xEnd, length=xCount), range(yStart, yEnd, length=yCount))
    output = cat(xGrid, yGrid, dims=3)
    return output
end

function createTrajectory(xStart, xEnd, count)
    xVals = range(xStart, xEnd, Int(count))'
    return cat(xVals, trajectory_y(xVals), dims=1)
end

function initDistanceArray(sensor, trajectory)
    return Array{Float64}(undef, (size(trajectory, 2), size(sensor, 1), size(sensor, 2), 3))
end

function calculateEField(electron_position, sensor, vecData, t)
    XDiff = electron_position[1] .- sensor[:, :, 1]
    YDiff = electron_position[2] .- sensor[:, :, 2]
    vecData[:, :, 1] .= XDiff
    vecData[:, :, 2] .= YDiff
    vecData[:, :, 3] = t .+ (XDiff .^ 2 .+ YDiff .^ 2) .^ 0.5 ./ c
    return vecData
end

function trajectory_y(xVals)
    return sin.(2 .* pi .* xVals .* 1 ./ 1e-4) * 1e-5
end

function calculateRadiationField(xSynced::Array, ySynced::Array, distSynced::Array, dt::Real)
  

end
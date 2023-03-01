function createGridArray(xStart, xEnd, xCount, yStart, yEnd, yCount)
    output = Array{Float64}(undef, (yCount, xCount, 2))
    xVals = range(xStart, xEnd, Int(xCount))
    yVals = range(yStart, yEnd, Int(yCount))

    @threads for j in eachindex(xVals)
        for i in eachindex(yVals)
            output[i, j, 2] = yVals[i]
            output[i, j, 1] = xVals[j]
        end
    end
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
    data[:, :, 1] .= distances
    data[:, :, 2] .= t .+ distances ./ c0
end
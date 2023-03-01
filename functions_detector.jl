function createGridArray(xStart, xEnd, xCount, yStart, yEnd, yCount)
    output = Array{Float64}(undef, (Int(yCount), Int(xCount), 2))
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
    xVals = range(xStart,xEnd,Int(count))'
    return cat(xVals,zeros(size(xVals)),dims=1)
end
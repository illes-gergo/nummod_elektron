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
    vX = myDiffTime(xSynced) / dt
    vY = myDiffTime(ySynced) / dt
    aX = myDiffTime(vX) / dt
    aY = myDiffTime(vY) / dt
    absV = (vX .^ 2 + vY .^ 2) .^ 0.5

    term11 = qe / 4 / pi / e0 * (1 .- (absV / c) .^ 2)
    term12xnum = xSynced .- distSynced .* vX / c
    term12ynum = ySynced .- distSynced .* vY / c
    termnom = (distSynced - myDotProduct(xSynced, ySynced, vX / c, vY / c)) .^ 3
    term21 = qe * mu0 / 4
    term2in = myCrossProduct(term12xnum, term12ynum, aX, aY)
    term2outX, term2outY = myCrossProduct(xSynced, ySynced, term2in)
    outX = @. (term11 * term12xnum + term21 * term2outX) / termnom
    outY = @. (term11 * term12ynum + term21 * term2outY) / termnom
    return outX, outY
end

function myDiffTime(timeData::Array)
    return cat(zeros(1, size(timeData)[2:end]...), diff(timeData, dims=1), dims=1)
end

function myDotProduct(X::Array, Y::Array, vX::Array, vY::Array)
    X .* vX .+ Y .* vY
end

function myCrossProduct(Ax::Array, Ay::Array, Bx::Array, By::Array)
    return myLength(Ax, Ay) .* myLength(Bx, By) .* sin.(myAngle(Ax, Ay, Bx, By))
end

function myCrossProduct(Ax::Array, Ay::Array, Bz::Array)
    outLen = myLength(Ax, Ay) .* Bz
    outAng = mySingleAngle(Ay, Ay) .- pi ./ 2
    outCPLX = @. outLen * exp(1im * outAng)
    outX = @. real(outCPLX)
    outY = @. imag(outCPLX)
end

function myAngle(Ax::Array, Ay::Array, Bx::Array, By::Array)
    aCPLX = @. Ax + 1im * Ay
    angA = @. angle(aCPLX)
    bCPLX = @. Bx + 1im * By
    angB = @. angle(bCPLX)
    return angB .- angA
end

function myLength(Ax::Array, Ay::Array)
    return @. (Ax^2 + Ay^2)^0.5
end

function mySingleAngle(Ax::Array, Ay::Array)
    vecCPLX = @. Ax + 1im * Ay
    return angle.(vecCPLX)
end
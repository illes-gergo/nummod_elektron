function createParticlesGauss(n, sigma1, sigma2, velocity)
    electrons = zeros(2, 3, n)
    electrons[1, 1, :] = randn(size(electrons[1, 1, :])) * sigma1
    electrons[2, 1, :] = randn(size(electrons[1, 1, :])) * sigma2
    electrons[1, 2, :] .= velocity
    #electrons[2, 2, :] .= -(electrons[2, 1, :]) .* velocity .* 1e5
    return electrons
end

function coulombForce(electrons)
    k = 8.9875517923e9
    q = 1.60217663e-19
    DistVecs = pairwise(Euclidean(), electrons[:, 1, :], dims=2)
    XComp = @spawn begin
        #SubVecsX = pairwise(-, electrons[1, 1, :])
        SubVecsX = pairwiseSubtract(electrons[1, 1, :])
        ForceX = -k .* q .^ 2 ./ DistVecs .^ 3 .* SubVecsX
        ForceX[I(size(electrons, 3))] .= 0
        sForceX = sum(ForceX, dims=1)
    end

    YComp = @spawn begin
        #SubVecsY = pairwise(-, electrons[2, 1, :])
        SubVecsY = pairwiseSubtract(electrons[2, 1, :])
        ForceY = -k .* q .^ 2 ./ DistVecs .^ 3 .* SubVecsY
        ForceY[I(size(electrons, 3))] .= 0
        sForceY = sum(ForceY, dims=1)
    end

    return cat(fetch(XComp), fetch(YComp), dims=1)
end

function lorentzForce(electrons)
    q = 1.60217663e-19
    magnetic = magneticField(electrons[1, 1, :], electrons[2, 1, :])
    ForceX = -1.0 .* electrons[2, 2, :] .* magnetic .* q
    ForceY = electrons[1, 2, :] .* magnetic .* q
    return cat(transpose(ForceX), transpose(ForceY), dims=1)
end

function magneticField(x, y)
    fieldStrength = 0.2
    period = 10e-6
    MagField = cos.(2 * pi * 1 ./ period .* x) .* fieldStrength
    return MagField
end

function stepInTime(electrons, dt, Forces)
    m_electron = 9.1093837e-31
    m_e = rel_mass(m_electron, electrons[:, 2, :])
    electrons[:, 3, :] .= Forces ./ m_e
    electrons[:, 2, :] .+= electrons[:, 3, :] .* dt
    electrons[:, 1, :] .+= electrons[:, 2, :] .* dt
    return electrons
end

function rel_mass(m, v)
    c = 3e8
    m_rel = m ./ sqrt.(1 .- normPair(v) .^ 2 ./ c .^ 2)
    return cat(m_rel, m_rel, dims=1)
end

function pairwiseSubtract(a)
    output = Array{Float64}(undef, length(a), length(a))
	@threads for i in eachindex(a)
	for j in eachindex(a)
        output[i, j] = a[i] - a[j]
end
    end
    return output
end
function normPair(v)
	res = Array{Float64}(undef,1,size(v,2))
	@threads for i = axes(v,2)
		res[i] = norm(v[:,i])
	end
	return res;
end

 

function createParticlesGauss(n, sigma1, sigma2, velocity)
    electrons = zeros(2, 3, n)
    electrons[1, 1, :] = randn(size(electrons[1, 1, :])) * sigma1
    electrons[2, 1, :] = randn(size(electrons[1, 1, :])) * sigma2
    electrons[1, 2, :] .= velocity
    return electrons
end

function coulombForce(electrons)
    k = 8.9875517923e9
    q = 1.60217663e-19
    @. SubVecs = pairwise(-, eachcol(electrons[:, 1, :]))
    @. DistVecs = pairwise(Euclidean(), electrons[:, 1, :], dims=2)
    @. Force = k .* q .^ 2 ./ DistVecs .^ 3 .* SubVecs
    return Force
end


function stepInTime(electrons, t1, t2)
    dt = t2 - t1
    electrons[:, 2, :] = electrons[:, 3, :] * dt
    electrons[:, 1, :] = electrons[:, 2, :] * dt
    electrons[:, 3, :] .= 0
    return electrons
end

function rel_mass(m, v)
    m_rel = m / sqrt(1 - norm(v)^2 / c^2)
    return m_rel
end
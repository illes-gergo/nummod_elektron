function createParticlesGauss(n, sigma1, sigma2, velocity)
    electrons = zeros(2, 3, n)
    electrons[1, 1, :] = randn(size(electrons[1, 1, :])) * sigma1
    electrons[2, 1, :] = randn(size(electrons[1, 1, :])) * sigma2
    electrons[1, 2, :] .= velocity
    return electrons
end

function coulombForce(electrons)
    @spawn for i = 1:axes(electrons,3)
        R_array = zeros(size(electrons,3));
        R_array[i] = norm(electrons[:,1,i]);
    end
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
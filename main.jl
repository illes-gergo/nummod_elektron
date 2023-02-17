using Plots, LinearAlgebra, Random, Base.Threads, StatsBase, Distances

include("functions.jl");

m_e = 9.1093837e-31;
q_e = 1.60217663e-13;

electrons = createParticlesGauss(200, 1.3e-6, 3e-6, 3e5);


t = 0;
dt = 1e-13;
t_end = 6 * 3e2 * dt;
let electrons_looped = electrons, t = t, K = 0
    global anim = @animate while t < t_end
        stepInTime(electrons_looped, dt, coulombForce(electrons_looped) .+ lorentzForce(electrons_looped))
        (scatter(electrons_looped[1, 1, :], electrons_looped[2, 1, :]))
        scatter!(size=(1600,900),aspect_ratio=:equal)
        #xlims!((-1e-5, 1e-4))
        ylims!((-0.1e-4, 0.1e-4))
        t += dt
        K += 1
        if K == 100
            K = 0
            println(t / t_end * 100)
        end
    end
end
gif(anim, "./temp.mp4", fps=30)
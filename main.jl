using Plots, LinearAlgebra, Random, Base.Threads, StatsBase, Distances

include("functions.jl");

m_e = 9.1093837e-31;
q_e = 1.60217663e-13;

electrons = createParticlesGauss(5000, 1.3e-3, 3e-3, 3e6);


t = 0;
dt = 1e-12;
t_end = 1e3 * dt
let electrons_looped = electrons, t = t, K = 0
    global anim = @animate while t < t_end
        @time begin
            magF =  @spawn coulombForce(electrons_looped)
            electF =  @spawn lorentzForce(electrons_looped)
            stepInTime(electrons_looped, dt, fetch(electF) .+ fetch(magF))
            (scatter(electrons_looped[1, 1, :], electrons_looped[2, 1, :]))
            scatter!(size=(1600, 900), aspect_ratio=:equal)
            #xlims!((-1e-5, 1e-4))
            ylims!((-5e-3, 5e-3))
            t += dt
            #println(t / t_end * 100)
        end
    end
end
gif(anim, "./temp.mp4", fps=30)
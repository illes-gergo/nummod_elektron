using Plots, LinearAlgebra, Random, Base.Threads, StatsBase, Distances

include("functions.jl");

m_e = 9.1093837e-31;
q_e = 1.60217663e-13;

electrons = createParticlesGauss(200, 1.3e-6, 1.3e-6, 0);


t = 0;
dt = 1e-13;
t_end = 2e3*dt;
let electrons_looped = electrons, t = t
global anim = @animate  while t < t_end
        electrons_looped = stepInTime(electrons_looped, dt, coulombForce(electrons_looped))
        (scatter(electrons_looped[1,1,:],electrons_looped[2,1,:]))
        t += dt
        #sleep(0.1)
        display(t/t_end*100)
    end
end
gif(anim, "./temp.mp4", fps=30)
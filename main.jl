using Plots, LinearAlgebra, Random, Base.Threads, StatsBase, Distances

include("functions.jl");

m_e = 9.1093837e-31;
q_e = 1.60217663e-13;

electrons = createParticlesGauss(2000, 1e-5, 1e-6, 1e2);

display(scatter(electrons[1, 1, :], electrons[2, 1, :], label="Electron"));

@time proba = coulombForce(electrons);

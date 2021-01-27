module TestODEModel
using CovidModels
using LabelledArrays
using OrdinaryDiffEq
import OrdinaryDiffEq: ODEProblem
using Petri
using Test
using Plots

function makeplots_seir(sol, prefix)
    mkpath(dirname(prefix))
    p1 = plot(
        sol,
        vars = [1, 2, 3, 4],
        xlabel = "",
        ylabel = "people",
        linewidth = 3,
        title = "Cities",
        legend = false,
    )
    p = plot(
        p1,
        layout = (1, 1),
        linewidth = 3,
        link = :both,
    )
    savefig(p, "$(prefix)simple.pdf")
    p
end

T = [
    # City 1 SEIR
    ([1, 2], [3, 2]),
    ([3], [2]), # E→I
    ([2], [4]), # I→R
]

# m = Petri.Model(1:12, T, missing, missing)
# nS = length(m.S)
# β = ones(Float64, length(T))
# @test fluxes(m)(zeros(Float64, nS), ones(Float64, nS), β, 1) |> length ==
#       length(m.S)
#
# tspan = (0.0, 60.0)
# prob = ODEProblem(fluxes(m), u₀(m, 1, 0), tspan, β)
# sol = OrdinaryDiffEq.solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-8)
# @test sol.u[end][end-3] > 0.90
# @test sol.u[end][end-2] == 0
# @test sol.u[end][end-1] == 0
# @test sol.u[end][end] == 0
#

observations(x,y,σ) = begin
    r = σ*randn(length(x))
    yhat = y.+(r.*sqrt.(y))
end

add_data!(x, y, p, σ=1) = begin
    yhat = observations(x,y,σ)
    scatter!(p, x, yhat)
    return yhat
end

add_data!(sol::ODESolution, p::Plots.Plot, i::Int, σ=1) = begin
    x = sol.t
    y = [u[i] for u in sol.u]
    add_data!(x,y,p,σ)
end

# add_data!(p::Plots.Plot, t::Vector, y::Vector) = scatter!(p, t, y)

seirs = Petri.Model([:S,:E,:I,:R],LVector(
                                    exp=(LVector(S=1.0,I=1.0), LVector(I=1.0,E=1.0)),
                                    inf=(LVector(E=1.0),     LVector(I=1.0)),
                                    rec=(LVector(I=1.0),     LVector(R=1.0)),
                                    wan=(LVector(R=1.0),     LVector(S=1.0)),
                                    ))
u0 = LVector(S=100.0, E=1.0, I=0.0, R=0.0)
p = (exp=0.0035, inf=0.05, rec=0.05, wan=0.03)
f = toODE(seirs)
prob = ODEProblem(f,u0,(0.0,365.0),p)
sol = OrdinaryDiffEq.solve(prob,Tsit5())
seirsplt = plot(sol, linewidth=3.0)
# add_data!(sol, seirsplt, 3)
that = sol.t
yhat = observations(sol.t, [u[3] for u in sol.u], 1)
scatter!(seirsplt, that, yhat, label="obs_infection")

sir = Petri.Model([:S,:I,:R],LVector(
                                    exp=(LVector(S=1.0,I=1.0), LVector(I=2.0)),
                                    rec=(LVector(I=1.0),     LVector(R=1.0))))
u0 = LVector(S=100.0, I=1.0, R=0.0)
p = (exp=0.0035, rec=0.05)
f = toODE(sir)
prob = ODEProblem(f,u0,(0.0,365.0),p)
sol = OrdinaryDiffEq.solve(prob,Tsit5())
sirplt = plot(sol, linewidth=3.0)
# add_data!(sol, sirplt, 2)
scatter!(sirplt, that, yhat, label="obs_infection")

seir = Petri.Model([:S,:E,:I,:R],LVector(
                                    exp=(LVector(S=1.0,I=1.0), LVector(I=1.0,E=1.0)),
                                    inf=(LVector(E=1.0),     LVector(I=1.0)),
                                    rec=(LVector(I=1.0),     LVector(R=1.0))))
u0 = LVector(S=100.0, E=1.0, I=0.0, R=0.0)
p = (exp=0.0035, inf=0.05, rec=0.05)
f = toODE(seir)
prob = ODEProblem(f,u0,(0.0,365.0),p)
sol = OrdinaryDiffEq.solve(prob,Tsit5())
seirplt = plot(sol, linewidth=3.0)
# add_data!(sol, seirplt, 3)
scatter!(seirplt, that, yhat, label="obs_infection")


seird = Petri.Model([:S,:E,:I,:R,:D],LVector(
                                    exp=(LVector(S=1.0,I=1.0), LVector(I=1.0,E=1.0)),
                                    inf=(LVector(E=1.0),     LVector(I=1.0)),
                                    rec=(LVector(I=1.0),     LVector(R=1.0)),
                                    dth=(LVector(I=1.0),     LVector(D=1.0))))
u0 = LVector(S=100.0, E=1.0, I=0.0, R=0.0, D=0.0)
p = (exp=0.0035, inf=0.05, rec=0.05, dth=0.03)
f = toODE(seird)
prob = ODEProblem(f,u0,(0.0,365.0),p)
sol = OrdinaryDiffEq.solve(prob,Tsit5())
seirdplt = plot(sol, linewidth=3.0)
scatter!(seirdplt, that, yhat, label="obs_infection")

plotlocs = [
    (sirplt, "simple.sir.pdf")
    (seirplt, "simple.seir.pdf")
    (seirsplt, "simple.seirs.pdf")
    (seirdplt, "simple.seird.pdf")
]
map(plotlocs) do p
    savefig(p...)
end

end #module

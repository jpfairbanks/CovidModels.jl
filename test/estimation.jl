# module TestODEModel
using CovidModels
using CovidModels.Data
using OrdinaryDiffEq
import OrdinaryDiffEq: ODEProblem
using Petri
using Test
using Plots

function makeplots_seir(sol, prefix)
    mkpath(dirname(prefix))
    p1 = plot(sol,vars=[1,2,3,4], xlabel="", ylabel="people", linewidth=3,title="Cities", legend=false)
    p2 = plot(sol,vars=[5,6,7,8], xlabel="", ylabel="people", linewidth=3, legend=false)
    p3 = plot(sol,vars=[9,10,11,12], xlabel="time", ylabel="people", linewidth=3, legend=false)
    p4 = plot(sol,vars=[2,6,10], xlabel="", linewidth=3, labels=["i1" "i2" "i3"], legend=true)
    p5 = plot(sol,vars=[3,7,11], xlabel="", linewidth=3,title="Populations", labels=["e1" "e2" "e3"], legend=true)
    p6 = plot(sol,vars=[4,8,12], xlabel="time", linewidth=3, labels=["r1" "r2" "r3"], legend=true)
    p = plot(p1, p5, p2, p4, p3, p6, layout=(3,2), linewidth=3, link=:both)
    savefig(p, "$(prefix)combined.pdf")
    p
end

function makeplots_seird(sol, prefix)
    mkpath(dirname(prefix))
    p1 = plot(sol,vars=[1,2,3,4,5], xlabel="", ylabel="people", linewidth=3,title="Cities", legend=false)
    p2 = plot(sol,vars=[6,7,8,9,10], xlabel="", ylabel="people", linewidth=3, legend=false)
    p3 = plot(sol,vars=[11,12,13,14,15], xlabel="time", ylabel="people", linewidth=3, legend=false)
    p4 = plot(sol,vars=[2,7,12], xlabel="", linewidth=3, labels=["i1" "i2" "i3"], legend=true)
    p5 = plot(sol,vars=[3,8,13], xlabel="", linewidth=3,title="Populations", labels=["e1" "e2" "e3"], legend=true)
    p6 = plot(sol,vars=[5,10,15], xlabel="time", linewidth=3, labels=["d1" "d2" "d3"], legend=true)
    p = plot(p1, p5, p2, p4, p3, p6, layout=(3,2), linewidth=3, link=:both)
    savefig(p, "$(prefix)combined.pdf")
    p
end


T = [
    # City 1 SEIR
    ([1, 2], [3, 2]),
    ([3], [2]), # E→I
    ([2], [4]), # I→R
    # outflow 1→2
    ([1], [5]), # S→S′
    ([2], [6]), # I→I′
    ([3], [7]), # E→E′
    # City 2 SEIR
    ([5, 6], [7, 6]),
    ([7], [6]), # E→I
    ([6], [8]), # I→R
    # outflow 2→3
    ([5], [9]), # S→S′
    ([6], [10]),# I→I′
    ([7], [11]),# E→E′
    # City 3 SEIR
    ([9, 10], [11, 10]),
    ([11], [10]), # E→I
    ([10], [12]) # I→R
]

m = Petri.Model(1:12, T, missing, missing)
nS = length(m.S)
β = ones(Float64, length(T))

tspan = (0.0,60.0)
prob = ODEProblem(fluxes(m), u₀(m, 1,0), tspan, β)
sol = OrdinaryDiffEq.solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)


u0 = zeros(Float64, length(m.S))
u0[1]  = 10000
u0[5]  = 10000
u0[9]  = 10000
u0[2]  = 1
u0
βseir = [10/sum(u0), 1/2, 1/5]
βtravel = [1/2, 1/2, 1/2]/1000
β = vcat(βseir, βtravel, βseir, βtravel, βseir)
@show β
prob = ODEProblem(fluxes(m), u0, tspan, β)
sol = OrdinaryDiffEq.solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)
@show sol.u[end]

using Plots
makeplots_seir(sol, "img/estimation/seir")



T = [
    # City 1 SEIₐIRD
    ([1, 2], [3, 2]), # S+Iₐ → E+Iₐ
    ([3], [2]),       # E  → Iₐ
    ([2], [4]),       # Iₐ → I
    ([2], [5]),       # Iₐ → R
    ([4], [5]),       # I  → R
    ([4], [6]),       # I  → D
]

m = Petri.Model(1:6, T, missing, missing)
nS = length(m.S)
β = ones(Float64, length(T))

tspan = (0.0,60.0)
prob = ODEProblem(fluxes(m), u₀(m, 1,0), tspan, β)
sol = OrdinaryDiffEq.solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)


u0 = zeros(Float64, length(m.S))
u0[1]  = 60000000
u0[2]  = 1
u0

βseir = [10/5sum(u0), 1/2, 1/5, 1/3, 1/4, 1/7]
# βtravel = [1/2, 1/2, 1/2]/1000
# β = vcat(βseir, βtravel, βseir, βtravel, βseir)
β = βseir
@show β
prob = ODEProblem(fluxes(m), u0, tspan, β)
sol = OrdinaryDiffEq.solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)
@show sol.u[end]

# makeplots_seir(sol, "img/estimation/seiird")

plot(sol, vars=[3, 4, 5, 6], label=["IA" "I" "R" "D"])
# end

using CSV
using Query
using DataFrames

df = CSV.read("../eucovid.csv")

italy = @from i in df begin
    @where i.Country == "Italy" && i.Cases > 0
    @select i
    @collect DataFrame
end

caseplot(italy)

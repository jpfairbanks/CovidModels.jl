# module TestODEModel
using CovidModels
using CovidModels.Data
using OrdinaryDiffEq
import OrdinaryDiffEq: ODEProblem
using Petri
using Test
using Plots
using StatsPlots
using Colors

function fitplot(sol, df)
    p = plot(sol, vars=[2, 4, 5, 6], label=["IA" "I" "R" "D"])
    @df df scatter!(p,
            0:length(italy.Deaths)-1,
            [:Cases, :Deaths, :Recoveries],
            label=["cases" "deaths" "recoveries"],
            marker=true,
            color=[RGB(1,0,0) RGB(0.5,0,0.5) RGB(0,1,0)]
            )
    return p
end

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

function lossfunc(df::DataFrame, m, u0)
    tspan = (0.0, length(df.Deaths)-1)
    mean(x) = sum(x)/length(x)
    t = tspan[1]:tspan[end]
    f(sol) = begin
        deaths = [u[end] for u in sol.(t)]
        recovr = [u[end-1] for u in sol.(t)]
        infsym = [u[end-2] for u in sol.(t)]
        dt = mean((deaths .- df.Deaths).^2) / mean(df.Deaths)
        rt = mean((recovr .- df.Recoveries).^2) / mean(df.Recoveries)
        it = mean((infsym .- df.Cases).^2) / mean(df.Cases)
        return dt + rt + it
    end

    g(θ) = begin
        prob = ODEProblem(fluxes(m), u0, tspan, θ)
        sol = OrdinaryDiffEq.solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)
        neg = sum(abs.(θ[θ .< 0]))
        λ = 10000000
        # if final deaths > final recoveries then penalty
        # deathpenalty = sol.u[end][end] > sol.u[end][end-1] ? 1 : 0
        # iapenalty = sum(tanh.([u[3] - u[2] for u in sol.(t)]).+1)
        iapenalty = 0
        lossval, regval = f(sol), λ*(neg + iapenalty)
        # lossval, regval = f(sol) , λ*(neg + deathpenalty)
        println("loss:$lossval\tregval:$regval")
        lossval+regval
    end
    return g
end

u0_italy = [60000000.0, 0, 1, 0, 0, 0]
ℓ = lossfunc(italy, m, u0_italy)

ℓ(β)

using Optim

function calibrate(k)
    β = [40/5sum(u0), 11/2, 1/1, 2/13, 1/18, 1/35]/k
    prob = ODEProblem(fluxes(m), u0_italy, (0.0,95), β)
    sol = OrdinaryDiffEq.solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)
    # plot(sol, vars=[2, 4, 5, 6, 3], label=["IA" "I" "R" "D" "E"])
    # p = fitplot(sol, italy)
    p = plot(sol, vars=[4, 5, 6], label=["I" "R" "D"],
            color=[RGB(1,0,0) RGB(0,1,0) RGB(0.5,0,0.5)])
    @df df scatter!(p,
            0:length(italy.Deaths)-1,
            [:Cases, :Deaths, :Recoveries],
            label=["cases" "deaths" "recoveries"],
            marker=true,
            color=[RGB(1,0,0) RGB(0.5,0,0.5) RGB(0,1,0)]
            )
    return β, sol, p
end

calibrate(15.5)[end]

β = [40/5sum(u0), 11/2, 1/1, 2/13, 1/18, 1/35]/15
# β = [60/5sum(u0), 11/2, 1/1, 2/13, 1/24, 1/39]/2
prob = ODEProblem(fluxes(m), u0_italy, (0.0,64), β)
sol = OrdinaryDiffEq.solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)
plot(sol, vars=[2, 4, 5, 6, 3], label=["IA" "I" "R" "D" "E"])
fitplot(sol, italy)

lb, ub = zeros(Float64, 6), 100*ones(Float64, 6)

res = optimize(ℓ, β)
# res = optimize(ℓ, lb, ub, β, Fminbox(NelderMead()))
@show res.minimizer
@show res.minimizer -β

# candidate minimizer
β̂ = [
 9.478916256341902e-9
 0.2915773795849941
 0.09434402982528367
 0.0061044262074341285
 0.016241417062038518
 0.02234216514573622
]

prob = ODEProblem(fluxes(m), u0_italy, (0.0,128), res.minimizer)
sol = OrdinaryDiffEq.solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)
fitplot(sol, italy)

prob = ODEProblem(fluxes(m), u0_italy, (0.0,70), res.minimizer)
sol = OrdinaryDiffEq.solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)
fitplot(sol, df)

# trying to add constraints to the parameter space
realparams(θ) = [θ[1:4]..., θ[5], θ[5]/5]
res = optimize(θ->ℓ(realparams(θ)), β[1:6])
@show res.minimizer
@show res.minimizer -β

prob = ODEProblem(fluxes(m), u0_italy, (0.0,128), realparams(res.minimizer))
sol = OrdinaryDiffEq.solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)
fitplot(sol, df)

using OrdinaryDiffEq
import OrdinaryDiffEq: ODEProblem
using Petri

export linreg, peakgap, paramsweep, fluxes, u₀, makeplots_seir, makeplots_seird, savedata

#this function was removed in Julia 0.7
"""    linreg(x,y)

perform a simple linear regression for scalar case. Coefficients are returned
as (intercept, scalar)
"""
linreg(x, y) = hcat(fill!(similar(x), 1), x) \ y

"""    peakgap(s1, s2)

identify the time delay between the peak in two time series. Time series should
be stored as (x(t), t) for consistency with the tuples(sol) from OrdinaryDiffEq.
"""
function peakgap(s1, s2)
    @show p1, i1 = findmax(first.(s1))
    @show p2, i2 = findmax(first.(s2))
    return s2[i2][end] - s1[i1][end]
end

"""    peakgap(sol, i::Int, j::Int)

identify the time delay between the peak in two dimensions of an ODE solution structure.
"""
function peakgap(sol, i::Int, j::Int)
    s1 = [(u[i],t) for (u,t) in tuples(sol)]
    s2 = [(u[j],t) for (u,t) in tuples(sol)]
    return peakgap(s1,s2)
end

"""    paramsweep(f::Function, m::Petri.Model, initials, tspan, params)

solve an ODEProblem for each parameter setting in params and apply the function f to the solution.

Usage:
    paramsweep(m,initials,tspan,[p1,p2,p3]) do sol
        return f(sol)
    end
"""
function paramsweep(f::Function, m::Petri.Model, initials, tspan, params)
    map(params) do p
        prob = ODEProblem(fluxes(m), initials, tspan, p)
        sol = OrdinaryDiffEq.solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)
        return f(sol)
    end
end

"""    fluxes(m::Petri.Model)

a PetriNet interpreter that computes the mass action kinetics of a petri net.

Usage: pass this to ODEProblem to set up an ODE for a given model.
"""
function fluxes(m::Petri.Model)
    S = m.S
    T = m.Δ
    nS = length(S)
    nT = length(T)
    ϕ = zeros(Float64, nT)
    f(du, u, p, t) = begin
        for (i, t) in enumerate(T)
            ins = t[1]
            # TODO: accomodate multiplicites here
            ϕ[i] = p[i]*prod(u[ins])
        end
        for i in S
            du[i] = 0
        end
        for (i, t) in enumerate(T)
            ins = t[1]
            out = t[2]
            for s in ins
                # TODO: accomodate multiplicites here
                du[s] -= ϕ[i]
            end
            for s in out
                # TODO: accomodate multiplicites here
                du[s] += ϕ[i]
            end
        end
        return du
    end
    return f
end

u₀(m::Petri.Model) = begin
    zeros(Float64, length(m.S))
end

u₀(m::Petri.Model, initialS, initialI=1) = begin
    u0=zeros(Float64, length(m.S))
    u0[1] = initialS
    u0[2] = initialI
    return u0
end

function savedata(sol::ODESolution, colnames, fname::String)
    open(fname, "w") do fp
        println(fp, "time, $colnames")
        map(tuples(sol)) do (u,t)
            print(fp, "$t")
            map(u) do x
                print(fp, ",$x")
            end
            println(fp, "")
        end
    end
end

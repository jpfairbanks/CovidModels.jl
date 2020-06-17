module TestValidation

using CovidModels
using OrdinaryDiffEq
import OrdinaryDiffEq: ODEProblem
using Petri
using Test
using SemanticModels.PetriModels

mkpath("outputs")

# states are [S, I, R]
Psir = PetriModel(Petri.Model(1:3, [
  ([1, 2], [2, 2]), # exposure
  ([2], [3]),     # recovery
]))

m = Psir.model

u0 = zeros(Float64, length(m.S))
u0[1] = 10000
u0[2] = 1

# probability of infection, recovery
β = [10 / sum(u0), 1 / 5]

tspan = (0, 100.0)
prob = ODEProblem(fluxes(m), u0, tspan, β)
sol = OrdinaryDiffEq.solve(prob, alg = Tsit5())

savedata(sol, "S,I,R", "outputs/sirdata.csv")

# states are [S, I, R, D]
Psird = PetriModel(Petri.Model(
  1:4,
  [
    ([1, 2], [2, 2]), # exposure
    ([2], [3]),     # recovery
    ([2], [4]),     # death
  ]
))

m = Psird.model

u0 = zeros(Float64, length(m.S))
u0[1] = 10000
u0[2] = 1

# probability of infection, recovery, death
β = [10 / sum(u0), 1 / 5, 1 / 10]

tspan = (0, 100.0)
prob = ODEProblem(fluxes(m), u0, tspan, β)
sol = OrdinaryDiffEq.solve(prob, alg = Tsit5())

savedata(sol, "S,I,R,D", "outputs/sirddata.csv")

# states are [S, I, E, R]
Pseir = PetriModel(Petri.Model(
  1:5,
  [
    ([1, 2], [3, 2]), # exposure
    ([3], [2]),     # onset
    ([2], [4]),     # recovery
  ],
))

m = Pseir.model

u0 = zeros(Float64, length(m.S))
u0[1] = 10000
u0[2] = 1

# probability of exposure, progression, recovery
β = [10 / sum(u0), 1 / 2, 1 / 5]

tspan = (0, 100.0)
prob = ODEProblem(fluxes(m), u0, tspan, β)
sol = OrdinaryDiffEq.solve(prob, alg = Tsit5())

savedata(sol, "S,I,E,R", "outputs/seirdata.csv")


# states are [S, I, E, R, D]
Pseird = PetriModel(Petri.Model(
  1:5,
  [
    ([1, 2], [3, 2]), # exposure
    ([3], [2]),     # onset
    ([2], [4]),     # recovery
    ([2], [5]),     # death
  ],
))

m = Pseird.model

u0 = zeros(Float64, length(m.S))
u0[1] = 10000
u0[2] = 1

seirdparams() = begin
  # probability of exposure, progression, recovery, death
  βseird = [10 / sum(u0), 1 / 2, 1 / 5, 1 / 16]
end

tspan = (0, 100.0)
prob = ODEProblem(fluxes(m), u0, tspan, seirdparams())
sol = OrdinaryDiffEq.solve(prob, alg = Tsit5())

savedata(sol, "S,I,E,R,D", "outputs/seirddata.csv")
end

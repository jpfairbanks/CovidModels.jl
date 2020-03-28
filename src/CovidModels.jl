module CovidModels

using SemanticModels
using OrdinaryDiffEq

include("odemodel.jl")
include("covid.jl")

end

module TestCovidModels

using CovidModels

println("running validation tests")
include("validation.jl")
println("running odemodel tests")
include("odemodel.jl")
println("running CT tests")
include("covid.jl")
end

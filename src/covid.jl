using Catlab
using Catlab.Doctrines
using Catlab.Graphics
using Catlab.WiringDiagrams
using Catlab.Programs
import Base.Multimedia: display
import Catlab.Graphics: to_graphviz, LeftToRight
import Base: (==), length, show
using Petri
using SemanticModels.ModelTools.CategoryTheory
import SemanticModels.ModelTools.CategoryTheory: undecorate, ⊔
using SemanticModels.ModelTools.PetriModels
using SemanticModels.ModelTools.PetriCospans
import SemanticModels.ModelTools.PetriCospans: otimes_ipm, compose_pushout
import SemanticModels.ModelTools: model

import Catlab.Doctrines:
  Ob, Hom, dom, codom, compose, ⋅, ∘, id, oplus, otimes, ⊗, ⊕, munit, mzero, braid,
  dagger, dunit, dcounit, mcopy, Δ, delete, ◊, mmerge, ∇, create, □,
  plus, zero, coplus, cozero, meet, top, join, bottom

export wd, draw, model, Disease

wd = to_wiring_diagram
draw(d::WiringDiagram) = to_graphviz(add_junctions(d), orientation=LeftToRight, labels=true)
draw(d::HomExpr) = draw(wd(d))

model(c::PetriCospan) = left(c.f).d[1].model

println("Making COVID Model")

@present Disease(FreeEpidemiology) begin
    S::Ob
    E::Ob
    I::Ob
    R::Ob
end

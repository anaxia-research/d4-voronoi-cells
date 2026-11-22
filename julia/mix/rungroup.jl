# usage: julia rungroup.jl cert.npz vectors.json pre|full
include(joinpath(@__DIR__, "..", "src", "D4Checks.jl"))
using .D4Checks, JSON
const C = D4Checks.Combo
base = C.load_case(joinpath(@__DIR__, "case_t1.json"))
vecs = [Int.(v) for v in JSON.parsefile(ARGS[2])]
case = C.Case(basename(ARGS[2]), base.M, base.types, base.trange, base.edges, base.bintype, vecs, base.raw)
C.run(ARGS[1], 10; case=case, mixtype=true, bb=(ARGS[3] == "full"), border=true, mg2=2e-6, mg3=5e-6, mgm=1e-6)

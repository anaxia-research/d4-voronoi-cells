include(joinpath(@__DIR__, "..", "src", "D4Checks.jl"))
using .D4Checks
R = D4Checks.Combo.run(ARGS[1], 10; case=ARGS[2], mixtype=true, bb=true, border=true, mg2=2e-6, mg3=5e-6, mgm=1e-6)

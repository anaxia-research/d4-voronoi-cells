# Build the columns of the mixing programme: every proved 27-centre certificate, as a function
# of the seven bins (one type per bin), sampled on pairs and triples of bins.
include(joinpath(@__DIR__, "..", "src", "D4Checks.jl"))
using .D4Checks
const C = D4Checks.Combo
using NPZ, JSON, Random, Printf, LinearAlgebra
const R = joinpath(@__DIR__, "..", "data")
const EDGES = [2.0, 2.0161, 2.05, 2.1, 2.15, 2.2, 2.35, sqrt(6.0)]
const D3 = 10
names = [l for l in readlines(joinpath(@__DIR__, "certs.txt")) if !isempty(l)]
struct Col; name; Af; z; t; m; groups; end
function loadcol(nm, cf)
    MT = JSON.parsefile(joinpath(@__DIR__, "certmt.json"))[nm]
    Z = npzread("$R/combo/$(nm)_d10.npz", MT ? ["A", "z", "t", "m", "x3g", "gtyp"] : ["A", "z", "t", "m", "x3"])
    J = JSON.parsefile("$R/combo/cases/$cf")
    if haskey(Z, "gtyp")
        GT = Int.(Z["gtyp"])
        groups = [C.Group(C.Blocks(k => C.psd_exact(Mx) for (k, Mx) in C.blocks_from_x3(Float64.(Z["x3g"][g, :]), D3, 3)), [GT[g, i] for i in 1:7]) for g in 1:size(GT, 1)]
        x3 = [Float64.(Z["x3g"][g, :]) for g in 1:size(GT, 1)]; gt = [GT[g, :] for g in 1:size(GT, 1)]
    else
        T = String.(J["types"]); tm = [findfirst(==(String(b[1])), T) - 1 for b in J["bins"]]
        groups = [C.Group(C.Blocks(k => C.psd_exact(Mx) for (k, Mx) in C.blocks_from_x3(Float64.(vec(Z["x3"])), D3, length(T))), tm)]
        x3 = [Float64.(vec(Z["x3"]))]; gt = [tm]
    end
    Af = Float64.(Z["A"]); (Col(nm, [Af[k, :, :] for k in 1:13], Float64.(Z["z"]), Float64(Z["t"][]), Float64.(vec(Z["m"])), groups), x3, gt)
end
cases = JSON.parsefile(joinpath(@__DIR__, "certcases.json"))
cols = []; X3 = []; GTs = []
for nm in names
    c, x3, gt = loadcol(nm, cases[nm]); push!(cols, c); push!(X3, x3); push!(GTs, gt)
end
println("columns: ", length(cols)); flush(stdout)
nc = length(cols)
# point terms and constants
PT = [Float64(C.g_point(c.groups, D3, b)) for c in cols, b in 1:7]
TT = [c.t for c in cols]; MM = [c.m[b] for c in cols, b in 1:7]
rng = Xoshiro(1)
top(p, q) = (p*p + q*q - 4)/(2p*q)
npair = parse(Int, get(ENV, "NPAIR", "1200")); ntri = parse(Int, get(ENV, "NTRI", "1200"))
prows = Vector{Vector{Float64}}(); pkind = Int[]; ppts = Vector{Vector{Float64}}()
pk = [(a, b) for a in 1:7 for b in a:7]
for (ki, (a, b)) in enumerate(pk)
    Pus = [Float64.(C.g_pair(c.groups, D3, a, b)) for c in cols]
    for i in 1:npair
        p = EDGES[a] + (EDGES[a+1] - EDGES[a])*rand(rng); q = EDGES[b] + (EDGES[b+1] - EDGES[b])*rand(rng)
        r = rand(rng); r < 0.1 && (p = EDGES[a]); 0.1 <= r < 0.2 && (q = EDGES[b])
        t = top(p, q); u = rand(rng) < 0.5 ? t - 0.02*rand(rng)^2 : -1 + (t + 1)*sqrt(rand(rng))
        pi = C.pair_float(p/2, q/2, u)
        push!(prows, [C.kfloat(c.Af, p, q, u) + C.polyval(Pus[j], u) - pi for (j, c) in enumerate(cols)])
        push!(pkind, ki); push!(ppts, [p, q, u])
    end
    @printf("pair kind %d (%d,%d) done\n", ki, a, b); flush(stdout)
end
tk = [(a, b, c) for a in 1:7 for b in a:7 for c in b:7]
trows = Vector{Vector{Float64}}(); tkind = Int[]; tpts = Vector{Vector{Float64}}()
for (ki, (a, b, c3)) in enumerate(tk)
    Ps = [C.FPoly(C.g_triple(c.groups, D3, a, b, c3)) for c in cols]
    T = (top(EDGES[a+1], EDGES[b+1]), top(EDGES[a+1], EDGES[c3+1]), top(EDGES[b+1], EDGES[c3+1]))
    got = 0
    while got < ntri
        X = randn(rng, 3, rand(rng, 2:4)); for i in 1:3; X[i, :] ./= norm(X[i, :]); end
        g = [dot(X[1, :], X[2, :]), dot(X[1, :], X[3, :]), dot(X[2, :], X[3, :])]
        # push toward the tops: random convex move
        if rand(rng) < 0.5
            g = [min(g[i], T[i]) for i in 1:3]; s = rand(rng)
            g2 = [T[i] - 0.05*rand(rng)^2 for i in 1:3]; g = s .* g .+ (1 - s) .* g2
        end
        (all(-1 <= g[i] <= T[i] for i in 1:3) && C.gramdet(g) >= 0) || continue
        push!(trows, [P(g) for P in Ps]); push!(tkind, ki); push!(tpts, g); got += 1
    end
    @printf("triple kind %d (%d,%d,%d) done\n", ki, a, b, c3); flush(stdout)
end
npzwrite(joinpath(@__DIR__, get(ENV, "OUT", "samples.npz")), Dict("P" => reduce(hcat, prows)', "pk" => pkind, "pp" => reduce(hcat, ppts)',
    "T" => reduce(hcat, trows)', "tk" => tkind, "tp" => reduce(hcat, tpts)', "PT" => PT, "TT" => TT, "MM" => MM))
println("saved")

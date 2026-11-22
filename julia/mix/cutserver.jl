# Build the columns of the mixing programme: every proved 27-centre certificate, as a function
# of the seven bins (one type per bin), sampled on pairs and triples of bins.
include("/home/user/d4-voronoi-cells/julia/src/D4Checks.jl")
using .D4Checks
const C = D4Checks.Combo
using NPZ, JSON, Random, Printf, LinearAlgebra
const R = "/home/user/d4-voronoi-cells/julia/data"
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

# Cutting planes: for a mixture lam of the columns, find the near-maximisers of each pair and triple
# function (as the exact check's threshold search does), and record every column's value there.
top(p, q) = (p*p + q*q - 4)/(2p*q)
pk = [(a, b) for a in 1:7 for b in a:7]
tk = [(a, b, c) for a in 1:7 for b in a:7 for c in b:7]
PU = [[Float64.(C.g_pair(c.groups, D3, a, b)) for c in cols] for (a, b) in pk]
PS0 = [[C.g_triple(c.groups, D3, a, b, c3) for c in cols] for (a, b, c3) in tk]
FS = [[C.FPoly(P) for P in Ps] for Ps in PS0]
println("polys ready"); flush(stdout)
function bestpts(f, feas, pts, repair; k=12)
    vals = [f(x) for x in pts]
    out = Vector{Float64}[]
    for i in partialsortperm(vals, 1:min(k, length(vals)); rev=true)
        x = copy(pts[i]); fx = f(x); h = 0.02
        dirs = [collect(d) for d in Iterators.product(ntuple(_ -> (-1, 0, 1), 3)...) if any(!=(0), d)]
        while h > 1e-9
            moved = false
            for d in dirs, y in repair(x .+ h .* d)
                feas(y) || continue
                fy = f(y); if fy > fx; x = y; fx = fy; moved = true; end
            end
            moved || (h /= 2)
        end
        push!(out, x); push!(out, pts[i])
    end
    out
end
rng = Xoshiro(parse(Int, get(ENV, "SEED", "5")))
NR = parse(Int, get(ENV, "NR", "20000"))
Q = joinpath(@__DIR__, "cutq"); mkpath(Q)
function docuts(lam)
    act = findall(>(1e-9), lam)
    prows = Vector{Vector{Float64}}(); pkind = Int[]; trows = Vector{Vector{Float64}}(); tkind = Int[]
    for (ki, (a, b)) in enumerate(pk)
        a0, a1, b0, b1 = EDGES[a], EDGES[a+1], EDGES[b], EDGES[b+1]
        Af = sum(lam[j] .* cols[j].Af for j in act); Pu = sum(lam[j] .* PU[ki][j] for j in act)
        f(x) = C.kfloat(Af, x[1], x[2], x[3]) + C.polyval(Pu, x[3]) - C.pair_float(x[1]/2, x[2]/2, x[3])
        feas(x) = a0 <= x[1] <= a1 && b0 <= x[2] <= b1 && -1 <= x[3] <= top(x[1], x[2])
        rep(y) = (q1 = clamp(y[1], a0, a1); q2 = clamp(y[2], b0, b1); ([q1, q2, clamp(y[3], -1.0, top(q1, q2))],))
        pts = Vector{Float64}[]
        for i in 1:NR
            p = a0 + (a1 - a0)*rand(rng); q = b0 + (b1 - b0)*rand(rng)
            r = rand(rng); r < 0.15 && (p = a0); 0.15 <= r < 0.3 && (q = b0); 0.3 <= r < 0.4 && (p = a1); 0.4 <= r < 0.5 && (q = b1)
            t = top(p, q); u = rand(rng) < 0.3 ? t - 3e-3*rand(rng) : -1 + (t + 1)*sqrt(rand(rng))
            push!(pts, [p, q, u])
        end
        for x in bestpts(f, feas, pts, rep)
            p, q, u = x; pi = C.pair_float(p/2, q/2, u)
            push!(prows, [C.kfloat(c.Af, p, q, u) + C.polyval(PU[ki][j], u) - pi for (j, c) in enumerate(cols)]); push!(pkind, ki)
        end
    end
    for (ki, (a, b, c3)) in enumerate(tk)
        T = (top(EDGES[a+1], EDGES[b+1]), top(EDGES[a+1], EDGES[c3+1]), top(EDGES[b+1], EDGES[c3+1]))
        Pd = Dict{NTuple{3,Int},Float64}(); for j in act, (e, co) in zip(FS[ki][j].exps, FS[ki][j].cf); Pd[e] = get(Pd, e, 0.0) + lam[j]*co; end
        it = sort(collect(Pd); by=first)
        f = C.FPoly([e for (e, _) in it], [c for (_, c) in it], ntuple(v -> maximum(e[v] for (e, _) in it; init=0), 3))
        feas(x) = all(-1 <= x[i] <= T[i] for i in 1:3) && C.gramdet(x) >= 0
        repair = C.gram_repair(-ones(3), collect(T))
        pts = Vector{Float64}[]
        for dim in (4, 3, 2)
            got = 0
            while got < NR
                X = randn(rng, 3, dim); for i in 1:3; X[i, :] ./= norm(X[i, :]); end
                g = [dot(X[1, :], X[2, :]), dot(X[1, :], X[3, :]), dot(X[2, :], X[3, :])]
                if rand(rng) < 0.5
                    g = [min(g[i], T[i]) for i in 1:3]; s = rand(rng); g2 = [T[i] - 0.05*rand(rng)^2 for i in 1:3]; g = s .* g .+ (1 - s) .* g2
                end
                (all(-1 <= g[i] <= T[i] for i in 1:3) && C.gramdet(g) >= 0) || continue
                push!(pts, g); got += 1
            end
        end
        ax = [vcat(-1.0, range(-1.0, T[i], length=30), T[i]) for i in 1:3]
        for x in ax[1], y in ax[2], z in ax[3]; g = [x, y, z]; C.gramdet(g) >= 0 && push!(pts, g); end
        for g in bestpts(f, feas, pts, repair)
            push!(trows, [F(g) for F in FS[ki]]); push!(tkind, ki)
        end
    end
    (prows, pkind, trows, tkind)
end
println("server ready"); flush(stdout)
while true
    fs = sort([f for f in readdir(Q) if endswith(f, ".json")])
    if isempty(fs); sleep(5); continue; end
    f = fs[1]; lam = Float64.(JSON.parsefile(joinpath(Q, f)))
    t0 = time(); prows, pkind, trows, tkind = docuts(lam)
    out = joinpath(Q, replace(f, ".json" => ".npz"))
    npzwrite(out[1:end-4] * "_tmp.npz", Dict("P" => reduce(hcat, prows)', "pk" => pkind, "T" => reduce(hcat, trows)', "tk" => tkind))
    mv(out[1:end-4] * "_tmp.npz", out; force=true); mv(joinpath(Q, f), joinpath(Q, f * ".done"); force=true)
    println("cuts ", f, " ", round(time() - t0, digits=1), " s"); flush(stdout)
end

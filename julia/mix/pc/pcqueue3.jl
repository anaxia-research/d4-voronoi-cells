# Runs the exact checks listed in a jobs file, N at a time, one julia process per job.
# usage (from the julia folder):  julia --project=. mix/pc/pcqueue3.jl JOBSFILE [N]
# Each line: certificate.npz vectors.json KIND   (KIND = PAIRS or a triple kind such as A1A2B4)
# The jobs file and mix/pc/N.txt (if present) are read again every loop, so jobs can be added
# and N changed without a restart. A job whose log shows its PASS line is done. A job whose log
# has no PASS or FAIL line but was written in the last 8 hours is taken to be running already
# (for example by an earlier queue process) and is left alone.
const PC = @__DIR__
const M = dirname(PC)
const J = dirname(M)
JF = joinpath(PC, ARGS[1]); N0 = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 5
logdir = joinpath(PC, "logs"); mkpath(logdir)
tag(j) = replace(j[2], "_vecs.json" => "") * "_" * j[3]
logf(j) = joinpath(logdir, tag(j) * ".log")
passline(j) = j[3] == "PAIRS" ? "PASS f <= m_b on every bin" : "PASS TRIPLE3_" * j[3] * " "
logtext(j) = isfile(logf(j)) ? read(logf(j), String) : ""
passed(j) = occursin(passline(j), logtext(j))
failed(j) = occursin(r"\bFAIL", logtext(j))
readjobs() = [String.(split(l)) for l in readlines(JF) if !isempty(strip(l))]
getN() = (f = joinpath(PC, "N.txt"); isfile(f) ? parse(Int, strip(read(f, String))) : N0)
running = Any[]; bad = Set{String}(); t0 = time()
ext(j) = isfile(logf(j)) && !passed(j) && !failed(j) && time() - mtime(logf(j)) < 8*3600
mine(j) = any(r -> r[2] == j, running)
while true
    for r in copy(running)
        p, j, fio = r
        if !process_running(p)
            close(fio); filter!(x -> x !== r, running)
            ok = passed(j) && !failed(j); ok || push!(bad, tag(j))
            println(ok ? "PASS " : "PROBLEM ", tag(j)); flush(stdout)
        end
    end
    jobs = readjobs()
    others = [j for j in jobs if !mine(j) && ext(j)]
    pending = [j for j in jobs if !mine(j) && !passed(j) && !(tag(j) in bad) && !ext(j)]
    while length(running) + length(others) < getN() && !isempty(pending)
        j = popfirst!(pending)
        args = j[3] == "PAIRS" ? ["pairs"] : (kf = joinpath(logdir, tag(j) * ".kind"); write(kf, j[3] * "\n"); ["only", kf])
        cmd = `$(Base.julia_cmd()) --project=$J $(joinpath(M, "rungroup2.jl")) $(joinpath(M, j[1])) $(joinpath(M, j[2])) $args`
        fio = open(logf(j), "w")
        push!(running, (run(pipeline(cmd; stdout=fio, stderr=fio); wait=false), j, fio))
        println("start ", tag(j)); flush(stdout)
    end
    np = count(passed, jobs)
    open(joinpath(PC, "status.txt"), "w") do io
        println(io, "elapsed $(round((time()-t0)/3600, digits=2)) h; passed $np of $(length(jobs)); running $(length(running)) (+$(length(others)) from an earlier queue); waiting $(length(pending))")
        for (p, j, fio) in running; println(io, "running ", tag(j)); end
        for j in others; println(io, "running (earlier queue) ", tag(j)); end
        for t in bad; println(io, "PROBLEM ", t); end
        np == length(jobs) && println(io, "ALL DONE")
    end
    sleep(30)
end

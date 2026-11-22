# Runs the exact checks listed in jobs.txt, N at a time, one julia process per job.
# usage (from the julia folder):  julia --project=. mix/pc/pcqueue.jl [N]
# Each line of jobs.txt: certificate.npz vectors.json KIND   (KIND = PAIRS or a triple kind such as A1A2B4)
# Resumable: a job whose log already shows its PASS line is skipped. Progress goes to mix/pc/status.txt.
const PC = @__DIR__
const M = dirname(PC)
const J = dirname(M)
N = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 12
jobs = [String.(split(l)) for l in readlines(joinpath(PC, "jobs.txt")) if !isempty(strip(l))]
logdir = joinpath(PC, "logs"); mkpath(logdir)
tag(j) = replace(j[2], "_vecs.json" => "") * "_" * j[3]
logf(j) = joinpath(logdir, tag(j) * ".log")
passline(j) = j[3] == "PAIRS" ? "PASS f <= m_b on every bin" : "PASS TRIPLE3_" * j[3] * " "
logtext(j) = isfile(logf(j)) ? read(logf(j), String) : ""
passed(j) = occursin(passline(j), logtext(j))
failed(j) = occursin(r"\bFAIL", logtext(j))
pending = filter(!passed, jobs)
running = Any[]
finished_bad = String[]
t0 = time()
function status()
    np = count(passed, jobs)
    open(joinpath(PC, "status.txt"), "w") do io
        println(io, "elapsed $(round((time()-t0)/3600, digits=2)) h; passed $np of $(length(jobs)); running $(length(running)); waiting $(length(pending))")
        for (p, j, fio) in running; println(io, "running ", tag(j)); end
        for t in finished_bad; println(io, "PROBLEM ", t); end
        isempty(pending) && isempty(running) && println(io, "ALL DONE")
    end
end
while !isempty(pending) || !isempty(running)
    for r in copy(running)
        p, j, fio = r
        if !process_running(p)
            close(fio); filter!(x -> x !== r, running)
            ok = passed(j) && !failed(j)
            ok || push!(finished_bad, tag(j))
            println(ok ? "PASS " : "PROBLEM ", tag(j))
        end
    end
    while length(running) < N && !isempty(pending)
        j = popfirst!(pending)
        args = if j[3] == "PAIRS"
            ["pairs"]
        else
            kf = joinpath(logdir, tag(j) * ".kind"); write(kf, j[3] * "\n"); ["only", kf]
        end
        cmd = `$(Base.julia_cmd()) --project=$J $(joinpath(M, "rungroup2.jl")) $(joinpath(M, j[1])) $(joinpath(M, j[2])) $args`
        get(ENV, "DRY", "") == "1" && (cmd = `$(Base.julia_cmd()) -e "println(\"PASS TRIPLE3_$(j[3]) dry\"); println(\"PASS f <= m_b on every bin\")"`)
        fio = open(logf(j), "w")
        p = run(pipeline(cmd; stdout=fio, stderr=fio); wait=false)
        push!(running, (p, j, fio))
        println("start ", tag(j))
    end
    status()
    sleep(get(ENV, "DRY", "") == "1" ? 1 : 30)
end
status()
println(isempty(finished_bad) ? "ALL PASS" : "PROBLEMS: " * join(finished_bad, ", "))

# Exact checks for the 27-centre groups (O1)

From the `julia` folder of this repository:

    julia --project=. -e "using Pkg; Pkg.instantiate()"
    julia --project=. mix/pc/pcqueue.jl 12

The number is how many checks run at once. `mix/pc/jobs.txt` lists the
checks (certificate, vectors, kind); `mix/pc/status.txt` shows progress and
`mix/pc/logs/` holds one log per check. The queue can be stopped and started
again: finished checks are skipped.

Together with `mix/cloudlogs/` (positivity, pair inequalities, bins and the
triple kinds finished earlier), a group is proved when every one of its
jobs shows its PASS line.

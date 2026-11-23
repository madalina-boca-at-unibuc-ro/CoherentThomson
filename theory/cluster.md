# Cluster (Slurm)

`theory/srun.txt` and `theory/slurm_head_node_report.md`, joined on 2026-11-23 with their text unchanged.

## 1. Running jobs (was `srun.txt`)

```text
srun --partition=debug --time=01:00:00 --exclusive --pty bash

or

srun --partition=debug --time=01:00:00 --cpus-per-task=64 --pty bash
(cpus-per-task can be up to 256)


or

srun --partition=interactive --time=04:00:00 --exclusive --pty bash


or 

srun --partition=interactive --time=04:00:00 --exclusive --exclude=hpc-db01,hpc-ctrl01 --pty bash


srun --partition=cpu-wide --time=24:00:00 --exclusive --exclude=hpc-db01,hpc-ctrl01 --pty bash




# or pick one:  --nodelist=hpc-node05

NOTE: every CPU partition (debug, interactive, cpu-*) also contains hpc-ctrl01 (head node, slurmctld) and
hpc-db01, with the same Weight as the compute nodes. Without --exclude, srun can land on hpc-ctrl01 itself:
the prompt looks unchanged ("still on ctrl01"), but you are inside a job with only 128 usable CPUs
(CoreSpecCount reserves half) on a loaded node. Check `hostname` / `squeue -u $USER` after srun.
Ending the job: Ctrl+D once. In the VS Code Remote-SSH terminal a second Ctrl+D only closes the terminal tab.
Reported to the admin: section 2 below

Partition time limits (sinfo, 2026-10-28): debug 1 h, interactive 4 h, cpu-wide 1 day, cpu-medium 1.5 days,
cpu-standard 3 days, cpu-narrow 5 days, cpu-long 12 days.


LONG RUNS THAT SURVIVE A DISCONNECT (tmux on the head node)
============================================================
A job lives only as long as the terminal holding its srun. A NetBird/SSH drop can kill that terminal (VS Code
sometimes keeps it alive, not guaranteed), and with it the job and the run. tmux runs on hpc-ctrl01 itself,
independent of the connection, so a drop only takes the view away.

1. On the head node (hpc-ctrl01), start a named session:
     tmux new -s r40
2. Inside tmux, get the node and start the run:
     srun --partition=cpu-wide --time=24:00:00 --exclusive --exclude=hpc-db01,hpc-ctrl01 --pty bash
     cd ~/work/c++/CoherentThomson
     python3 py_scripts/compile_and_run.py <config_file>
3. Detach (the run keeps going): Ctrl-b, release, then d. You are back in the head-node shell.
   Do NOT type exit inside tmux: it closes the srun shell and ends the job.
4. Reattach later (after a drop, or to look):
     tmux attach -t r40
   Scroll inside tmux: Ctrl-b [ , arrows/PgUp, q to leave scroll mode.

RELEASE THE NODE WHEN THE RUN ENDS (unattended runs, e.g. overnight)
--------------------------------------------------------------------
With `srun ... --pty bash` the job is the shell, not the solver: when the run finishes the node stays allocated
(and accounted) until the shell ends or the time limit. Instead, make the run itself the job, inside tmux:
     srun --partition=cpu-wide --time=12:00:00 --exclusive --exclude=hpc-db01,hpc-ctrl01 \
          bash -c "cd ~/work/c++/CoherentThomson && python3 py_scripts/compile_and_run.py config/work/<my_copy>.cfg"
It builds on the node (-march=native), runs, and the job ends as soon as the command does (success or error),
after radiation_field.dat and run_log.txt are written. The output stays visible in the tmux window.
Alternatives in an interactive srun shell: end the command with `; exit`, or, while the solver runs, type `exit` +
Enter (bash runs it when the solver returns).

Checking from the head node without attaching:
  tmux ls                     # "r40: 1 windows ..." (no "(attached)"): alive and detached
  squeue -u $USER             # the job is listed while it is alive (TIME_LEFT column: -o "%.8i %.10M %.10L %R")
  srun --jobid=<jobid> --overlap -n1 ps -u $USER -o etime,%cpu,args | grep [c]oherent
                              # solver on the compute node: elapsed time, ~25000 %CPU while computing
  ls ~/output/<run>           # radiation_field.dat appears only when the run ends; run_log.txt at the start,
                              # marked "STATUS: IN PROGRESS", rewritten with the timing at the end

Still ends the run: the job's time limit, or a reboot of hpc-ctrl01.
After reconnecting, check squeue -u $USER for forgotten jobs (an orphaned interactive job keeps its whole node
idle until its time limit) and free them with scancel <jobid>.
```

## 2. Slurm: jobs are scheduled on the controller node `hpc-ctrl01` (was `slurm_head_node_report.md`)

Observed 2026-10-24, Slurm 26.05.5.

### What happened

An ordinary interactive request, with no node constraint,

```
srun --partition=debug --time=01:00:00 --exclusive --pty bash
```

was given `hpc-ctrl01`, the node that runs `slurmctld` (`SlurmctldHost = hpc-ctrl01.acc-ub.local`) and that users
log in to. It was not given one of the compute nodes, although `hpc-node[06-20]` were idle at the time:

```
JOBID  PARTITION  NAME  STATE    TIME  TIME_LIMIT  NODES  NODELIST
30255  debug      bash  RUNNING  2:39  1:00:00     1      hpc-ctrl01
```

Because of `--exclusive`, the job held all of the node's allocatable resources: `CPUAlloc=128`,
`AllocMem=1003235M`.

### Why

Every CPU partition lists the controller node and the accounting node (`AccountingStorageHost =
hpc-db01.acc-ub.local`) alongside the compute nodes:

| Partition    | MaxTime     | Nodes                                      |
|--------------|-------------|--------------------------------------------|
| cpu-long     | 12-00:00:00 | `hpc-ctrl01,hpc-db01,hpc-node[01-20]`       |
| cpu-narrow   | 5-00:00:00  | same                                       |
| cpu-standard | 3-00:00:00  | same                                       |
| cpu-medium   | 1-12:00:00  | same                                       |
| cpu-wide     | 1-00:00:00  | same                                       |
| debug        | 01:00:00    | same                                       |
| interactive  | 04:00:00    | same                                       |
| gpu          | 08:00:00    | `hpc-nodegpu01`                            |

Both nodes have the same `Weight=1` as the compute nodes, so the scheduler does not prefer the compute nodes and
can pick `hpc-ctrl01` or `hpc-db01` even when compute nodes are free.

Both nodes do reserve resources for system use, which suggests that running jobs there is intended:

```
hpc-ctrl01 / hpc-db01:  CPUTot=256  CPUEfctv=128  CoreSpecCount=64  CPUSpecList=0-127  MemSpecLimit=1048576
```

### Effects seen from the user side

- A job intended for a compute node silently gets half the CPUs (128 instead of 256), on a node that also carries
  `slurmctld`, user logins and VS Code Remote servers. The background `CPULoad` on `hpc-ctrl01` was 41–66 during
  this session, so timings there are not comparable with a compute node.
- An `--exclusive` job there blocks other jobs from using the controller's allocatable half for its whole time
  limit (up to 12 days in `cpu-long`).
- From inside the job you can't tell it apart from the login shell (same hostname, same prompt), so it looks as if
  `srun` did nothing.

### Questions / possible changes

1. Is it intended that `hpc-ctrl01` and `hpc-db01` take jobs in every CPU partition?
2. If yes: could they get a higher `Weight` (e.g. `Weight=100`), so the scheduler uses them only after the compute
   nodes are full?
3. If no: could they be removed from the partitions' `Nodes=` lists (or put in a separate partition)?

Current workaround on our side:

```
srun --partition=interactive --time=04:00:00 --exclusive --exclude=hpc-ctrl01,hpc-db01 --pty bash
```

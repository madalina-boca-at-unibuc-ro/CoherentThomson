# Slurm: jobs are scheduled on the controller node `hpc-ctrl01`

Observed 2026-10-24, Slurm 26.05.5.

## What happened

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

## Why

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

## Effects seen from the user side

- A job intended for a compute node silently gets half the CPUs (128 instead of 256), on a node that also carries
  `slurmctld`, user logins and VS Code Remote servers. The background `CPULoad` on `hpc-ctrl01` was 41–66 during
  this session, so timings there are not comparable with a compute node.
- An `--exclusive` job there blocks other jobs from using the controller's allocatable half for its whole time
  limit (up to 12 days in `cpu-long`).
- From inside the job you can't tell it apart from the login shell (same hostname, same prompt), so it looks as if
  `srun` did nothing.

## Questions / possible changes

1. Is it intended that `hpc-ctrl01` and `hpc-db01` take jobs in every CPU partition?
2. If yes: could they get a higher `Weight` (e.g. `Weight=100`), so the scheduler uses them only after the compute
   nodes are full?
3. If no: could they be removed from the partitions' `Nodes=` lists (or put in a separate partition)?

Current workaround on our side:

```
srun --partition=interactive --time=04:00:00 --exclusive --exclude=hpc-ctrl01,hpc-db01 --pty bash
```

# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`Flow-Otter` is a Python/Julia workflow framework for large-scale materials-simulation pipelines. It chains molecular dynamics (LAMMPS) → electronic Hamiltonian construction (Hamster ML or empirical tight-binding) → optoelectronic properties (band gap/DOS, COHP, PDOS, conductivity). Orchestration is done with [PerQueue](https://gitlab.com/asm-dtu/perqueue), which sits on top of the `myqueue` job scheduler and requires a SLURM/PBS/LSF cluster to actually run jobs.

A vendored copy of PerQueue lives in `src/perqueue/` (installed editable via `requirements/requirements_flow_otter.txt`).

## Running

There is no build step and the test files under `tasks/*/test/` are pipeline *stages*, not a unit-test suite — there is no `pytest`/CI harness to invoke locally.

```bash
# one-time per working directory
pq init

# launch a workflow with a config file (pick the entrypoint, see below)
python /path/to/flow.py /path/to/otter.yaml

# inspect / drive the running workflow
pq ls                       # list tasks + status
pq resubmit -i ID           # resubmit by job id (also -m myqueue-id, -n name, -s status)
pq modify r -i ID           # change resource string of a job
```

Almost every parameter is re-read from the YAML config at task runtime, so editing the config and `pq resubmit`-ing a task picks up the change. Exceptions (baked in by PerQueue at submit time): all `resources_*` strings, `num_simulations`, the config path itself, and the generated `branch_config.yaml` files. To re-edit branch configs without restarting everything, resubmit the first task `check_simulation_type`; set `just_update_branch_config_file: true` to update branch configs and stop.

## Entrypoints (top-level `flow_*.py`)

Each top-level script builds one PerQueue `Workflow` graph and submits it. They differ only in which stages they wire together:
- `flow.py` — full pipeline: check → MD → Hamiltonian → optoelectronic property (+ plotting).
- `flow_MD.py` — MD only.
- `flow_MD_DFT.py` — MD + DFT (`tasks/DFT/`).

## Architecture

**Workflow = graph of Tasks built from PerQueue group primitives.** The entrypoint files assemble the graph; the logic lives in per-stage task scripts. Key PerQueue primitives:
- `Task(script_path, kwargs_or_None, resources_string, name=...)` — one job; runs `script_path`'s `main()`.
- `SwitchGroup({option: task_or_subgraph, ...})` — picks a branch at runtime based on the `SWITCHGROUP_KEY` returned by the preceding task. This is how config options like `MD_type`, `H_type`, `optoelec_type` select which path runs.
- `StaticWidthGroup(subgraph, width=N)` — runs N copies in parallel ("sweep").
- `CyclicalGroup(subgraph, N)` — runs N copies sequentially, each restarting from the previous ("cascade").

**Task contract.** Every task script exposes:
```python
def main(path_configWF="workflow_config.yaml", num_simulations=1, **kwargs) -> Tuple[bool, dict]:
    ...
    return True, {"path_configWF": ..., "num_simulations": ..., SWITCHGROUP_KEY: <branch>}
```
The returned dict is passed as kwargs into dependent tasks — this is how state (config path, branch selection, computed values) flows down the graph. Returning `SWITCHGROUP_KEY` (from `perqueue.constants`) selects the next `SwitchGroup` branch. PerQueue injects `pq_index` (sweep position) and `pq_iteration` (cascade iteration) via `kwargs`; under multi-simulation these arrive as dicts/lists, so tasks unwrap them (see `tasks/MD/check_MD_type.py`).

**`tasks/` vs `codes/`.** `tasks/` holds the orchestration glue — `check_*` (decide branch), `prep_*` (write inputs), `run_*`/`<engine>_*` (launch the external code), `skip_*` (no-op branch), `*_post*`/`postprocessing_*` (parse outputs), and `test/` stages that sanity-check a stage's results. `codes/` holds the actual compute kernels (Python and Julia, e.g. `codes/MD/`, `codes/optoelec/`) that the tasks invoke. Pipeline stage prefixes: `check_simulation_type` (0) → MD (1) → H (2) → optoelec (3).

**Multi-simulation branching.** When `num_simulations > 1`, `check_simulation_type.py` writes one `branch_config.yaml` per branch into `dir_project/<param_to_vary>_<value>/`, varying up to six `param_to_vary*` parameters across branches (optionally nested under a `param_group_for_vary*`). A `pre_branch_config.yaml`, if present, is deep-merged in. Downstream tasks locate their branch's config from `pq_index`/`pq_iteration`.

**Config is the API.** Behavior is driven almost entirely by the YAML config (`examples/*.yaml`). `docs/parameters.md` is the authoritative parameter reference (global params + per-section groups `lammps:`, `gap+dos:`, `cohp:`, `pdos:`, `conductivity:`). When adding a parameter, read it via `configWF.get(...)` in the relevant task and document it in `docs/parameters.md`.

## External software & Julia

External engines are optional and selected by config: LAMMPS (+ MACE / VASP ML force fields) for MD, Hamster + Vampires (Julia) for ML Hamiltonians, VASP for DFT. Julia deps come from `requirements/Project.toml`/`Manifest.toml` or `requirements/install_julia_requirements.jl`; optional packages (Vampires, Hamster) are added via `requirements/add_packages/add_*.jl` after the external code is installed. `preambles/preamble_*.sh` load environment-specific modules per stage. Because engines often need incompatible Python/Julia environments, a failed task can be safely resubmitted from the correct environment.

See `docs/setup.md` for full setup and `examples/examples.md` for worked configurations.
## Cluster information
The HPC cluster that I am working are accessible under the ssh aliases fau_fritz and fau_alex, which are cpu and gpu clusters, respectively. My home directory on these clusters is /home/hpc/b299bb/b299bb25. The folder where our group's system administrators have stored important codes are /home/atuin/b299bb/TheoFEM/. On the cluster, i have a version of flow-otter that you can rsync to. You have to run module use /home/atuin/b299bb/TheoFEM/Modules beforehand.

The password for the ssh connection is stored in the secret.txt file. I know this is unprofessional, but I added it to the .gitignore file.

### VASP on alex
A working VASP environment is set up via `/home/hpc/b299bb/b299bb25/flow-otter/load_env.sh` — `source` it to put `vasp_std`/`vasp_gam`/`vasp_ncl` on `$PATH` (used by `tasks/DFT/run_DFT.py`, which runs `srun vasp_std`). Notes:
- The group VASP modulefiles under `/home/atuin/b299bb/TheoFEM/Modules/VASP` currently fail to `module load` — they pin dependency versions no longer on alex (e.g. `mkl/2023.2.0`, `intel/2023.2.1`, `cuda/12.1.1`). `load_env.sh` tries the `VASP/6.6.0_intel` module first, then falls back to loading `intel/2024.2.1 intelmpi/2021.17.0 mkl/2024.2.2` and adding `/home/atuin/b299bb/TheoFEM/Codes/VASP/vasp.6.6.0_intel/bin` to `$PATH`.
- This is the **intel CPU build** (all libs resolve and `vasp_std` runs). The GPU `nvhpc_erlangen` build is not usable yet — its `bin/` is empty (not compiled).

## DFT pipeline (`flow_MD_DFT.py`, `tasks/DFT/`)

Graph: `check_simulation_type` → MD switch (`sweep`/`cascade`, currently `skip_MD` for the CuTaN2 runs) → `swg1_DFT` → `postprocessing_DFT`.

- `prep_DFT.py` — obtains the trajectory, then runs `vamp supercell sample --N N_snapshots` to create one `config_*` snapshot dir each, then writes `INCAR` (from the `vasp:` config group), `KPOINTS` (from `kpoints:`) and copies the `potcar_file` into every snapshot dir. Two ways to get the trajectory into `<branch>/2-DFT/XDATCAR`: (1) if `trajectory_file` is set, that external VASP XDATCAR is copied in (the original CuTaN2 workflow); (2) otherwise the **on-the-fly** path runs — `lammps_trajectory_to_xdatcar()` reads the branch's MD dump (`1-MD/position.lammpstrj`, name overridable via `md_trajectory_file`) with ASE, regroups atoms by species in `lammps: elements` order (must match POTCAR concatenation order), and writes the XDATCAR. So `flow_MD_DFT.py` with `MD_type: lammps` runs MD per branch and feeds it straight into DFT; `examples/otter_MD_DFT_LJ.yaml` is a worked end-to-end test using the LJ placeholder potential. (`flow_MD_DFT.py`'s graph already had both MD switch branches — only `prep_DFT.py` needed the lammps-trajectory support.)
- `run_DFT.py` — runs **one** snapshot. It is wrapped in a nested `StaticWidthGroup(width=N_snapshots)` inside the per-branch `StaticWidthGroup(width=num_simulations)`, so snapshots run as independent parallel jobs (total `num_simulations × N_snapshots`) instead of serially in one job. PerQueue delivers `pq_index = [snapshot_idx, branch_idx]`; the snapshot dir is picked from `sorted_dft_dirs()` (numeric sort on the `config_*` suffix). Each job: `srun vasp_std` → `vamp eigenval read --par bandgap` (`bandgap.log`) → `vamp doscar read` (`dos.h5`).
- `postprocessing_DFT.py` — parses `bandgap.log` from every snapshot into `postprocessing_results.json`, and writes band-gap + DOS plots to `<dir_project>/analysis/` (`bandgap_distribution.png`, `bandgap_vs_parameter.png`, `dos_per_branch.png`, `dos_overlay.png`). Notes:
  - Band gaps are pulled with a regex (`bandgap is:\s*(...)`) — the log line ends in a period (`...is: 0.69133.`), so the previous `float(line.split()[-1])` raised `ValueError` and silently dropped every snapshot (`postprocessing_results.json` came out `[]`).
  - DOS is read straight from each snapshot's `DOSCAR` (full spectrum + Fermi energy), **not** from `dos.h5`: the `dos.h5` written by `vamp doscar read` is broken — it stores only the first three rows of the DOSCAR, mislabelled as `energy`/`total_dos`/`integrated_dos`. DOS curves are shifted so E_F = 0 and ensemble-averaged per branch (branch label/sort-key derived from the `param_to_vary*` values).
  - Plotting needs `matplotlib`/`numpy` in the task's environment; if absent it is skipped (prints `Plotting skipped...`) and the rest of postprocessing still succeeds.
  - There is also a standalone `analyze_dft_results.py` at the repo root that does the same analysis offline against an rsync'd copy of the results (defaults to `./dft_analysis_data`).

INCAR/KPOINTS guidance for the large CuTaN2 supercell (~232 atoms): `NCORE` must divide `cores/KPAR` (e.g. 48 cores → `KPAR 4 / NCORE 4`); set `LWAVE false` for gap/DOS runs to avoid huge WAVECAR I/O; use `type: Gamma` so the zone centre is sampled; `LREAL Auto` for the large cell. Expect ~1–3 h per snapshot SCF on ~48 CPU cores. `resources_DFT` walltime now applies per snapshot.

### Band-gap accuracy (PBE underestimation)

Plain PBE (`GGA: PE`, the setup in `examples/otter_DFT.yaml`) consistently returns CuTaN2 band gaps that are **too small** — the usual semilocal-GGA gap underestimation, dominated here by self-interaction error on the localized Cu 3d / Ta 5d states. `examples/otter_DFT_2.yaml` addresses this with a rotationally-invariant **Dudarev DFT+U** (`LDAUTYPE: 2`) on the metal d-states (Cu `U_eff ≈ 4 eV`, Ta `≈ 2 eV`), at negligible extra cost over PBE. It is otherwise identical to `otter_DFT.yaml` (same k-points/ENCUT/snapshots/temperatures) and writes to a separate `dir_project` (`flow-otter-test-DFTU/`) for a clean PBE-vs-PBE+U comparison.

- **`LDAU*` arrays are per-species, ordered to match the POTCAR concatenation order.** For `/home/hpc/b299bb/b299bb25/POTCAR` that order is **Ta, Cu, N** (verified from the POTCAR `TITEL` lines and the XDATCAR header `Ta Cu N` / `48 48 96`) — so `LDAUL: [2, 2, -1]`, `LDAUU: [2.0, 4.0, 0.0]`. Re-order if the POTCAR is ever rebuilt. `LMAXMIX: 4` is required for correct l=2 occupancy mixing.
- The `U` values are tunable knobs, not constants — if the gap is still off, scan `LDAUU` across branches (the sweep machinery) to calibrate against a target gap.
- For gold-standard gaps a hybrid (HSE06) is the next step, but it is ~50–100× PBE cost — likely infeasible for ~232 atoms × 450 snapshots within the 2 h/snapshot walltime; reserve it for a few representative snapshots.

Because the INCAR is written generically by `write_incar()` (every `vasp:` key → `TAG = value`, booleans → `.TRUE.`/`.FALSE.`, lists → space-separated), switching functionals/methods (+U, meta-GGA, hybrids) is purely a YAML edit — no task-code changes needed.

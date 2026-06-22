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

⚠️ `flow_MD_DFT.py` currently runs `os.system("scancel -u $USER")` at import and has a `delete_everything_in_dir` helper — both marked `# TODO: REMOVE`. Do not treat these as intended behavior; flag/remove before relying on this entrypoint.

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

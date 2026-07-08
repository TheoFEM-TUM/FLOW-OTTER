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

### LAMMPS+MACE on alex

The group's Apptainer-based LAMMPS+MACE (KOKKOS GPU) environment (`lammps_modules/lammps-mace_nhvhpc`, backed by `/home/atuin/b299bb/TheoFEM/Environments/mace_lammps-apptainer`) is broken — its `.def` file never actually compiles anything at build time (compilation was deferred to first container run and apparently never completed). We built our own LAMMPS+MACE from source instead.

**Build.** `~/build_lammps_mace.sbatch` on alex. Key points, each the fix for a real failure hit along the way:
- Must clone **`https://github.com/ACEsuit/lammps.git` branch `mace`** (the ACEsuit fork with native `pair_style mace`/`mace/kk` in `src/ML-MACE` + `src/KOKKOS/pair_mace_kokkos.cpp`), *not* mainline `lammps/lammps`. Enable it with `-D PKG_ML-MACE=ON` (not `PKG_ML-MLIAP`, which isn't a real package name and is silently ignored by CMake).
- Needs LibTorch downloaded and passed via `-D CMAKE_PREFIX_PATH=<libtorch dir>` (`find_package(Torch REQUIRED)` in `cmake/Modules/Packages/ML-MACE.cmake`). LibTorch 2.7.0 cu128 (`https://download.pytorch.org/libtorch/cu128/...`) works fine even though the system CUDA module is 13.2 — LibTorch bundles its own CUDA runtime; only the driver needs to be new enough.
- Also needs `module load mkl/2024.2.2` — LibTorch's own `Caffe2Config.cmake` calls `find_package(MKL)` and bakes the result into the `torch` imported target's include dirs even though LAMMPS+MACE doesn't use MKL itself; without it, configure fails with `MKL_INCLUDE_DIR-NOTFOUND`.
- Also needs a small injected CMake shim (`-D CMAKE_PROJECT_INCLUDE=<path>/nvtx_shim.cmake` defining a no-op `CUDA::nvToolsExt` interface target) because CUDA 13.2 ships neither the legacy NVTX lib nor a discoverable `nvtx3` header set that LibTorch's `cuda.cmake` expects; NVTX is profiling-only and irrelevant here.
- `sbatch` runs a non-login shell, so `module` is undefined unless you `source /etc/profile.d/zz-rrze-local.sh` first.
- The a100 **compute nodes cannot reach `github.com`** (times out; confirmed via a dedicated diagnostic job), even though the login node can and other CDNs (e.g. `download.pytorch.org`, backed by CloudFront) work fine from compute nodes. Fix: `git clone` the LAMMPS repo from the **login node** first, then let the batch job build from the already-staged copy.
- `/home/atuin` is **group-shared** storage with a group-wide inode/file-count quota (was 492k/500k when hit) — a full LAMMPS clone (~13k files) can blow past it. Build under `$HOME` instead (per-user quota, plenty of headroom); the group's other big consumers include the broken Apptainer sandbox's ~67k extracted files.
- The final `lmp` binary's RPATH hardcodes an absolute path to LibTorch's `lib/` dir (not `$ORIGIN`-relative) — safe to copy the binary anywhere as long as the LibTorch directory it points at isn't moved.

Working binary: `/home/atuin/b299bb/b299bb25/lammps-mace/lmp` (LibTorch stays under `$HOME/software/lammps_mace_build/`). Loadable via `module load lammps-mace` (custom modulefile at `~/modulefiles/lammps-mace/1.0`, wired into `~/enable_term_alex.sh`) — loads `openmpi/5.0.8-nvhpc26.3-cuda13.2` + `mkl/2024.2.2` and puts the binary on `PATH`. `~/enable_term_alex.sh` (sourced by `preambles/preamble_alex_lammps_mace.sh`, the MD stage's PerQueue `activation_script`) loads `python/3.12-base` → `julia` → `lammps-mace` → activates `.venv-alex` → sources `~/.bashrc`; re-verified this order does not conflict (`lammps-mace`'s nested `openmpi`/`mkl` loads succeed fine after `julia`). Two loose ends noted but not yet shown to matter: `~/.bashrc` unconditionally re-exports `THEOFEM_DIR` to a different path (`.../b299bb11/TheoFEM`, not `.../TheoFEM/Modules`) after the preamble already used the correct value for `module use`, and its final `source .bashrc` uses a path relative to the job's cwd (`~/flow-otter`, not `$HOME`) so it may silently no-op there.

**Fixed bug: `libnvrtc.so.12` dlopen failure (TorchScript JIT), distinct from the Kokkos crash below.** Symptom: `lammps_MD.py`'s `srun lmp -k on g 1 -sf kk -pk kokkos ...` dies immediately (before any force evaluation, before `log.lammps` gets a single line) with `RuntimeError: Error in dlopen for library libnvrtc.so.12and libnvrtc-XXXXXXXX.so.12` from inside the TorchScript interpreter. Root cause: `$HOME/software/lammps_mace_build/libtorch/lib/` only contained the hash-named `libnvrtc-c8f577df.so.12` — the plain-named `libnvrtc.so.12` symlink that a normal `nvidia-cuda-nvrtc-cu12` pip wheel ships (and that LibTorch's CUDA JIT fuser `dlopen()`s by exact fixed name at runtime) was missing, apparently lost when this LibTorch 2.7.0/cu128 archive was manually downloaded/extracted rather than pip-installed. Fix: `ln -s libnvrtc-c8f577df.so.12 libnvrtc.so.12` in that directory. Confirmed this gets a run past the dlopen error (it then proceeds to actually launch LAMMPS).

Separately, one resubmitted job (after the above fix) failed at the `srun lmp ...` step itself with `execve(): lmp: No such file or directory`, even though `lmp` resolves fine on `PATH` both interactively and via a fresh `source ~/enable_term_alex.sh`, and two earlier jobs with the identical preamble had successfully found and executed the same binary on the same node. Not reproduced on retry-by-hand; treat as a possibly transient/flaky compute-node issue (e.g. an `/home/atuin` NFS-mount race right after node (re)allocation) rather than a real preamble bug, unless it recurs.

**Model export gotcha.** MACE's `mace_create_lammps_model` CLI has two mutually incompatible output modes: default (`--format libtorch`) actually JIT-compiles the model into a real TorchScript archive for `pair_style mace`/`mace/kk`; `--format mliap` just pickles it (`torch.save`, no JIT) for LAMMPS's separate, unrelated `pair_style mliap unified` (mainline ML-IAP) interface. Pointing `pair_style mace` at an mliap-format export fails with `PytorchStreamReader failed locating file constants.pkl: file not found` — this is exactly the bug that cost significant debugging time when a supervisor-provided `model.pt` turned out to be mliap-format. Use `codes/MD/export_mace_lammps_model.sh` (wraps the CLI correctly) and see `docs/parameters.md` under `path_FF_MD`.

**Open bug: `pair_style mace/kk` (Kokkos GPU) crashes on every run.** Minimal repro: a plain Cu FCC box (any size tested: 10.8/14.5/32.5 Å), `pair_style mace no_domain_decomposition`, run via `lmp -k on g 1 -sf kk -pk kokkos newton on neigh half` (the exact flags `tasks/MD/lammps_MD.py` uses). The model loads and reports itself correctly (`r_max`, layer count, species mapping all sane), then the very first force evaluation crashes with:
  ```
  /pytorch/aten/src/ATen/native/cuda/ScatterGatherKernel.cu:144: operator(): ... Assertion `idx_dim >= 0 && idx_dim < index_size && "index out of bounds"` failed.
  cudaFree(arg_alloc_ptr) error( cudaErrorAssert): device-side assert triggered .../lib/kokkos/core/src/Cuda/Kokkos_CudaSpace.cpp:379
  ```
  Ruled out as the cause, each independently confirmed:
  - **Box size** — crashes identically from 10.8 Å up to 32.5 Å (well past `2 × ghost_atom_cutoff`, which LAMMPS itself reports as 14 Å for this 2-layer model — not just `r_max`, since multi-layer message passing extends the effective interaction range).
  - **Foundation model checkpoint** — identical crash with both the original Dec-2023 `2023-12-10-mace-128-L0_energy_epoch-249.model` and the current default `mace-mpa-0-medium.model`.
  - **Precision** — identical crash (same block/thread numbers) with both `--dtype float64` (default) and `--dtype float32`.

  This points to a structural bug in `src/KOKKOS/pair_mace_kokkos.cpp` for this exact build (LAMMPS `mace` branch @ `4d222cb`, mace-torch 0.3.16, LibTorch 2.7.0/cu128) — likely a fixed-size buffer or index mapping issue in the Kokkos neighbor-list marshaling, independent of model content. That commit is the *only* commit that's ever touched `pair_mace_kokkos.cpp` on the `mace` branch, so there's no newer upstream fix to pull in.

  **Plain (non-Kokkos) `pair_style mace`** (no `-sf kk`/`-pk kokkos`, still GPU-capable via LibTorch's own CUDA dispatch) does *not* crash, but hangs indefinitely (10+ min, no output past module loading — stuck in setup, not just slow) — not yet root-caused.

  Untried next steps: an older LibTorch version (the Kokkos glue code may predate LibTorch 2.7's ABI/kernel changes), or filing an upstream issue with this reproduction.

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

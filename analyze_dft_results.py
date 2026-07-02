#!/usr/bin/env python3
"""Analyze Flow-Otter DFT results (band gaps + DOS) and make plots.

The DFT stage of the pipeline produces, for every MD snapshot, a directory

    <results>/<branch>/2-DFT/config_<n>/

containing a ``bandgap.log`` (parsed band gap) and a VASP ``DOSCAR`` (total
density of states).  ``<branch>`` is whatever ``check_simulation_type`` named
the branch -- for the CuTaN2 runs this is ``temperature_<T>``.

This script walks that tree, collects the band gaps and densities of states per
branch, and writes a few summary plots plus a CSV of the band gaps.

It also reads each snapshot's ``POSCAR`` (if present and ASE is installed) and
computes first-shell interatomic-distance descriptors -- mean/disorder of the
Cu-N and Ta-N bonds -- then correlates them with the band gap (pooled and
within each temperature) to quantify how thermal bond-length fluctuations move
the gap.  See ``gap_vs_distance.png`` / ``gap_distance_correlation.png``.

Note on ``dos.h5``: the file written by ``vamp doscar read`` in this dataset is
broken -- it stores only the first three rows of the DOSCAR, mislabelled as
``energy``/``total_dos``/``integrated_dos``.  We therefore read the full
spectrum straight from ``DOSCAR`` instead.

Usage
-----
    python analyze_dft_results.py [RESULTS_DIR] [-o OUTDIR]

``RESULTS_DIR`` defaults to ``./dft_analysis_data`` (the rsync'd copy of the
cluster results).  Plots and the CSV are written to ``OUTDIR``
(default ``<RESULTS_DIR>/analysis``).
"""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless: just write files
import matplotlib.pyplot as plt
import numpy as np

try:  # structural analysis is optional; degrade gracefully if ASE is absent
    from ase.io import read as _ase_read
except Exception:  # pragma: no cover
    _ase_read = None

try:  # scipy gives us Spearman + clean p-values; fall back to numpy otherwise
    from scipy import stats as _scipy_stats
except Exception:  # pragma: no cover
    _scipy_stats = None

# matches e.g. "The value for bandgap is: 0.69133."  (trailing period is fine)
_BANDGAP_RE = re.compile(r"bandgap is:\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)")

# --------------------------------------------------------------------------- #
# Interatomic-distance descriptors
# --------------------------------------------------------------------------- #
# CuTaN2 first coordination shell (verified from the POSCARs): every Cu has two
# short Cu–N bonds (~1.84 A, linear Cu(I)), every Ta has six Ta–N bonds (~2.0 A,
# octahedral), with a clear gap before the next shell.  We therefore define the
# bonds by a *fixed coordination number* (the N nearest N atoms of each metal)
# rather than a distance cutoff, which makes the descriptors insensitive to how
# much thermal disorder stretches the bonds.
COORDINATION = {"Cu": ("N", 2), "Ta": ("N", 6)}  # metal -> (ligand, CN)

# (key, human-readable label) for every descriptor we compute per snapshot.
DESCRIPTORS = [
    ("d_CuN_mean", "mean Cu–N bond (Å)"),
    ("d_CuN_std", "Cu–N bond disorder σ (Å)"),
    ("d_CuN_min", "shortest Cu–N bond (Å)"),
    ("d_TaN_mean", "mean Ta–N bond (Å)"),
    ("d_TaN_std", "Ta–N bond disorder σ (Å)"),
    ("d_TaN_min", "shortest Ta–N bond (Å)"),
    ("d_metalN_mean", "mean metal–N bond (Å)"),
    ("d_metalN_std", "metal–N bond disorder σ (Å)"),
]


@dataclass
class Snapshot:
    branch: str
    name: str  # e.g. "config_1"
    index: int  # numeric suffix for sorting
    bandgap: float | None = None
    energy: np.ndarray | None = None  # DOS energy grid, shifted so E_Fermi = 0
    total_dos: np.ndarray | None = None
    e_fermi: float | None = None
    desc: dict = field(default_factory=dict)  # interatomic-distance descriptors


@dataclass
class Branch:
    name: str
    label: str  # human-friendly (e.g. "300 K" if a temperature branch)
    sort_key: float
    snapshots: list[Snapshot] = field(default_factory=list)

    @property
    def gaps(self) -> np.ndarray:
        return np.array([s.bandgap for s in self.snapshots if s.bandgap is not None])


# --------------------------------------------------------------------------- #
# Parsing
# --------------------------------------------------------------------------- #
def parse_bandgap(log_path: Path) -> float | None:
    """Extract the band gap (eV) from a bandgap.log, robust to trailing punctuation."""
    try:
        text = log_path.read_text()
    except OSError:
        return None
    matches = _BANDGAP_RE.findall(text)
    if not matches:
        return None
    try:
        return float(matches[-1])
    except ValueError:
        return None


def parse_doscar(doscar_path: Path):
    """Return (energy, total_dos, e_fermi) from a VASP DOSCAR (total DOS only).

    Handles both spin-unpolarised (3 columns: E, DOS, intDOS) and
    spin-polarised (5 columns: E, DOS_up, DOS_dn, intDOS_up, intDOS_dn) files;
    in the spin case the two channels are summed into a single total DOS.
    Returns (None, None, None) on any parsing problem.
    """
    try:
        lines = doscar_path.read_text().splitlines()
    except OSError:
        return None, None, None
    if len(lines) < 7:
        return None, None, None

    header = lines[5].split()
    try:
        nedos = int(float(header[2]))
        e_fermi = float(header[3])
    except (IndexError, ValueError):
        return None, None, None

    energy, total = [], []
    for line in lines[6 : 6 + nedos]:
        cols = line.split()
        if len(cols) < 2:
            continue
        try:
            vals = [float(c) for c in cols]
        except ValueError:
            continue
        energy.append(vals[0])
        if len(vals) >= 5:  # spin-polarised: up + down
            total.append(vals[1] + vals[2])
        else:
            total.append(vals[1])

    if not energy:
        return None, None, None
    return np.array(energy), np.array(total), e_fermi


def compute_descriptors(poscar_path: Path) -> dict:
    """Per-snapshot interatomic-distance descriptors from a VASP POSCAR.

    For each metal we take the ``CN`` nearest N atoms (minimum-image distances,
    so the periodic supercell is handled correctly) and summarise the resulting
    bond lengths.  Returns ``{}`` if ASE is unavailable or the file can't be read.
    """
    if _ase_read is None:
        return {}
    try:
        atoms = _ase_read(str(poscar_path), format="vasp")
    except Exception:
        return {}

    symbols = np.array(atoms.get_chemical_symbols())
    dmat = atoms.get_all_distances(mic=True)  # minimum-image, general cell

    desc: dict = {}
    all_bonds = []
    for metal, (ligand, cn) in COORDINATION.items():
        im = np.where(symbols == metal)[0]
        il = np.where(symbols == ligand)[0]
        if len(im) == 0 or len(il) < cn:
            continue
        sub = np.sort(dmat[np.ix_(im, il)], axis=1)  # sort ligands by distance
        bonds = sub[:, :cn].ravel()  # the cn nearest ligands of every metal
        desc[f"d_{metal}{ligand}_mean"] = float(bonds.mean())
        desc[f"d_{metal}{ligand}_std"] = float(bonds.std())
        desc[f"d_{metal}{ligand}_min"] = float(bonds.min())
        all_bonds.append(bonds)

    if all_bonds:
        allb = np.concatenate(all_bonds)
        desc["d_metalN_mean"] = float(allb.mean())
        desc["d_metalN_std"] = float(allb.std())
    return desc


def _branch_label(name: str) -> tuple[str, float]:
    """Pretty label + numeric sort key for a branch directory name."""
    m = re.match(r"temperature_([-+]?\d*\.?\d+)", name)
    if m:
        t = float(m.group(1))
        return f"{t:g} K", t
    m = re.search(r"([-+]?\d*\.?\d+)", name)
    return name, (float(m.group(1)) if m else 0.0)


def collect(results_dir: Path) -> list[Branch]:
    """Walk RESULTS_DIR/<branch>/2-DFT/config_*/ and gather all snapshots."""
    branches: list[Branch] = []
    branch_dirs = sorted(
        d for d in results_dir.iterdir() if d.is_dir() and (d / "2-DFT").is_dir()
    )
    # also support a flat single-simulation layout (config_* directly under 2-DFT)
    if not branch_dirs and (results_dir / "2-DFT").is_dir():
        branch_dirs = [results_dir]

    for bdir in branch_dirs:
        label, sort_key = _branch_label(bdir.name)
        branch = Branch(name=bdir.name, label=label, sort_key=sort_key)
        config_dirs = sorted(
            (bdir / "2-DFT").glob("config_*"),
            key=lambda p: int(re.search(r"(\d+)$", p.name).group(1)),
        )
        for cdir in config_dirs:
            idx = int(re.search(r"(\d+)$", cdir.name).group(1))
            snap = Snapshot(branch=bdir.name, name=cdir.name, index=idx)
            snap.bandgap = parse_bandgap(cdir / "bandgap.log")
            snap.desc = compute_descriptors(cdir / "POSCAR")
            energy, total, ef = parse_doscar(cdir / "DOSCAR")
            if energy is not None:
                snap.energy = energy - ef  # align to Fermi level (E - E_F)
                snap.total_dos = total
                snap.e_fermi = ef
            branch.snapshots.append(snap)
        branches.append(branch)

    branches.sort(key=lambda b: b.sort_key)
    return branches


# --------------------------------------------------------------------------- #
# Plotting
# --------------------------------------------------------------------------- #
def plot_bandgap_distribution(branches: list[Branch], out: Path) -> None:
    """Box + jittered-scatter of band gaps per branch."""
    data = [b.gaps for b in branches]
    labels = [b.label for b in branches]
    rng = np.random.default_rng(0)

    fig, ax = plt.subplots(figsize=(7, 5))
    # 'tick_labels' (mpl>=3.9) replaced 'labels'; fall back for older mpl
    try:
        ax.boxplot(data, tick_labels=labels, showfliers=False, widths=0.5,
                   medianprops=dict(color="black"))
    except TypeError:
        ax.boxplot(data, labels=labels, showfliers=False, widths=0.5,
                   medianprops=dict(color="black"))
    for i, gaps in enumerate(data, start=1):
        if len(gaps) == 0:
            continue
        x = rng.normal(i, 0.05, size=len(gaps))
        ax.scatter(x, gaps, alpha=0.7, color="tab:blue", zorder=3, s=30)
    ax.set_xlabel("MD temperature")
    ax.set_ylabel("Band gap (eV)")
    ax.set_title("DFT band gap distribution per snapshot ensemble")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def plot_bandgap_vs_temperature(branches: list[Branch], out: Path) -> None:
    """Mean band gap +/- std vs branch sort key (temperature)."""
    xs, means, stds = [], [], []
    for b in branches:
        g = b.gaps
        if len(g) == 0:
            continue
        xs.append(b.sort_key)
        means.append(g.mean())
        stds.append(g.std(ddof=1) if len(g) > 1 else 0.0)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.errorbar(xs, means, yerr=stds, marker="o", capsize=4,
                color="tab:red", lw=1.5)
    ax.set_xlabel("MD temperature (K)")
    ax.set_ylabel("Mean band gap (eV)")
    ax.set_title("Mean band gap vs. temperature (error bars = std over snapshots)")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def plot_dos(branches: list[Branch], out: Path) -> None:
    """One panel per branch: individual snapshot DOS (thin) + ensemble mean (thick)."""
    branches_with_dos = [
        b for b in branches if any(s.total_dos is not None for s in b.snapshots)
    ]
    if not branches_with_dos:
        return

    n = len(branches_with_dos)
    fig, axes = plt.subplots(n, 1, figsize=(8, 3 * n), sharex=True, squeeze=False)
    axes = axes[:, 0]

    for ax, b in zip(axes, branches_with_dos):
        snaps = [s for s in b.snapshots if s.total_dos is not None]
        # interpolate every snapshot onto a common (E - E_F) grid before averaging
        emin = max(s.energy.min() for s in snaps)
        emax = min(s.energy.max() for s in snaps)
        grid = np.linspace(emin, emax, 600)
        stack = []
        for s in snaps:
            interp = np.interp(grid, s.energy, s.total_dos)
            ax.plot(grid, interp, color="tab:gray", alpha=0.25, lw=0.8)
            stack.append(interp)
        stack = np.array(stack)
        mean = stack.mean(axis=0)
        ax.plot(grid, mean, color="tab:blue", lw=2, label="ensemble mean")
        ax.fill_between(grid, mean - stack.std(axis=0), mean + stack.std(axis=0),
                        color="tab:blue", alpha=0.2, label=r"$\pm$ std")
        ax.axvline(0.0, color="k", ls="--", lw=1, label=r"$E_\mathrm{F}$")
        ax.set_ylabel("DOS (states/eV)")
        ax.set_title(f"Total DOS — {b.label}  (n={len(snaps)})")
        ax.legend(loc="upper left", fontsize=8)
        ax.grid(alpha=0.3)

    axes[-1].set_xlabel(r"$E - E_\mathrm{F}$ (eV)")
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def plot_dos_overlay(branches: list[Branch], out: Path) -> None:
    """Ensemble-mean DOS of every branch overlaid in one axes for comparison."""
    branches_with_dos = [
        b for b in branches if any(s.total_dos is not None for s in b.snapshots)
    ]
    if not branches_with_dos:
        return

    fig, ax = plt.subplots(figsize=(8, 5))
    cmap = plt.get_cmap("viridis")
    for i, b in enumerate(branches_with_dos):
        snaps = [s for s in b.snapshots if s.total_dos is not None]
        emin = max(s.energy.min() for s in snaps)
        emax = min(s.energy.max() for s in snaps)
        grid = np.linspace(emin, emax, 600)
        stack = np.array([np.interp(grid, s.energy, s.total_dos) for s in snaps])
        color = cmap(i / max(1, len(branches_with_dos) - 1))
        ax.plot(grid, stack.mean(axis=0), color=color, lw=2, label=b.label)

    ax.axvline(0.0, color="k", ls="--", lw=1)
    ax.set_xlabel(r"$E - E_\mathrm{F}$ (eV)")
    ax.set_ylabel("Mean DOS (states/eV)")
    ax.set_title("Ensemble-mean total DOS by temperature")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Band gap vs. interatomic distance
# --------------------------------------------------------------------------- #
def _pearson(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Pearson r and two-sided p-value (p is NaN if scipy is unavailable)."""
    if len(x) < 3:
        return float("nan"), float("nan")
    if _scipy_stats is not None:
        r, p = _scipy_stats.pearsonr(x, y)
        return float(r), float(p)
    return float(np.corrcoef(x, y)[0, 1]), float("nan")


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    """Spearman rank correlation (rank with numpy if scipy is absent)."""
    if len(x) < 3:
        return float("nan")
    if _scipy_stats is not None:
        return float(_scipy_stats.spearmanr(x, y).correlation)
    rx = np.argsort(np.argsort(x))
    ry = np.argsort(np.argsort(y))
    return float(np.corrcoef(rx, ry)[0, 1])


def _paired(branches: list[Branch], key: str) -> tuple[np.ndarray, np.ndarray]:
    """All (descriptor, band gap) pairs across every branch for one descriptor."""
    xs, ys = [], []
    for b in branches:
        for s in b.snapshots:
            if s.bandgap is not None and key in s.desc:
                xs.append(s.desc[key])
                ys.append(s.bandgap)
    return np.array(xs), np.array(ys)


def correlation_stats(branches: list[Branch]) -> list[dict]:
    """For every descriptor: pooled Pearson/Spearman + slope, and per-branch r.

    The pooled correlation mixes the temperature trend (hotter -> longer, more
    disordered bonds *and* smaller gap) with the structural sensitivity, so we
    also report the within-branch Pearson r at each fixed temperature, which
    isolates "how does the gap respond to a bond fluctuation at constant T".
    """
    rows = []
    for key, label in DESCRIPTORS:
        x, y = _paired(branches, key)
        if len(x) < 3:
            continue
        r, p = _pearson(x, y)
        rho = _spearman(x, y)
        slope = float(np.polyfit(x, y, 1)[0])  # dE_gap / d(descriptor), eV/Å
        per_branch = {}
        for b in branches:
            bx = np.array([s.desc[key] for s in b.snapshots
                           if s.bandgap is not None and key in s.desc])
            by = np.array([s.bandgap for s in b.snapshots
                           if s.bandgap is not None and key in s.desc])
            per_branch[b.label] = _pearson(bx, by)[0] if len(bx) >= 3 else float("nan")
        within = np.array([v for v in per_branch.values() if not np.isnan(v)])
        rows.append({
            "key": key, "label": label, "n": len(x),
            "r": r, "p": p, "rho": rho, "slope": slope,
            "per_branch": per_branch,
            "within_mean": float(within.mean()) if len(within) else float("nan"),
        })
    return rows


def print_correlations(rows: list[dict], branches: list[Branch]) -> None:
    if not rows:
        print("\nNo structural descriptors available "
              "(ASE missing or POSCARs not found) — skipping distance analysis.")
        return
    blabels = [b.label for b in branches]
    head = (f"\n{'descriptor':<26}{'n':>5}{'pooled r':>10}{'p':>10}"
            f"{'spearman':>10}{'slope':>11}{'within-T r':>12}")
    print(head)
    print("-" * len(head))
    for row in rows:
        pstr = "  n/a" if np.isnan(row["p"]) else f"{row['p']:.1e}"
        print(f"{row['label']:<26}{row['n']:>5}{row['r']:>10.3f}{pstr:>10}"
              f"{row['rho']:>10.3f}{row['slope']:>11.3f}{row['within_mean']:>12.3f}")
    print(f"\n  slope = dE_gap/d(descriptor) in eV/Å (or eV per unit of the σ "
          f"descriptors).\n  per-branch (within-T) Pearson r: " +
          ", ".join(blabels))
    for row in rows:
        vals = ", ".join(f"{row['per_branch'][l]:+.3f}" for l in blabels)
        print(f"    {row['label']:<26} {vals}")
    print()


def plot_gap_vs_distance(branches: list[Branch], rows: list[dict], out: Path,
                         keys: list[str]) -> None:
    """Scatter of band gap vs. each chosen descriptor, coloured by branch,
    with a pooled linear fit overlaid."""
    keys = [k for k in keys if any(k in s.desc for b in branches for s in b.snapshots)]
    if not keys:
        return
    by_key = {row["key"]: row for row in rows}
    ncol = 2
    nrow = (len(keys) + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(6.2 * ncol, 4.6 * nrow),
                             squeeze=False)
    axes = axes.ravel()
    cmap = plt.get_cmap("viridis")
    nb = max(1, len(branches) - 1)

    for ax, key in zip(axes, keys):
        label = dict(DESCRIPTORS).get(key, key)
        for i, b in enumerate(branches):
            x = np.array([s.desc[key] for s in b.snapshots
                          if s.bandgap is not None and key in s.desc])
            y = np.array([s.bandgap for s in b.snapshots
                          if s.bandgap is not None and key in s.desc])
            if len(x):
                ax.scatter(x, y, s=22, alpha=0.6, color=cmap(i / nb), label=b.label)
        x_all, y_all = _paired(branches, key)
        if len(x_all) >= 2:
            xs = np.linspace(x_all.min(), x_all.max(), 50)
            m, c = np.polyfit(x_all, y_all, 1)
            ax.plot(xs, m * xs + c, "k--", lw=1.5)
            row = by_key.get(key, {})
            ax.set_title(f"{label}\npooled r = {row.get('r', float('nan')):.2f}, "
                         f"within-T r = {row.get('within_mean', float('nan')):.2f}")
        ax.set_xlabel(label)
        ax.set_ylabel("Band gap (eV)")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)

    for ax in axes[len(keys):]:
        ax.set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def plot_correlation_bars(rows: list[dict], out: Path) -> None:
    """Bar chart comparing pooled vs. within-temperature Pearson r per descriptor."""
    if not rows:
        return
    labels = [r["label"] for r in rows]
    pooled = [r["r"] for r in rows]
    within = [r["within_mean"] for r in rows]
    y = np.arange(len(labels))
    h = 0.38

    fig, ax = plt.subplots(figsize=(8, 0.7 * len(labels) + 2))
    ax.barh(y + h / 2, pooled, height=h, color="tab:purple", label="pooled (all snapshots)")
    ax.barh(y - h / 2, within, height=h, color="tab:orange", label="within-T mean")
    ax.axvline(0, color="k", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Pearson correlation with band gap")
    ax.set_title("Band-gap sensitivity to interatomic-distance descriptors")
    ax.legend()
    ax.grid(axis="x", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


# --------------------------------------------------------------------------- #
def write_csv(branches: list[Branch], out: Path) -> None:
    desc_keys = [k for k, _ in DESCRIPTORS]
    with out.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["branch", "label", "config", "bandgap_eV", "e_fermi_eV"] + desc_keys)
        for b in branches:
            for s in b.snapshots:
                row = [b.name, b.label, s.name,
                       "" if s.bandgap is None else f"{s.bandgap:.6f}",
                       "" if s.e_fermi is None else f"{s.e_fermi:.6f}"]
                row += ["" if k not in s.desc else f"{s.desc[k]:.6f}" for k in desc_keys]
                w.writerow(row)


def print_summary(branches: list[Branch]) -> None:
    print(f"\n{'branch':<20}{'n':>4}{'mean gap':>12}{'std':>10}"
          f"{'min':>10}{'max':>10}  (eV)")
    print("-" * 76)
    for b in branches:
        g = b.gaps
        if len(g) == 0:
            print(f"{b.label:<20}{0:>4}   (no band gaps parsed)")
            continue
        std = g.std(ddof=1) if len(g) > 1 else 0.0
        print(f"{b.label:<20}{len(g):>4}{g.mean():>12.4f}{std:>10.4f}"
              f"{g.min():>10.4f}{g.max():>10.4f}")
    print()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("results_dir", nargs="?", default="dft_analysis_data",
                    help="directory holding <branch>/2-DFT/config_*/ (default: ./dft_analysis_data)")
    ap.add_argument("-o", "--outdir", default=None,
                    help="output directory for plots/CSV (default: <results_dir>/analysis)")
    args = ap.parse_args()

    results_dir = Path(args.results_dir).expanduser().resolve()
    if not results_dir.is_dir():
        ap.error(f"results directory not found: {results_dir}")
    outdir = Path(args.outdir).expanduser().resolve() if args.outdir \
        else results_dir / "analysis"
    outdir.mkdir(parents=True, exist_ok=True)

    branches = collect(results_dir)
    n_snaps = sum(len(b.snapshots) for b in branches)
    n_gaps = sum(len(b.gaps) for b in branches)
    print(f"Found {len(branches)} branch(es), {n_snaps} snapshot(s), "
          f"{n_gaps} band gap(s).")
    print_summary(branches)

    write_csv(branches, outdir / "bandgaps.csv")
    plot_bandgap_distribution(branches, outdir / "bandgap_distribution.png")
    plot_bandgap_vs_temperature(branches, outdir / "bandgap_vs_temperature.png")
    plot_dos(branches, outdir / "dos_per_branch.png")
    plot_dos_overlay(branches, outdir / "dos_overlay.png")

    # band gap vs. interatomic distance
    n_desc = sum(1 for b in branches for s in b.snapshots if s.desc)
    if n_desc:
        print(f"Structural descriptors computed for {n_desc} snapshot(s).")
        rows = correlation_stats(branches)
        print_correlations(rows, branches)
        plot_gap_vs_distance(branches, rows, outdir / "gap_vs_distance.png",
                             keys=["d_metalN_mean", "d_metalN_std",
                                   "d_CuN_mean", "d_TaN_mean"])
        plot_correlation_bars(rows, outdir / "gap_distance_correlation.png")
    elif _ase_read is None:
        print("\nASE not available — skipping interatomic-distance analysis "
              "(`pip install ase` to enable).")
    else:
        print("\nNo POSCARs found next to the snapshots — "
              "skipping interatomic-distance analysis.")

    print(f"Wrote results to {outdir}:")
    for p in sorted(outdir.iterdir()):
        print(f"  {p.name}")


if __name__ == "__main__":
    main()

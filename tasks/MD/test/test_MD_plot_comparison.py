from typing import Tuple
import yaml
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import scipy.integrate as integrate
import subprocess


def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:

    print("Start task: test_MD_plot_comparison", flush=True)

    # Read in global configurations
    with open(path_configWF, "r") as f:
        configWF = yaml.safe_load(f)

    run_test_MD = configWF.get("run_test_MD", True)

    if run_test_MD:

        dir_project = Path(configWF.get("dir_project", "./"))

        # determine units for plotting
        units_type = configWF["lammps"].get("units")

        if units_type == "real":
            units = ["A", "(A/fs)", "((kcal/mol)/A)", "PHz", "fs"]
        elif units_type == "metal":
            units = ["A", "(A/ps)", "(eV/A)", "THz", "ps"]
        elif units_type == "si":
            units = ["m", "(m/s)", "N", "Hz", "s"]
        elif units_type == "cgs":
            units = ["cm", "(cm/s)", "dynes", "Hz", "s"]
        elif units_type == "electron":
            units = ["Bohr", "(Bohr/atomic time units)", "(Hartrees/Bohr)", "PHz", "fs"]
        elif units_type == "micro": 
            units = ["μm", "(m/s)", "nN", "MHz", "μs"]
        elif units_type == "nano":
            units = ["nm", "(m/s)", "pN", "GHz", "ns"]
        else:
            if configWF["lammps"].get("units_array") is not None:
                units = configWF["lammps"]["units_array"][5:]
            else:
                print("Unknown units type. Using no units.", flush=True)
                units = ["", "", "", "", ""]


        # Handle multiple simulations (branching)
        if num_simulations > 1:

            dir_plots = dir_project / "plots/"
            dir_plots.mkdir(parents=True, exist_ok=True)

            # determine correct branch config file
            param_to_vary = configWF["param_to_vary"]
            array_to_vary = configWF["array_to_vary"]

            # read in element names for each atom type
            elements = configWF["lammps"]["elements"]

            # get maximal frequency for plotting if specified
            omega_max = configWF.get("vdos_omega_max", None)

            dir_plot_vdos = dir_plots / "MD/vdos/"
            dir_plot_vdos.mkdir(parents=True, exist_ok=True)

            # plot VDOS comparison for each atom type and total VDOS
            for e in ["total", *elements]:

                fig1, (ax1) = plt.subplots()
                plt.title(f"VDOS ({e}) comparison between different MD trajectories")

                # plot vdos for each branch
                for i in range(len(array_to_vary)):

                    # read in branch configuration
                    dir_project_i = dir_project / f"{param_to_vary}_{array_to_vary[i]}/"

                    with open(str(dir_project_i / 'branch_config.yaml'), 'r') as f:
                        configWF_i = yaml.safe_load(f)

                    dir_MD = Path(configWF_i.get("dir_MD", str(dir_project_i / "1-MD/")))

                    # read in VDOS
                    freq, vdos = np.loadtxt(str(dir_MD / f"test_MD/vdos/vdos_{e}.txt"), unpack=True, skiprows=1)

                    # plot VDOS
                    label = f"{param_to_vary} {array_to_vary[i]}"
                    plt.plot(freq, vdos, label=label)

                # finalize plot
                plt.xlabel(f"Frequency/{units[3]}")
                plt.ylabel(f"VDOS (arb. units)")
                if omega_max is not None:
                    plt.xlim(0, omega_max)
                plt.legend()

                outfile = dir_plot_vdos / f"vdos_{e}_comparison.pdf"
                plt.tight_layout()
                plt.savefig(outfile)
                plt.close(fig1)
                print(f"✅ Saved VDOS ({e}) comparison plot → {outfile}", flush=True)

            dir_plot_dis = dir_plots / "MD/displacement_distribution"
            dir_plot_dis.mkdir(parents=True, exist_ok=True)

            # plot distribution comparison for each atom type
            for e in elements: 

                fig2, axes2 = plt.subplots(3, 1)
                plt.suptitle(f"Displacement distribution ({e}) comparison between different MD trajectories")

                # plot distributions for each branch
                for i in range(len(array_to_vary)):

                    # read in branch configuration
                    dir_project_i = dir_project / f"{param_to_vary}_{array_to_vary[i]}/"

                    with open(str(dir_project_i / 'branch_config.yaml'), 'r') as f:
                        configWF_i = yaml.safe_load(f)

                    dir_MD = Path(configWF_i.get("dir_MD", str(dir_project_i / "1-MD/")))

                    # read in distributions
                    disp_x, dens_x = np.loadtxt(str(dir_MD / f"test_MD/histograms/displacements_histogram/displacements_histogram_{e}_Δx.txt"), unpack=True, skiprows=1)
                    disp_y, dens_y = np.loadtxt(str(dir_MD / f"test_MD/histograms/displacements_histogram/displacements_histogram_{e}_Δy.txt"), unpack=True, skiprows=1)
                    disp_z, dens_z = np.loadtxt(str(dir_MD / f"test_MD/histograms/displacements_histogram/displacements_histogram_{e}_Δz.txt"), unpack=True, skiprows=1)

                    # plot distributions
                    label = f"{param_to_vary} {array_to_vary[i]}"
                    axes2[0].plot(disp_x, dens_x, label=label)
                    axes2[1].plot(disp_y, dens_y, label=label)
                    axes2[2].plot(disp_z, dens_z, label=label)

                # finalize plot
                axes2[0].set_xlabel(f"Δx/{units[0]}")
                axes2[0].set_ylabel("Probability density")
                axes2[1].set_xlabel(f"Δy/{units[0]}")
                axes2[1].set_ylabel("Probability density")
                axes2[2].set_xlabel(f"Δz/{units[0]}")
                axes2[2].set_ylabel("Probability density")

                plt.legend()
                outfile = dir_plot_dis / f"displacement_distribution_{e}_comparison.pdf"
                plt.tight_layout()
                plt.savefig(outfile)
                plt.close(fig2)
                print(f"✅ Saved Displacement distribution ({e}) comparison plot → {outfile}", flush=True)


            if configWF_i["lammps"].get("compute_msd", True):
                dir_plot_msd = dir_plots / "MD/msd/"
                dir_plot_msd.mkdir(parents=True, exist_ok=True)

                # convert step to time
                plot_MD_time = configWF.get("plot_MD_time", True)
                if plot_MD_time:
                    dt = configWF["lammps"]["dt"] * configWF["lammps"].get("thermo_output_step_size", 100)

                for e in ["all", *elements]:

                    # plot MSD comparison
                    fig3, (ax3) = plt.subplots()
                    plt.title(f"MSD ({e}) comparison between different MD trajectories")

                    # plot MSD for each branch
                    for i in range(len(array_to_vary)):

                        # read in branch configuration
                        dir_project_i = dir_project / f"{param_to_vary}_{array_to_vary[i]}/"

                        with open(str(dir_project_i / 'branch_config.yaml'), 'r') as f:
                            configWF_i = yaml.safe_load(f)

                        dir_MD = Path(configWF_i.get("dir_MD", str(dir_project_i / "1-MD/")))

                        # read in MSD
                        step, msd_x, msd_y, msd_z, msd_tot = np.loadtxt(str(dir_MD / f"msd/msd_{e}.txt"), unpack=True, skiprows=2)

                        step -= step[0]

                        # convert step to time
                        if plot_MD_time:
                            step *= dt

                        # plot MSD
                        label = f"{param_to_vary} {array_to_vary[i]}"
                        plt.plot(step, msd_tot, label=label)

                    # finalize plot
                    if plot_MD_time:
                        plt.xlabel(f"Time/{units[4]}")
                    else:
                        plt.xlabel("Step")
                    plt.ylabel(f"MSD/{units[0]}^2")
                    plt.legend()

                    outfile = dir_plot_msd / f"msd_{e}_comparison.pdf"
                    plt.tight_layout()
                    plt.savefig(outfile)
                    plt.close(fig3)
                    print(f"✅ Saved MSD comparison plot → {outfile}", flush=True)


            if configWF_i.get("compute_rdf", True):
                
                dir_plot_rdf = dir_plots / "MD/rdf/"
                dir_plot_rdf.mkdir(parents=True, exist_ok=True)

                for e in ["all", *elements]:

                    # plot RDF comparison
                    fig4, (ax4) = plt.subplots()
                    plt.title(f"RDF ({e}) comparison between different MD trajectories")

                    # plot RDF for each branch
                    for i in range(len(array_to_vary)):

                        # read in branch configuration
                        dir_project_i = dir_project / f"{param_to_vary}_{array_to_vary[i]}/"

                        with open(str(dir_project_i / 'branch_config.yaml'), 'r') as f:
                            configWF_i = yaml.safe_load(f)

                        dir_MD = Path(configWF_i.get("dir_MD", str(dir_project_i / "1-MD/")))

                        # read in RDF
                        if e == "all":
                            _, r, rdf, _ = np.loadtxt(str(dir_MD / f"rdf/rdf_all.txt"), unpack=True, skiprows=4)
                        else:
                            _, r, rdf, _ = np.loadtxt(str(dir_MD / f"rdf/rdf_{e}-{e}.txt"), unpack=True, skiprows=4)

                        # plot RDF
                        label = f"{param_to_vary} {array_to_vary[i]}"
                        plt.plot(r, rdf, label=label)

                    # finalize plot
                    plt.xlabel(f"r/{units[0]}")
                    plt.ylabel("g(r)")
                    plt.legend()

                    outfile = dir_plot_rdf / f"rdf_{e}_comparison.pdf"
                    plt.tight_layout()
                    plt.savefig(outfile)
                    plt.close(fig4)
                    print(f"✅ Saved RDF comparison plot → {outfile}", flush=True)


        else:
            print("No comparison of different MDs trajectories needed since only one MD trajectory was calculated.", flush=True)

    else:
        print("Skip test_MD_plot_comparison.", flush=True)

    print("Finish task: test_MD_plot_comparison", flush=True)

    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations}

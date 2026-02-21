from typing import Tuple
import yaml
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from perqueue.constants import CYCLICALGROUP_KEY

def moving_average(a, w):
    cumsum = np.cumsum(np.insert(a, 0, 0)) 
    return (cumsum[w:] - cumsum[:-w]) / w


def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:

    print("Start task: test_MD_equilibration", flush=True)

    # Read in global configurations
    with open(path_configWF, "r") as f:
        configWF = yaml.safe_load(f)

    dir_project = Path(configWF.get("dir_project", "./"))
    simulation_type = configWF.get("simulation_type", "sweep")
    cg_criteria = False


    # Handle multiple simulations (branching)
    if num_simulations > 1:

        if simulation_type == "cascade":
            i = kwargs['pq_iteration'][0]
        else:
            i = kwargs['pq_index'][0]

        # determine correct branch config file
        param_to_vary = configWF["param_to_vary"]
        array_to_vary = configWF["array_to_vary"]
    
        if i == (len(array_to_vary)-1):
            cg_criteria = True
            
        dir_project_i = dir_project / f"{param_to_vary}_{array_to_vary[i]}/"

        with open(str(dir_project_i / 'branch_config.yaml'), 'r') as f:
            configWF_i = yaml.safe_load(f)


    else:
        configWF_i = configWF.copy()
        dir_project_i = dir_project

    dir_MD = Path(configWF_i.get("dir_MD", str(dir_project_i / "1-MD/")))
    equilibrate = configWF_i.get("equilibrate", False)
    npt_equilibrate = configWF_i.get("npt_equilibrate", True)

    if equilibrate:

        # read in thermodynamic configurations
        if "temperature" in configWF_i:
            T_set = configWF_i["temperature"]  
        else:
            T_set = configWF_i["lammps"]["T"]

        p_set = configWF_i["lammps"].get("P", 0.0)
        T_start = configWF_i["lammps"].get("T_start", T_set)
        p_start = configWF_i["lammps"].get("P_start", p_set)

        # get indeces where to split the trajectory for different ensembles
        if T_start != T_set:
            n1 = configWF_i["lammps"]['eqsteps_nvt_heating']
        else:
            n1 = 0

        n2 = configWF_i["lammps"]['eqsteps_nvt']

        if npt_equilibrate:

            if p_start != p_set:
                n3 = configWF_i["lammps"]['eqsteps_npt_expansion'] 
            else: 
                n3 = 0

            n4 = configWF_i["lammps"]['eqsteps_npt']

        else:
            n3 = 0
            n4 = 0


        # Read in LAMMPS dump data
        dir_test = dir_MD / "test_equilibrate/"
        dir_test.mkdir(parents=True, exist_ok=True)

        N_atoms = np.loadtxt(str(dir_MD / "position.lammpstrj"), unpack=True, skiprows=3, max_rows=1)

        data_MD = np.loadtxt(str(dir_MD / "thermo_output.txt"), unpack=True, skiprows=1)

        step = data_MD[0,:]
        T = data_MD[1,:]
        E = data_MD[2,:]
        L = data_MD[3:6,:]
        V = data_MD[6,:]
        p = data_MD[7,:]

        n_tot = step[-1]
        x1 = n1/n_tot
        x2 = (n1 + n2)/n_tot
        x3 = (n1 + n2 + n3)/n_tot
        x4 = (n1 + n2 + n3 + n4)/n_tot

        w = configWF_i.get("window_size", 20)  # window size for moving average

        # convert step to time
        plot_MD_time = configWF_i.get("plot_MD_time", True)
        if plot_MD_time:
            dt = configWF_i["lammps"]["dt"] 
            step *= dt  
            n1 *= dt
            n2 *= dt
            n3 *= dt
            n4 *= dt

        # determine units for plotting
        units_type = configWF_i["lammps"].get("units")

        if units_type == "real":
            units = ["K", "kcal/mol", "A", "A^3", "atm", "fs"]
        elif units_type == "metal":
            units = ["K", "eV", "A", "A^3", "bar", "ps"]
        elif units_type == "si":
            units = ["K", "J", "m", "m^3", "Pa", "s"]
        elif units_type == "cgs":
            units = ["K", "erg", "cm", "cm^3", "barye", "s"]
        elif units_type == "electron":
            units = ["K", "Hartrees", "Bohr", "Bohr^3", "Pa", "fs"]
        elif units_type == "micro": 
            units = ["K", "fJ", "μm", "μm^3", "kPa", "μs"]
        elif units_type == "nano":
            units = ["K", "zJ", "nm", "nm^3", "MPa", "ns"]
        else:
            if configWF_i["lammps"].get("units_array") is not None:
                units = configWF_i["lammps"]["units_array"]
            else:
                print("WARNING: Unknown units type. Using no units.", flush=True)
                units = ["", "", "", "", "", ""]


        ### ### ### ### ### ### ### 
        ### temperature
        ### ### ### ### ### ### ### 

        # calculate statistical temperature values for temperature equilibration criterion
        T_eq_avg_nvt = np.mean(T[(n1 < step) & (step <= n1 + n2)][-2*w:])
        T_eq_std_nvt = np.std(T[(n1 < step) & (step <= n1 + n2)][-2*w:])
        T_avg_nvt = np.mean(T[(n1 < step) & (step <= n1 + n2)])
        T_std_nvt = np.std(T[(n1 < step) & (step <= n1 + n2)])
        if npt_equilibrate:
            T_eq_avg_npt = np.mean(T[(n1 + n2 + n3 < step) & (step <= n1 + n2 + n3 + n4)][-2*w:])
            T_eq_std_npt = np.std(T[(n1 + n2 + n3 < step) & (step <= n1 + n2 + n3 + n4)][-2*w:])
            T_avg_npt = np.mean(T[(n1 + n2 + n3 < step) & (step <= n1 + n2 + n3 + n4)])
            T_std_npt = np.std(T[(n1 + n2 + n3 < step) & (step <= n1 + n2 + n3 + n4)])
        T_avg_prod = np.mean(T[n1 + n2 + n3 + n4 < step])
        T_std_prod = np.std(T[n1 + n2 + n3 + n4 < step])

        # temperature equilibration criterion
        dT_nvt = T_avg_nvt * np.sqrt(1.5 / N_atoms)
        if npt_equilibrate:
            dT_npt = T_avg_npt * np.sqrt(1.5 / N_atoms)

        with open(str(dir_MD / "test_equilibrate/thermo_avg_std.txt"), "w") as f:
            f.write("TEMPERATURE: \n\n")
            f.write(f"  Equilibrated average temperature (NVT): {T_eq_avg_nvt} +/- {T_eq_std_nvt} {units[0]}\n")
            f.write(f"  Overall average temperature (NVT): {T_avg_nvt} +/- {T_std_nvt} {units[0]}\n")
            f.write(f"  Expected thermal fluctuation (NVT): {dT_nvt} {units[0]}\n\n")
            if npt_equilibrate:
                f.write(f"  Equilibrated average temperature (NPT): {T_eq_avg_npt} +/- {T_eq_std_npt} {units[0]}\n")
                f.write(f"  Overall average temperature (NPT): {T_avg_npt} +/- {T_std_npt} {units[0]}\n")
                f.write(f"  Expected thermal fluctuation (NPT): {dT_npt} {units[0]}\n\n")
            f.write(f"  Average temperature (production): {T_avg_prod} +/- {T_std_prod} {units[0]}\n\n\n")


        #if "temperature" in configWF:
        #    T_set = configWF_i["temperature"]  
        #else:
        #    T_set = configWF_i["lammps"]["T"]

        error = False

        # check for temperature equilibration criterion
        if (dT_nvt < T_eq_std_nvt):
            print("dT_nvt < T_eq_std_nvt: ", dT_nvt, f" {units[0]}", " < ", T_eq_std_nvt, f" {units[0]}", flush=True)
            error = True
        if npt_equilibrate:
            if (dT_npt < T_eq_std_npt):
                print("dT_npt < T_eq_std_npt: ", dT_npt, f" {units[0]}", " < ", T_eq_std_npt, f" {units[0]}", flush=True)
                error = True
        if (abs(T_eq_avg_nvt - T_set) > 2 * T_eq_std_nvt):
            print("abs(T_eq_avg_nvt - T_set) > 2 * T_eq_std_nvt: ", abs(T_eq_avg_nvt - T_set), f" {units[0]}", " > ", 2 * T_eq_std_nvt, f" {units[0]}", flush=True)
            print("T_eq_avg_nvt: ", T_eq_avg_nvt, f" {units[0]}", flush=True)
            print("T_set: ", T_set, f" {units[0]}", flush=True)
            error = True
        if npt_equilibrate:
            if (abs(T_eq_avg_npt - T_set) > 2 * T_eq_std_npt):
                print("abs(T_eq_avg_npt - T_set) > 2 * T_eq_std_npt: ", abs(T_eq_avg_npt - T_set), f" {units[0]}", " > ", 2 * T_eq_std_npt, f" {units[0]}", flush=True)
                print("T_eq_avg_npt: ", T_eq_avg_npt, f" {units[0]}", flush=True)
                print("T_set: ", T_set, f" {units[0]}", flush=True)
                error = True

        # plot temperature data
        fig1, (ax1) =  plt.subplots()
        fig1.suptitle("Temperature")
        if plot_MD_time:
            plt.xlabel("Time/" + units[5])
        else:
            plt.xlabel("Step")
        plt.ylabel(f"Temperature/{units[0]}")
        plt.plot(step, T, color='blue', label="Temperature")
        plt.plot(moving_average(step, w), moving_average(T, w), color='orange', label="moving average")
        plt.axhline(y=T_avg_nvt, xmin=x1, xmax=x2, color='red', label="Average")
        plt.axhline(y=T_avg_nvt+T_std_nvt, xmin=x1, xmax=x2, linestyle='--', color='red', label="Std")
        plt.axhline(y=T_avg_nvt-T_std_nvt, xmin=x1, xmax=x2, linestyle='--', color='red')
        if npt_equilibrate:
            plt.axhline(y=T_avg_npt, xmin=x3, xmax=x4, color='red')
            plt.axhline(y=T_avg_npt+T_std_npt, xmin=x3, xmax=x4, linestyle='--', color='red')
            plt.axhline(y=T_avg_npt-T_std_npt, xmin=x3, xmax=x4, linestyle='--', color='red')
        plt.axhline(y=T_avg_prod, xmin=x4, xmax=1.0, color='red')
        plt.axhline(y=T_avg_prod+T_std_prod, xmin=x4, xmax=1.0, linestyle='--', color='red')
        plt.axhline(y=T_avg_prod-T_std_prod, xmin=x4, xmax=1.0, linestyle='--', color='red')

        if n1 != 0:
            plt.axvline(x=n1, linestyle='--', color="black")
        plt.axvline(x=n1+n2, linestyle='--', color="black")
        if n3 != 0:
            plt.axvline(x=n1+n2+n3, linestyle='--', color="black")
        if n4 != 0:
            plt.axvline(x=n1+n2+n3+n4, linestyle='--', color="black")

        plt.xlim(0, step[-1])
        plt.legend()
        plt.savefig(str(dir_MD / "test_equilibrate/T.pdf"))
        plt.close()


        ### ### ### ### ### ### ### 
        ### energy
        ### ### ### ### ### ### ### 

        # calculate statistical energy values
        E_eq_avg_nvt = np.mean(E[(n1 < step) & (step <= n1 + n2)][-2*w:])
        E_eq_std_nvt = np.std(E[(n1 < step) & (step <= n1 + n2)][-2*w:])
        E_avg_nvt = np.mean(E[(n1 < step) & (step <= n1 + n2)])
        E_std_nvt = np.std(E[(n1 < step) & (step <= n1 + n2)])
        if npt_equilibrate:
            E_eq_avg_npt = np.mean(E[(n1 + n2 + n3 < step) & (step <= n1 + n2 + n3 + n4)][-2*w:])
            E_eq_std_npt = np.std(E[(n1 + n2 + n3 < step) & (step <= n1 + n2 + n3 +n4)][-2*w:])
            E_avg_npt = np.mean(E[(n1 + n2 + n3 < step) & (step <= n1 + n2 + n3 + n4)])
            E_std_npt = np.std(E[(n1 + n2 + n3 < step) & (step <= n1 + n2 + n3 +n4)])
        E_avg_prod = np.mean(E[n1 + n2 + n3 + n4 < step])
        E_std_prod = np.std(E[n1 + n2 + n3 + n4 < step])

        with open(str(dir_MD / "test_equilibrate/thermo_avg_std.txt"), "a") as f:
            f.write("ENERGY PER ATOM: \n\n")
            f.write(f"  Equilibrated average energy per atom (NVT): {E_eq_avg_nvt} +/- {E_eq_std_nvt} {units[1]}\n")
            f.write(f"  Overall average energy per atom (NVT): {E_avg_nvt} +/- {E_std_nvt} {units[1]}\n\n")
            if npt_equilibrate:
                f.write(f"  Equilibrated average energy per atom (NPT): {E_eq_avg_npt} +/- {E_eq_std_npt} {units[1]}\n")
                f.write(f"  Overall average energy per atom (NPT): {E_avg_npt} +/- {E_std_npt} {units[1]}\n\n")
            f.write(f"  Average energy per atom (production): {E_avg_prod} +/- {E_std_prod} {units[1]}\n\n\n")

        # plot energy data
        fig1, (ax1) =  plt.subplots()
        fig1.suptitle("Energy per atom")
        if plot_MD_time:
            plt.xlabel("Time/" + units[5])
        else:
            plt.xlabel("Step")
        plt.ylabel(f"Energy per atom/{units[1]}")
        plt.plot(step, E/N_atoms, color='blue', label="energy")
        plt.plot(moving_average(step, w), moving_average(E, w)/N_atoms, color='orange', label="moving average")
        plt.axhline(y=E_avg_nvt/N_atoms, xmin=x1, xmax=x2, color='red', label="Average")
        plt.axhline(y=(E_avg_nvt+E_std_nvt)/N_atoms, xmin=x1, xmax=x2, linestyle='--', color='red', label="Std")
        plt.axhline(y=(E_avg_nvt-E_std_nvt)/N_atoms, xmin=x1, xmax=x2, linestyle='--', color='red')
        if npt_equilibrate:
            plt.axhline(y=E_avg_npt/N_atoms, xmin=x3, xmax=x4, color='red')
            plt.axhline(y=(E_avg_npt+E_std_npt)/N_atoms, xmin=x3, xmax=x4, linestyle='--', color='red')
            plt.axhline(y=(E_avg_npt-E_std_npt)/N_atoms, xmin=x3, xmax=x4, linestyle='--', color='red')
        plt.axhline(y=E_avg_prod/N_atoms, xmin=x4, xmax=1.0, color='red')
        plt.axhline(y=(E_avg_prod+E_std_prod)/N_atoms, xmin=x4, xmax=1.0, linestyle='--', color='red')
        plt.axhline(y=(E_avg_prod-E_std_prod)/N_atoms, xmin=x4, xmax=1.0, linestyle='--', color='red')

        if n1 != 0:
            plt.axvline(x=n1, linestyle='--', color="black")
        plt.axvline(x=n1+n2, linestyle='--', color="black")
        if n3 != 0:
            plt.axvline(x=n1+n2+n3, linestyle='--', color="black")
        if n4 != 0:
            plt.axvline(x=n1+n2+n3+n4, linestyle='--', color="black")

        plt.xlim(0, step[-1])
        plt.legend()
        plt.savefig(str(dir_MD / "test_equilibrate/E.pdf"))
        plt.close()


        ### ### ### ### ### ### ### 
        ### volume
        ### ### ### ### ### ### ### 

        # calculate statistical volume values
        with open(str(dir_MD / "test_equilibrate/thermo_avg_std.txt"), "a") as f:
            if npt_equilibrate:
                V_eq_avg_npt = np.mean(V[(n1 + n2 + n3 < step) & (step <= n1 + n2 + n3 + n4)][-2*w:])
                V_eq_std_npt = np.std(V[(n1 + n2 + n3 < step) & (step <= n1 + n2 + n3 + n4)][-2*w:])
                V_avg_npt = np.mean(V[(n1 + n2 + n3 < step) & (step <= n1 + n2 + n3 + n4)])
                V_std_npt = np.std(V[(n1 + n2 + n3 < step) & (step <= n1 + n2 + n3 + n4)])

                f.write("VOLUME: \n\n")
                f.write(f"  Equilibrated average volume (NPT): {V_eq_avg_npt} +/- {V_eq_std_npt} {units[3]}\n")
                f.write(f"  Overall average volume (NPT): {V_avg_npt} +/- {V_std_npt} {units[3]}\n\n")

            f.write(f"  Volume (production): {V[-1]} {units[3]}\n\n\n")

        # plot volume data
        fig1, (ax1) =  plt.subplots()
        fig1.suptitle("Volume")
        if plot_MD_time:
            plt.xlabel("Time/" + units[5])
        else:
            plt.xlabel("Step")
        plt.ylabel(f"Volume/{units[3]}")
        plt.plot(step, V, color='blue', label="Volume")
        plt.plot(moving_average(step, w), moving_average(V, w), color='orange', label="moving average")
        if npt_equilibrate:
            plt.axhline(y=V_avg_npt, xmin=x3, xmax=x4, color='red', label="Average")
            plt.axhline(y=V_avg_npt+V_std_npt, xmin=x3, xmax=x4, linestyle='--', color='red', label="Std")
            plt.axhline(y=V_avg_npt-V_std_npt, xmin=x3, xmax=x4, linestyle='--', color='red')

        if n1 != 0:
            plt.axvline(x=n1, linestyle='--', color="black")
        plt.axvline(x=n1+n2, linestyle='--', color="black")
        if n3 != 0:
            plt.axvline(x=n1+n2+n3, linestyle='--', color="black")
        if n4 != 0:
            plt.axvline(x=n1+n2+n3+n4, linestyle='--', color="black")

        plt.xlim(0, step[-1])
        if min(V) == max(V):
            plt.ylim(min(V)*0.99, max(V)*1.01)

        plt.legend()
        plt.savefig(str(dir_MD / "test_equilibrate/V.pdf"))
        plt.close()


        ### ### ### ### ### ### ### 
        ### lattice constants
        ### ### ### ### ### ### ### 

        # calculate statistical lattice constants values
        if npt_equilibrate:
            L_eq_avg_npt = np.mean(L[:, (n1 + n2 + n3 < step) & (step <= n1 + n2 + n3 + n4)][-2*w:], axis=1)
            L_eq_std_npt = np.std(L[:, (n1 + n2 + n3 < step) & (step <= n1 + n2 + n3 + n4)][-2*w:], axis=1)
            L_avg_npt = np.mean(L[:, (n1 + n2 + n3 < step) & (step <= n1 + n2 + n3 + n4)], axis=1)
            L_std_npt = np.std(L[:, (n1 + n2 + n3 < step) & (step <= n1 + n2 + n3 + n4)], axis=1)


        with open(str(dir_MD / "test_equilibrate/thermo_avg_std.txt"), "a") as f:
            f.write("LATTICE CONSTANTS: \n\n")

        for i in range(3):
            if npt_equilibrate:
                with open(str(dir_MD / "test_equilibrate/thermo_avg_std.txt"), "a") as f:
                    f.write(f"  Equilibrated average lattice constant {i+1} (NPT): {L_eq_avg_npt[i]} +/- {L_eq_std_npt[i]} {units[2]}\n")
                    f.write(f"  Overall average lattice constant {i+1} (NPT): {L_avg_npt[i]} +/- {L_std_npt[i]} {units[2]}\n")

            with open(str(dir_MD / "test_equilibrate/thermo_avg_std.txt"), "a") as f:
                f.write(f"  Lattice constant {i+1} (production): {L[i, -1]} {units[2]}\n\n")

            # plot lattice constants data
            fig1, (ax1) =  plt.subplots()
            fig1.suptitle(f"Lattice constant {i+1}")
            if plot_MD_time:
                plt.xlabel("Time/" + units[5])
            else:
                plt.xlabel("Step")
            plt.ylabel(f"Lattice constant/{units[2]}")
            plt.plot(step, L[i], color="blue", label=f"lattice constant {i+1}")
            plt.plot(moving_average(step, w), moving_average(L[i], w), color='orange', label="moving average")
            if npt_equilibrate:
                plt.axhline(y=L_avg_npt[i], xmin=x3, xmax=x4, color="red", label="Average")
                plt.axhline(y=L_avg_npt[i]+L_std_npt[i], xmin=x3, xmax=x4, linestyle='--', color="red", label="Std")
                plt.axhline(y=L_avg_npt[i]-L_std_npt[i], xmin=x3, xmax=x4, linestyle='--', color="red")

            if n1 != 0:
                plt.axvline(x=n1, linestyle='--', color="black")
            plt.axvline(x=n1+n2, linestyle='--', color="black")
            if n3 != 0:
                plt.axvline(x=n1+n2+n3, linestyle='--', color="black")
            if n4 != 0:
                plt.axvline(x=n1+n2+n3+n4, linestyle='--', color="black")

            plt.xlim(0, step[-1])
            if min(L[i]) == max(L[i]):
                plt.ylim(min(L[i])*0.99, max(L[i])*1.01)
                
            plt.legend()
            plt.savefig(str(dir_MD / f"test_equilibrate/L{i+1}.pdf"))
            plt.close()


        if npt_equilibrate:
            with open(str(dir_MD / "test_equilibrate/thermo_avg_std.txt"), "a") as f:
                f.write(f"  Equilibrated average lattice constant averaged over all dimensions (NPT): {np.mean(L_eq_avg_npt)} +/- {np.mean(L_eq_std_npt)} {units[2]}\n")
                f.write(f"  Overall average lattice constant averaged over all dimensions (NPT): {np.mean(L_avg_npt)} +/- {np.mean(L_std_npt)} {units[2]}\n")

        with open(str(dir_MD / "test_equilibrate/thermo_avg_std.txt"), "a") as f:
            f.write(f"  Lattice constant averaged over all dimensions (production): {np.mean(L[:, -1])} {units[2]}\n\n\n")


        ### ### ### ### ### ### ### 
        ### pressure
        ### ### ### ### ### ### ### 

        # calculate statistical pressure values
        p_eq_avg_nvt = np.mean(p[(n1 < step) & (step <= n1 + n2)][-2*w:])
        p_eq_std_nvt = np.std(p[(n1 < step) & (step <= n1 + n2)][-2*w:])
        p_avg_nvt = np.mean(p[(n1 < step) & (step <= n1 + n2)])
        p_std_nvt = np.std(p[(n1 < step) & (step <= n1 + n2)])
        if npt_equilibrate:
            p_eq_avg_npt = np.mean(p[(n1 + n2 + n3 < step) & (step <= n1 + n2 + n3 + n4)][-2*w:])
            p_eq_std_npt = np.std(p[(n1 + n2 + n3 < step) & (step <= n1 + n2 + n3 + n4)][-2*w:])
            p_avg_npt = np.mean(p[(n1 + n2 + n3 < step) & (step <= n1 + n2 + n3 + n4)])
            p_std_npt = np.std(p[(n1 + n2 + n3 < step) & (step <= n1 + n2 + n3 + n4)])
        p_avg_prod = np.mean(p[n1 + n2 + n3 + n4 < step])
        p_std_prod = np.std(p[n1 + n2 + n3 + n4 < step])

        with open(str(dir_MD / "test_equilibrate/thermo_avg_std.txt"), "a") as f:
            f.write("PRESSURE: \n\n")
            f.write(f"  Equilibrated average pressure (NVT): {p_eq_avg_nvt} +/- {p_eq_std_nvt} {units[4]}\n")
            f.write(f"  Overall average pressure (NVT): {p_avg_nvt} +/- {p_std_nvt} {units[4]}\n\n")
            if npt_equilibrate:
                f.write(f"  Equilibrated average pressure (NPT): {p_eq_avg_npt} +/- {p_eq_std_npt} {units[4]}\n")
                f.write(f"  Overall average pressure (NPT): {p_avg_npt} +/- {p_std_npt} {units[4]}\n\n")
            f.write(f"  Average pressure (production): {p_avg_prod} +/- {p_std_prod} {units[4]}\n\n\n")

        # plot pressure data
        fig1, (ax1) =  plt.subplots()
        fig1.suptitle("Pressure")
        if plot_MD_time:
            plt.xlabel("Time/" + units[5])
        else:
            plt.xlabel("Step")
        plt.ylabel(f"Pressure/{units[4]}")
        plt.plot(step, p, color='blue', label="Pressure")
        plt.plot(moving_average(step, w), moving_average(p, w), color='orange', label="moving average")
        plt.axhline(y=p_avg_nvt, xmin=x1, xmax=x2, color='red', label="Average")
        plt.axhline(y=p_avg_nvt+p_std_nvt, xmin=x1, xmax=x2, linestyle='--', color='red', label="Std")
        plt.axhline(y=p_avg_nvt-p_std_nvt, xmin=x1, xmax=x2, linestyle='--', color='red')
        if npt_equilibrate:
            plt.axhline(y=p_avg_npt, xmin=x3, xmax=x4, color='red')
            plt.axhline(y=p_avg_npt+p_std_npt, xmin=x3, xmax=x4, linestyle='--', color='red')
            plt.axhline(y=p_avg_npt-p_std_npt, xmin=x3, xmax=x4, linestyle='--', color='red')
        plt.axhline(y=p_avg_prod, xmin=x4, xmax=1.0, color='red')
        plt.axhline(y=p_avg_prod+p_std_prod, xmin=x4, xmax=1.0, linestyle='--', color='red')
        plt.axhline(y=p_avg_prod-p_std_prod, xmin=x4, xmax=1.0, linestyle='--', color='red')

        if n1 != 0:
            plt.axvline(x=n1, linestyle='--', color="black")
        plt.axvline(x=n1+n2, linestyle='--', color="black")
        if n3 != 0:
            plt.axvline(x=n1+n2+n3, linestyle='--', color="black")
        if n4 != 0:
            plt.axvline(x=n1+n2+n3+n4, linestyle='--', color="black")

        plt.xlim(0, step[-1])
        plt.legend()
        plt.savefig(str(dir_MD / "test_equilibrate/p.pdf"))
        plt.close()

        #p_set = configWF_i["lammps"].get("P", 0.0)

        if error:    
            raise Exception("Trajectory is NOT equilibrated with respect to temperature.")
        else:
            print("Trajectory is equilibrated with respect to temperature.", flush=True)

        # check if pressure is equilibrated
        if npt_equilibrate:
            if (abs(p_avg_npt - p_set) > p_std_npt):
                print("abs(p_avg_npt - p_set) > p_std_npt: ", abs(p_avg_npt - p_set), f" {units[4]}", " > ", p_std_npt, f" {units[4]}")
                print("p_avg_npt: ", p_avg_npt, f" {units[4]}", flush=True)
                print("p_set: ", p_set, f" {units[4]}", flush=True)
                error = True
        #if not np.allclose(moving_average(p, w)[-2*w:], p_set, atol=abs(p_set + p_std_npt)/p_set):
        #    print("Pressure is fluctuating too much around the set point.")
        #    error = True

            if error:    
                raise Exception("Trajectory is NOT equilibrated with respect to pressure.")
            else:
                print("Trajectory is equilibrated with respect to pressure.", flush=True)

        # human_in_loop flag allows for human check before continuing
        if configWF_i.get("human_in_loop", False):
            print("human_in_loop is set to True. Stopping workflow after each test.", flush=True)
            raise Exception("Equilibration test done. Please check the plots in " + str(dir_MD / "test_equilibrate/") + " and continue the workflow manually.")

    else:
        print("No equilibration was performed.", flush=True)

    print("Finished task: test_MD_equilibration", flush=True)

    return True, {CYCLICALGROUP_KEY: cg_criteria, "path_configWF": path_configWF, "num_simulations": num_simulations}
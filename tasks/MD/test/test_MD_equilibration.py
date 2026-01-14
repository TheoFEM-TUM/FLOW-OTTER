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


    if equilibrate:

        # read in temperature
        if "temperature" in configWF_i:
            T = configWF_i["temperature"]  
        else:
            T = configWF_i["lammps"]["T"]

        T_start = configWF_i["lammps"].get("T_start", T)

        # get indeces where to split the trajectory for different ensembles
        if T_start != T:
            n1 = configWF_i["lammps"]['eqsteps_nvt_heating']
        else:
            n1 = 0

        n2 = configWF_i["lammps"]['eqsteps_nvt']
        n3 = configWF_i["lammps"]['eqsteps_npt']


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

        w = configWF_i["lammps"].get("window_size", 20)  # window size for moving average

        # determine units for plotting
        units_type = configWF_i["lammps"].get("units")

        if units_type == "real":
            units = ["Temperature (K)", "Energy per atom (kcal/mol)", "Lattice constant (A)", "Volume (A^3)", "Pressure (atm)"]
        elif units_type == "metal":
            units = ["Temperature (K)", "Energy per atom (eV)", "Lattice constant (A)", "Volume (A^3)", "Pressure (bars)"]
        elif units_type == "si":
            units = ["Temperature (K)", "Energy per atom (J)", "Lattice constant (m)", "Volume (m^3)", "Pressure (Pa)"]
        elif units_type == "cgs":
            units = ["Temperature (K)", "Energy per atom (erg)", "Lattice constant (cm)", "Volume (cm^3)", "Pressure (barye)"]
        elif units_type == "electron":
            units = ["Temperature (K)", "Energy per atom (Hartrees)", "Lattice constant (Bohr)", "Volume (Bohr^3)", "Pressure (Pa)"]
        elif units_type == "micro": 
            units = ["Temperature (K)", "Energy per atom (fJ)", "Lattice constant (μm)", "Volume (μm^3)", "Pressure (kPa)"]
        elif units_type == "nano":
            units = ["Temperature (K)", "Energy per atom (zJ)", "Lattice constant (nm)", "Volume (nm^3)", "Pressure (MPa)"]
        else:
            if configWF_i["lammps"].get("units_array") is not None:
                units = configWF_i["lammps"]["units_array"]
            else:
                print("Unknown units type. Using no units.")
                units = ["Temperature", "Energy", "Lattice constant", "Volume", "Pressure"]


        ### ### ### ### ### ### ### 
        ### temperature
        ### ### ### ### ### ### ### 

        # calculate statistical temperature values for temperature equilibration criterion
        T_avg_nvt = np.mean(T[(n1 < step) & (step <= n1 + n2)][-w:])
        T_std_nvt = np.std(T[(n1 < step) & (step <= n1 + n2)][-w:])
        T_avg_npt = np.mean(T[(n1 + n2 < step) & (step <= n1 + n2 + n3)][-w:])
        T_std_npt = np.std(T[(n1 + n2 < step) & (step <= n1 + n2 + n3)][-w:])

        # temperature equilibration criterion
        dT_nvt = T_avg_nvt * np.sqrt(1.5 / N_atoms)
        dT_npt = T_avg_npt * np.sqrt(1.5 / N_atoms)

        print("NVT:")
        print("Average temperature: ", T_avg_nvt, " +/- ", T_std_nvt)
        print("Thermal fluctuation: ", dT_nvt)

        print("NPT:")
        print("Average temperature: ", T_avg_npt, " +/- ", T_std_npt)
        print("Thermal fluctuation: ", dT_npt)
    
        if "temperature" in configWF:
            T_set = configWF_i["temperature"]  
        else:
            T_set = configWF_i["lammps"]["T"]

        error = False

        # check for temperature equilibration criterion
        if (dT_nvt < T_std_nvt):
            print("dT_nvt < T_std_nvt: ", dT_nvt, " < ", T_std_nvt)
            error = true
        if (dT_npt < T_std_npt):
            print("dT_npt < T_std_npt: ", dT_npt, " < ", T_std_npt)
            error = true
        if (abs(T_avg_nvt - T_set) > T_std_nvt):
            print("abs(T_avg_nvt - T_set) > T_std_nvt: ", abs(T_avg_nvt - T_set), " > ", T_std_nvt)
            error = true
        if (abs(T_avg_npt - T_set) > T_std_npt):
            print("abs(T_avg_npt - T_set) > T_std_npt: ", abs(T_avg_npt - T_set), " > ", T_std_npt)
            error = true

        if error:    
            raise Exception("Trajectory is NOT equilibrated with respect to temperature.")
        else:
            print("Trajectory is equilibrated with respect to temperature.")

        # plot temperature data
        fig1, (ax1) =  plt.subplots()
        fig1.suptitle("Temparture")
        plt.xlabel("Step")
        plt.ylabel(units[0])
        plt.plot(step, T, color='blue', label="Temperature")
        plt.plot(moving_average(step, w), moving_average(T, w), color='orange', label="moving average")
        plt.axhline(y=T_avg_nvt, xmin=x1, xmax=x2, color='red', label="Average")
        plt.axhline(y=T_avg_nvt+T_std_nvt, xmin=x1, xmax=x2, linestyle='--', color='red', label="Std")
        plt.axhline(y=T_avg_nvt-T_std_nvt, xmin=x1, xmax=x2, linestyle='--', color='red')
        plt.axhline(y=T_avg_npt, xmin=x2, xmax=x3, color='red')
        plt.axhline(y=T_avg_npt+T_std_npt, xmin=x2, xmax=x3, linestyle='--', color='red')
        plt.axhline(y=T_avg_npt-T_std_npt, xmin=x2, xmax=x3, linestyle='--', color='red')
        if n1 != 0:
            plt.axvline(x=n1, linestyle='--', color="black")
        plt.axvline(x=n1+n2, linestyle='--', color="black")
        plt.axvline(x=n1+n2+n3, linestyle='--', color="black")

        plt.xlim(0, step[-1])
        plt.legend()
        plt.savefig(str(dir_MD / "test_equilibrate/T.pdf"))
        plt.close()


        ### ### ### ### ### ### ### 
        ### energy
        ### ### ### ### ### ### ### 

        # calculate statistical energy values
        E_avg_nvt = np.mean(E[(n1 < step) & (step <= n1 + n2)][-w:])
        E_std_nvt = np.std(E[(n1 < step) & (step <= n1 + n2)][-w:])
        E_avg_npt = np.mean(E[(n1 + n2 < step) & (step <= n1 + n2 + n3)][-w:])
        E_std_npt = np.std(E[(n1 + n2 < step) & (step <= n1 + n2 + n3)][-w:])

        print("\n")
        print("Average energy (NVT): ", E_avg_nvt, " +/- ", E_std_nvt)
        print("Average energy (NPT): ", E_avg_npt, " +/- ", E_std_npt)

        # plot energy data
        fig1, (ax1) =  plt.subplots()
        fig1.suptitle("Energy per atom")
        plt.xlabel("Step")
        plt.ylabel(units[1])
        plt.plot(step, E/N_atoms, color='blue', label="energy")
        plt.plot(moving_average(step, w), moving_average(E, w)/N_atoms, color='orange', label="moving average")
        plt.axhline(y=E_avg_nvt/N_atoms, xmin=x1, xmax=x2, color='red', label="Average")
        plt.axhline(y=(E_avg_nvt+E_std_nvt)/N_atoms, xmin=x1, xmax=x2, linestyle='--', color='red', label="Std")
        plt.axhline(y=(E_avg_nvt-E_std_nvt)/N_atoms, xmin=x1, xmax=x2, linestyle='--', color='red')
        plt.axhline(y=E_avg_npt/N_atoms, xmin=x2, xmax=x3, color='red')
        plt.axhline(y=(E_avg_npt+E_std_npt)/N_atoms, xmin=x2, xmax=x3, linestyle='--', color='red')
        plt.axhline(y=(E_avg_npt-E_std_npt)/N_atoms, xmin=x2, xmax=x3, linestyle='--', color='red')
        if n1 != 0:
            plt.axvline(x=n1, linestyle='--', color="black")
        plt.axvline(x=n1+n2, linestyle='--', color="black")
        plt.axvline(x=n1+n2+n3, linestyle='--', color="black")

        plt.xlim(0, step[-1])
        plt.legend()
        plt.savefig(str(dir_MD / "test_equilibrate/E.pdf"))
        plt.close()


        ### ### ### ### ### ### ### 
        ### volume
        ### ### ### ### ### ### ### 

        # calculate statistical volume values
        V_avg_npt = np.mean(V[(n1 + n2 < step) & (step <= n1 + n2 + n3)][-w:])
        V_std_npt = np.std(V[(n1 + n2 < step) & (step <= n1 + n2 + n3)][-w:])

        print("\n")
        print("Average volume (NPT): ", V_avg_npt, " +/- ", V_std_npt)
   
        # plot volume data
        fig1, (ax1) =  plt.subplots()
        fig1.suptitle("Volume")
        plt.xlabel("Step")
        plt.ylabel(units[3])
        plt.plot(step, V, color='blue', label="Volume")
        plt.plot(moving_average(step, w), moving_average(V, w), color='orange', label="moving average")
        plt.axhline(y=V_avg_npt, xmin=x2, xmax=x3, color='red', label="Average")
        plt.axhline(y=V_avg_npt+V_std_npt, xmin=x2, xmax=x3, linestyle='--', color='red', label="Std")
        plt.axhline(y=V_avg_npt-V_std_npt, xmin=x2, xmax=x3, linestyle='--', color='red')
        if n1 != 0:
            plt.axvline(x=n1, linestyle='--', color="black")
        plt.axvline(x=n1+n2, linestyle='--', color="black")
        plt.axvline(x=n1+n2+n3, linestyle='--', color="black")

        plt.xlim(0, step[-1])
        plt.legend()
        plt.savefig(str(dir_MD / "test_equilibrate/V.pdf"))
        plt.close()


        ### ### ### ### ### ### ### 
        ### lattice constants
        ### ### ### ### ### ### ### 

        # calculate statistical lattice constants values
        L_avg_npt = np.mean(L[:, (n1 + n2 < step) & (step <= n1 + n2 + n3)][-w:], axis=1)
        L_std_npt = np.std(L[:, (n1 + n2 < step) & (step <= n1 + n2 + n3)][-w:], axis=1)

        print("\n")
        for i in range(3):
            print(f"Average lattice constant {i+1} (NPT): ", L_avg_npt[i], " +/- ", L_std_npt[i])

            # plot lattice constants data
            fig1, (ax1) =  plt.subplots()
            fig1.suptitle(f"Lattice constant {i+1}")
            plt.xlabel("Step")
            plt.ylabel(units[2])
            plt.plot(step, L[i], color="blue", label=f"lattice constant {i+1}")
            plt.plot(moving_average(step, w), moving_average(L[i], w), color='orange', label="moving average")
            plt.axhline(y=L_avg_npt[i], xmin=x2, xmax=x3, color="red", label="Average")
            plt.axhline(y=L_avg_npt[i]+L_std_npt[i], xmin=x2, xmax=x3, linestyle='--', color="red", label="Std")
            plt.axhline(y=L_avg_npt[i]-L_std_npt[i], xmin=x2, xmax=x3, linestyle='--', color="red")
            if n1 != 0:
                plt.axvline(x=n1, linestyle='--', color="black")
            plt.axvline(x=n1+n2, linestyle='--', color="black")
            plt.axvline(x=n1+n2+n3, linestyle='--', color="black")

            plt.xlim(0, step[-1])
            plt.legend()
            plt.savefig(str(dir_MD / f"test_equilibrate/L{i+1}.pdf"))
            plt.close()


        ### ### ### ### ### ### ### 
        ### pressure
        ### ### ### ### ### ### ### 

        # calculate statistical pressure values
        p_avg_nvt = np.mean(p[(n1 < step) & (step <= n1 + n2)][-w:])
        p_std_nvt = np.std(p[(n1 < step) & (step <= n1 + n2)][-w:])
        p_avg_npt = np.mean(p[(n1 + n2 < step) & (step <= n1 + n2 + n3)][-w:])
        p_std_npt = np.std(p[(n1 + n2 < step) & (step <= n1 + n2 + n3)][-w:])

        print("\n")
        print("Average pressure (NVT): ", p_avg_nvt, " +/- ", p_std_nvt)
        print("Average pressure (NPT): ", p_avg_npt, " +/- ", p_std_npt)

        # plot pressure data
        fig1, (ax1) =  plt.subplots()
        fig1.suptitle("Pressure")
        plt.xlabel("Step")
        plt.ylabel(units[4])
        plt.plot(step, p, color='blue', label="Pressure")
        plt.plot(moving_average(step, w), moving_average(p, w), color='orange', label="moving average")
        plt.axhline(y=p_avg_nvt, xmin=x1, xmax=x2, color='red', label="Average")
        plt.axhline(y=p_avg_nvt+p_std_nvt, xmin=x1, xmax=x2, linestyle='--', color='red', label="Std")
        plt.axhline(y=p_avg_nvt-p_std_nvt, xmin=x1, xmax=x2, linestyle='--', color='red')
        plt.axhline(y=p_avg_npt, xmin=x2, xmax=x3, color='red')
        plt.axhline(y=p_avg_npt+p_std_npt, xmin=x2, xmax=x3, linestyle='--', color='red')
        plt.axhline(y=p_avg_npt-p_std_npt, xmin=x2, xmax=x3, linestyle='--', color='red')
        if n1 != 0:
            plt.axvline(x=n1, linestyle='--', color="black")
        plt.axvline(x=n1+n2, linestyle='--', color="black")
        plt.axvline(x=n1+n2+n3, linestyle='--', color="black")

        plt.xlim(0, step[-1])
        plt.legend()
        plt.savefig(str(dir_MD / "test_equilibrate/p.pdf"))
        plt.close()

        p_set = configWF_i["lammps"].get("P", 0.0)

        # check if pressure is equilibrated
        if (abs(p_avg_npt - p_set) > p_std_npt):
            print("abs(T_avg_npt - T_set) > T_std_npt: ", abs(p_avg_npt - p_set), " > ", p_std_npt)
            error = true
        #if not np.allclose(moving_average(p, w)[-w:], p_set, atol=abs(p_set + p_std_npt)/p_set):
        #    print("Pressure is fluctuating too much around the set point.")
        #    error = true

        if error:    
            raise Exception("Trajectory is NOT equilibrated with respect to pressure.")
        else:
            print("Trajectory is equilibrated with respect to pressure.")

        # human_in_loop flag allows for human check before continuing
        if configWF_i.get("human_in_loop", False):
            print("human_in_loop is set to True. Stopping workflow after each test.")
            raise Exception("Equilibration test done. Please check the plots in " + str(dir_MD / "test_equilibrate/") + " and continue the workflow manually.")

    else:
        print("No equilibration was performed.")


    return True, {CYCLICALGROUP_KEY: cg_criteria, "path_configWF": path_configWF, "num_simulations": num_simulations}
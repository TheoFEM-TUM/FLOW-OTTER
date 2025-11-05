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


    with open(path_configWF, "r") as f:
        configWF = yaml.safe_load(f)

    dir_project = Path(configWF.get("dir_project", "./"))
    simulation_type = configWF.get("simulation_type", "sweep")
    cg_criteria = False

    if num_simulations > 1:

        if simulation_type == "cascade":
            i = kwargs['pq_iteration'][0]
        else:
            i = kwargs['pq_index'][0]

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

    dir_MD = dir_project_i / configWF_i.get("dir_MD", "1-MD/")
    equilibrate = configWF_i.get("equilibrate", False)


    if equilibrate:

        if "temperature" in configWF_i:
            T = configWF_i["temperature"]  
        else:
            T = configWF_i["lammps"]["T"]

        T_start = configWF_i["lammps"].get("T_start", T)


        if T_start != T:
            n1 = configWF_i["lammps"]['eqsteps_nvt_heating']
        else:
            n1 = 0

        n2 = configWF_i["lammps"]['eqsteps_nvt']
        n3 = configWF_i["lammps"]['eqsteps_npt']


        dir_test = dir_MD / "test_equilibrate/"
        dir_test.mkdir(parents=True, exist_ok=True)
        #os.system(f"mkdir -p {dir_MD}/test_equilibrate/")

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

        w = configWF_i["lammps"].get("window_size", 50)  # window size for moving average

        # temperature
        T_avg_nvt = np.mean(T[(n1 < step) & (step <= n1 + n2)][-w:])
        T_std_nvt = np.std(T[(n1 < step) & (step <= n1 + n2)][-w:])
        T_avg_npt = np.mean(T[(n1 + n2 < step) & (step <= n1 + n2 + n3)][-w:])
        T_std_npt = np.std(T[(n1 + n2 < step) & (step <= n1 + n2 + n3)][-w:])

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

        fig1, (ax1) =  plt.subplots()
        fig1.suptitle("Temparture")
        plt.xlabel("Step")
        plt.ylabel("Temperature (K)")
        plt.plot(step, T, color='blue', label="Temperature")
        plt.axhline(y=T_avg_nvt, xmin=x1, xmax=x2, color='red', label="Average")
        plt.plot(moving_average(step, w), moving_average(T, w), color='orange', label="moving average")
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


        # energy
        E_avg_nvt = np.mean(E[(n1 < step) & (step <= n1 + n2)][-w:])
        E_std_nvt = np.std(E[(n1 < step) & (step <= n1 + n2)][-w:])
        E_avg_npt = np.mean(E[(n1 + n2 < step) & (step <= n1 + n2 + n3)][-w:])
        E_std_npt = np.std(E[(n1 + n2 < step) & (step <= n1 + n2 + n3)][-w:])

        print("\n")
        print("Average energy (NVT): ", E_avg_nvt, " +/- ", E_std_nvt)
        print("Average energy (NPT): ", E_avg_npt, " +/- ", E_std_npt)
   
        fig1, (ax1) =  plt.subplots()
        fig1.suptitle("Energy per atom")
        plt.xlabel("Step")
        plt.ylabel("Energy per atom")
        plt.plot(step, E/N_atoms, color='blue', label="energy")
        plt.plot(moving_average(step, w), moving_average(E, w), color='orange', label="moving average")
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


        # volume
        V_avg_npt = np.mean(V[(n1 + n2 < step) & (step <= n1 + n2 + n3)][-w:])
        V_std_npt = np.std(V[(n1 + n2 < step) & (step <= n1 + n2 + n3)][-w:])

        print("\n")
        print("Average volume (NPT): ", V_avg_npt, " +/- ", V_std_npt)
   
        fig1, (ax1) =  plt.subplots()
        fig1.suptitle("Volume")
        plt.xlabel("Step")
        plt.ylabel("Volume")
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

        # lattice constants
        L_avg_npt = np.mean(L[:, (n1 + n2 < step) & (step <= n1 + n2 + n3)][-w:], axis=1)
        L_std_npt = np.std(L[:, (n1 + n2 < step) & (step <= n1 + n2 + n3)][-w:], axis=1)

        print("\n")
        for i in range(3):
            print(f"Average lattice constant {i+1} (NPT): ", L_avg_npt[i], " +/- ", L_std_npt[i])
   
            fig1, (ax1) =  plt.subplots()
            fig1.suptitle(f"Lattice constant {i+1}")
            plt.xlabel("Step")
            plt.ylabel(f"Lattice constant {i+1}")
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


        # pressure
        p_avg_nvt = np.mean(p[(n1 < step) & (step <= n1 + n2)][-w:])
        p_std_nvt = np.std(p[(n1 < step) & (step <= n1 + n2)][-w:])
        p_avg_npt = np.mean(p[(n1 + n2 < step) & (step <= n1 + n2 + n3)][-w:])
        p_std_npt = np.std(p[(n1 + n2 < step) & (step <= n1 + n2 + n3)][-w:])

        print("\n")
        print("Average pressure (NVT): ", p_avg_nvt, " +/- ", p_std_nvt)
        print("Average pressure (NPT): ", p_avg_npt, " +/- ", p_std_npt)
   
        fig1, (ax1) =  plt.subplots()
        fig1.suptitle("Pressure")
        plt.xlabel("Step")
        plt.ylabel("Pressure")
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

        if configWF_i.get("human_in_loop", False):
            print("human_in_loop is set to True. Stopping workflow after each test.")
            raise Exception("Equilibration test done. Please check the plots in " + str(dir_MD / "test_equilibrate/") + " and continue the workflow manually.")

    else:
        print("No equilibration was performed.")


    return True, {CYCLICALGROUP_KEY: cg_criteria, "path_configWF": path_configWF, "num_simulations": num_simulations}
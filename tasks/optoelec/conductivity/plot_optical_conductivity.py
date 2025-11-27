import subprocess
from typing import Tuple
from pathlib import Path


def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:

    with open(path_configWF, "r") as f:
        configWF = yaml.safe_load(f)

    dir_project = Path(configWF.get("dir_project", "./"))
    dir_conductivity = Path(configWF["dir_conductivity"])

    if num_simulations > 1:

        param_to_vary = configWF["param_to_vary"]
        array_to_vary = configWF["array_to_vary"]
            
        fig1, (ax1) = plt.subplots()
        plt.title(f"Average density of states")

        outputs = []

        for i in range(len(array_to_vary)):

            dir_project_i = dir_project / f"{param_to_vary}_{array_to_vary[i]}/"

            with open(str(dir_project_i / 'branch_config.yaml'), 'r') as f:
                configWF_i = yaml.safe_load(f)

            dir_config = dir_project / configWF.get("dir_config", "3-conductivity/")
            dir_output = Path(configWF_i.get("dir_output", dir_config))  

            outputs.append(dir_output)


        result = subprocess.run([str(dir_conductivity / "plotting_scripts/multi_plot/multi_plot.sh"), outputs], check=True)

    else:
        print("No comparison of different conductivity output needed since only one calculation was performed.")


    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations}
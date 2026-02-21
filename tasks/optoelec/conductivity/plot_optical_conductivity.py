import subprocess
from typing import Tuple
from pathlib import Path


def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:

    print("Start task: plot_optical_conductivity", flush=True)

    with open(path_configWF, "r") as f:
        configWF = yaml.safe_load(f)

    dir_project = Path(configWF.get("dir_project", "./"))
    dir_conductivity = Path(configWF["dir_conductivity"])

    if num_simulations > 1:

        param_to_vary = configWF["param_to_vary"]
        array_to_vary = configWF["array_to_vary"]

        outputs = []
        labels = []
        dir_plots = dir_project / "plots/optoelec/conductivity/"

        dir_plots.mkdir(parents=True, exist_ok=True)

        for i in range(len(array_to_vary)):

            dir_project_i = dir_project / f"{param_to_vary}_{array_to_vary[i]}/"

            with open(str(dir_project_i / 'branch_config.yaml'), 'r') as f:
                configWF_i = yaml.safe_load(f)

            dir_config = Path(configWF_i.get("dir_config", str(dir_project_i / "3-conductivity/")))
            dir_output = Path(configWF_i.get("dir_output", str(dir_project_i /"3-conductivity/")))

            unit_to_vary = configWF_i.get("unit_to_vary", "")

            outputs.append(dir_output)
            labels.append(f"{param_to_vary} = {array_to_vary[i]} {unit_to_vary}")


        result = subprocess.run([
            str(dir_conductivity / "plotting_scripts/cluster_plot.sh"), 
            outputs,
            "--labels",
            *labels,
            "--outfolder",
            str(dir_plots)
            ], check=True)

    else:
        print("No comparison of different conductivity output needed since only one calculation was performed.", flush=True)

    print("Finish task: plot_optical_conductivity", flush=True)

    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations}
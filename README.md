# MD_TB_PQ_workflow



Here you find the definition of parameters that can be used within the workflow manager:

example:
**param** (type): meaning (*"default_option"*, "option1", "option2", ...)


global:
**dir_project** (str): directory of project output (*"./"*)
**dir_code** (str): directory of this repo 

**resources_MD** (str): myqueue string specifying resources for MD jobs 
**resources_H**  (str): myqueue string specifying resources for H jobs 
**resources_optoelec** (str): myqueue string specifying resources for optoelec jobs 
**resources_instant** (str): myqueue string specifying resources for jobs which should be finished instantaneously (*resources of resources_H with walltime = 5m)
**resources_short** (str): myqueue string specifying resources for jobs which should only run shortly (*resources of resources_H with walltime = 2h)
**resources_long** (str): myqueue string specifying resources for jobs which should run very long (*resources of resources_H with walltime = 1d)

**simulation_type** (str):  (*"sweep"*, "cascade")
**num_simulations** (int): (*1*)
**param_to_vary** (str):
**param_group_for_vary** (str): 
**array_to_vary** (array of str): 
**simulation_index** (int):  (internally set by PQ)





**N_avg** (int): (*1)

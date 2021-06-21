from io import StringIO
import pathlib
import numpy as np
from os import path
import os, os.path, sys
import subprocess
from itertools import count, product

from tqdm import tqdm 
from bsx_data import study_ranges, normal_ranges, folder_param_dict


# Edit this
MODEL = '/Users/mirunaserian/Desktop/PHD_Rotation2/BCMD/build/bsx.model'
VALUE_STEPS = 50
increase_range = 50
save_folder = f"/Users/mirunaserian/Desktop/PHD_Rotation2/haemodynamics-metabolism/bsx/new14jun_comb_x{increase_range}"

def create_input(init_params, changing_params, folder_name, param_type, count=None):
    folder_name = folder_name.replace("increase", "i")
    folder_name = folder_name.replace("decrease", "d")
    folder_name = folder_name.replace("constant_modified", "c")

    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
    change_in_params = {k:(changing_params[k] - v) for k, v in init_params.items()}
    if param_type == "variable":
        init_params = {k: init_params[k]+v for k, v in change_in_params.items()}
        init_params["u"] = 1.25
        change_in_params = {k: 0 if k!="u" else 1.5-init_params["u"] \
            for k in init_params.keys() }

    f_out = "@ 201 \n"
    f_out = f_out + f': {len(init_params.keys())} ' + ' '.join([*init_params]) + '\n'
    f_out = f_out + f">>> {len(output)} " + ' '.join(output) + '\n'
    f_out = f_out + "!!!"+'\n'
    f_out = f_out + f"= 0 0 " +' '.join(str(v) for k, v in init_params.items())+ '\n' 
    f_out = f_out + "!!!"+'\n'
    f_out = f_out + f"* 120 1 " + " ".join(["0"] *len(init_params.keys())) +'\n'
    f_out = f_out + f"* 1 1 " + ' '.join(str(v) for k, v in change_in_params.items())+ '\n'
    f_out = f_out + f"* 29 1 " + " ".join(["0"] *len(init_params.keys())) +'\n'
    f_out = f_out + f"* 1 1 " + ' '.join(str(-v) for k, v in change_in_params.items())+ '\n'
    f_out = f_out + f"* 49 1 " + " ".join(["0"] *len(init_params.keys())) +'\n'
    f_out = f_out + "!!!"+'\n'
    
    file_name = count


    if path.exists(f"{folder_name}/bsx_{file_name}.input") is False:
        text_file = open(f"{folder_name}/bsx_{file_name}.input", "wt")
        text_file.write(f_out)
        text_file.close()
        return file_name


output =  \
    ["t",  "P_a", "Pa_CO2", "SaO2sup", "u", "r", "H", "O", "a", "CCO", \
        "Vmca", "CBF", "TOI", "DHbT", "HbO2", "HHb", "CMRO2", *normal_ranges.keys()]

changed_values = {
    "u": 1.5
}
def sc_type(study_ranges, action_type, p):
    if "i_" in action_type or  "d_" in action_type:
        x = action_type.split("_")[-1]
        x = float(x)/100
        study_ranges =  {
        "SaO2sup": (0.9, 0.96),
        **{k: ((1-x)*normal_ranges[k], (1+x)*normal_ranges[k]) \
            for k in normal_ranges.keys() if k is not "SaO2sup"}
        } 

    if "i_" in action_type or "increase" in action_type:
        return [normal_ranges[p]+ x for x in np.arange(0, (study_ranges[p][1]-normal_ranges[p])+1, (study_ranges[p][1]-normal_ranges[p])/VALUE_STEPS)]
    elif "d_" in action_type or "decrease" in action_type:
        return [normal_ranges[p]- x for x in np.arange(0, (normal_ranges[p]-study_ranges[p][0])+1, (normal_ranges[p]-study_ranges[p][0])/VALUE_STEPS)]
    elif action_type == "constant_modified":
        return [changed_values[p] for i in range(VALUE_STEPS)]
    elif action_type=="constant":
        return[normal_ranges[p] for i in range(VALUE_STEPS)]


param_list = ["phi", "kCV", "k_aut", "P_a", "Pa_CO2"]
action_list = ["i_25", "i_50", "d_25", "d_50"]


# generating combinations
temp = product(action_list, repeat = len(param_list))

# constructing dicts using combinations
action_type = [{key : val for (key , val) in zip(param_list, ele)} for ele in temp]


def get_scenarios(action_type, param_type, increase_range, sa_params):
    if param_type == "variable" and sa_params.any():
        changing_params = [
            {k: sa_params[i][list(action_type.keys()).index(k)-1] if k!="u" else 1.5 \
                    for k, v in action_type.items()} \
                    for i in range(0, sa_params.shape[0])
                    ]
    else:
        changing_params = [
            {k: sc_type(study_ranges, v, k)[i] for k, v in action_type.items()} \
                    for i in range(0, VALUE_STEPS)
        ]
    return {
        f"{('_').join(action_type.keys())}": {
            "init_params": {k: normal_ranges[k] for k in action_type.keys()},
            "changing_params": changing_params
        } }  
        

# Depending on the simulation type, set action_type_list to one of the lists below
action_type_params_u_constant = [*[{k: "increase", "u": "constant_modified"} for k in normal_ranges.keys()],
                     *[{k: "decrease", "u": "constant_modified"} for k in normal_ranges.keys()]]

action_type_params = [*[{k: "increase"} for k in normal_ranges.keys() if k is not "SaO2sup"],
                     *[{k: "decrease"} for k in normal_ranges.keys()]]

action_type_combinations = [{**d, "u": "constant_modified"} for d in action_type]
def run_model(MODEL, input_file, output_file):
    subprocess.call(
        [MODEL, 
        '-i', input_file, 
        '-o', output_file]
        )


def run_sims(action_type_list):
    param_type = "variable"

    for action in action_type_list:
        for scenario_name, v in tqdm(get_scenarios(action, param_type, increase_range,  np.empty(shape=[0, 0])).items()):
            action_folder = "-".join(f"{k}_{val}" for k, val in action.items() if k in v["init_params"].keys())
            action_folder = action_folder.replace("increase", "i")
            action_folder = action_folder.replace("decrease", "d")
            action_folder=action_folder.replace("constant_modified", "c")
            action_folder = save_folder + "/" + scenario_name + "/"  + action_folder
            if not os.path.exists(action_folder):
                count = 0
                for c_p in v["changing_params"]:
                    file_name = create_input(v["init_params"],c_p, action_folder, param_type, count)
                    run_model(MODEL,os.path.join(scenario_name, f"{action_folder}/bsx_{file_name}.input"),
                    os.path.join(scenario_name, f"{action_folder}/bsx_{file_name}.csv" ))
                    count +=1 
                
if __name__ == '__main__':
    action_type_list = action_type_combinations
    run_sims(action_type_list)
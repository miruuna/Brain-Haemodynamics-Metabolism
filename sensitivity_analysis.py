import SALib.sample.morris
import SALib.analyze.morris
import SALib.analyze.sobol

import SALib.sample.fast_sampler
import SALib.analyze.fast
import SALib.sample.saltelli

import sys
import json 

import os

from scipy import stats
from bsx_data import folder_param_dict, study_ranges, normal_ranges
import numpy as np
import pandas as pd
from os import path
from create_input import create_input, get_scenarios, run_model
import matplotlib.pyplot as plt
from tqdm import tqdm 
from save_datasets import save_indiv_data
import math

FOLDER_BASE = f"/Users/mirunaserian/Desktop/PHD_Rotation2/haemodynamics-metabolism/bsx"
range_type = 35


def truncate(number, decimals=0):
    """
    Return a value truncated to a specific number of decimal places.
    """
    if not isinstance(decimals, int):
        raise TypeError("decimal places must be an integer.")
    elif decimals < 0:
        raise ValueError("decimal places has to be 0 or more.")
    elif decimals == 0:
        return math.trunc(number)

    factor = 10.0 ** decimals
    return math.trunc(number * factor) / factor

def plot_sobol_indices(Si, problem, folder_base, folder_name, title=None, factor_names=None):

    # Get indices
    factors = problem['names']
    s1 = Si['S1']
    s1_conf = Si['S1_conf']
    st = Si['ST']
    st_conf = Si['ST_conf']

    # Set up plot
    fig, ax = plt.subplots(nrows=1, ncols=1)
    
    x_shift = 0.08
    x_vals_1 = np.arange(0 - x_shift, len(s1) - x_shift, 1)
    x_vals_2 = np.arange(0 + x_shift, len(st) + x_shift, 1)

    # Plot indices
    ax.errorbar(x_vals_1, s1, yerr=s1_conf, marker='o', markersize=9,
        linewidth=0.8, linestyle='None', color='k', capsize=0,
        label='first order', markerfacecolor='w', markeredgewidth=0.8)
    ax.errorbar(x_vals_2, st, yerr=st_conf, marker='s', markersize=9,
        linewidth=0.8, linestyle='None', color='k', capsize=0,
        label='total order', markeredgewidth=0.8)

    if title != None:
        ax.set_title(title)

    ax.legend(handletextpad=0.001)

    # Set x axis.

    if factor_names == None:
        xlabels = factors
    else:
        xlabels = factor_names
    
    ax.set_xticks(range(len(factors)))
    ax.set_xticklabels(xlabels)
    ax.set_xlabel('factors')

    # Set y axis.

    ax.set_ylim([-0.1, 1.1])
    ax.set_yticks(np.linspace(0, 1, 3))
    ax.set_yticks([0.25, 0.75], minor=True)
    plt.setp(ax.get_yminorticklabels(), visible=False)
    ax.set_ylabel('Sobol index')

    # Set horizontal lines.

    for y in np.linspace(0, 1,5):
        ax.axhline(y=y, linestyle=':', color='0.8', linewidth=0.5)

    plt.tight_layout()
    fig.savefig(f'{folder_base}/{folder_name}.jpg', dpi=fig.dpi)


def get_si(df, params, problem, sample, sa_type, y_type):
    """
    Return numpy array with sensitivity indexes.

    """
    Y = np.zeros([sample.shape[0]])
    for i, X in enumerate(sample):
        x=X[0]
        if not df[df[params[0]].astype(str).str.contains(str(truncate(x, 7)))].empty:
            if y_type in ["impaired_hhb", "impared_hbo2"]:
                signal_type = y_type.split("_")[-1]
                if df[df[params[0]].astype(str).str.contains(str(truncate(x, 7)))]["rPWR_hhb"].tolist()[0] == 0:
                    Y[i] = df[df[params[0]].astype(str).str.contains(str(truncate(x, 7)))]["rCST_hhb"].tolist()[0]
                else:
                    Y[i] = df[df[params[0]].astype(str).str.contains(str(truncate(x, 7)))]["rCST_hhb"].tolist()[0]/df[df[params[0]].astype(str).str.contains(str(truncate(x, 7)))]["rPWR_hhb"].tolist()[0]
            elif y_type == "impaired_hbo2":
                if df[df[params[0]].astype(str).str.contains(str(truncate(x, 7)))]["rPWR_hbo2"].tolist()[0] == 0:
                    Y[i] = df[df[params[0]].astype(str).str.contains(str(truncate(x, 7)))]["rCST_hbo2"].tolist()[0]
                else:
                    Y[i] = df[df[params[0]].astype(str).str.contains(str(truncate(x, 7)))]["rCST_hbo2"].tolist()[0]/df[df[params[0]].astype(str).str.contains(str(truncate(x, 7)))]["rPWR_hbo2"].tolist()[0]
            else:
                if len(df[df[params[0]].astype(str).str.contains(str(truncate(x, 7)))][y_type].tolist()) !=0:
                    Y[i] = df[df[params[0]].astype(str).str.contains(str(truncate(x, 7)))][y_type].tolist()[0]
        else:
            continue
        

    if sa_type == "morris":
        Si = SALib.analyze.morris.analyze(problem, sample, Y, num_levels=10)
    elif sa_type == "sobol":
        Si = SALib.analyze.sobol.analyze(problem, Y)
    elif sa_type == "fast":
        Si = SALib.analyze.fast.analyze(problem, Y)
    return Si

def get_problem(params):
    return {
        'num_vars': 1,
        'names': params,
        'groups': None,
        'bounds': [[study_ranges[p][0], study_ranges[p][1]] for p in params]
    }

MODEL = '/Users/mirunaserian/Desktop/PHD_Rotation2/BCMD/build/bsx.model'

def run_Sa_analysis(sa_type):
    # set the folder where data will be stored
    folder_name = f"new17jun_comb_sensitivity_analysis_{sa_type}_x{range_type}"

    save_folder = \
        f"/Users/mirunaserian/Desktop/PHD_Rotation2/haemodynamics-metabolism/bsx/{folder_name}"
    param_list = ["P_v", "r_n", "R_autu","C_im", "R_auto","phi", "kCV", "k_aut", "P_a", "Pa_CO2", "SaO2sup"]
    folder_f= ("_").join(param_list)

    folder_param_dict = {p:[p] for p in param_list}
    sensitivity_dict = {}
    for f, p in folder_param_dict.items():
        sensitivity_dict = {param: {} for param in p}

        problem = {
            'num_vars': len(p),
            'names': p,
            'groups': None,
            'bounds': [[study_ranges[param][0], study_ranges[param][1]] for param in p]
        }
        directory= f
        file_path = f"{FOLDER_BASE}/{folder_name}/{directory}/{directory}.csv"
        action_type = {"u": "constant", **{param:"modify" for param in p}}

        sa_params = None
        if sa_type == "fast":
            sa_params = SALib.sample.fast_sampler.sample(problem, 100)
        elif sa_type == "morris":
            sa_params = SALib.sample.morris.sample(problem, 100, num_levels=100)
        elif sa_type == "sobol":
            sa_params = SALib.sample.saltelli.sample(problem, 100)

        scenarios = get_scenarios(action_type, "variable", [], sa_params)
        for scenario_name, v in tqdm(scenarios.items()):
            action_folder = "_".join(f"{k}_{val}" for k, val in action_type.items() if k in v["init_params"].keys())
            action_folder = save_folder + "/" + scenario_name + "/"  + action_folder

            count = 0

            for c_p in v["changing_params"]:
                file_name = create_input(v["init_params"],c_p, f"{save_folder}/{f}", "variable", count)
                run_model(MODEL, f"{save_folder}/{f}/bsx_{file_name}.input",
                os.path.join(scenario_name, f"{save_folder}/{f}/bsx_{file_name}.csv" ))
                count +=1 

        df = save_indiv_data(f, p, FOLDER_BASE, folder_name)
        if df is not None:
            for y_type in ["rPWR_hhb", "rPWR_hbo2", "rCST_hhb", "rCST_hbo2", ]:
                Si = get_si(df, p, problem, sa_params, sa_type, y_type)    
                Si_dict = {k: v if isinstance(v, list) else v.tolist() for k,v in Si.items()}
                with open(f"{FOLDER_BASE}/{folder_name}/{directory}.json", 'w') as outfile:
                    json.dump(Si_dict, outfile)
                for i in range(len(p)):
                    if sa_type == "sobol":
                        plot_sobol_indices(Si, problem, f"{FOLDER_BASE}/{folder_name}", directory, title=f"{f} - {y_type}")
                        sensitivity_dict[p[i]][y_type] = {"S1": Si_dict["S1"][i], "ST": Si_dict["ST"][i]}
                    elif sa_type == "morris":
                        sensitivity_dict[p[i]][y_type] = {
                            "mu": Si_dict["mu"][i], 
                            "mu_star": Si_dict["mu_star"][i],
                            "mu_star_conf": Si_dict["mu_star_conf"][i],
                            "sigma": Si_dict["sigma"][i]}
                    elif sa_type == "fast":
                        sensitivity_dict[p[i]][y_type] = {"S1": Si_dict["S1"][i], "ST": Si_dict["ST"][i]}
        else:
        df_melt = pd.json_normalize(sensitivity_dict, max_level=2)

        # df_all = pd.DataFrame.from_dict(sensitivity_dict, orient="columns")
        df_melt.to_csv(f"{save_folder}/{f}/all_sensitivities_{('_').join(p)}.csv")


if __name__ == "__main__":
   run_Sa_analysis(sys.argv[1])
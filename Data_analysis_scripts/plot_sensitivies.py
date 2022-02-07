from typing import no_type_check
from bsx_data import folder_param_dict
import pandas as pd
import numpy as no_type_check
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import json
plt.style.use('ggplot')



sa_type = "morris"
folder_base = \
    f"/Users/mirunaserian/Desktop/PHD_Rotation2/haemodynamics-metabolism/bsx/newsensitivity_analysis_{sa_type}_x50"

def unpack(df, column, fillna=None):
    ret = None
    if fillna is None:
        tmp = pd.DataFrame((d for idx, d in df[column].iteritems()))
        ret = pd.concat([df.drop(column,axis=1), tmp], axis=1)
    else:
        tmp = pd.DataFrame((d for idx, d in 
        df[column].iteritems())).fillna(fillna)
        ret = pd.concat([df.drop(column,axis=1), tmp], axis=1)
    return ret

def plot_all_indices_sep(sa_type):
    #  Categorical Data
    a = 6  # number of rows
    b = 4 # number of columns
    c = 1  # initialize plot counter
    fig = plt.figure(figsize=(18,18))
    si_name = "mu_star" if sa_type == "morris" else "S1"
    for f, p in folder_param_dict.items():
        csv_file_name = f"{folder_base}/{f}/all_sensitivities_{('_').join(p)}.csv"
        if len(p) == 1:
            df = pd.read_csv(csv_file_name)
            cols_to_plot = ["rPWR_hhb", "rCST_hhb", "rPWR_hbo2", "rCST_hbo2"]
            vals_to_plot = [df[f"{p[0]}.{col}.{si_name}"].values[0] for col in cols_to_plot]
            if not all(v == 0.000 for v in vals_to_plot):

                plt.subplot(a, b, c)
                plt.title('{}'.format(p[0]))

                x_pos = [i for i, _ in enumerate(cols_to_plot)]
                plt.bar(x_pos, vals_to_plot)
                plt.xticks(x_pos, cols_to_plot)

                c = c + 1

    # df= df.apply(lambda row : row[0], axis=1)

def plot_all_indices(sa_type):
    #  Categorical Data
    a = 6  # number of rows
    b = 4 # number of columns
    c = 1  # initialize plot counter
    d = {}
    si_name = "sigma" if sa_type == "morris" else "S1"
    for f, p in folder_param_dict.items():
        csv_file_name = f"{folder_base}/{f}/all_sensitivities_{('_').join(p)}.csv"
        if len(p) == 1:
            df = pd.read_csv(csv_file_name)
            cols_to_plot = ["rPWR_hhb", "rCST_hhb", "rPWR_hbo2", "rCST_hbo2"]
            vals_to_plot = [df[f"{p[0]}.{col}.{si_name}"].values[0] for col in cols_to_plot]
            if not all(v == 0.000 for v in vals_to_plot):
                d[p[0]] = {}
                for col in cols_to_plot:
                    d[p[0]][col] = df[f"{p[0]}.{col}.{si_name}"].values[0]
    df = pd.DataFrame.from_dict(d)
    fig = plt.figure(figsize=(25,25))
    df.T.plot.barh(figsize=(10,10))
    # df= df.apply(lambda row : row[0], axis=1)


plot_all_indices(sa_type)
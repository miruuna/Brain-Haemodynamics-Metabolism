
from scipy.signal import find_peaks
import numpy as np
from scipy import stats
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec



def rot(angle):
        return np.array([[np.cos(angle), np.sin(angle)],
                    [-np.sin(angle), np.cos(angle)]])

def get_max(data_hhb, data_hbo2):
    # get min of both columns
    min_V = min(np.amin(np.concatenate([data_hhb, data_hbo2], axis=1)),
                np.amin(np.concatenate([data_hhb, data_hbo2], axis=0)))
    max_V = max(np.amax(np.concatenate([data_hhb, data_hbo2], axis=1)), 
                np.amax(np.concatenate([data_hhb, data_hbo2], axis=0)))
    if min_V < 0 :
        min_V = -min_V
    max_val = max([min_V, max_V])
    return max_val

def get_rPWR(df, normalise):
    # start at t=1 so perecentage increases can be calculated
    time = df["t"]
    first_row = df.iloc[[0]].values[0]
    df = df.apply(lambda row: row - first_row, axis=1)

    df = df[1:]
    df["t"] = time[1:]


#   normalise using a baseline of 10 seconds before burst
    df_baseline = df[(df.t.isin(range(114,119)))]
    df_sub = df.sub(df_baseline.mean())
    z_scored_baseline = df_sub.div(df_baseline.std())
    z_scored_baseline.t = time

# use the first 10 seconds after a 5 second buffer
    df = df[(125 <= df.t) & (df.t <= 135)]
    cco = df.CCO
    hhb = df.HHb
    hbo2 = df.HbO2
    if normalise == True:
        z_cco = z_scored_baseline.CCO
        z_hhb = z_scored_baseline.HHb
        z_hbo2 = z_scored_baseline.HbO2

        mean_z_cco = np.mean(cco)
        mean_z_hhb = np.mean(hhb)
        mean_z_hbo2 = np.mean(hbo2)
    else:
        mean_z_cco = np.mean(cco)
        mean_z_hhb = np.mean(hhb)
        mean_z_hbo2 = np.mean(hbo2)
    
    data_hhb = np.array([[mean_z_hhb], [mean_z_cco]])
    data_hbo2 = np.array([[mean_z_hbo2], [mean_z_cco]])
    data_hhb = np.nan_to_num(data_hhb)
    data_hbo2 = np.nan_to_num(data_hbo2)


    rotated_data_hhb = np.matmul(rot((np.pi)/4), data_hhb)
    rotated_data_hbo2 = np.matmul(rot((np.pi)/4), data_hbo2)
 
    rPWR_hhb = rotated_data_hhb[0,:]
    rCST_hhb = rotated_data_hhb[1,:]

    rPWR_hbo2 = rotated_data_hbo2[0,:]
    rCST_hbo2 = rotated_data_hbo2[1,:]

    return {
        "mean_z_cco": mean_z_cco, 
        "mean_z_hhb": mean_z_hhb,
        "mean_z_hbo2": mean_z_hbo2,
        "rPWR_hhb": rPWR_hhb,
        "rPWR_hbo2": rPWR_hbo2,
        "rCST_hhb": rCST_hhb,
        "rCST_hbo2": rCST_hbo2}

def get_data_set(df_list, param, normalise):
    mean_z_cco_list, mean_z_hhb_list, mean_z_hbo2_list, rPWR_hhb_list, \
        rCST_hhb_list, p_val, rPWR_hbo2_list, rCST_hbo2_list = \
            [], [], [], [], [], [], [], []
    
    param_list = param
    param_arr = np.empty((0,len(param_list)), dtype=np.float64)

    for df in df_list:
        non_delta_df = df
        res_dict = get_rPWR(df, normalise)

        # use the first 10 seconds after a 5 second buffer
        df = df[(125 <= df.t) & (df.t <= 135)]

        mean_z_cco_list.append(res_dict["mean_z_cco"])
        mean_z_hhb_list.append(res_dict["mean_z_hhb"])
        mean_z_hbo2_list.append(res_dict["mean_z_hbo2"])
        rCST_hhb_list.append(-res_dict["rCST_hhb"][0])
        rPWR_hhb_list.append(-res_dict["rPWR_hhb"][0])
        rCST_hbo2_list.append(res_dict["rCST_hbo2"][0])
        rPWR_hbo2_list.append(res_dict["rPWR_hbo2"][0])  
        
        
        param_arr = np.append(param_arr, np.array([[df[p].max() for p in param_list]], dtype=np.float64), axis=0)

    df_new = pd.DataFrame({"mean_z_cco": mean_z_cco_list,
                        "mean_z_hhb": mean_z_hhb_list,
                        "mean_z_hbo2": mean_z_hbo2_list,
                        "rCST_hhb": rCST_hhb_list,
                        "rCST_hbo2": rCST_hbo2_list,
                        "rPWR_hhb": rPWR_hhb_list,
                        "rPWR_hbo2": rPWR_hbo2_list,
                        **{param_list[x]: param_arr[:, x] for x in range(len(param_list))}})
    return df_new


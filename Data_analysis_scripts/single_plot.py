from scipy.signal import find_peaks
import numpy as np
from scipy import stats
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

def get_rPWR(df, normalise):
    # start at t=1 so perecentage increases can be calculated
    time = df.t
    df = df.diff()
    df = df[1:]
    df.t = time[1:]


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
    cbf = df.CBF


    # z-score the values
    # z_cco = stats.zscore(cco)
    # z_hhb = stats.zscore(hhb)
    # z_hbo2 = stats.zscore(hbo2)

    if normalise == True:
        z_cco = z_scored_baseline.CCO
        z_hhb = z_scored_baseline.HHb
        z_hbo2 = z_scored_baseline.HbO2
        z_cbf = z_scored_baseline.CBF

        mean_z_cco = np.mean(cco)
        mean_z_hhb = np.mean(hhb)
        mean_z_hbo2 = np.mean(hbo2)
        mean_z_cbf = np.mean(cbf)
    else:
        mean_z_cco = np.mean(cco)
        mean_z_hhb = np.mean(hhb)
        mean_z_hbo2 = np.mean(hbo2)
        mean_z_cbf= np.mean(cbf)
    
    data_hhb = np.array([[mean_z_hhb], [mean_z_cco]])
    data_hbo2 = np.array([[mean_z_hbo2], [mean_z_cco]])

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
        "mean_z_cbf": mean_z_cbf,
        "rPWR_hhb": rPWR_hhb,
        "rPWR_hbo2": rPWR_hbo2,
        "rCST_hhb": rCST_hhb,
        "rCST_hbo2": rCST_hbo2}


def get_PWR_plots(df_list, param, delta, normalise, plot_type, save_fig, folder_name, study_ranges=None):
    mean_z_cco_list, mean_z_hhb_list, mean_z_hbo2_list, rPWR_hhb_list, rCST_hhb_list, p_val, rPWR_hbo2_list, rCST_hbo2_list = [], [], [], [], [], [], [], []
    
    mean_z_cbf_list = []
    p_val_2 = []
    # Plot figure 1
    if plot_type == "all":
        fig = plt.figure(constrained_layout=True)
    al = 0
    al_list = []
    for df in df_list:
        al += 1
        al_list.append(al)
        if study_ranges:
            df_max = df[param].max() if df[param].max()!=normal_ranges[param] else df[param].min()
            if df_max < study_ranges[param][0] and df_max > study_ranges[param][1]:
                continue 
        non_delta_df = df[(125 <= df.t) & (df.t <= 135)]
        non_delta_time = df.t
        

        res_dict = get_rPWR(df, normalise)

        # use the first 10 seconds after a 5 second buffer
        df = df[(125 <= df.t) & (df.t <= 135)]
        
        time = df.t

        mean_z_cco_list.append(res_dict["mean_z_cco"])
        mean_z_hhb_list.append(res_dict["mean_z_hhb"])
        mean_z_hbo2_list.append(res_dict["mean_z_hbo2"])
        mean_z_cbf_list.append(res_dict["mean_z_cbf"])
        rCST_hhb_list.append(-res_dict["rCST_hhb"])
        rPWR_hhb_list.append(-res_dict["rPWR_hhb"])
        rCST_hbo2_list.append(res_dict["rCST_hbo2"])
        rPWR_hbo2_list.append(res_dict["rPWR_hbo2"])  
    delta_sign = "Δ"

    if plot_type == "all":
        plt.plot(al_list, mean_z_hhb_list)
        plt.ylabel(f"{delta_sign}{param}")
        plt.xlabel('time')


plot_type = "all"
save_fig = False
folder_name =  "try30secs_x1"
for f, p in {"only_rautc": "u"}.items():
    directory= f
    time=0
    df_list = []
    folder_path = f"/Users/mirunaserian/Desktop/PHD_Rotation2/haemodynamics-metabolism/bsx/{folder_name}/"
    for root, dirs, files in os.walk(f"{folder_path}/{directory}"):
        for subdir in dirs:
            sub_dir = os.path.join(root, subdir)
            print("SUBDIR", sub_dir)
            for file in os.listdir(sub_dir):
                filename = os.fsdecode(file)
                if filename.endswith(f".csv"):
                    filename = os.path.join(sub_dir, file)
                    df2 = pd.read_csv(filename, delimiter="\t")
                    df_list.append(df2)
    for delta in ["incremental-diff"]: # "normal", "percentage-change", 
        get_PWR_plots(df_list, p, delta,True, plot_type, True, sub_dir, study_ranges=None)    

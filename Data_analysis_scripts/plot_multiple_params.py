# %%
from scipy.signal import find_peaks
import numpy as np
from scipy import stats
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns

normal_ranges = {
    "u": 1,
    "P_a": 100,
    "Pa_CO2": 40
    }

def rot(angle):
        return np.array([[np.cos(angle), np.sin(angle)],
                    [-np.sin(angle), np.cos(angle)]])

def get_rotated_axes(rotated_data_hhb, rotated_data_hbo2, ax1, ax2, delta_sign):
    # get min of both columns
    min_V = min(np.amin(np.concatenate([rotated_data_hhb, rotated_data_hbo2], axis=1)),
                np.amin(np.concatenate([rotated_data_hhb, rotated_data_hbo2], axis=0)))
    max_V = max(np.amax(np.concatenate([rotated_data_hhb, rotated_data_hbo2], axis=1)), 
                np.amax(np.concatenate([rotated_data_hhb, rotated_data_hbo2], axis=0)))
    if min_V < 0 :
        min_V = -min_V
    max_val = max([min_V, max_V])

    axes = np.array([[-max_val, max_val], [-max_val, max_val]])
    rotated_axes = np.matmul(rot((np.pi)/2), axes)
    for a in [ax1, ax2]:
        a.axhline(0, color='black')
        a.axvline(0, color='black')
        a.annotate(s='', xy=(-max_val, -max_val), xytext=(max_val, max_val),  arrowprops=dict(arrowstyle='<-',linestyle='--', color="b"))
        a.text(max_val-2, max_val-0.5, "rPWR")
        a.annotate(s='', xy=(-max_val, max_val), xytext=(max_val, -max_val), arrowprops=dict(arrowstyle='->',linestyle='--', color="r"))
        a.text(-max_val+0.5, max_val-0.5, "rCST")
        a.plot([-max_val, max_val], [-max_val, max_val], color = 'black', linestyle="--", linewidth=0.2)
        a.set_ylabel(f'z({delta_sign}oxCCO)')
        a.axis('scaled')

def get_rPWR(df, normalise, delta=None):
    # start at t=1 so perecentage increases can be calculated
    time = df.t
    if delta == "change-from-first":
        df = df - df.iloc[0]
        df = df[1:]
    elif delta == "incremental-diff":
        df = df.diff()
        df = df[1:]
    elif delta == "percentage-change":
        df = df.pct_change()
        df = df[1:]

    # use the first 10 seconds after a 5 second buffer

    df.t = time
#   normalise

    df_baseline = df[(df.t.isin(range(114,119)))]
    df_sub = df.sub(df_baseline.mean())
    z_scored_baseline = df_sub.div(df_baseline.std())
    z_scored_baseline.t = time

    df = df[(125 <= df.t) & (df.t <= 135)]
    cco = df.CCO
    hhb = df.HHb
    hbo2 = df.HbO2


    # z-score the values
    # z_cco = stats.zscore(cco)
    # z_hhb = stats.zscore(hhb)
    # z_hbo2 = stats.zscore(hbo2)

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
        "rPWR_hhb": rPWR_hhb[0],
        "rPWR_hbo2": rPWR_hbo2[0],
        "rCST_hhb": rCST_hhb[0],
        "rCST_hbo2": rCST_hbo2[0]}


def get_PWR_plots(df_list, param, delta, normalise, plot_type, save_fig, folder_name, study_ranges=None):
    
    mean_z_cco_list, mean_z_hhb_list, mean_z_hbo2_list, rPWR_hhb_list, rCST_hhb_list, p_val, rPWR_hbo2_list, rCST_hbo2_list = [], [], [], [], [], [], [], []
    
    # Plot figure 1
    if plot_type == "all":
        fig = plt.figure(constrained_layout=True)
        fig.tight_layout()
        fig.set_size_inches(11.69, 8.27, forward=True)
        fig.subplots_adjust(hspace=0.5, wspace=0.5)
        fig.suptitle(f'Incremental changes of {get_title(folder_name)}', fontsize=16)


        gs = GridSpec(4, 3)
        
        ax1 = fig.add_subplot(gs[:2, 0])
        ax1.set_ylabel('z({delta_sign}HbO2)')
        ax2 = fig.add_subplot(gs[2:, 0])
        ax2.set_ylabel('z({delta_sign}HHb)')
        ax3 = fig.add_subplot(gs[0, 1])
        ax4 = fig.add_subplot(gs[1, 1])
        ax5 = fig.add_subplot(gs[2, 1])
        ax6 = fig.add_subplot(gs[3, 1])
        ax7 = fig.add_subplot(gs[0, 2])
        ax8 = fig.add_subplot(gs[1, 2])
        ax9 = fig.add_subplot(gs[2, 2])
        ax10 = fig.add_subplot(gs[3, 2])
    else:
        fig = plt.figure()
        fig.tight_layout()
        fig.set_size_inches( 12, 12, forward=True)
        fig.subplots_adjust(hspace=0.5, wspace=0.5)
        fig.suptitle(f'Incremental changes of {get_title(folder_name)}', fontsize=16)

        gs3 = GridSpec(2, 2)
        ax3_1 = fig.add_subplot(gs3[0, 0])
        ax3_2 = fig.add_subplot(gs3[0, 1])
        ax3_3 = fig.add_subplot(gs3[1, 0])
        ax3_4 = fig.add_subplot(gs3[1, 1])

    if len(param) == 2:
        param1, param2 = param
    elif len(param) == 3:
        param1, param2, param3 = param
    elif len(param) == 4:
        param1, param2, param3, param4 = param
    count_sets = 0
    for df in df_list:
        if study_ranges:
            df_max = df[param1].max() if df[param1].max()!=normal_ranges[param1] else df[param1].min()
            if df_max < study_ranges[param1][0] and df_max > study_ranges[param1][1]:
                continue 
        non_delta_param1 = df[param1]
        non_delta_param2 = df[param2]
        non_delta_time = df.t
        

        res_dict = get_rPWR(df, normalise, delta)
        # df = df[(df.t.isin(range(124,136)))]
        df = df[(125 <= df.t) & (df.t <= 135)]

        count_sets += 1
        param_changes = df[param1].mean()

        delta_sign = ""
        time = df.t
        if delta == "change-from-first":
            df = df - df.iloc[0]
            df = df[1:]
            time = time[1:]
            delta_sign = "Δ"
        elif delta == "incremental-diff":
            df = df.diff()
            df = df[1:]
            time = time[1:]
            delta_sign = "Δ"
        elif delta == "percentage-change":
            df = df.pct_change()
            df = df[1:]
            time = time[1:]
            delta_sign = "Δ"
        

        mean_z_cco_list.append(res_dict["mean_z_cco"])
        mean_z_hhb_list.append(res_dict["mean_z_hhb"])
        mean_z_hbo2_list.append(res_dict["mean_z_hbo2"])
        rCST_hhb_list.append(-res_dict["rCST_hhb"])
        rPWR_hhb_list.append(-res_dict["rPWR_hhb"])
        rCST_hbo2_list.append(res_dict["rCST_hbo2"])
        rPWR_hbo2_list.append(res_dict["rPWR_hbo2"])  
        p_val.append(non_delta_param2.mean())


        if plot_type == "all":
            ax3.plot(non_delta_time, non_delta_param1)
            ax3.set_ylabel(f"{delta_sign}{param1}")
            ax3.set_xlabel('time')

            ax4.scatter(param_changes, df.CCO.mean())
            ax4.set_ylabel(f'Mean {delta_sign}CCO')
            ax4.set_xlabel(f"max {param}")

            ax5.scatter(param_changes, df.HHb.mean())
            ax5.set_ylabel(f'Mean {delta_sign}HHb')
            ax5.set_xlabel(f"max {param}")

            ax6.scatter(param_changes, df.HbO2.mean())
            ax6.set_ylabel(f'Mean {delta_sign}HbO2')
            ax6.set_xlabel(f"max {param}")

            ax7.scatter(non_delta_time, non_delta_param2)
            ax7.set_ylabel(param2)
            ax7.set_xlabel(f"time")

            ax8.plot(time, df.CCO)
            ax8.set_ylabel(f'{delta_sign}CCO')
            ax8.set_xlabel(f"time")

            ax9.plot(time, df.HHb)
            ax9.set_ylabel(f'{delta_sign}HHb')
            ax9.set_xlabel(f"time")

            ax10.plot(time, df.HbO2)
            ax10.set_ylabel(f'{delta_sign}HbO2')
            ax10.set_xlabel(f"time")

        if "mean_only" in plot_type:
            ax3_1.scatter(param_changes, df.CCO.mean())
            ax3_1.set_ylabel(f'Mean {delta_sign}CCO')
            ax3_1.set_xlabel(f"max {param}")

            ax3_2.scatter(param_changes, df.HHb.mean())
            ax3_2.set_ylabel(f'Mean {delta_sign}HHb')
            ax3_2.set_xlabel(f"max {param}")

            ax3_3.scatter(param_changes, df.HbO2.mean())
            ax3_3.set_ylabel(f'Mean {delta_sign}HbO2')
            ax3_3.set_xlabel(f"max {param}")
            
        if "normal_only" in plot_type:
            ax3_1.scatter(time, df[param])
            ax3_1.set_ylabel(f"{delta_sign}{param}")
            ax3_1.set_xlabel(f"time")

            ax3_2.plot(time, df.CCO)
            ax3_2.set_ylabel(f'{delta_sign}CCO')
            ax3_2.set_xlabel(f"time")

            ax3_3.plot(time, df.HHb)
            ax3_3.set_ylabel(f'{delta_sign}HHb')
            ax3_3.set_xlabel(f"time")

            ax3_4.plot(time, df.HbO2)
            ax3_4.set_ylabel(f'{delta_sign}HbO2')
            ax3_4.set_xlabel(f"time")

            
    if plot_type == "all":
        # invert hhb z scores
        data_hhb = np.array([stats.zscore(mean_z_hhb_list), stats.zscore(mean_z_cco_list)])
        data_hbo2 = np.array([stats.zscore(mean_z_hbo2_list), stats.zscore(mean_z_cco_list)])


        rotated_data_hhb = np.matmul(rot((np.pi)/4), data_hhb)
        rotated_data_hbo2 = np.matmul(rot((np.pi)/4), data_hbo2)

        rPWR_hhb = rotated_data_hhb[0,:]
        rCST_hhb = rotated_data_hhb[1,:]

        rPWR_hbo2 = rotated_data_hbo2[0,:]
        rCST_hbo2 = rotated_data_hbo2[1,:]

        ax1.scatter(stats.zscore([-x for x in mean_z_hhb_list]), stats.zscore(mean_z_cco_list), s=50, c=p_val, cmap='Greens')
        ax2.set_title("HbO2")
        ax2.set_xlabel(f'z(ΔHbO2)')
        ax2.scatter(stats.zscore(mean_z_hbo2_list), stats.zscore(mean_z_cco_list),s=50, c=p_val, cmap='Greens')
        ax1.set_title("Hhb")
        ax1.set_xlabel(f'z(ΔHHb)')

        # for i,  txt in enumerate(p_val):
        #     ax1.annotate("{:.1f}".format(txt), (stats.zscore(mean_z_hhb_list)[i], stats.zscore(mean_z_cco_list)[i]))
        #     ax2.annotate("{:.1f}".format(txt), (stats.zscore(mean_z_hbo2_list)[i], stats.zscore(mean_z_cco_list)[i]))
        # # plot figure 2
        # fig2 = plt.figure(constrained_layout=True)
        # fig2.tight_layout()
        # fig2.set_size_inches( 11.69, 8.27, forward=True)
        # fig2.subplots_adjust(hspace=0.4, wspace=0.2)

        # gs2 = GridSpec(2, 2)
        # ax12 = fig2.add_subplot(gs2[0, 0])
        # ax22 = fig2.add_subplot(gs2[1, 0])
        # ax32 = fig2.add_subplot(gs2[0, 1])
        # ax42 = fig2.add_subplot(gs2[1, 1])

        # ax12.scatter(rotated_data_hbo2[0,:], rotated_data_hbo2[1,:], s=50, c=p_val, cmap='Greens')
        # ax12.set_title("HbO2")
        # ax22.scatter(rotated_data_hhb[0,:], rotated_data_hhb[1,:],s=50, c=p_val, cmap='Greens')
        # ax22.set_title("Hhb")
        # ax12.scatter(stats.zscore(mean_z_hbo2_list), \
        #     stats.zscore(mean_z_cco_list), s=50, c=p_val, cmap='Blues')
        # ax22.scatter(stats.zscore(mean_z_hhb_list), \
        #     stats.zscore(mean_z_cco_list), s=50, c=p_val, cmap='Blues')
        # ax32.scatter(rPWR_hbo2_list, rCST_hbo2_list, \
        #     s=50, c=p_val, cmap='Greens')
        # ax42.scatter(rPWR_hhb_list, stats.zscore(rCST_hhb_list), \
        #     s=50, c=p_val, cmap='Greens')
    

    get_rotated_axes(rotated_data_hhb, rotated_data_hbo2, ax1, ax2, delta_sign)
    # Add rPWR and r CST axes to plots

        # ax1.text(1, 1, '\n'.join(('▲{delta_sign}HbO2', "▼{delta_sign}oxCCO")), transform=a.transAxes, fontsize=14,
        # verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        # ax1.text(0.05, 0.95, '\n'.join(('▲ΔHbO2', "▲ΔoxCCO")), transform=a.transAxes, fontsize=14,
        # verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    if save_fig == True:
        if plot_type ==  "all":
            fig.savefig(f'{folder_name}/normal_varying_{"_".join([k for k in param])}_{delta}.jpg', dpi=fig.dpi)
        else:
            fig.savefig(f'{folder_name}/{plot_type}_{"_".join([k for k in param])}_{delta}.jpg', dpi=fig.dpi)

        
def plot_only_rotation(df_list, param, normalise, plot_type, save_fig, folder_name, study_ranges=None):
    mean_z_cco_list, mean_z_hhb_list, mean_z_hbo2_list, rPWR_hhb_list, rCST_hhb_list, p_val, rPWR_hbo2_list, rCST_hbo2_list = [], [], [], [], [], [], [], []
    
    symbol_dict = {"increase": '▲',
                    "decrease": '▼'}
    fig = plt.figure(constrained_layout=True)
    fig.tight_layout()
    fig.set_size_inches(11.69, 8.27, forward=True)
    fig.subplots_adjust(hspace=0.5, wspace=0.5)
    title = folder_name.split("/")[-1]
    title_with_symbols = title.replace("increase", symbol_dict["increase"])
    title_with_symbols = title_with_symbols.replace("decrease", symbol_dict["decrease"])

    fig.suptitle(f'Incremental changes of {title_with_symbols.replace("_", " ")}', fontsize=16)


    gs = GridSpec(2, 1)
    
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_ylabel('z({delta_sign}HbO2)')
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.set_ylabel('z({delta_sign}HHb)')
    
    if len(param) == 2:
        param1, param2 = param
    elif len(param) == 3:
        param1, param3, param2 = param
    else:
        param1, param3, param2, param4 = param
    count_sets = 0
    for df in df_list:
        non_delta_param1 = df[param1]
        non_delta_param2 = df[param2]
        non_delta_time = df.t
        

        res_dict = get_rPWR(df, normalise, delta)
        # df = df[(df.t.isin(range(124,136)))]
        df = df[(125 <= df.t) & (df.t <= 135)]

        count_sets += 1
        param_changes = df[param1].mean()

        delta_sign = ""
        time = df.t
        df = df - df.iloc[0]
        df = df[1:]
        df.t = time[1:]
        time = time[1:]
        delta_sign = "Δ"

        

        mean_z_cco_list.append(res_dict["mean_z_cco"])
        mean_z_hhb_list.append(res_dict["mean_z_hhb"])
        mean_z_hbo2_list.append(res_dict["mean_z_hbo2"])
        rCST_hhb_list.append(-res_dict["rCST_hhb"])
        rPWR_hhb_list.append(-res_dict["rPWR_hhb"])
        rCST_hbo2_list.append(res_dict["rCST_hbo2"])
        rPWR_hbo2_list.append(res_dict["rPWR_hbo2"])  
        p_val.append(non_delta_param1.mean())

    data_hhb = np.array([stats.zscore(mean_z_hhb_list), stats.zscore(mean_z_cco_list)])
    data_hbo2 = np.array([stats.zscore(mean_z_hbo2_list), stats.zscore(mean_z_cco_list)])
    rotated_data_hhb = np.matmul(rot((np.pi)/4), data_hhb)
    rotated_data_hbo2 = np.matmul(rot((np.pi)/4), data_hbo2)

    ax1.scatter(stats.zscore([-x for x in mean_z_hhb_list]), stats.zscore(mean_z_cco_list), s=50, c=p_val, cmap='Greens')
    ax2.set_title("HbO2")
    ax2.set_xlabel(f'z(ΔHbO2)')
    ax2.scatter(stats.zscore(mean_z_hbo2_list), stats.zscore(mean_z_cco_list),s=50, c=p_val, cmap='Greens')
    ax1.set_title("Hhb")
    ax1.set_xlabel(f'z(ΔHHb)')

    # Add rPWR and r CST axes to plots
    get_rotated_axes(rotated_data_hhb, rotated_data_hbo2, ax1, ax2, delta_sign)

        # ax1.text(1, 1, '\n'.join(('▲{delta_sign}HbO2', "▼{delta_sign}oxCCO")), transform=a.transAxes, fontsize=14,
        # verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        # ax1.text(0.05, 0.95, '\n'.join(('▲ΔHbO2', "▲ΔoxCCO")), transform=a.transAxes, fontsize=14,
        # verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    print(f"{folder_name} - {get_title(folder_name)}")

    if save_fig == True:
        if plot_type ==  "all":
            fig.savefig(f'{folder_name}/only_rot_{plot_type}_{"_".join([k for k in param])}_{delta}.jpg', dpi=fig.dpi)


def plot_3d(df_list, param, delta, normalise, plot_type, save_fig, folder_name, study_ranges=None):
    mean_z_cco_list, mean_z_hhb_list, mean_z_hbo2_list, rPWR_hhb_list, rCST_hhb_list, p_val, rPWR_hbo2_list, rCST_hbo2_list = [], [], [], [], [], [], [], []
    
    if len(param) == 2:
        param1, param2 = param
    elif len(param) == 3:
        param1, param3, param2 = param
    else:
        param1, param3, param2, param4 = param
    count_sets = 0
    for df in df_list:
        non_delta_param1 = df[param1]
        non_delta_param2 = df[param2]
        non_delta_time = df.t
        

        res_dict = get_rPWR(df, normalise, delta)
        df = df[(125 <= df.t) & (df.t <= 135)]

        count_sets += 1
        param_changes = df[param1].mean()

        delta_sign = ""
        time = df.t
        if delta == "change-from-first":
            df = df - df.iloc[0]
            df = df[1:]
            time = time[1:]
            delta_sign = "Δ"
        elif delta == "incremental-diff":
            df = df.diff()
            df = df[1:]
            time = time[1:]
            delta_sign = "Δ"
        elif delta == "percentage-change":
            df = df.pct_change()
            df = df[1:]
            time = time[1:]
            delta_sign = "Δ"
        

        mean_z_cco_list.append(res_dict["mean_z_cco"])
        mean_z_hhb_list.append(res_dict["mean_z_hhb"])
        mean_z_hbo2_list.append(res_dict["mean_z_hbo2"])
        rCST_hhb_list.append(-res_dict["rCST_hhb"])
        rPWR_hhb_list.append(-res_dict["rPWR_hhb"])
        rCST_hbo2_list.append(res_dict["rCST_hbo2"])
        rPWR_hbo2_list.append(res_dict["rPWR_hbo2"])  
        p_val.append(non_delta_param1.mean())

    data_hhb = np.array([stats.zscore(mean_z_hhb_list), stats.zscore(mean_z_cco_list)])
    data_hbo2 = np.array([stats.zscore(mean_z_hbo2_list), stats.zscore(mean_z_cco_list)])
  
    
    fig = plt.figure(figsize=(8, 6))
    fig.suptitle(f'Incremental changes of {folder_name.split("/")[-1]}', fontsize=16)

    ax = fig.add_subplot(111, projection='3d')

    xs = p_val
    ys =  stats.zscore(mean_z_hbo2_list)
    zs = stats.zscore(mean_z_cco_list)
    ax.scatter(xs, ys, zs, s=50, alpha=0.6, edgecolors='w')

    ax.set_ylabel('hbo2')
    ax.set_xlabel('u')
    ax.set_zlabel('cco')

    max_val = get_max(data_hhb, data_hbo2)
    


def plot_pca(df_list, param, delta, normalise, plot_type, save_fig, folder_name, study_ranges=None):
    mean_z_cco_list, mean_z_hhb_list, mean_z_hbo2_list, rPWR_hhb_list, rCST_hhb_list, p_val, rPWR_hbo2_list, rCST_hbo2_list = [], [], [], [], [], [], [], []
    
    if len(param) == 2:
        param1, param2 = param
    elif len(param) == 3:
        param1, param3, param2 = param
    else:
        param1, param3, param2, param4 = param
    count_sets = 0
    for df in df_list:
        non_delta_df = df
        non_delta_param1 = df[param1]
        non_delta_param2 = df[param2]
        non_delta_time = df.t
        

        res_dict = get_rPWR(df, normalise, delta)
        df = df[(125 <= df.t) & (df.t <= 135)]

        count_sets += 1
        param_changes = df[param1].mean()

        delta_sign = ""
        time = df.t
        

        mean_z_cco_list.append(res_dict["mean_z_cco"])
        mean_z_hhb_list.append(res_dict["mean_z_hhb"])
        mean_z_hbo2_list.append(res_dict["mean_z_hbo2"])
        rCST_hhb_list.append(-res_dict["rCST_hhb"])
        rPWR_hhb_list.append(-res_dict["rPWR_hhb"])
        rCST_hbo2_list.append(res_dict["rCST_hbo2"])
        rPWR_hbo2_list.append(res_dict["rPWR_hbo2"])  
        p_val.append(non_delta_param1.mean())


    df_new = pd.DataFrame({"mean_z_cco": mean_z_cco_list,
                       "mean_z_hhb": mean_z_hhb_list,
                       "mean_z_hbo2": mean_z_hbo2_list,
                       "rCST_hhb": rCST_hhb_list,
                       "rCST_hbo2": rCST_hbo2_list,
                       "rPWR_hhb": rPWR_hhb_list,
                       "rPWR_hbo2": rPWR_hbo2_list,
                       **{p:non_delta_df[p].mean() for p in param}})
    import plotly.express as px
    from sklearn.preprocessing import StandardScaler    
    from sklearn.decomposition import PCA
    import seaborn as sns

    features = [c for c in list(df_new.columns) if c != "u"]

    pca =PCA(n_components=3)
    pca.fit((df_new))
    n2 = pca.transform(df_new)
    components = pd.DataFrame(data = n2, columns = ['Principal component 1', 'Principal component 2','Principal component 3'])
    x = components["Principal component 1"]
    y = components["Principal component 2"]
    # fig, ax = plt.subplots(2, 2, figsize=(12,8))
    # sns.kdeplot(components["Principal component 1"], components["Principal component 2"], cmap='Blues',
    #         shade=True, shade_lowest=False, palette="crest",
    # alpha=.6)

    fig2, ax = plt.subplots(1)
    plt.matshow(pca.components_,cmap='viridis')
    plt.yticks([0,1,2],['1st Comp','2nd Comp','3rd Comp'],fontsize=10)
    plt.colorbar()
    plt.xticks(range(len(features)),features,rotation=65,ha='left')
    plt.tight_layout()
    fig2.savefig(f'{folder_name}/score_pca_{folder_name.split("/")[-1]}')
    # fig1.save_fig(f'{folder_name}/pca_{folder_name.split("/")[-1]}')


def get_data(df_list, param_list, normalise):
    mean_z_cco_list, mean_z_hhb_list, mean_z_hbo2_list, rPWR_hhb_list, rCST_hhb_list, p_val, rPWR_hbo2_list, rCST_hbo2_list = [], [], [], [], [], [], [], []

    param_arr = np.empty((0,len(param_list)), float)
    mean_z_arr = np.empty((0,3), float)
    rotated_hhb_arr  = np.empty((0,2), float)
    rotated_hbo2_arr  = np.empty((0,2), float)

    for df in df_list:
        non_delta_df = df
        non_delta_time = df.t
        

        res_dict = get_rPWR(df, normalise, delta)
        df = df[(125 <= df.t) & (df.t <= 135)]

        time = df.t

        mean_z_cco_list.append(res_dict["mean_z_cco"])
        mean_z_hhb_list.append(res_dict["mean_z_hhb"])
        mean_z_hbo2_list.append(res_dict["mean_z_hbo2"])
        rCST_hhb_list.append(-res_dict["rCST_hhb"])
        rPWR_hhb_list.append(-res_dict["rPWR_hhb"])
        rCST_hbo2_list.append(res_dict["rCST_hbo2"])
        rPWR_hbo2_list.append(res_dict["rPWR_hbo2"])  

        param_arr = np.append(param_arr, np.array([[non_delta_df[p].mean() for p in param_list]]), axis=0)

    df_new = pd.DataFrame({"mean_z_cco": mean_z_cco_list,
                        "mean_z_hhb": mean_z_hhb_list,
                        "mean_z_hbo2": mean_z_hbo2_list,
                        "rCST_hhb": rCST_hhb_list,
                        "rCST_hbo2": rCST_hbo2_list,
                        "rPWR_hhb": rPWR_hhb_list,
                        "rPWR_hbo2": rPWR_hbo2_list,
                        **{param_list[x]: param_arr[:, x] for x in range(len(param_list))}})
    return df_new

def get_title(folder_name):
    symbol_dict = {"increase": '▲',
                    "decrease": '▼'}
    title = folder_name.split("/")[-1]
    title_with_symbols = title.replace("increase", symbol_dict["increase"])
    return \
        title_with_symbols.replace("decrease", symbol_dict["decrease"])

def get_corr_matrix(df_list, param_list, normalise, folder_name):
    df = get_data(df_list, param_list, normalise)
    df = df[[n for n in list(df.columns) if n not in ["mean_z_cco", "mean_z_hhb", "mean_z_hbo2"]]]
    corr = df.corr(method='spearman')

    # Generate a mask for the upper triangle
    mask = np.zeros_like(corr, dtype=np.bool)
    mask[np.triu_indices_from(mask)] = True

    # Set up the matplotlib figure
    fig, ax = plt.subplots(figsize=(6, 5))

    # Generate a custom diverging colormap
    cmap = sns.diverging_palette(220, 10, as_cmap=True, sep=100)

    # Draw the heatmap with the mask and correct aspect ratio
    sns.heatmap(corr, mask=mask, cmap=cmap, vmin=-1, vmax=1, center=0, linewidths=.5).set_title(f'Correlation matrix of features -  {get_title(folder_name)}', fontsize=15)

    fig.tight_layout()
    fig.savefig(f'{folder_name}/corr_{get_title(folder_name)}.jpg', dpi=fig.dpi)

def heatmap_rcst(df_list, param_lisr, normalise, folder_name):
    df = get_data(df_list, param_lisr, normalise)
    df[["mean_z_hhb", "mean_z_hbo2", "mean_z_cco"]].apply(stats.zscore)
    df["rcst_pos_hhb"] =  df.apply(lambda row : 1 if row["mean_z_hhb"]<0 and row["mean_z_cco"]>0 \
        else -1 if row["mean_z_hhb"]>0 and row["mean_z_cco"]<0 else 0, axis=1)
    df["rcst_pos_hbo2"] =  df.apply(lambda row : 1 if row["mean_z_hbo2"]<0 and row["mean_z_cco"]>0 \
        else -1 if row["mean_z_hbo2"]>0 and row["mean_z_cco"]<0 else 0, axis=1)
    # ax = sns.heatmap(df[["u", "rcst_pos_hbo2"]]).set_title(get_title(folder_name))

    df_wrong = df.loc[(df["rcst_pos_hhb"] != 0) & (df["rcst_pos_hbo2"] != 0)]
    print(get_title(folder_name))
    print(len(df_wrong))
    # print(df_wrong[[n for n in list(df.columns) \
    #     if "mean" not in n and "rCST" not in n and "rPWR" not in n]])

study_ranges = {
    "u": (0.9, 1.1),
    "P_a": (80, 120),
    "Pa_CO2": (38, 42),
    "SaO2": [90, 100]
} 
plot_type = "all"
save_fig = False
folder_name = "new_try30secs_x2"
# for f, p in {"u_and_Pa": ["u", "P_a"],
#             "u_and_Pa_CO2": ["u", "Pa_CO2"],
#             "u_and_SaO2sup": ["u", "SaO2sup"],
#             "u_and_Pa_and_PaCO2": ["u", "P_a", "Pa_CO2"],
#             "u_and_Pa_and_SaO2sup": ["u", "P_a", "SaO2sup"],
#             "u_and_SaO2sup_and_PaCO2": ["u", "SaO2sup", "Pa_CO2"],
for f, p in {"u_and_Pa_and_PaCO2_and_SaO2sup": ["u", "P_a", "SaO2sup", "Pa_CO2"]}.items():
    directory= f
    time=0
    
    folder_path = f"/Users/mirunaserian/Desktop/PHD_Rotation2/haemodynamics-metabolism/bsx/{folder_name}/"
    if not os.path.exists(folder_path):
        print(f"{folder_path} DOESN'T EXIST")
        continue
    # for file1 in os.listdir(f"{folder_path}/{directory}"):
    for root, dirs, files in os.walk(f"{folder_path}/{directory}"):
        for subdir in dirs:
            df_list = []
            sub_dir = os.path.join(root, subdir)
            for file in os.listdir(sub_dir):
                filename = os.fsdecode(file)
                if filename.endswith(f".csv"):
                    filename = os.path.join(sub_dir, file)
                    df2 = pd.read_csv(filename, delimiter="\t")
                    df_list.append(df2)
            # get_PWR_plots(df_list, p, "change-from-first" , True, "all", True, sub_dir, study_ranges=None)

            # get_PWR_plots(df_list, p, True, plot_type,True, sub_dir, study_ranges=None)      
            get_corr_matrix(df_list, p, True, sub_dir)
# %%

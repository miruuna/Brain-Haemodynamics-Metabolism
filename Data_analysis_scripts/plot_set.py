# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%

from scipy.signal import find_peaks
import numpy as np

from scipy import stats
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from bsx_data import study_ranges, normal_ranges, folder_param_dict
from utils_sims_bsx import get_data_set, get_max, get_rPWR, rot
from mpl_toolkits.mplot3d import Axes3D
plt.style.use('tableau-colorblind10')
print("styles", plt.style.available)
folder_base = ""

def get_PWR_plots(df_list, param, delta, normalise, plot_type, save_fig, folder_name, study_ranges=None):
    mean_z_cco_list, mean_z_hhb_list, mean_z_hbo2_list, rPWR_hhb_list, rCST_hhb_list, p_val, rPWR_hbo2_list, rCST_hbo2_list = [], [], [], [], [], [], [], []
    
    p_val_2 = []
    # Plot figure 1
    if plot_type == "all":
        fig = plt.figure(constrained_layout=True)
        fig.tight_layout()
        fig.set_size_inches(11.69, 8.27, forward=True)
        fig.subplots_adjust(hspace=0.5, wspace=0.5)
        fig.suptitle(f'Incremental changes of {param}', fontsize=16)


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
        fig = plt.figure(constrained_layout=True)
        fig.tight_layout()
        fig.set_size_inches( 12, 12, forward=True)
        fig.subplots_adjust(hspace=0.5, wspace=0.5)
        fig.suptitle(f'Incremental changes of {param}', fontsize=16)

        gs3 = GridSpec(2, 2)
        ax3_1 = fig.add_subplot(gs3[0, 0])
        ax3_2 = fig.add_subplot(gs3[0, 1])
        ax3_3 = fig.add_subplot(gs3[1, 0])
        ax3_4 = fig.add_subplot(gs3[1, 1])

    for df in df_list:
        if study_ranges:
            df_max = df[param].max() if df[param].max()!=normal_ranges[param] else df[param].min()
            if df_max < study_ranges[param][0] and df_max > study_ranges[param][1]:
                continue 
        non_delta_df = df
        non_delta_param = df[param]
        non_delta_time = df.t
        non_delta_u = df["u"]
        

        res_dict = get_rPWR(df, normalise)

        # use the first 10 seconds after a 5 second buffer
        df = df[(125 <= df.t) & (df.t <= 135)]
        param_changes = df[param].mean()
        delta_sign = "Δ"
        time = df.t


        mean_z_cco_list.append(res_dict["mean_z_cco"])
        mean_z_hhb_list.append(res_dict["mean_z_hhb"])
        mean_z_hbo2_list.append(res_dict["mean_z_hbo2"])
        rCST_hhb_list.append(-res_dict["rCST_hhb"])
        rPWR_hhb_list.append(-res_dict["rPWR_hhb"])
        rCST_hbo2_list.append(res_dict["rCST_hbo2"])
        rPWR_hbo2_list.append(res_dict["rPWR_hbo2"])  
        p_val.append(param_changes)


        if plot_type == "all":
            ax3.plot(non_delta_time, non_delta_param)
            ax3.set_ylabel(f"{delta_sign}{param}")
            ax3.set_xlabel('time')

            ax4.scatter(param_changes, res_dict["mean_z_cco"])
            ax4.set_ylabel(f'Mean {delta_sign}CCO')
            ax4.set_xlabel(f"max {param}")

            ax5.scatter(param_changes, df.HHb.mean())
            ax5.set_ylabel(f'Mean {delta_sign}HHb')
            ax5.set_xlabel(f"max {param}")

            ax6.scatter(param_changes, df.HbO2.mean())
            ax6.set_ylabel(f'Mean {delta_sign}HbO2')
            ax6.set_xlabel(f"max {param}")


            ax7.plot(non_delta_time, non_delta_u)
            ax7.set_ylabel("u")
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

        ax1.scatter(stats.zscore([-x for x in mean_z_hhb_list]), stats.zscore(mean_z_cco_list), s=50, c=p_val, cmap='Greens')
        ax1.set_title("Hhb")
        ax1.set_xlabel(f'z(ΔHHb)')

        ax2.scatter(stats.zscore(mean_z_hbo2_list), stats.zscore(mean_z_cco_list),s=50, c=p_val, cmap='Greens')
        ax2.set_title("HbO2")
        ax2.set_xlabel(f'z(ΔHbO2)')
        for i,  txt in enumerate(p_val):
            if i==0 or i==len(p_val)-1:
                ax1.annotate("{:.1f}".format(txt), (stats.zscore(mean_z_hhb_list)[i], stats.zscore(mean_z_cco_list)[i]))
                ax2.annotate("{:.1f}".format(txt), (stats.zscore(mean_z_hbo2_list)[i], stats.zscore(mean_z_cco_list)[i]))



    # Add rPWR and r CST axes to plots
        max_val = get_max(data_hhb, data_hbo2)
        
        axes = np.array([[-max_val, max_val], [-max_val, max_val]])
        rotated_axes = np.matmul(rot((np.pi)/2), axes)
        for a in [ax1, ax2]:
            a.axhline(0, color='black')
            a.axvline(0, color='black')
            # a.annotate(s='', xy=(-max_val, -max_val), xytext=(max_val, max_val),  arrowprops=dict(arrowstyle='<-',linestyle='--', color="b"))
            a.text(0.8, 2-0.5, "rPWR")
            # a.annotate(s='', xy=(-max_val, max_val), xytext=(max_val, -max_val), arrowprops=dict(arrowstyle='->',linestyle='--', color="r"))
            a.text(-2+0.5, 2-0.5, "rCST")
            a.plot([-max_val, max_val], [-max_val, max_val], color = 'black', linestyle="--", linewidth=0.2)
            a.set_ylabel(f'z({delta_sign}oxCCO)')
            a.axis('scaled')

    if save_fig == True:
        if plot_type ==  "all":
            fig.savefig(f'{folder_base}/{folder_name}/normal_varying_{param[0]}.jpg', dpi=fig.dpi)
        else:
            fig.savefig(f'{folder_name}/{plot_type}_{param}_{delta}.jpg', dpi=fig.dpi)

        

def get_PWR_plots_paper(df_list, param, delta, normalise, plot_type, save_fig, folder_name, study_ranges=None):
    mean_z_cco_list, mean_z_hhb_list, mean_z_hbo2_list, rPWR_hhb_list, rCST_hhb_list, p_val, rPWR_hbo2_list, rCST_hbo2_list = [], [], [], [], [], [], [], []
    
    p_val_2 = []
    # Plot figure 1
    fig = plt.figure(constrained_layout=True)
    fig.set_size_inches(11.69, 8.27, forward=True)
    fig.subplots_adjust(hspace=0.5, wspace=0.5)
    if param[0] == "u":
        t = fig.suptitle(f'Variations in demand (u)', fontsize=16)
    else:
        t = fig.suptitle(f'Changes of {param}', fontsize=16)

    t.set_y(1.02)

    gs = GridSpec(4, 2)
    
    ax1 = fig.add_subplot(gs[:2, 0])
    ax1.set_ylabel('z({delta_sign}HbO2)')
    ax2 = fig.add_subplot(gs[2:, 0])
    ax2.set_ylabel('z({delta_sign}HHb)')
    ax3 = fig.add_subplot(gs[0, 1])

    ax4 = fig.add_subplot(gs[1, 1])
    ax5 = fig.add_subplot(gs[2, 1])
    ax6 = fig.add_subplot(gs[3, 1])
 

    for df in df_list:
        if study_ranges:
            df_max = df[param].max() if df[param].max()!=normal_ranges[param] else df[param].min()
            if df_max < study_ranges[param][0] and df_max > study_ranges[param][1]:
                continue 
        non_delta_df = df
        non_delta_param = df[param]
        non_delta_time = df.t
        non_delta_u = df["u"]
        

        res_dict = get_rPWR(df, normalise)

        # use the first 10 seconds after a 5 second buffer
        df = df[(125 <= df.t) & (df.t <= 135)]

        param_changes = df[param].mean()

        delta_sign = "Δ"
        time = df.t


        mean_z_cco_list.append(res_dict["mean_z_cco"])
        mean_z_hhb_list.append(res_dict["mean_z_hhb"])
        mean_z_hbo2_list.append(res_dict["mean_z_hbo2"])
        rCST_hhb_list.append(-res_dict["rCST_hhb"])
        rPWR_hhb_list.append(-res_dict["rPWR_hhb"])
        rCST_hbo2_list.append(res_dict["rCST_hbo2"])
        rPWR_hbo2_list.append(res_dict["rPWR_hbo2"])  
        p_val.append(param_changes)

        cmap = "Pastel1"
        ax3.plot(non_delta_time, non_delta_param)
        ax3.set_ylabel(f"{delta_sign}{param[0]}")
        ax3.set_xlabel('Time')

        ax4.scatter(param_changes, df.CCO.mean(),cmap=cmap)
        ax4.set_ylabel(f'Mean {delta_sign}CCO')
        ax4.set_xlabel(f"{param[0]}")

        ax5.scatter(param_changes, df.HHb.mean(), cmap=cmap)
        ax5.set_ylabel(f'Mean {delta_sign}HHb')
        ax5.set_xlabel(f"{param[0]}")

        ax6.scatter(param_changes, df.HbO2.mean(), cmap=cmap)
        ax6.set_ylabel(f'Mean {delta_sign}HbO2')
        ax6.set_xlabel(f"{param[0]}")

        # invert hhb z scores
    data_hhb = np.array([stats.zscore(mean_z_hhb_list), stats.zscore(mean_z_cco_list)])
    data_hbo2 = np.array([stats.zscore(mean_z_hbo2_list), stats.zscore(mean_z_cco_list)])

    scatter1 = ax1.scatter(stats.zscore([-x for x in mean_z_hhb_list]), stats.zscore(mean_z_cco_list), s=50, c=p_val, cmap='Greens')  
    ax1.set_title("HHb")
    ax1.set_xlabel(f'z(ΔHHb)')
  


    ax2.scatter(stats.zscore(mean_z_hbo2_list), stats.zscore(mean_z_cco_list),s=50, c=p_val, cmap='Greens')
    ax2.set_title("HbO2")
    ax2.set_xlabel(f'z(ΔHbO2)')

    max_val = get_max(data_hhb, data_hbo2)

    cbar = plt.colorbar(scatter1, ax=[ax1, ax2], fraction=0.09, pad=0.03)
    cbar.ax.set_ylabel(f'{param[0]}', rotation=90)
    cbar.ax.get_yaxis().labelpad = 0.03

    for a in [ax1, ax2]:
        a.axhline(0, color='black')
        a.axvline(0, color='black')
        a.annotate(s='', xy=(-max_val, -max_val), xytext=(max_val, max_val),  arrowprops=dict(arrowstyle='<-',linestyle='--', color="b"))
        a.text(0.8, 2-0.5, "rPWR")
        a.annotate(s='', xy=(-max_val, max_val), xytext=(max_val, -max_val), arrowprops=dict(arrowstyle='->',linestyle='--', color="r"))
        a.text(-2+0.5, 2-0.5, "rCST")
        a.plot([-max_val, max_val], [-max_val, max_val], color = 'black', linestyle="--", linewidth=0.2)
        a.set_ylabel(f'z({delta_sign}oxCCO)',fontsize=12)
        a.axis('scaled')
    for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
        ax.xaxis.set_tick_params(labelsize=12)
        ax.yaxis.set_tick_params(labelsize=12)
        ax.xaxis.label.set_size(12)
        ax.yaxis.label.set_size(12)


    # ax1.text(-0.05, 1.2, '\n'.join(('▼HbO2', "▲oxCCO")), transform=ax1.transAxes, fontsize=12,
    # verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
    # ax1.text(1, 1.2, '\n'.join(('▲HbO2', "▲oxCCO")), transform=ax1.transAxes, fontsize=12,
    # verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
    fig.tight_layout()


    if save_fig == True:
        if plot_type ==  "all":
            fig.savefig(f'{folder_name}/normal_varying_{param}_{delta}.jpg', dpi=fig.dpi)
        else:
            fig.savefig(f'{folder_name}/{plot_type}_{param}_{delta}.jpg', dpi=fig.dpi)

        
def get_PWR_plots_paper_3drot(df_list, param, delta, normalise, plot_type, save_fig, folder_name, study_ranges=None):
    mean_z_cco_list, mean_z_hhb_list, mean_z_hbo2_list, rPWR_hhb_list, rCST_hhb_list, p_val, rPWR_hbo2_list, rCST_hbo2_list = [], [], [], [], [], [], [], []
    
    p_val_2 = []
    # Plot figure 1
    fig = plt.figure(constrained_layout=True)
    fig.tight_layout()
    fig.set_size_inches(11.69, 8.27, forward=True)
    fig.subplots_adjust(hspace=0.5, wspace=0.5)
    fig.suptitle(f'Incremental changes of {param}', fontsize=16)


    gs = GridSpec(2, 2)
    
    ax1 = fig.add_subplot(gs[0, 0], projection='3d')
    ax1.set_ylabel('z({delta_sign}HbO2)')
    ax2 = fig.add_subplot(gs[1, 0], projection='3d')
    ax2.set_ylabel('z({delta_sign}HHb)')
    ax3 = fig.add_subplot(gs[0, 1])

    ax4 = fig.add_subplot(gs[1, 1])
 

    for df in df_list:
        if study_ranges:
            df_max = df[param].max() if df[param].max()!=normal_ranges[param] else df[param].min()
            if df_max < study_ranges[param][0] and df_max > study_ranges[param][1]:
                continue 
        non_delta_df = df
        non_delta_param = df[param]
        non_delta_time = df.t
        non_delta_u = df["u"]
        

        res_dict = get_rPWR(df, normalise)

        # use the first 10 seconds after a 5 second buffer
        df = df[(125 <= df.t) & (df.t <= 135)]

        param_changes = df[param].mean()

        delta_sign = "Δ"

        mean_z_cco_list.append(res_dict["mean_z_cco"])
        mean_z_hhb_list.append(res_dict["mean_z_hhb"])
        mean_z_hbo2_list.append(res_dict["mean_z_hbo2"])
        rCST_hhb_list.append(-res_dict["rCST_hhb"])
        rPWR_hhb_list.append(-res_dict["rPWR_hhb"])
        rCST_hbo2_list.append(res_dict["rCST_hbo2"])
        rPWR_hbo2_list.append(res_dict["rPWR_hbo2"])  
        p_val.append(param_changes)



        # invert hhb z scores
    data_hhb = np.array([stats.zscore(mean_z_hhb_list), stats.zscore(mean_z_cco_list)])
    data_hbo2 = np.array([stats.zscore(mean_z_hbo2_list), stats.zscore(mean_z_cco_list)])

    x1 = stats.zscore(mean_z_hbo2_list)
    y1 = p_val
    z1 = stats.zscore(mean_z_cco_list)
    ax2.scatter(x1, y1, z1, \
        s=50, alpha=0.6, edgecolors='w')

    ax2.set_xlabel(f'z(ΔHbO2)')
    ax2.set_ylabel('u')
    ax2.set_zlabel(f'z(ΔoxCCO)')
    # ax1.azim =  45
    # ax1.dist = 10
    # ax1.elev = 30


    ax1.scatter(stats.zscore([-x for x in mean_z_hhb_list]), p_val, stats.zscore(mean_z_cco_list), \
    s=50, alpha=0.6, edgecolors='w')

    ax1.set_xlabel(f'z(ΔHHb)')
    ax1.set_ylabel('u')
    ax1.set_zlabel(f'z(ΔoxCCO)')
                
    

    # Add rPWR and r CST axes to plots
    max_val = get_max(data_hhb, data_hbo2)
    
    axes = np.array([[-max_val, max_val], [-max_val, max_val]])
    rotated_axes = np.matmul(rot((np.pi)/2), axes)
    for a in [ax1, ax2]:
        # a.axhline(0, color='black')
        # a.axvline(0, color='black')
        # a.annotate(s='', xy=(-max_val, -max_val), xytext=(max_val, max_val),  arrowprops=dict(arrowstyle='<-',linestyle='--', color="b"))
        # a.text(0.8, 2-0.5, "rPWR")
        # a.annotate(s='', xy=(-max_val, max_val), xytext=(max_val, -max_val), arrowprops=dict(arrowstyle='->',linestyle='--', color="r"))
        # a.text(-2+0.5, 2-0.5, "rCST")
        a.plot([-max_val, max_val], [-max_val, max_val], [-max_val, max_val], color = 'green', linestyle="--", linewidth=0.9)
        # a.set_ylabel(f'z({delta_sign}oxCCO)')
        # a.axis('scaled')
    
    ax3.scatter(stats.zscore([-x for x in mean_z_hhb_list]), stats.zscore(mean_z_cco_list), s=50, c=p_val, cmap='Greens')  
    ax3.set_title("Hhb")
    ax3.set_xlabel(f'z(ΔHHb)')

    ax4.scatter(stats.zscore(mean_z_hbo2_list), stats.zscore(mean_z_cco_list),s=50, c=p_val, cmap='Greens')
    ax4.set_title("HbO2")
    ax4.set_xlabel(f'z(ΔHbO2)')
    for a in [ax3, ax4]:
        a.axhline(0, color='black')
        a.axvline(0, color='black')
        a.annotate(s='', xy=(-max_val, -max_val), xytext=(max_val, max_val),  arrowprops=dict(arrowstyle='<-',linestyle='--', color="b"))
        a.text(0.8, 2-0.5, "rPWR")
        a.annotate(s='', xy=(-max_val, max_val), xytext=(max_val, -max_val), arrowprops=dict(arrowstyle='->',linestyle='--', color="r"))
        a.text(-2+0.5, 2-0.5, "rCST")
        a.plot([-max_val, max_val], [-max_val, max_val], color = 'black', linestyle="--", linewidth=0.2)
        a.set_ylabel(f'z({delta_sign}oxCCO)')
        a.axis('scaled')

    if save_fig == True:
        if plot_type ==  "all":
            fig.savefig(f'{folder_name}/normal_varying_{param}_{delta}.jpg', dpi=fig.dpi)
        else:
            fig.savefig(f'{folder_name}/{plot_type}_{param}_{delta}.jpg', dpi=fig.dpi)


def plot_all_plots(folder_name, plot_type=None):
    save_fig = False
    for f, p in folder_param_dict.items():
        directory= f
        time=0
        df_list = []
        folder_path = f"{folder_base}/{folder_name}/"
        for root, dirs, files in os.walk(f"{folder_path}/{directory}"):
            for subdir in dirs:
                sub_dir = os.path.join(root, subdir)
                for file in os.listdir(sub_dir):
                    filename = os.fsdecode(file)
                    if filename.endswith(f".csv"):
                        filename = os.path.join(sub_dir, file)

                        df2 = pd.read_csv(filename, delimiter="\t")
                        df_list.append(df2)
        if df_list:
            if plot_type == "paper":
                get_PWR_plots_paper(df_list, p, "incremental-diff", True, "all", False, folder_name, study_ranges=None)    
            elif plot_type == "three_d":
                get_PWR_plots_paper_3drot(df_list, p, "incremental-diff", True, "all", False, folder_name, study_ranges=None)    
            else:
                get_PWR_plots(df_list, p, "incremental-diff", True, "all", True, folder_name, study_ranges=None)    




def plot_individual_rotations(df_list, param, normalise, fig, gs, plot_column=None):
    mean_z_cco_list, mean_z_hhb_list, mean_z_hbo2_list, rPWR_hhb_list, \
        rCST_hhb_list, p_val, rPWR_hbo2_list, rCST_hbo2_list = \
            [], [], [], [], [], [], [], []
    
    dataset  = get_data_set(df_list, param, normalise)

    mean_z_hhb_list = dataset["mean_z_hhb"]
    mean_z_hbo2_list = dataset["mean_z_hbo2"]
    mean_z_cco_list = dataset["mean_z_cco"]
    p_val = dataset[param]

    # invert hhb z scores
    data_hhb = np.array([stats.zscore(mean_z_hhb_list), stats.zscore(mean_z_cco_list)])
    data_hbo2 = np.array([stats.zscore(mean_z_hbo2_list), stats.zscore(mean_z_cco_list)])
    data_hhb = np.nan_to_num(data_hhb)
    data_hbo2 = np.nan_to_num(data_hbo2)

    # rotated_data_hhb = np.matmul(rot((np.pi)/4), data_hhb)
    # rotated_data_hbo2 = np.matmul(rot((np.pi)/4), data_hbo2)

    # rPWR_hhb = rotated_data_hhb[0,:]
    # rCST_hhb = rotated_data_hhb[1,:]

    # rPWR_hbo2 = rotated_data_hbo2[0,:]
    # rCST_hbo2 = rotated_data_hbo2[1,:]

    # Plot
    fig.suptitle(f'Single parameter changes', fontsize=16)
    delta_sign = "Δ"
    ax1 = fig.add_subplot(gs[0, plot_column])
    ax1.set_ylabel('z({delta_sign}HbO2)', labelpad=0.1)
    ax2 = fig.add_subplot(gs[1, plot_column])
    ax2.set_ylabel('z({delta_sign}HHb)', labelpad=0.1)

    # This will be the title of each column
    ax1.set_title(param, fontsize=14, fontweight='bold')

    ax1.set_xlabel(f'z(ΔHHb)', labelpad=0.1) 
    ax2.set_xlabel(f'z(ΔHbO2)', labelpad=0.1)

    scatter = ax1.scatter(stats.zscore([-x for x in mean_z_hhb_list]), stats.zscore([-x for x in mean_z_cco_list]), s=50, c=p_val, cmap='Greens')

    cbar = plt.colorbar(scatter, ax=[ax1, ax2], fraction=0.09, pad=0.04)
    cbar.ax.set_ylabel(f'Max {param}', rotation=90)
    cbar.ax.get_yaxis().labelpad = 0.9

    ax2.scatter(stats.zscore(mean_z_hbo2_list), stats.zscore(mean_z_cco_list),s=50, c=p_val, cmap='Greens')


# Add rPWR and rCST axes to plots

    max_val = get_max(data_hhb, data_hbo2)
    axes = np.array([[-max_val, max_val], [-max_val, max_val]])
    rotated_axes = np.matmul(rot((np.pi)/2), axes)
    for a in [ax1, ax2]:
        a.axhline(0, color='black')
        a.axvline(0, color='black')
        a.annotate(s='', xy=(-max_val, -max_val), xytext=(max_val, max_val),  arrowprops=dict(arrowstyle='<-',linestyle='--', color="b"))
        a.text(0.3, max_val-0.5, "rPWR")
        a.annotate(s='', xy=(-max_val, max_val), xytext=(max_val, -max_val), arrowprops=dict(arrowstyle='->',linestyle='--', color="r"))
        a.text(-2, max_val-0.5, "rCST")
        a.plot([-max_val, max_val], [-max_val, max_val], color = 'black', linestyle="--", linewidth=0.2)
        a.set_ylabel(f'z({delta_sign}oxCCO)')
        a.axis('scaled')


def plot_all_rotations(folder_name, normalise, study_ranges):
    fig = plt.figure(constrained_layout=False)
    # fig.tight_layout()
    fig.set_size_inches(11.69, 5.27, forward=True)
    fig.subplots_adjust(hspace=0.01, wspace=0.6)

    gs = GridSpec(2, 4)

    column_count = 0 
    folder_base = f"/Users/mirunaserian/Desktop/PHD_Rotation2/haemodynamics-metabolism/bsx/"

    for f, p in folder_param_dict.items():
        directory= f
        time=0
        df_list = []
        folder_path = folder_base + f"/{folder_name}"
        for root, dirs, files in os.walk(f"{folder_path}/{directory}"):
            for subdir in dirs:
                sub_dir = os.path.join(root, subdir)
                for file in os.listdir(sub_dir):
                    filename = os.fsdecode(file)
                    if filename.endswith(f".csv"):
                        filename = os.path.join(sub_dir, file)
                        df2 = pd.read_csv(filename, delimiter="\t")
                        df_list.append(df2)
        if df_list:
            plot_individual_rotations(df_list, p, normalise, fig, gs, column_count)
            column_count += 1
            fig.savefig(f'{folder_base}/all_rotations_{folder_name}.jpg', dpi=fig.dpi)



folder_name = "tryjune2_x50"


# study_ranges =  {
#     "u": (0.9, 1.1),
#     "P_a": (80, 120),
#     "Pa_CO2": (38, 42),
#     "SaO2sup": (90, 96),
#     "CBFn": (0.5*0.0125, 1.5*0.0125),
#     "CMRO2": (0.5*0.033953, 1.5*0.033953),
#     "a_n": (0.5*0.6567, 1.5*0.6567)
#     } 

# plot_all_rotations(folder_name, True, study_ranges)
plot_all_plots(folder_name)

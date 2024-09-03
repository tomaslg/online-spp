# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 12:39:35 2022

@author: tol28
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from TS import format_distrib_name
from matplotlib.font_manager import FontProperties



def plot_regret_t(prior_Spatial_name,regret_NormalIG,regret_Spatial_,title_,pp,plt_):
    plt_.rcParams['text.usetex'] = True
    plt_.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    plt_.plot((list(range(1,len(regret_NormalIG)+1))), 
             regret_NormalIG, 
             label=r'Independent Gaussian',
             color='k')
    colors_ = ['r', "c", "g", 'b', "m"]
    style_ = [":", "-.", "--", "-"]
    avg_start_O_1_t = max_start_O_1_t = regret_NormalIG[0]
    cnt_ = 1
    for _style,_color, regret_Spatial, name, _title_ in zip(
            style_, colors_, regret_Spatial_, prior_Spatial_name,
            title_):
        avg_start_O_1_t = avg_start_O_1_t * cnt_ + regret_Spatial[0]
        cnt_ +=1
        avg_start_O_1_t /= cnt_
        max_start_O_1_t = max(max_start_O_1_t, regret_Spatial[0])
        plt_.plot((list(range(1,len(regret_Spatial)+1))), 
                 regret_Spatial, 
                 linestyle=_style,
                 label= name + f" ({_title_})",
                 color=_color)
    plt_.plot((list(range(1,len(regret_NormalIG)+1))), 
            [avg_start_O_1_t / ((t + 1) ** .5) for t in range(
                len(regret_NormalIG))], 
            label= r"$O(1/\sqrt{t})$",
            color=colors_[-1])
    plt_.legend(fontsize='large')
    font = FontProperties()
    font.set_family('serif')
    font_dict = {'fontsize': 'large', 'fontproperties': font}
    plt_.xlabel('t', fontdict=font_dict)
    plt_.ylabel(r'$\frac{\mathcal{R}_t}{t z^*}$ '
                , rotation=0, fontdict=font_dict)
    
    # Tweak spacing to prevent clipping of ylabel
    plt_.subplots_adjust(left=0.15)
    
    if pp:
        plt_.savefig(pp, format='pdf')
    else:
        plt_.show()
    plt_.close()



output_dir = os.path.join(os.path.dirname(
        os.getcwd()), "output_beta3_alpha1_20_150")

T_last_iteration = 150

pp = PdfPages(os.path.join(output_dir,
        f"regret_artificial_indp_vs_spatial.pdf"))

kernel_type = "exponential"
files = [os.path.join(output_dir, fl_) for fl_ in os.listdir(output_dir)
    if fl_.split(".")[-1]== "csv"# and "norm1" in fl_ 
    and kernel_type in fl_]


theta_parameters = [
                    (4/3., 8/3.), (2., 4.), 
                    (2., 8/3.), (2., 4.)
    ]
expert_objective = [0.7694762998838899,#one_border 1
                    202.49349548634692,#snake 3
                    1.4431074645596709,#two_border 4
                    0.12495338590266658]#with_obstacle 2
norm_spatial = [(2, 1), (1, 2), (2, 2), (2, 2)]
_columns = ['benchmark', 'distribution', 'theta', 'norm', 'T', 't']
instances_names = ["one_border",
                   "snake",
                   "two_border",
                   "with_obstacle"]
for k,thetas in enumerate(theta_parameters):
    file_PB = [fl_ for fl_ in files
                      if instances_names[k] in fl_
                      and f"norm{norm_spatial[k][0]}" in fl_][0]
    regret_data_PB = pd.read_csv(file_PB)
    file_PA = [fl_ for fl_ in files
                      if instances_names[k] in fl_
                      and f"norm{norm_spatial[k][1]}" in fl_][0]
    regret_data_PA = pd.read_csv(file_PA)
    data_for_plots_PB = regret_data_PB[(
        regret_data_PB["benchmark"]==kernel_type) & (
        regret_data_PB["theta"].between(thetas[0] - .1**6, thetas[0] + .1**6)
        )]
    data_for_plots_PA = regret_data_PA[(
        regret_data_PA["benchmark"]==kernel_type) & (
        regret_data_PA["theta"].between(thetas[1] - .1**6, thetas[1] + .1**6)
        )]
    
    avg_regret_df_PB = data_for_plots_PB.groupby(_columns).mean()
    avg_regret_df_PA = data_for_plots_PA.groupby(_columns).mean()
    
    Indpendent = {}
    PB = {}
    PA = {}
    for (kernel_, distrib_, theta_, norm_, last_iter, _t),(
            n_sim_,regret_) in  avg_regret_df_PB.iterrows():
        if "IGM" in distrib_:
            Indpendent[_t] = regret_ / expert_objective[k] / 2
        elif "PB" in distrib_:
            PB[_t] = regret_ / expert_objective[k]
    for (kernel_, distrib_, theta_, norm_, last_iter, _t),(
            n_sim_,regret_) in  avg_regret_df_PA.iterrows():
        if "IGM" in distrib_:
            Indpendent[_t] += regret_ / expert_objective[k] / 2
        elif "PA" in distrib_:
            PA[_t] = regret_ / expert_objective[k]
    
    title_ = [
        (f"$\\varphi = {np.round(thetas[i_],2)}$" + r" norm $=\|\cdot\|_{" +
        f"{norm_spatial[k][i_]}" + "}$")
        for i_ in [1, 0]
        ]
    
    plot_regret_t(["PA", "PB"], [Indpendent[_t] for _t in range(
        1, T_last_iteration + 1)],
        [[PA[_t] for _t in range(
            1,T_last_iteration + 1)], [PB[_t] for _t in range(
                1,T_last_iteration + 1)]],
        title_, pp, plt)
    print(Indpendent[10], Indpendent[T_last_iteration])
    
if pp: pp.close()

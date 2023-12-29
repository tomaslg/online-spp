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
             label=r'Independent Gaussian',#+f'{np.round(regret_NormalIG[1][-1],2)}', $\frac{\mathcal{R}_t}{t z^*}$ 
             color='k')
    colors_ = ['r', "c", "g", 'b', "m"]#'grey', 'lightgray', 
    style_ = [":", "-.", "--", "-"]
    avg_start_O_1_t = max_start_O_1_t = regret_NormalIG[0]
    cnt_ = 1
    for _style,_color, regret_Spatial, name in zip(
            style_, colors_, regret_Spatial_, prior_Spatial_name):
        avg_start_O_1_t = avg_start_O_1_t * cnt_ + regret_Spatial[0]
        cnt_ +=1
        avg_start_O_1_t /= cnt_
        max_start_O_1_t = max(max_start_O_1_t, regret_Spatial[0])
        plt_.plot((list(range(1,len(regret_Spatial)+1))), 
                 regret_Spatial, 
                 linestyle=_style,
                 label= name + f" ({title_})",#' = {np.round(regret_Spatial[1][-1],2)}', 
                 color=_color)
    plt_.plot((list(range(1,len(regret_NormalIG)+1))), 
            [max_start_O_1_t / ((t + 1) ** .5) for t in range(
                len(regret_NormalIG))], 
            label= r"$O(1/\sqrt{t})$",
            color=colors_[-1])
    plt_.legend(fontsize='large')
    font = FontProperties()
    font.set_family('serif')
    font_dict = {'fontsize': 'large', 'fontproperties': font}
    # plt_.title(title_, fontdict=font_dict)
    plt_.xlabel('t', fontdict=font_dict)
    plt_.ylabel(r'$\frac{\mathcal{R}_t}{t z^*}$ '#'$z^*$'
                , rotation=0, fontdict=font_dict)
    
    # Tweak spacing to prevent clipping of ylabel
    plt_.subplots_adjust(left=0.15)
    
    if pp:
        plt_.savefig(pp, format='pdf')
    else:
        plt_.show()
    plt_.close()



# output_dir = r"C:\Users\tol28\Dropbox\OTDSPP\codes\output9"
output_dir = os.path.join(os.path.dirname(
        os.getcwd()), "output_beta3_alpha1")

T_last_iteration = 50

pp = PdfPages(os.path.join(output_dir,
        f"regret_artificial_indp_vs_spatial.pdf"))


files = [os.path.join(output_dir, fl_) for fl_ in os.listdir(output_dir)
    if fl_.split(".")[-1]== "csv" and "norm1" in fl_ and "exponential" in fl_]


theta_parameters = [#8/3, 10/3, 2/3, 8/3
    #4., 4., 2/3, 8/3 #
                    # 4., 4., 
                    2., 4.
    # [4.,4.],[4.,4.],[4/3,4/3],[4.,4.]
    ]
expert_objective = [0.7694762998838899,#one_border 1
                    202.49349548634692,#snake 3
                    1.4431074645596709,#two_border 4
                    0.12495338590266658]#with_obstacle 2
# ['benchmark', 'distribution', 'theta', 'norm', 'n_sim', 'T', 't', 'regret']
_columns = ['benchmark', 'distribution', 'theta', 'norm', 'T', 't']

for k,thetas in enumerate(theta_parameters):
    regret = pd.read_csv(files[k])
    kernel_type = files[k].split("\\")[-1].split(".")[0].split("_")
    kernel_type,norm_str,instance = kernel_type[2], kernel_type[3], kernel_type[4:]
    
    data_for_plots = regret[(regret["benchmark"]==kernel_type) & (
        regret["theta"].between(thetas - .1**6, thetas + .1**6)#[0]
        )]
    
    avg_regret_df = data_for_plots.groupby(_columns).mean()
    
    Indpendent = {}
    PB = {}
    PA = {}
    for (kernel_, distrib_, theta_, norm_, last_iter, _t),(
            n_sim_,regret_) in  avg_regret_df.iterrows():
        if "IGM" in distrib_:
            Indpendent[_t] = regret_ / expert_objective[k]
        elif "PA" in distrib_:
            PA[_t] = regret_ / expert_objective[k]
        elif "PB" in distrib_:
            PB[_t] = regret_ / expert_objective[k]
    
    title_ = (
        # ("Case 1"if instance[0]=="one" else (
        # "Case 2" if instance[0]=="with" else (
        #     "Case 3" if instance[0]=="snake" else (
        #         "Case 4" if instance[0]=="two" else "")))
        # ) +
        # "Kernel {kernel_type} {norm_str} "
        f"$\\varphi = {np.round(thetas,2)}$"
        # f"$\\varphi = {np.round(thetas[0],2)}$"
        )
    
    plot_regret_t(["PA", "PB"], [Indpendent[_t] for _t in range(
        1, T_last_iteration + 1)],
        [[PA[_t] for _t in range(
            1,T_last_iteration + 1)], [PB[_t] for _t in range(
                1,T_last_iteration + 1)]],
        title_, pp, plt)
    print(Indpendent[10], Indpendent[T_last_iteration])
    
if pp: pp.close()

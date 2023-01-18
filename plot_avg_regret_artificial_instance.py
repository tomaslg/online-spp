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
             label=r'$\frac{\mathcal{R}_t}{t z^*}$ Independent Gaussian',#+f'{np.round(regret_NormalIG[1][-1],2)}', 
             color='k')
    colors_ = ['grey', 'lightgray', 'r', "c", "g", 'b', "m"]
    for _color, regret_Spatial, name in zip(colors_, regret_Spatial_, prior_Spatial_name): 
        plt_.plot((list(range(1,len(regret_Spatial)+1))), 
                 regret_Spatial, 
                 label=r'$\frac{\mathcal{R}_t}{t z^*}$ ' + name,#' = {np.round(regret_Spatial[1][-1],2)}', 
                 color=_color)
    plt_.legend(fontsize='large')
    font = FontProperties()
    font.set_family('serif')
    font_dict = {'fontsize': 'large', 'fontproperties': font}
    plt_.title(title_, fontdict=font_dict)
    plt_.xlabel('t', fontdict=font_dict)
    plt_.ylabel('$z^*$', rotation=0, fontdict=font_dict)
    
    # Tweak spacing to prevent clipping of ylabel
    plt_.subplots_adjust(left=0.15)
    
    if pp:
        plt_.savefig(pp, format='pdf')
    else:
        plt_.show()
    plt_.close()



output_dir = r"C:\Users\tol28\Dropbox\OTDSPP\codes\output9"


pp = PdfPages(os.path.join(output_dir,
        f"regret_artificial_indp_vs_spatial.pdf"))


files = [os.path.join(output_dir, fl_) for fl_ in os.listdir(output_dir)
    if fl_.split(".")[-1]== "csv" and "norm2" in fl_ and "exponential" in fl_]


theta_parameters = [[4.,4.],[4.,4.],[4/3,4/3],[4.,4.]]
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
    
    data_for_plots = regret[(regret["benchmark"]==kernel_type) & (regret[
        "theta"]==thetas[0])]
    
    avg_regret_df = data_for_plots.groupby(_columns).mean()
    
    Indpendent = {}
    PB = {}
    PA = {}
    for (kernel_, distrib_, theta_, norm_, last_iter, _t),(
            n_sim_,regret_) in  avg_regret_df.iterrows():
        if "NormalIG" in distrib_:
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
        f"$\\varphi = {np.round(thetas[0],2)}$")
    
    plot_regret_t(["PA", "PB"], [Indpendent[_t] for _t in range(1,51)],
        [[PA[_t] for _t in range(1,51)], [PB[_t] for _t in range(1,51)]],
        title_, pp, plt)
    print(Indpendent[10], Indpendent[50])
    
if pp: pp.close()

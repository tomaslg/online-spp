# -*- coding: utf-8 -*-
"""
Created on Sun Nov  6 21:08:57 2022

@author: tol28
"""


import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from TS import format_distrib_name
from matplotlib.font_manager import FontProperties

def plot_regret_t(prior_Spatial_name,regret_NormalIG,regret_Spatial_,pp,plt_):
#     # fig, ax1 = plt_.subplots()
#     # ax2 = ax1.twinx()
# # ax1.plot(x, y1, 'g-')
# # ax2.plot(x, y2, 'b-')

# # ax1.set_xlabel('X data')
# # ax1.set_ylabel('Y1 data', color='g')
#     font = {'family': 'serif',
#         # 'color':  'darkred',
#         'weight': 'normal',
#         'size': 8,
#         }
    
    plt_.plot((list(range(1,len(regret_NormalIG[1])+1))), 
             regret_NormalIG[1], 
             label=r'$\mathcal{R}_t/t$ Independent Gaussian',#+f'{np.round(regret_NormalIG[1][-1],2)}', 
             color='k')
    plt_.plot((list(range(1,len(regret_NormalIG[0])+1))), 
             regret_NormalIG[0], linestyle=':', linewidth=.2,
             label=r'$\Delta_t \quad$ Independent Gaussian',
             # f'{np.round(regret_NormalIG[0][-1],2)}', 
             color='k')
    colors_ = ['r', 'b', "g", "m", "c"]
    for _color, regret_Spatial, name in zip(colors_, regret_Spatial_, prior_Spatial_name): 
        plt_.plot((list(range(1,len(regret_Spatial[1])+1))), 
                 regret_Spatial[1], 
                 label=r'$\mathcal{R}_t/t$ '+f'{format_distrib_name(name)}',#' = {np.round(regret_Spatial[1][-1],2)}', 
                 color=_color)
        plt_.plot((list(range(1,len(regret_Spatial[0])+1))), 
                 regret_Spatial[0], linestyle=':', linewidth=.2,
                  label=r'$\Delta_t \quad$ '+f'{format_distrib_name(name)}',
                 # f'{name})/T={np.round(regret_Spatial[1][-1],2)}', 
                 color=_color)
    plt_.legend(fontsize='x-small')
    font = FontProperties()
    font.set_family('serif')
    font_dict = {'fontsize': 'large', 'fontproperties': font}
    # plt_.title(r'$\frac{\mathcal{R}_t}{t}$  $\Delta_t$', fontdict=font)
    plt_.xlabel('t', fontdict=font_dict)
    plt_.ylabel('[min]', fontdict=font_dict)
    
    # Tweak spacing to prevent clipping of ylabel
    plt_.subplots_adjust(left=0.15)
    
    if pp:
        plt_.savefig(pp, format='pdf')
    else:
        plt_.show()
    plt_.close()



def plot_regret_t_naive(prior_Spatial_name,regret_Spatial_,pp,plt_):
    colors_ = ['k', 'b','r', "g", "m", "c"]
    for _color, regret_Spatial, name in zip(colors_, regret_Spatial_, prior_Spatial_name): 
        plt_.plot((list(range(1,len(regret_Spatial[1])+1))), 
                 regret_Spatial[1], 
                 label=r'$\mathcal{R}_t/t$ '+f'{format_distrib_name(name)}',#' = {np.round(regret_Spatial[1][-1],2)}', 
                 color=_color)
        plt_.plot((list(range(1,len(regret_Spatial[0])+1))), 
                 regret_Spatial[0], linestyle=':', linewidth=.2,
                  label=r'$\Delta_t \quad$ '+f'{format_distrib_name(name)}',
                 # f'{name})/T={np.round(regret_Spatial[1][-1],2)}', 
                 color=_color)
    plt_.legend(fontsize='x-small')
    font = FontProperties()
    font.set_family('serif')
    font_dict = {'fontsize': 'large', 'fontproperties': font}
    # plt_.title(r'$\frac{\mathcal{R}_t}{t}$  $\Delta_t$', fontdict=font)
    plt_.xlabel('t', fontdict=font_dict)
    plt_.ylabel('[min]', fontdict=font_dict)
    
    # Tweak spacing to prevent clipping of ylabel
    plt_.subplots_adjust(left=0.15)
    
    if pp:
        plt_.savefig(pp, format='pdf')
    else:
        plt_.show()
    plt_.close()

        
output_dir = os.path.join(os.path.join(os.path.dirname(
    os.getcwd()), "data"), "new_large")
# r"C:\Users\tol28\Dropbox\OTDSPP\codes\data\new_large"
T = 150
instances_naive = 5
instances_shown_start = 1
instances_shown_end = 1
files = [os.path.join(output_dir, fl_) for fl_ in os.listdir(output_dir)
    if f"SIMULATION_Regret_{T}_Naive" in fl_  and fl_.split(".")[-1]== "csv"]
pp = PdfPages(os.path.join(output_dir,
        f"regret_indp_vs_naive_T{T}_N{len(files)}.pdf"))
# pp = False

regret = pd.read_csv(files[0])


for fl_ in files[1:]:
    regret += pd.read_csv(files[0])


plot_regret_t(regret.columns[instances_shown_start:(instances_shown_end + 1)],
    [np.array(regret[regret.columns[instances_naive + 1]])/len(
    files)/60, np.array(regret[regret.columns[0]])/len(files)/60], [
        [np.array(regret[regret.columns[j_]])/len(files)/60,
        np.array(regret[regret.columns[i_]])/len(files)/60]
    for i_,j_ in zip(range(instances_shown_start,(instances_shown_end + 1)),
            range((instances_naive + instances_shown_start + 1), (
                instances_naive + instances_shown_end + 2)))
    ], pp, plt)

plot_regret_t_naive(regret.columns[1:(instances_naive + 1)],
        [[np.array(regret[regret.columns[j_]])/len(files)/60,
        np.array(regret[regret.columns[i_]])/len(files)/60]
    for i_,j_ in zip(range(1,(instances_naive + 1)),
            range((instances_naive + 2), (
                instances_naive + instances_naive + 2)))
    ], pp, plt)

if pp: pp.close()


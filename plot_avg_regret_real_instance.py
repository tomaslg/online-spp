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
from collections import defaultdict

greedy_name = r"$\epsilon$-greedy"

def plot_regret_t(prior_Spatial_name,regret_NormalIG,regret_Spatial_,pp,plt_):
    
    plt_.plot((list(range(1,len(regret_NormalIG[1])+1))), 
             regret_NormalIG[1], 
             label=r'$\mathcal{R}_t/t$ Independent Gaussian',
             color='k')
    plt_.plot((list(range(1,len(regret_NormalIG[0])+1))), 
             regret_NormalIG[0], linestyle=':', linewidth=.8,
             label=r'$\Delta_t \quad$ Independent Gaussian',
             color='k')
    colors_ = ['r', 'b', "g", "m", "c"]
    for _color, regret_Spatial, name in zip(colors_, regret_Spatial_, prior_Spatial_name): 
        plt_.plot((list(range(1,len(regret_Spatial[1])+1))), 
                 regret_Spatial[1], 
                 label=r'$\mathcal{R}_t/t$ ' +
                 f'{format_distrib_name(name.replace("Naive", greedy_name))}',
                 color=_color)
        plt_.plot((list(range(1,len(regret_Spatial[0])+1))), 
                 regret_Spatial[0], linestyle=':', linewidth=.8,
                  label=r'$\Delta_t \quad$ ' +
                  f'{format_distrib_name(name.replace("Naive", greedy_name))}',
                 color=_color)
    plt_.legend(fontsize='x-small')
    font = FontProperties()
    font.set_family('serif')
    font_dict = {'fontsize': 'large', 'fontproperties': font}
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
                 label=r'$\mathcal{R}_t/t$ ' +
                 f'{format_distrib_name(name.replace("Naive", greedy_name))}',
                 color=_color)
    plt_.legend(fontsize='x-small')
    font = FontProperties()
    font.set_family('serif')
    font_dict = {'fontsize': 'large', 'fontproperties': font}
    plt_.xlabel('t', fontdict=font_dict)
    plt_.ylabel('[min]', fontdict=font_dict)
    
    # Tweak spacing to prevent clipping of ylabel
    plt_.subplots_adjust(left=0.15)
    
    if pp:
        plt_.savefig(pp, format='pdf')
    else:
        plt_.show()
    plt_.close()

if __name__ == "__main__":
    output_dir = os.path.join(os.path.join(os.path.dirname(
        os.getcwd()), "data"), "new_large")
    T = 150
    instances_naive = 5
    instances_shown_start = 1
    instances_shown_end = 1
    files = [os.path.join(output_dir, fl_) for fl_ in os.listdir(output_dir)
        if f"SIMULATION_Regret_{T}_Naive" in fl_  and fl_.split(".")[-1]== "csv"]
    pp = PdfPages(os.path.join(output_dir,
            f"regret_indp_vs_naive_T{T}_N{len(files)}.pdf"))
    
    regret = pd.read_csv(files[0])
    
    for fl_ in files[1:]:
        regret += pd.read_csv(fl_)
    
    approaches_keys = [c_ for c_ in regret.columns if "_delta" not in c_]
    stopping_time = {}
    
    for i, fl_ in enumerate(files):
        regret_realization = pd.read_csv(fl_)
        for c_ in approaches_keys:
            stopping_time[c_, i + 1] = -1
        for index, row in regret_realization.iterrows():
            regret_realization.loc[index, "period"] = index + 1
            regret_realization.loc[index, "experiment"] = i + 1
            for c_ in approaches_keys:
                if stopping_time[c_, i + 1]==-1 and regret_realization.loc[
                        index, c_ + "_delta"] <= .1 ** 6:
                    stopping_time[c_, i + 1] = index + 1
        if i==0:
            regret_all_realizations = regret_realization
        else:
            regret_all_realizations = pd.concat([
                regret_all_realizations, regret_realization], axis=0,
                ignore_index=True)
    regret_all_realizations.reset_index(inplace=True, drop=True)
    
    experiments_with_the_same_initial_regret = {regret_val : indices
     for regret_val, indices in 
     regret_all_realizations.groupby(
        "Naive_0.1").indices.items(
            ) if len(indices)>=10}
    _min_total_regret_per_run_index = {index : min([
                regret_all_realizations.loc[index + T - 1, c_]
                for c_ in approaches_keys])
        for indices in experiments_with_the_same_initial_regret.values()
        for index in indices}
    approaches_instance_convergence = defaultdict(dict)
    approaches_instance_not_converged = defaultdict(dict)
    approaches_instance_sum_total_regret = defaultdict(dict)
    approaches_instance_sum_pseudo_regret = defaultdict(dict)
    approaches_instance_num_wins = defaultdict(dict)
    for regret_val, indices in experiments_with_the_same_initial_regret.items():
        for c_ in approaches_keys:
            approaches_instance_sum_total_regret[
                regret_val][c_] = sum(regret_all_realizations.loc[
                    [index + T - 1 for index in indices], c_])
            approaches_instance_sum_pseudo_regret[
                regret_val][c_] = sum(regret_all_realizations.loc[
                    [index + T - 1 for index in indices], c_ + "_delta"])
            approaches_instance_num_wins[regret_val][c_] = 0
            for index in indices:
                approaches_instance_num_wins[regret_val][c_] += 1 if (
                    regret_all_realizations.loc[index + T - 1, c_
                        ] <= _min_total_regret_per_run_index[index]) else 0
        for index in indices:
            period, experiment = regret_all_realizations.loc[
                index, ["period", "experiment"]]
            if period>1: continue
            period, experiment = int(period), int(experiment)
            for c_ in approaches_keys:
                if c_ not in approaches_instance_convergence[regret_val]:
                    approaches_instance_convergence[regret_val][c_] = 0
                    approaches_instance_not_converged[regret_val][c_] = 0
                if stopping_time[c_, experiment] > 0:
                    approaches_instance_convergence[regret_val][c_] += 1
                else:
                    approaches_instance_not_converged[regret_val][c_] += 1
    print("regret val", end=" & ")
    for c_ in approaches_keys:
        print(c_, end=(" & " if c_!=approaches_keys[-1] else "\n"))
    for regret_val in approaches_instance_convergence.keys():
        print(f"{regret_val}", end=" & ")
        for c_ in approaches_keys:
            print(approaches_instance_convergence[regret_val][c_],
                  ",", approaches_instance_sum_pseudo_regret[regret_val][c_],
                  ",", approaches_instance_sum_total_regret[regret_val][c_],
                  end=(" & " if c_!=approaches_keys[-1] else "\n"))
        
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
    
    

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 08:42:23 2022

@author: tomas
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import colorConverter as cc
from matplotlib.backends.backend_pdf import PdfPages

from distrib import stats

def fix_long_decimal(arr):
    new_arr = []
    for a in arr:
        if type(a) is str:
            new_arr.append(a.replace("66666666666666","").replace("33333333333333",""))
        else:
            new_arr.append(a)
    return new_arr

ik=2;
def plot_mean_and_CI(plt_,mean, lb, ub, color_mean=None, color_shading=None):
    # plot the shaded range of the confidence intervals
    xx=[i for i in range(1,mean.shape[0]+1)]
    plt_.fill_between(xx, ub, lb,
                     color=color_shading, alpha=.5)
    # plot the mean on top
    plt_.plot(xx,mean, color_mean)

 
class LegendObject(object):
    def __init__(self, facecolor='red', edgecolor='white', dashed=False):
        self.facecolor = facecolor
        self.edgecolor = edgecolor
        self.dashed = dashed
 
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        patch = mpatches.Rectangle(
            # create a rectangle that is filled with color
            [x0, y0], width, height, facecolor=self.facecolor,
            # and whose edges are the faded color
            edgecolor=self.edgecolor, lw=3)
        handlebox.add_artist(patch)
 
        # if we're creating the legend for a dashed line,
        # manually add the dash in to our rectangle
        if self.dashed:
            patch1 = mpatches.Rectangle(
                [x0 + 2*width/5, y0], width/5, height, facecolor=self.edgecolor,
                transform=handlebox.get_transform())
            handlebox.add_artist(patch1)
 
        return patch

Output_dir = os.path.dirname(os.getcwd())
Output_dir = [os.path.join(Output_dir,"output_beta3_alpha1_100_50")
              ]
Name_out_file = ["path_based.pdf","path_aggregated.pdf"
                 ]

T_iter_set=range(50,51)
Theta__={"hyperbolic" : [1000.], 
              "gaussian" : [2,5/3,4/3,1], 
            "exponential" : [2,5/3,4/3,1] 
        }
T_snapshots_print = [10,
                     50]
all_instances_ = [#reversed,one-border,start-get-out-of-the-slow-zone-and-come-back,
                      #avoid-an-obstacle,snake##,chessboard-non-smooth-no-theta
        [True,True,False,False,False],#15(2)-1--15(3)-2
        [True,False,False,False,True],#15(2)-1--15(4)-2
        [True,False,False,True,False],#15(5)-1--11(5)-2
        [True,False,True,False,False],#15(5)-1--11(6)-2# 
        [True,False,False,False,False],
        [False,True,False,False,False],#11(7)-1--37(1)-2
        [False,False,False,False,True],#11(7)-1--37(2)-2
        [False,False,False,True,False],#11(8)-1--37(3)-2
        [False,False,True,False,False],#11(8)-1--37(4)-2#,
        [False,False,False,False,False]
        ]
_columns = ['benchmark', 'distribution', 'theta', 'norm', 'T', 't']

    
dataframe_all_results = {}
_columns = ['benchmark', 'theta', 'norm', 'T', 't', 'n_sim']

for output_dir in Output_dir:
    for _col in _columns:
        dataframe_all_results[_col] = []
    dataframe_all_results["instance"] = []
    dataframe_all_results["regret"+"_PA"] = []
    dataframe_all_results["regret"+"_PB"] = []
    dataframe_all_results["regret_indpendent"] = []
    dataframe_all_results["wins"+"_PA"] = []
    dataframe_all_results["wins"+"_PB"] = []
    st_ = ""
    for (_reversed__,one_border,borders_faster,obstacle_slower,snake_opt_path
         ) in all_instances_:
        for norm_ in range(1,3):
            instance___ = ("_two_border" if borders_faster else (
                                "_one_border" if one_border else (
                                "_with_obstacle" if obstacle_slower else (
                                "_chessboard" if not snake_opt_path else 
                                "_snake"))) 
                                )
            for kernel_name,kernel_thetas in Theta__.items():
                try:
                    registered_results = pd.read_csv(
                        os.path.join(
                            output_dir,
                        f"SIMULATION_Regret_{kernel_name}_"+
                        f"norm{norm_}{'_reversed' if _reversed__ else ''}"+
                            instance___+".csv"))
                    registered_results.drop(
                        registered_results[
                            registered_results["distribution"]=="Naive "
                            ].index, inplace=True)
                    registered_results.reset_index(drop=True, inplace=True)
                    for key_exp,index in registered_results.groupby(
                            _columns).groups.items():
                        for ind_col,_col in enumerate(_columns):
                            dataframe_all_results[_col
                                ].append(key_exp[ind_col])
                        instance_name_label_sep = instance___.split(
                            "_")[1:]
                        instance_name_label = ""
                        for segment in instance_name_label_sep:
                            instance_name_label += " " + segment
                        dataframe_all_results["instance"].append(
                            instance_name_label)
                        dataframe_all_results["regret_indpendent"].append(
                              registered_results.loc[min(index),"regret"])
                        #
                        ratio_regret_PA_INdp = (
                            registered_results.loc[max(index),"regret"]/
                            registered_results.loc[min(index),"regret"])
                        ind_PB = 0
                        for i_ in index:
                            ind_PB = (ind_PB if i_==min(index) or i_==max(
                                index) else i_)
                        ratio_regret_PB_INdp = (
                            registered_results.loc[ind_PB,"regret"]/
                            registered_results.loc[min(index),"regret"])
                        dataframe_all_results["regret"+"_PA"].append(
                            ratio_regret_PA_INdp)
                        dataframe_all_results["regret"+"_PB"].append(
                            ratio_regret_PB_INdp)
                        dataframe_all_results["wins"+"_PA"].append(
                            0 if ratio_regret_PA_INdp>=1. else (
                                0 if ratio_regret_PA_INdp>ratio_regret_PB_INdp 
                                else 1))
                        dataframe_all_results["wins"+"_PB"].append(
                            0 if ratio_regret_PB_INdp>=1. else (
                                0 if ratio_regret_PB_INdp>ratio_regret_PA_INdp 
                                else 1))
                except FileNotFoundError:
                    continue


dataframe_all_results = pd.DataFrame(dataframe_all_results)  
_columns = ['instance','benchmark', 'theta', 'norm', 'T', 't']
registered_results = dataframe_all_results[[col_ 
                        for col_ in dataframe_all_results.columns 
                        if not "n_sim" in col_]]
avg_regret_df = registered_results.groupby(_columns).mean()
std_regret_df = registered_results.groupby(_columns).std()
avg_regret_group_indp_df = registered_results[[col_ 
                                    for col_ in registered_results.columns
                                    if col_!="theta" and
                                    col_!='norm' and
                                    col_!='benchmark'
                                    ]].groupby([col_ 
                                    for col_ in _columns
                                    if col_!="theta" and 
                                    col_!='norm' and 
                                    col_!='benchmark'
                                    ]).mean()




st_ = ""
non_repeating_keys = ["","","","",""]
for key_exp,curve_indices in avg_regret_df.groupby([
        _col for _col in _columns if "t" != _col]
        ).groups.items():
    first_ = True
    curve = avg_regret_df.loc[curve_indices]
    print(key_exp[0])
    for key_comp_index,key_comp  in enumerate(key_exp[1:2]):
        st_ += ("" if non_repeating_keys[key_comp_index]==key_comp else ( 
            r"\multirow{xx}{*}{"+("" if not first_ else (
            key_comp if type(key_comp) is str else (
            str(key_comp) if type(key_comp) is int else 
            str(np.round(key_comp,2))))) + "}")) + " & "
        non_repeating_keys[key_comp_index] = key_comp 
    for key_comp  in key_exp[2:-1]:
        st_ += (key_comp if type(key_comp) is str else (
            str(key_comp) if type(key_comp) is int else 
            str(np.round(key_comp,2)))) + " & "
    for t_snapshots in T_snapshots_print:
        # st_ += ""
        
        st_ += str(np.round(curve["regret"
                        +"_PA"][t_snapshots-1],2))
        st_ += " & "
        st_ += str(np.round(curve["regret"
                        +"_PB"][t_snapshots-1],2))
        st_ += " & "
    for t_snapshots in T_snapshots_print:
        st_ += str(np.round(curve["wins"+"_PA"][
                                t_snapshots-1],2))
        st_ += " & "
        st_ += str(np.round(curve["wins"+"_PB"][
                                t_snapshots-1],2))
        if t_snapshots!=T_snapshots_print[-1]:
            st_ += " & "
        else:
            st_ += r" \\"
    st_ += "\n"
        

print(r"\multirow{3}{*}{\bf Kernel} & \multirow{3}{*}{\bf $\dpar$} & \multirow{3}{*}{\bf Norm} & \multicolumn{4}{c|}{\bf Pseudo-regret} & \multicolumn{4}{c}{\bf \% of wins} \\")
print(r"&  &  & \multicolumn{2}{c}{\bf $t=10$} & \multicolumn{2}{c|}{\bf $t=50$} & \multicolumn{2}{c}{\bf $t=10$} & \multicolumn{2}{c}{\bf $t=50$}  \\")
print(r" &  &  & \textbf{PA} & \textbf{PB} & \textbf{PA} & \textbf{PB} & \textbf{PA} & \textbf{PB} & \textbf{PA} & \textbf{PB} \\\hline")



print(st_)


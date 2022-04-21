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
from matplotlib.colors import colorConverter as cc # import matplotlib.colors as mcolors
from matplotlib.backends.backend_pdf import PdfPages

from distrib import stats
# from plot_util import plot_distribution,plot_histogram_sigma,get_plot

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

output_dir=os.path.join(#os.path.dirname(os.path.dirname(
    os.path.dirname(os.getcwd()),#)),
    "output6")

T_iter_set=range(50,51)
N=20
Theta__={"hyperbolic" : [1000.], 
              "gaussian" : [2,5/3,4/3,1], 
            "exponential" : [2,5/3,4/3,1] 
        }

pp=False
pp=True
if pp:pp=PdfPages(os.path.join(output_dir,"all_results.pdf"))
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
for (_reversed__,one_border,borders_faster,obstacle_slower,snake_opt_path
     ) in all_instances_:
    # for instance___ in ["_two_border","_one_border","_with_obstacle","_snake"]:#for _reversed__ in [True,False]:
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
                    f"SIMULATION_Regret_{kernel_name}_norm{norm_}{'_reversed' if _reversed__ else ''}"+
                        instance___+".csv")
                    
                    )
                registered_results = registered_results[
                    [col_ for col_ in registered_results.columns 
                     if not "mu_(" in col_ and col_!="n_sim"
                     and col_!="theta" and col_!="norm"]]
                # registered_paths = pd.read_csv(os.path.join(output_dir,f"SIMULATION_Observed_Paths_{kernel_name}_norm{norm_}.csv"))
                # Average Curve plot
                _columns = ['benchmark', 'distribution',# 'theta', 'norm', 
                            'T', 't']
                avg_regret_df = registered_results.groupby(_columns).mean()
                for keys_,curve_indices in avg_regret_df.groupby(
                        [_columns[0]]+_columns[2:-1]).groups.items():
                    curve = avg_regret_df.loc[curve_indices]
                    
                    fig, ax1 = plt.subplots()
                    spatial_curve = curve[keys_[-1]:].groupby(
                            [_columns[0]]+_columns[-2:]).mean()["regret"]
                    
                    ax1.plot(range(1,keys_[-1]+1),curve[:keys_[-1]], color="r")
                    ax1.plot(range(1,keys_[-1]+1),spatial_curve, color="b")
                    ax1.legend([curve_indices[0][1] , 
                                curve_indices[keys_[-1]][0]])
                    
                    # ax2 = ax1.twinx()
                    
                    # ax2.plot(range(1,keys_[-1]+1),(np.array(curve[:keys_[-1]])-
                    #                                      np.array(spatial_curve)), 
                    #          color="g")
                    # ax2.legend(["Difference"])
                    plt.xlabel(_columns[-1])
                    plt.ylabel("Average Regret")
                    
                    plt.suptitle(f"{[_columns[0]]+_columns[2:-2]}={fix_long_decimal(keys_[:-1])}"+
                                 f"_{'_reversed' if _reversed__ else ''}{instance___}",
                                 fontsize=8)
                    # plt.title(f"{[_columns[0]]+_columns[2:-2]}={keys_[:-1]}")
                    if pp==False:
                        plt.show()
                    else:
                        plt.savefig(pp, format='pdf')
        
                    plt.close()
                for keys_,curve_indices in registered_results.groupby(
                        [_columns[0]#,_columns[2]
                         ]).groups.items():
                    partition_results_by_approach_all = registered_results.loc[curve_indices]
                    partition_results_by_approach = partition_results_by_approach_all[
                        partition_results_by_approach_all[
                            "distribution"]=="NormalIG"].copy()
                    partition_results_by_approach.loc[:,"regret_NormalIG"] = list(
                        partition_results_by_approach["regret"])
                    partition_results_by_approach.loc[:,"regret_Spatial"] = list(
                        partition_results_by_approach_all[
                            partition_results_by_approach_all[
                                "distribution"]!="NormalIG"]["regret"])
                    partition_results_by_approach.loc[:,"distribution"] = list(
                        partition_results_by_approach_all[partition_results_by_approach_all[
                            "distribution"]!="NormalIG"]["distribution"])
                    partition_results_by_approach.loc[:,"Spatial_wins"] = (
                        .9 if False else 1
                        )*partition_results_by_approach[
                            "regret_NormalIG"]>=partition_results_by_approach["regret_Spatial"]
                    partition_results_by_approach.loc[:,"NormalIG_wins"] = (
                        .9 if False else 1
                        )*partition_results_by_approach[
                            "regret_Spatial"]>=partition_results_by_approach["regret_NormalIG"]
                    partition_results_by_approach.loc[:,"ties"] = (np.ones(
                        len(partition_results_by_approach["Spatial_wins"])) - 
                        partition_results_by_approach["Spatial_wins"] - 
                        partition_results_by_approach["NormalIG_wins"])
                    
                    partition_results_by_approach_groupby_t = partition_results_by_approach.groupby(["t"])
                    column_bars = ['Spatial_wins', 'NormalIG_wins'#,'ties'
                                   ]
                    # ax = 
                    if False:
                        partition_results_by_approach_groupby_t.sum(
                            column_bars)[column_bars].plot.bar(capsize=2,#stacked=True,
                                color=['b','orange','g'],error_kw={'ecolor': '0.3'},
                                yerr=partition_results_by_approach_groupby_t.std(
                                    )[column_bars]*stats.norm().ppf(.975)*500**.5)#,rot=0)#(x="t",y="NormalIG_wins")
                        plt.legend(["NormalIG" , f"Spatial_{fix_long_decimal(keys_)}", "Ties" ])
                    else:
                        colors___ = ["b","r","g","k"]
                        fig = plt.figure(1, figsize=(15, 5))
                        for ind_colr,colr in enumerate(colors___):
                            if len(column_bars)<=ind_colr: continue
                            err__ = partition_results_by_approach_groupby_t.std(
                                )[column_bars[ind_colr]]*stats.norm(
                                    ).ppf(.975)/500**.5
                            mean__ = partition_results_by_approach_groupby_t.mean(
                                )[column_bars[ind_colr]]
                            plot_mean_and_CI(plt,mean__, mean__-err__,
                                             mean__+err__, color_mean=colr, 
                                             color_shading=colr)
                        # plt.xscale('log',basex=4)
                        # plt.xlim(4**ik,4**(k-1))
                        # plt.ylim(0.975, 1.002);
                        bg = np.array([1, 1, 1])  # background of the legend is white
                        colors_faded = [(np.array(cc.to_rgb(colr)) + bg
                                         ) / 2.0 for colr in colors___]
                        plt.legend(list(range(min(len(column_bars)
                                                  ,len(colors___)))), 
                                   column_bars[:min(len(column_bars)
                                                  ,len(colors___))],
                            handler_map={__i_ : LegendObject(
                                colors___[__i_], colors_faded[__i_]) for 
                                __i_ in range(min(len(column_bars)
                                                  ,len(colors___)))})
                    # fig = plt.figure(1, figsize=(15*k/7, 5))
                    # if adp_b:
                    #     plot_mean_and_CI(plt,mean1, ub1, lb1, color_mean='b', color_shading='b')
                    # plot_mean_and_CI(plt,mean0, ub0, lb0, color_mean='k', color_shading='k')
                    # if rh_b:
                    #     plot_mean_and_CI(plt,mean2, ub2, lb2, color_mean='g', color_shading='g')
                    # plot_mean_and_CI(plt,mean3, ub3, lb3, color_mean='r', color_shading='r')
                    # plt.xscale('log',basex=4)
                    # plt.xlim(4**ik,4**(k-1))
                    # plt.ylim(0.975, 1.002);
                    
                    # bg = np.array([1, 1, 1])  # background of the legend is white
                    # if adp_b:
                    #     colors = ['black', 'blue', 'green', 'red']
                    # elif rh_b:
                    #     colors = ['black', 'green', 'red']
                    # else:
                    #     colors = ['black', 'green', 'red']
                    # # with alpha = .5, the faded color is the average of the background and color
                    # colors_faded = [(np.array(cc.to_rgb(color)) + bg) / 2.0 for color in colors]
                     
                    # if adp_b:
                    #     plt.legend([0, 1, 2, 3], ['2S', 'ADP', 'RH', 'PK'],
                    #            handler_map={
                    #                0: LegendObject(colors[0], colors_faded[0]),
                    #                1: LegendObject(colors[1], colors_faded[1]),
                    #                2: LegendObject(colors[2], colors_faded[2]),
                    #                3: LegendObject(colors[3], colors_faded[3])
                    #             })
                    # elif rh_b:
                    #     plt.legend([0, 2, 3], ['2S', 'RH', 'PK'],
                    #            handler_map={
                    #                0: LegendObject(colors[0], colors_faded[0]),
                    #                2: LegendObject(colors[1], colors_faded[1]),
                    #                3: LegendObject(colors[2], colors_faded[2])
                    #             })
                    # else:
                    #     plt.legend([0, 3], ['2S', 'PK'],
                    #            handler_map={
                    #                0: LegendObject(colors[0], colors_faded[0]),
                    #                3: LegendObject(colors[2], colors_faded[2])
                    #             })
                     
                    # plt.xlabel("#drill-holes");
                    # plt.ylabel("gap_wr_PK");
                    # plt.title('mean and confidence interval plot')
                    # #plt.show()
                    
                    # plt.savefig(pp, format='pdf')
                    # pp.close()                                       
                            
                    
                    
# means.plot.bar(yerr=errors, ax=ax, capsize=4, rot=0);
                    # partition_results_by_approach[["NormalIG_wins", "Spatial_wins",'ties']].plot.box(#boxplot(column=, 
                    #     by=["t"])
                    # ax =partition_results_by_approach.groupby(["t"]).sum(
                    #     ['Spatial_wins', 'NormalIG_wins'])[].plot.bar()#(x="t",y="Spatial_wins")
                    # ax =partition_results_by_approach.groupby(["t"]).sum(
                    #     ['Spatial_wins', 'NormalIG_wins',"ties"]
                    #     )['NormalIG_wins'].plot(x="t",y="NormalIG_wins")
                    # ax =partition_results_by_approach.groupby(["t"]).sum(
                    #     ['Spatial_wins', 'NormalIG_wins',"ties"]
                    #     )["Spatial_wins"].plot(x="t",y="Spatial_wins")
                    # ax =partition_results_by_approach.groupby(["t"]).sum(
                    #     ['Spatial_wins', 'NormalIG_wins',"ties"]
                    #     )["ties"].plot(x="t",y="ties")
                    
                    plt.suptitle(f"{_columns[1:-2]}={fix_long_decimal(keys_)}{norm_}_"+
                                 f"{'_reversed' if _reversed__ else ''}{instance___}",
                                 fontsize=12)
                    plt.ylabel("Number_wins")
                    plt.tight_layout()
                    plt.grid()
                    if pp==False:
                        plt.show()
                    else:
                        plt.savefig(pp, format='pdf')
        
                    plt.close()
                    # partition_results_by_approach.merge(
                    #     partition_results_by_approach_all[partition_results_by_approach_all["distribution"]!="NormalIG"],
                    #     how="right",on=_columns[:1]+_columns[2:])
                       #left_on=_columns[:1]+_columns[2:],right_on=_columns[:1]+_columns[2:],suffixes=('_NormalIG', '_Spatial'))
                    # print("")
            except FileNotFoundError:
                continue



if pp:pp.close()

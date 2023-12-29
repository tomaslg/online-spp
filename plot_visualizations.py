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

# os.path.join(#os.path.dirname(os.path.dirname(
Output_dir = os.path.dirname(os.getcwd())#)),
    # "output5")
Output_dir = [os.path.join(Output_dir,"output_beta3_alpha1")#,os.path.join(Output_dir,"output7")
              ]
Name_out_file = ["path_based.pdf","path_aggregated.pdf"
                 ]

T_iter_set=range(50,51)
Theta__={"hyperbolic" : [1000.], 
              "gaussian" : [2,5/3,4/3,1], 
            "exponential" : [2,5/3,4/3,1] 
        }
T_snapshots_print = [10,#20,30,40,
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

# for output_dir,name_out_file in zip(Output_dir,Name_out_file):
#     pp=False
#     pp=True
#     if pp:pp=PdfPages(os.path.join(output_dir,name_out_file))
#     for (_reversed__,one_border,borders_faster,obstacle_slower,snake_opt_path
#          ) in all_instances_:
#         # for instance___ in ["_two_border","_one_border","_with_obstacle","_snake"]:#for _reversed__ in [True,False]:
#         for norm_ in range(1,3):
#             instance___ = ("_two_border" if borders_faster else (
#                                 "_one_border" if one_border else (
#                                 "_with_obstacle" if obstacle_slower else (
#                                 "_chessboard" if not snake_opt_path else 
#                                 "_snake"))) 
#                                 )
#             for kernel_name,kernel_thetas in Theta__.items():
#                 try:
#                     registered_results = pd.read_csv(
#                         os.path.join(
#                             output_dir,
#                         f"SIMULATION_Regret_{kernel_name}_norm{norm_}{'_reversed' if _reversed__ else ''}"+
#                             instance___+".csv")
                        
#                         )
#                     registered_results = registered_results[[col_ for col_ in registered_results.columns if not "mu_(" in col_ and col_!="n_sim"]]
#                     # registered_paths = pd.read_csv(os.path.join(output_dir,f"SIMULATION_Observed_Paths_{kernel_name}_norm{norm_}.csv"))
#                     # Average Curve plot
#                     # _columns = ['benchmark', 'distribution', 'theta', 'norm', 'T', 't']
#                     avg_regret_df = registered_results.groupby(_columns).mean()
#                     for keys_,curve_indices in avg_regret_df.groupby(
#                             [_columns[0]]+_columns[2:-1]).groups.items():
#                         curve=avg_regret_df.loc[curve_indices]
                        
#                         fig, ax1 = plt.subplots()
                        
#                         ax1.plot(range(1,keys_[-1]+1),curve[:keys_[-1]], color="r")
#                         ax1.plot(range(1,keys_[-1]+1),curve[keys_[-1]:], color="b")
#                         ax1.legend([curve_indices[0][1] , 
#                                     curve_indices[keys_[-1]][1]])
                        
#                         ax2 = ax1.twinx()
                        
#                         ax2.plot(range(1,keys_[-1]+1),(np.array(curve[:keys_[-1]])-
#                                                              np.array(curve[keys_[-1]:])), 
#                                  color="g")
#                         ax2.legend(["Difference"])
#                         plt.xlabel(_columns[-1])
#                         plt.ylabel("Average Regret")
                        
#                         plt.suptitle(f"{[_columns[0]]+_columns[2:-2]}={fix_long_decimal(keys_[:-1])}"+
#                                      f"_{'_reversed' if _reversed__ else ''}{instance___}",
#                                      fontsize=8)
#                         # plt.title(f"{[_columns[0]]+_columns[2:-2]}={keys_[:-1]}")
#                         if pp==False:
#                             plt.show()
#                         else:
#                             plt.savefig(pp, format='pdf')
            
#                         plt.close()
#                     for keys_,curve_indices in registered_results.groupby(
#                             [_columns[0],_columns[2]]).groups.items():
#                         partition_results_by_approach_all = registered_results.loc[curve_indices]
#                         partition_results_by_approach = partition_results_by_approach_all[
#                             partition_results_by_approach_all[
#                                 "distribution"]=="NormalIG"].copy()
#                         partition_results_by_approach.loc[:,"regret_NormalIG"] = list(
#                             partition_results_by_approach["regret"])
#                         partition_results_by_approach.loc[:,"regret_Spatial"] = list(
#                             partition_results_by_approach_all[
#                                 partition_results_by_approach_all[
#                                     "distribution"]!="NormalIG"]["regret"])
#                         partition_results_by_approach.loc[:,"distribution"] = list(
#                             partition_results_by_approach_all[partition_results_by_approach_all[
#                                 "distribution"]!="NormalIG"]["distribution"])
#                         partition_results_by_approach.loc[:,"Spatial_wins"] = (
#                             .9 if False else 1
#                             )*partition_results_by_approach[
#                                 "regret_NormalIG"]>=partition_results_by_approach["regret_Spatial"]
#                         partition_results_by_approach.loc[:,"NormalIG_wins"] = (
#                             .9 if False else 1
#                             )*partition_results_by_approach[
#                                 "regret_Spatial"]>=partition_results_by_approach["regret_NormalIG"]
#                         partition_results_by_approach.loc[:,"ties"] = (np.ones(
#                             len(partition_results_by_approach["Spatial_wins"])) - 
#                             partition_results_by_approach["Spatial_wins"] - 
#                             partition_results_by_approach["NormalIG_wins"])
                        
#                         partition_results_by_approach_groupby_t = partition_results_by_approach.groupby(["t"])
#                         column_bars = ['Spatial_wins', 'NormalIG_wins'#,'ties'
#                                        ]
#                         # ax = 
#                         if False:
#                             partition_results_by_approach_groupby_t.sum(
#                                 column_bars)[column_bars].plot.bar(capsize=2,#stacked=True,
#                                     color=['b','orange','g'],error_kw={'ecolor': '0.3'},
#                                     yerr=partition_results_by_approach_groupby_t.std(
#                                         )[column_bars]*stats.norm().ppf(.975)*500**.5)#,rot=0)#(x="t",y="NormalIG_wins")
#                             plt.legend(["NormalIG" , f"Spatial_{fix_long_decimal(keys_)}", "Ties" ])
#                         else:
#                             colors___ = ["b","r","g","k"]
#                             fig = plt.figure(1, figsize=(15, 5))
#                             for ind_colr,colr in enumerate(colors___):
#                                 if len(column_bars)<=ind_colr: continue
#                                 err__ = partition_results_by_approach_groupby_t.std(
#                                     )[column_bars[ind_colr]]*stats.norm(
#                                         ).ppf(.975)/100**.5
#                                 mean__ = partition_results_by_approach_groupby_t.mean(
#                                     )[column_bars[ind_colr]]
#                                 plot_mean_and_CI(plt,mean__, mean__-err__,
#                                                  mean__+err__, color_mean=colr, 
#                                                  color_shading=colr)
#                             # plt.xscale('log',basex=4)
#                             # plt.xlim(4**ik,4**(k-1))
#                             # plt.ylim(0.975, 1.002);
#                             bg = np.array([1, 1, 1])  # background of the legend is white
#                             colors_faded = [(np.array(cc.to_rgb(colr)) + bg
#                                              ) / 2.0 for colr in colors___]
#                             plt.legend(list(range(min(len(column_bars)
#                                                       ,len(colors___)))), 
#                                        column_bars[:min(len(column_bars)
#                                                       ,len(colors___))],
#                                 handler_map={__i_ : LegendObject(
#                                     colors___[__i_], colors_faded[__i_]) for 
#                                     __i_ in range(min(len(column_bars)
#                                                       ,len(colors___)))})                              
#     # means.plot.bar(yerr=errors, ax=ax, capsize=4, rot=0);
#                         # partition_results_by_approach[["NormalIG_wins", "Spatial_wins",'ties']].plot.box(#boxplot(column=, 
#                         #     by=["t"])
#                         # ax =partition_results_by_approach.groupby(["t"]).sum(
#                         #     ['Spatial_wins', 'NormalIG_wins'])[].plot.bar()#(x="t",y="Spatial_wins")
#                         # ax =partition_results_by_approach.groupby(["t"]).sum(
#                         #     ['Spatial_wins', 'NormalIG_wins',"ties"]
#                         #     )['NormalIG_wins'].plot(x="t",y="NormalIG_wins")
#                         # ax =partition_results_by_approach.groupby(["t"]).sum(
#                         #     ['Spatial_wins', 'NormalIG_wins',"ties"]
#                         #     )["Spatial_wins"].plot(x="t",y="Spatial_wins")
#                         # ax =partition_results_by_approach.groupby(["t"]).sum(
#                         #     ['Spatial_wins', 'NormalIG_wins',"ties"]
#                         #     )["ties"].plot(x="t",y="ties")
                        
#                         plt.suptitle(f"{_columns[1:-2]}={fix_long_decimal(keys_)}{norm_}_"+
#                                      f"{'_reversed' if _reversed__ else ''}{instance___}",
#                                      fontsize=12)
#                         plt.ylabel("Number_wins")
#                         plt.tight_layout()
#                         plt.grid()
#                         if pp==False:
#                             plt.show()
#                         else:
#                             plt.savefig(pp, format='pdf')
            
#                         plt.close()
#                         # partition_results_by_approach.merge(
#                         #     partition_results_by_approach_all[partition_results_by_approach_all["distribution"]!="NormalIG"],
#                         #     how="right",on=_columns[:1]+_columns[2:])
#                            #left_on=_columns[:1]+_columns[2:],right_on=_columns[:1]+_columns[2:],suffixes=('_NormalIG', '_Spatial'))
#                         # print("")
#                 except FileNotFoundError:
#                     continue
#     if pp!=False:pp.close()
    
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
        # for instance___ in ["_two_border","_one_border","_with_obstacle","_snake"]:#for _reversed__ in [True,False]:
        # attempt_successful = False
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
                                # +"_"+name_out_file.replace(".pdf", "")
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
# first_indicator = True
# for output_dir,name_out_file in zip(Output_dir,Name_out_file):
#     if first_indicator:
#         for _col in _columns:
#             dataframe_all_results[_col] = []
#         dataframe_all_results["instance"] = []
#     dataframe_all_results["regret"+"_"+name_out_file.replace(".pdf", "")] = []
#     dataframe_all_results["wins"+"_"+name_out_file.replace(".pdf", "")] = []
#     dataframe_all_results["regret_indpendent"+
#                           "_"+name_out_file.replace(".pdf", "")] = []
#     st_ = ""
#     for (_reversed__,one_border,borders_faster,obstacle_slower,snake_opt_path
#          ) in all_instances_:
#         # for instance___ in ["_two_border","_one_border","_with_obstacle","_snake"]:#for _reversed__ in [True,False]:
#         # attempt_successful = False
#         for norm_ in range(1,3):
#             instance___ = ("_two_border" if borders_faster else (
#                                 "_one_border" if one_border else (
#                                 "_with_obstacle" if obstacle_slower else (
#                                 "_chessboard" if not snake_opt_path else 
#                                 "_snake"))) 
#                                 )
#             for kernel_name,kernel_thetas in Theta__.items():
#                 try:
#                     registered_results = pd.read_csv(
#                         os.path.join(
#                             output_dir,
#                         f"SIMULATION_Regret_{kernel_name}_"+
#                         f"norm{norm_}{'_reversed' if _reversed__ else ''}"+
#                             instance___+".csv"))
#                     for key_exp,index in registered_results.groupby(
#                             _columns).groups.items():
#                         if first_indicator:
#                             for ind_col,_col in enumerate(_columns):
#                                 dataframe_all_results[_col
#                                     # +"_"+name_out_file.replace(".pdf", "")
#                                     ].append(key_exp[ind_col])
#                             instance_name_label_sep = instance___.split(
#                                 "_")[1:]
#                             instance_name_label = ""
#                             for segment in instance_name_label_sep:
#                                 instance_name_label += " " + segment
#                             dataframe_all_results["instance"].append(
#                                 instance_name_label)
#                         dataframe_all_results["regret_indpendent"+
#                           "_"+name_out_file.replace(".pdf", "")].append(
#                               registered_results.loc[min(index),"regret"])
#                         ratio_regret = (
#                             registered_results.loc[max(index),"regret"]/
#                             registered_results.loc[min(index),"regret"])
#                         dataframe_all_results["regret"
#                             +"_"+name_out_file.replace(".pdf", "")].append(
#                                 ratio_regret)
#                         dataframe_all_results["wins"
#                             +"_"+name_out_file.replace(".pdf", "")].append(
#                                 1 if ratio_regret<1. else 0)
#                 except FileNotFoundError:
#                     continue
#     first_indicator = False

dataframe_all_results = pd.DataFrame(dataframe_all_results)  
_columns = ['instance','benchmark', 'theta', 'norm', 'T', 't']
# for name_out_file in Name_out_file:
#     for _col in :
#         _columns.append(_col+"_"+name_out_file.replace(".pdf", ""))    
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


# st_ = ""
# non_repeating_keys = ["","","","",""]
# for key_exp,curve_indices in avg_regret_df.groupby([
#         _col for _col in _columns if "t" != _col]
#         ).groups.items():
#     curve = avg_regret_df.loc[curve_indices]
#     curve_err = std_regret_df.loc[curve_indices]/500**.5
#     for key_comp_index,key_comp  in enumerate(key_exp[:2]):
#         st_ += ("" if non_repeating_keys[key_comp_index]==key_comp else ( 
#             r"\multirow{xx}{*}{"+(
#             key_comp if type(key_comp) is str else (
#             str(key_comp) if type(key_comp) is int else 
#             str(np.round(key_comp,2)))) + "}")) + " & "
#         non_repeating_keys[key_comp_index] = key_comp 
#     for key_comp  in key_exp[2:-1]:
#         st_ += (key_comp if type(key_comp) is str else (
#             str(key_comp) if type(key_comp) is int else 
#             str(np.round(key_comp,2)))) + " & "
#     for t_snapshots in T_snapshots_print:
#         # st_ += ""
        
#         st_ += str(np.round(curve["regret"
#                         +"_PA"][t_snapshots-1],2))
#         st_ += " & "
#         st_ += str(np.round(curve["regret"
#                         +"_PB"][t_snapshots-1],2))
#         st_ += " & "
#         st_ += str(np.round(avg_regret_group_indp_df.loc[
#                 key_exp[:1]+key_exp[-1:]+(t_snapshots,)][
#                     "regret_indpendent"],2))+" & "
#     for t_snapshots in T_snapshots_print:
#         # st_ += ""
#         st_ += str(np.round(curve["wins"+"_PA"][
#                                 t_snapshots-1],2))
#         st_ += " & "
#         st_ += str(np.round(curve["wins"+"_PB"][
#                                 t_snapshots-1],2))
#         if t_snapshots!=T_snapshots_print[-1]:
#             st_ += " & "
#         else:
#             st_ += r" \\\hline"
#     st_ += "\n"
        
# print(r"\multirow{3}{*}{Instance} & \multirow{3}{*}{Kernel} & "+
#       r"\multirow{3}{*}{$\theta$} & \multirow{3}{*}{Norm} & "+
#       r"\multicolumn{6}{|c|}{"+"Pseudo-regret} & "+
#       r"\multicolumn{4}{|c|}{"+r"\% of wins}"+
#       r" \\\cline{5-14}")
# # (PB/Ind,PA/Ind)
# print(" &  &  &  & "+
#       r"\multicolumn{3}{|c|}{"+"$t=10$} & "+
#       r"\multicolumn{3}{|c|}{"+"$t=50$} & "+
#       r"\multicolumn{2}{|c|}{"+"$t=10$} & "+
#       r"\multicolumn{2}{|c|}{"+"$t=50$} "+
#       r" \\\cline{5-14}")
# print(" &  &  &  & PA & PB & Ind. & PA & PB & Ind. & PA & PB & PA & PB"+
#       r" \\\hline")
# # print(" &  &  &  & $t=10$ & $t=20$ & $t=30$ & $t=40$ & $t=50$"+
# #       r" \\\hline")


st_ = ""
non_repeating_keys = ["","","","",""]
for key_exp,curve_indices in avg_regret_df.groupby([
        _col for _col in _columns if "t" != _col]
        ).groups.items():
    first_ = True
    curve = avg_regret_df.loc[curve_indices]
    # curve_err = std_regret_df.loc[curve_indices]/500**.5
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
        # st_ += str(np.round(avg_regret_group_indp_df.loc[
        #         key_exp[:1]+key_exp[-1:]+(t_snapshots,)][
        #             "regret_indpendent"],2))+" & "
    for t_snapshots in T_snapshots_print:
        # st_ += ""
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
# \multirow{12}{*}{exponential} & \multirow{2}{*}{0.67} & 1 & 0.82 & 0.86 & 0.68 & 0.72 & 0.48 & 0.41 & 0.5 & 0.43 \\ 
# &  & 2 & 0.7 & 0.78 & 0.59 & 0.65 & 0.55 & 0.38 & 0.58 & 0.4 \\ \cline{2-11}
# & \multirow{2}{*}{1.0} & 1 & 0.67 & 0.72 &  0.57 & 0.61 &  0.52 & 0.41 & 0.54 & 0.43 \\
# &  & 2 & 0.6 & 0.69 &  0.5 & 0.57 &  0.57 & 0.4 & 0.61 & 0.38 \\ \cline{2-11}
# & \multirow{2}{*}{1.33} & 1 & 0.58 & 0.63 & 0.49 & 0.54 &  0.52 & 0.45 & 0.54 & 0.44 \\
# &  & 2 & 0.57 & 0.63 &  0.47 & 0.54 & 0.55 & 0.42 & 0.57 & 0.41 \\ \cline{2-11}
# & \multirow{2}{*}{1.67} & 1 & 0.57 & 0.63 &  0.47 & 0.53 &  0.57 & 0.4 & 0.61 & 0.38 \\
# &  & 2 & 0.55 & 0.59 & 0.45 & 0.5 & 0.54 & 0.43 & 0.6 & 0.39 \\ \cline{2-11}
# & \multirow{2}{*}{2.0} & 1 & 0.55 & 0.59 &  0.46 & 0.52 & 0.52 & 0.46 & 0.58 & 0.41 \\
# &  & 2 & 0.52 & 0.58 & 0.43 & 0.5 & 0.58 & 0.39 & 0.63 & 0.36 \\ \cline{2-11}
# & \multirow{2}{*}{2.33} & 1 & 0.53 & 0.63 & 0.43 & 0.54 & 0.59 & 0.39 & 0.64 & 0.35 \\
# &  & 2 & 0.5 & 0.54 & 0.41 & 0.48 & 0.56 & 0.42 & 0.63 & 0.36 \\ \hline 
# \multirow{6}{*}{gaussian} & 0.67 & \multirow{6}{*}{2} & 0.9 & 0.97  & 0.77 & 0.82 & 0.43 & 0.36 & 0.49 & 0.41 \\
# & 1.0 &  & 0.82 & 0.85 & 0.66 & 0.68 & 0.48 & 0.41 & 0.52 & 0.44 \\
# & 1.33 &  & 0.75 & 0.74 & 0.69 & 0.62 & 0.43 & 0.48 & 0.39 & 0.58 \\
# & 1.67 &  & 0.92 & 0.72 & 1.26 & 0.6 & 0.33 & 0.57 & 0.03 & 0.9 \\
# & 2.0 &  & 1.24 & 0.7 & 2.01 & 0.58 & 0.13 & 0.73 & 0.0 & 0.89 \\
# & 2.33 &  & 1.49 & 0.7 & 2.35 & 0.58 & 0.09 & 0.74 & 0.0 & 0.9 \\\hline


print(st_)


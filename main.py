#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 16:02:14 2021

@author: tomas
"""


# import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats
from scipy.sparse import dok_matrix

from time import time
from _meta_util import timer_func
from distrib import (NormalIG,Spatial_Metric,Spatial_Metric_2,Spatial_Metric3,
                     np)
from TS import get_nodes_and_edges,TS_Stationary,nx
    
from plot_util import (plot_distribution,plot_histogram_sigma,get_plot,
                       get_graph_instance_1,get_graph_instance_2,plot_regret_t)


# def shperical(rho):
#     return (1- 1.5*rho + .5*rho**3 )

def hyperbolic(theta,rho):
    return 1/(1 + rho/theta)

def gaussian(theta,rho):
    return np.exp(- (rho/theta)**2)

def exponential(theta,rho):
    return np.exp(-rho/theta)
Kernels={"hyperbolic" : hyperbolic , 
          "gaussian" : gaussian , 
        "exponential" : exponential
         }
def general_multiplication(__theta_,dist__):
    if type(__theta_) is str:
        __theta_ = [float(i) for i in __theta_.split("_")]
    if type(__theta_) is list:
        st_=""
        for i in range(len(__theta_)):
            __theta_[i] = __theta_[i]*dist__
            st_ += ("_" if i>0 else "") + str(__theta_[i]) 
        return st_#__theta_
    else:
        return __theta_*dist__


output_dir=os.path.join(os.path.dirname(os.getcwd()),"output7")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

norm_range_ = range(1,3)
T_iter_set = range(50,51)
N = 20
nruns = 100
output_data_file = True
___theta____ = [2/3,1,4/3,
                5/3,2,7/3]
all_instances_ = [#reversed,one-border,start-get-out-of-the-slow-zone-and-come-back,
                  #avoid-an-obstacle,snake##,chessboard-non-smooth-no-theta
    # [True,True,False,False,False],#37(1)
    # [True,False,False,False,True],#37(1)
    # [True,False,False,True,False],#37(2)
    # [True,False,True,False,False],#37(2) 
    # [True,False,False,False,False],#37(1)
    # [False,True,False,False,False],#11(1)
    # [False,False,False,False,True],#11(2)
    # [False,False,False,True,False],#11(3)
    [False,False,True,False,False],#11(4)
    # [False,False,False,False,False]#11(1)
    ]

# _reversed__ = True
# borders_faster = True
# obstacle_slower = False
# chessboard_instance = False
for (_reversed__,one_border,borders_faster,obstacle_slower,snake_opt_path
     ) in all_instances_:
        
    
    if _reversed__:
        if one_border:
            #one border
            re_do_list = [
                # (1,"exponential",'0.6666666666666666_2'),
                #   (1,"exponential",'1_1.6666666666666667'),
                #   (1,"exponential",'1_2'),
                #   (1,"exponential",'2.3333333333333335_2.3333333333333335'),
                #   (2,"gaussian",'0.6666666666666666_0.6666666666666666'),
                #   (2,"gaussian",'2.3333333333333335_2.3333333333333335'),
                #   # (2,"gaussian",'1.6666666666666667_1.3333333333333333'),
                #   # (2,"gaussian",'2.3333333333333335_2'),
                #   (2,"exponential",'1.6666666666666667_1.3333333333333333'),
                #   (2,"exponential",'2.3333333333333335_1.3333333333333333')
                #   # (2,"exponential",'1_0.6666666666666666')
                  ]
        # elif borders_faster or chessboard_instance_bool:
            # re_do_list = []
        elif snake_opt_path:
            re_do_list = [
                # (1,"exponential",'0.6666666666666666_2.3333333333333335'),
                #           (1,"exponential",'1_2.3333333333333335'),
                #           (1,"exponential",'2.3333333333333335_1.3333333333333333'),
                #           (2,"gaussian",'1.3333333333333333_0.6666666666666666'),
                #           (2,"gaussian",'1_0.6666666666666666'),
                #           (2,"gaussian",'2.3333333333333335_1.3333333333333333'),
                #           (2,"exponential",'1.3333333333333333_0.6666666666666666'),
                #           (2,"exponential",'1_2.3333333333333335'),
                ]
        elif obstacle_slower:
            #w/ obstacle reversed
            re_do_list = [
                # #(1,"exponential",'0.6666666666666666_0.6666666666666666'),
                #           # (1,"exponential",'0.6666666666666666_2'),
                #           (1,"exponential",'0.6666666666666666_2.3333333333333335'),
                #           # (1,"exponential",'1.3333333333333333_0.6666666666666666'),
                #           (1,"exponential",'1_1.6666666666666667'),
                #           # (1,"exponential",'1.6666666666666667_0.6666666666666666'),
                #           # (1,"exponential",'1_2'),
                #           # (1,"exponential",'2.3333333333333335_0.6666666666666666'),
                #           # (1,"exponential",'2.3333333333333335_1'),
                #           # (1,"exponential",'2.3333333333333335_2.3333333333333335'),
                #           # (1,"exponential",'2.3333333333333335_1.3333333333333333'),
                #           # (1,"exponential",'2.3333333333333335_1.6666666666666667'),
                #           # (1,"exponential",'2_0.6666666666666666'),
                #           # (1,"exponential",'2_2'),
                #           (2,"gaussian",'0.6666666666666666_2.3333333333333335'),
                #           # (2,"gaussian",'0.6666666666666666_1.3333333333333333'),
                #           (2,"gaussian",'1.3333333333333333_0.6666666666666666'),
                #           # (2,"gaussian",'1.3333333333333333_1.6666666666666667'),
                #           # (2,"gaussian",'1.3333333333333333_2'),
                #           # (2,"gaussian",'1.6666666666666667_1.6666666666666667'),
                #           (2,"gaussian",'1_2'),
                #           (2,"gaussian",'2.3333333333333335_1.6666666666666667'),
                #           (2,"gaussian",'2.3333333333333335_2.3333333333333335'),
                #           # (2,"gaussian",'2_1'),
                #           # (2,"exponential",'1.6666666666666667_1.3333333333333333'),
                #           # (2,"exponential",'0.6666666666666666_1.3333333333333333'),
                #           (2,"exponential",'0.6666666666666666_2'),
                #           # (2,"exponential",'1_0.6666666666666666'),
                #           (2,"exponential",'1_2'),
                #           # (2,"exponential",'1_23333333333333335'),
                #           (2,"exponential",'23333333333333335_23333333333333335')
                  ]
        elif borders_faster:
            re_do_list = [
                # (1,"exponential",'0.6666666666666666_0.6666666666666666'),
                #           (1,"exponential",'0.6666666666666666_2'),
                #           (1,"exponential",'1_0.6666666666666666'),
                #           (1,"exponential",'1.3333333333333333_0.6666666666666666'),
                #           # (1,"exponential",'1_1.6666666666666667'),
                #           (1,"exponential",'1_2'),
                #           # (1,"exponential",'1.3333333333333333_1.6666666666666667'),
                #           # (1,"exponential",'1.6666666666666667_0.6666666666666666'),
                #           # (1,"exponential",'1_2.3333333333333335'),
                #           # (1,"exponential",'2.3333333333333335_0.6666666666666666'),
                #           # (1,"exponential",'2.3333333333333335_1'),
                #           # (1,"exponential",'2.3333333333333335_1.3333333333333333'),
                #           (1,"exponential",'2.3333333333333335_2.3333333333333335'),
                #           # (1,"exponential",'2_0.6666666666666666'),
                #           # (1,"exponential",'2_2'),
                #           # (2,"gaussian",'0.6666666666666666_0.6666666666666666'),
                #           # (2,"gaussian",'0.6666666666666666_1.3333333333333333'),
                #           (2,"gaussian",'1_0.6666666666666666'),
                #           # (2,"gaussian",'1.3333333333333333_0.6666666666666666'),
                #           # (2,"gaussian",'1.3333333333333333_1.6666666666666667'),
                #           # (2,"gaussian",'1.3333333333333333_2'),
                #            (2,"gaussian",'2.3333333333333335_1.6666666666666667'),
                #           # (2,"gaussian",'1_2'),
                #           # (2,"gaussian",'2.3333333333333335_0.6666666666666666'),
                #           (2,"gaussian",'2.3333333333333335_2.3333333333333335'),
                #           # (2,"gaussian",'2_1'),
                #           (2,"exponential",'0.6666666666666666_0.6666666666666666'),
                #           # (2,"exponential",'1.6666666666666667_1.3333333333333333'),
                #           (2,"exponential",'2.3333333333333335_1.6666666666666667'),
                #           # (2,"exponential",'0.6666666666666666_2'),
                #           (2,"exponential",'1_0.6666666666666666'),
                #           # (2,"exponential",'1_2.3333333333333335')
                #           (2,"exponential",'1_2')
                  ]
        else:
            re_do_list = []
    else:
        if one_border:
            #one border sequential
            re_do_list = [
                # #(1,"exponential",'1.3333333333333333_0.6666666666666666'),
                #   # (1,"exponential",'1.3333333333333333_1'),
                #   # (1,"exponential",'1.3333333333333333_1.3333333333333333'),
                #   # (1,"exponential",'1.3333333333333333_2'),
                #   # (1,"exponential",'1.6666666666666667_0.6666666666666666'),
                #   # (1,"exponential",'1.6666666666666667_1.3333333333333333'),
                #   # (1,"exponential",'1.6666666666666667_2'),
                #   # (1,"exponential",'1_0.6666666666666666'),
                #   # (1,"exponential",'1_1.6666666666666667'),
                #   (1,"exponential",'1_2'),
                #   # (1,"exponential",'1_2.3333333333333335'),
                #   # (1,"exponential",'2.3333333333333335_1.6666666666666667'),
                #   # (1,"exponential",'2_0.6666666666666666'),
                #   # (2,"gaussian",'0.6666666666666666_1'),
                #   (2,"gaussian",'1_0.6666666666666666'),
                #   # (2,"gaussian",'1.6666666666666667_2'),
                #   # (2,"gaussian",'2_1.6666666666666667'),
                #   # (2,"exponential",'0.6666666666666666_2.3333333333333335'),
                #   # (2,"exponential",'1.3333333333333333_0.6666666666666666'),
                #   # (2,"exponential",'1.3333333333333333_1'),
                #   # (2,"exponential",'1.3333333333333333_1.3333333333333333'),
                #   # (2,"exponential",'1.6666666666666667_0.6666666666666666'),
                #   (2,"exponential",'0.6666666666666666_2'),
                #   # (2,"exponential",'1_1.3333333333333333'),
                #   # (2,"exponential",'1_2')
                  ]
        elif snake_opt_path:
            re_do_list = [
                # (1,"exponential",'0.6666666666666666_2'),
                #           (1,"exponential",'1.3333333333333333_0.6666666666666666'),
                #           (1,"exponential",'1_0.6666666666666666'),
                #           (1,"exponential",'1_2'),
                #           (1,"exponential",'1_2.3333333333333335'),
                #           (1,"exponential",'2.3333333333333335_1.3333333333333333'),
                #           (2,"gaussian",'0.6666666666666666_0.6666666666666666'),
                #           (2,"gaussian",'0.6666666666666666_2.3333333333333335'),
                #           (2,"gaussian",'1_1.6666666666666667'),
                #           (2,"gaussian",'2.3333333333333335_1.6666666666666667'),
                #           (2,"gaussian",'2.3333333333333335_2.3333333333333335'),
                #           (2,"exponential",'1_1.6666666666666667'),
                #           (2,"exponential",'2.3333333333333335_1.6666666666666667')
                ]
        elif obstacle_slower:
            re_do_list = [
                # (1,"exponential",'0.6666666666666666_0.6666666666666666'),
                #           (1,"exponential",'0.6666666666666666_2.3333333333333335'),
                #           (1,"exponential",'2.3333333333333335_1.6666666666666667'),
                #           (1,"exponential",'2.3333333333333335_2.3333333333333335'),
                #           # (1,"exponential",'1.3333333333333333_0.6666666666666666'),
                #           # (1,"exponential",'1.6666666666666667_0.6666666666666666'),
                #           # (1,"exponential",'1.6666666666666667_1'),
                #           # (1,"exponential",'1.6666666666666667_1.3333333333333333'),
                #           # (1,"exponential",'1_1.6666666666666667'),
                #           # (2,"gaussian",'1.3333333333333333_0.6666666666666666'),
                #           # (2,"gaussian",'1.3333333333333333_1'),
                #           # (2,"gaussian",'1.3333333333333333_2'),
                #           # (2,"gaussian",'2.3333333333333335_2.3333333333333335'),
                #           # (2,"gaussian",'2_1.3333333333333333'),
                #           # (2,"gaussian",'2_2.3333333333333335'),
                #           (2,"exponential",'0.6666666666666666_2.3333333333333335'),
                #           (2,"exponential",'2.3333333333333335_1.6666666666666667'),
                #           # (2,"exponential",'0.6666666666666666_1.3333333333333333'),
                #           # (2,"exponential",'1.3333333333333333_0.6666666666666666'),
                #           # (2,"exponential",'1_1.3333333333333333'),
                #           # (2,"exponential",'1_2.3333333333333335')
                  ]
        elif borders_faster:
            re_do_list = [
                # (1,"exponential",'0.6666666666666666_2'),
                #           (1,"exponential",'0.6666666666666666_2.3333333333333335'),
                #           (1,"exponential",'1_0.6666666666666666'),
                #           # (1,"exponential",'1_2.3333333333333335'),
                #           # (1,"exponential",'2.3333333333333335_1.3333333333333333'),
                #           (2,"gaussian",'0.6666666666666666_2'),
                #           (2,"gaussian",'1.6666666666666667_1.3333333333333333'),
                #           (2,"gaussian",'2.3333333333333335_1.3333333333333333'),
                #           (2,"gaussian",'1.3333333333333333_0.6666666666666666'),
                #           # (2,"gaussian",'1_0.6666666666666666'),
                #           (2,"exponential",'1.3333333333333333_0.6666666666666666'),
                #           (2,"exponential",'1_1.6666666666666667'),
                #           (2,"exponential",'1_2.3333333333333335'),
                          ]
        else:
            re_do_list = []
    # re_do_list = []
    # G_1=get_graph_instance_1(N)
    G_2 = get_graph_instance_2(N)
    V,A,centers_ = get_nodes_and_edges(G_2,True)
    
    multipliers_theta = (___theta____#[f"{theta1}_{theta2}"
        #for theta1 in ___theta____ for theta2 in ___theta____
        # if (theta1,theta2)==(2/3,7/3) or
        # (theta1,theta2)==(5/3,4/3) or
        # (theta1,theta2)==(1,2/3) or
        # (theta1,theta2)==(1,5/3) or
        # (theta1,theta2)==(1,2) or
        # (theta1,theta2)==(2/3,2/3) or
        # (theta1,theta2)==(2/3,2) or
        # (theta1,theta2)==(1,7/3) or
        # (theta1,theta2)==(7/3,4/3) or
        # (theta1,theta2)==(7/3,5/3) or
        # (theta1,theta2)==(7/3,7/3) or
        # (theta1,theta2)==(4/3,2/3)] 
        if any([one_border,borders_faster,obstacle_slower,snake_opt_path]
                 ) else [1000.])
    
    for norm_ in norm_range_:
        Theta__ = {kernel_name : multipliers_theta 
                   for kernel_name in Kernels.keys()
            #"hyperbolic" : multipliers_theta, 
                #   "gaussian" : multipliers_theta, 
                # "exponential" : multipliers_theta
            }
        # def cardinal_sin(rho):
        #     return np.sin(rho)/rho
          
        for kernel_name,kernel in Kernels.items():
            if any([one_border,borders_faster,obstacle_slower,snake_opt_path]
                   ) and kernel_name=="hyperbolic": continue
            elif norm_==1 and kernel_name=="gaussian": continue
            elif norm_!=1 and kernel_name=="hyperbolic": continue
            elif not any([one_border,borders_faster,obstacle_slower,
                          snake_opt_path]) and kernel_name!="hyperbolic": continue
            # if norm_==2 and kernel_name=="exponential": continue
            registered_results={}
            registered_results["benchmark"] = []
            registered_results["distribution"] = []
            registered_results["theta"] = []
            registered_results["norm"] = []
            registered_results["n_sim"] = []
            registered_results["T"] = []
            registered_results["t"] = []
            registered_results["regret"] = []
            # registered_paths = {}
            # registered_paths["distribution"] = []
            # registered_paths["theta"] = []
            # registered_paths["norm"] = []
            # registered_paths["n_sim"] = []
            # registered_paths["T"] = []
            # registered_paths["t"] = []
            # registered_paths["index_arc_in_path"] = []
            # registered_paths["arc"] = []
            # registered_paths["observed_value"] = []
            
            
            # for i in range(len(A)):
            #     registered_results[f"mu_{A[i]}"] = []
    
            for __theta_p in Theta__[kernel_name]:
                # if (norm_,kernel_name,__theta_p)!=(1,"exponential",'1_2'):
                #     continue
                if (re_do_list != [] and 
                    not (norm_,kernel_name,__theta_p) in re_do_list): continue
                
                local_theta_neighborhood = [
                    __theta_p,
                    (f"{1/6. + float(__theta_p.split('_')[0])}"
                    +f"_{__theta_p.split('_')[1]}"),
                    (f"{-1/6. + float(__theta_p.split('_')[0])}"
                    +f"_{__theta_p.split('_')[1]}"),
                    (f"{__theta_p.split('_')[0]}"
                    +f"_{1/6. + float(__theta_p.split('_')[1])}"),
                    (f"{__theta_p.split('_')[0]}"
                    +f"_{-1/6. + float(__theta_p.split('_')[1])}")] if any(
                        [one_border,borders_faster,
                         obstacle_slower,snake_opt_path]
                        ) and False else [__theta_p]
                        
                iteration_experiment_set = [(local_point,it)  
                                            for it in range(nruns)
                                for local_point in local_theta_neighborhood]
                
                for __theta_,it in iteration_experiment_set:
                    
                    for T_iter in T_iter_set:
                    
                #T_iter=2
                    # N=3
                        
                        pp = False
                        # pp = it==0
                        # pp = True
                        # pp = it in [49,99,149,199]
                        
                        
                        
                        M=3#3 for regular chessboard
                        Delta= N / (2**M)
                        
                        _alpha__ = 11.
                        _beta__ = 13.
                        Plot_Instance_1 = False
                        Plot_Instance_2 = not Plot_Instance_1
                        
                        
                        
                        
                        max_x_y= (  max( [centers_[i][0] 
                                          for i in range(len(A))] ), 
                                  max( [centers_[i][1] 
                                        for i in range(len(A))] )) 
                        min_x_y= (  min( [centers_[i][0] 
                                          for i in range(len(A))] ), 
                                  min( [centers_[i][1] 
                                        for i in range(len(A))] )) 
                        
                        
                        
                        
                        v_real=50#km/h
                        dist__=1.5
                        d_all_arcs={a : dist__ for a in A}
                        
                        if pp:pp=PdfPages(os.path.join(output_dir,
                            f"test_N{N}_ntest{it}_T{T_iter}_norm{norm_}"+
                            f"{'_reversed' if _reversed__ else ''}"+
                            "_{kernel_name}"+
                            f"_theta{str(__theta_).replace('.','-')}"+
                            ("_two_border" if borders_faster else (
                            "_one_border" if one_border else (
                            "_with_obstacle" if obstacle_slower else (
                            "_chessboard" if not snake_opt_path else 
                            "_snake"))) 
                            )+".pdf"))
                        
                        if Plot_Instance_1:
                            edge_color_instance_1=['r' if int(
                                centers_[i][0]/Delta )%2==int( 
                                    centers_[i][1]/Delta )%2 else 'b' 
                                        for (i,a) in enumerate(A)]
                            # get_plot(G=G_2,edge_color_=edge_color_instance_1,pp=pp,plt_=plt,Title="Instance")
                            chessboard_truth= np.random.random(
                                len(A))*.1 + np.array([ np.log(v_real*(
                                .1 if int( centers_[i][0]/Delta )%2 == int( 
                                    centers_[i][1]/Delta )%2 
                                else .9
                                )) for (i,a) in enumerate(A) ])
                            instance__= chessboard_truth
                        
                        if Plot_Instance_2:
                            V_arcs_dict={}
                            # _edge_color___={}
                        #    Sampled_color_indices=[]
                            color_keys = list(mcolors.CSS4_COLORS.keys())
                        #    color_index = color_keys[np.random.randint(len(mcolors.CSS4_COLORS))]
                            if borders_faster:
                                for i in range(N+1):
                                    for i_a,a in enumerate(A):
                                        if (((centers_[i_a][0]==i*.5 
                                              or centers_[i_a][0]==(N-1- i*.5)) 
                                             and i*.5<=centers_[i_a][1]<=(
                                                 N-1- i*.5)) or 
                                            ((centers_[i_a][1]==i*.5 or 
                                              centers_[i_a][1]==(N-1- i*.5))  
                                             and i*.5<=centers_[i_a][0]<=(
                                                 N-1- i*.5))):
                                            V_arcs_dict[a]=np.log(v_real * (
                                                2/(1 + 2**.5 ))**(i-2) )
                                            # _edge_color___[a]=mcolors.CSS4_COLORS[color_keys[14+i]]#color_index]
                            elif one_border:
                                for i in range(2*N - 1):
                            #        while color_index in Sampled_color_indices or all(
                            #                [ coord_>=.8 for coord_ in 
                            #                 mcolors.to_rgba(mcolors.CSS4_COLORS[color_index]) ]):
                            #            color_index = color_keys[np.random.randint(len(mcolors.CSS4_COLORS))]
                            #        Sampled_color_indices.append(color_index)
                                    for i_a,a in enumerate(A):
                                        if (centers_[i_a][0]==i*.5 and 
                                            N-1-centers_[i_a][0]>=
                                            centers_[i_a][1]) or (
                                                centers_[i_a][1]==(2*N - 2 - i
                                                                   )*.5 
                                                and centers_[i_a][0]>=
                                                N-1-centers_[i_a][1]):
                                            V_arcs_dict[a]=np.log(
                                                v_real * (2/(1 + 2**.5 )
                                                          )**(i-2) )
                            elif obstacle_slower:
                                slow_gravity_center = ((N-1)/2. , 0)
                                orientation_effect = (10. , 3.5)
                                for i_a,a in enumerate(A):
                                #     if (centers_[i_a][0]==i*.5 and N-1-centers_[i_a][0]>=centers_[i_a][1]
                                #         ) or (
                                #             centers_[i_a][1]== (2*N - 2 - i)*.5 
                                #             and centers_[i_a][0]>=N-1-centers_[i_a][1]):
                                    V_arcs_dict[a] = np.log(
                                        v_real * ((1 + 2**.5 )/2)**(
                                            max(orientation_effect[0]*max(
                                                centers_[i_a][0]-
                                                slow_gravity_center[0],
                                                slow_gravity_center[0]-
                                                centers_[i_a][0]),
                                                orientation_effect[1]*max(
                                                centers_[i_a][1]-
                                                slow_gravity_center[1],
                                                slow_gravity_center[1]-
                                                centers_[i_a][1]) 
                                                )-60))
                                            # _edge_color___[a]=mcolors.CSS4_COLORS[color_keys[65+i]]#mcolors.CSS4_COLORS[color_index]
                            elif snake_opt_path:
                                _slow_gravity_center = [((N-1)/3. , (N-1) /3 ),
                                                       (2*(N-1)/3. , 2*(N-1)/3
                                                        ) ]
                                _orientation_effect = [(10. , 3.5),
                                                      (10. , 3.5)]
                                for i_a,a in enumerate(A):
                                #     if (centers_[i_a][0]==i*.5 and N-1-centers_[i_a][0]>=centers_[i_a][1]
                                #         ) or (
                                #             centers_[i_a][1]== (2*N - 2 - i)*.5 
                                #             and centers_[i_a][0]>=N-1-centers_[i_a][1]):
                                    V_arcs_dict[a] = min([ np.log(
                                        v_real * ((1 + 2**.5 )/2)**(
                                            max(orientation_effect[0]*max(
                                                centers_[i_a][0]-
                                                slow_gravity_center[0],
                                                slow_gravity_center[0]-
                                                centers_[i_a][0]),
                                                orientation_effect[1]*max(
                                                centers_[i_a][1]-
                                                slow_gravity_center[1],
                                                slow_gravity_center[1]-
                                                centers_[i_a][1]) 
                                                )-60)) for 
                                        slow_gravity_center,orientation_effect
                                        in zip(_slow_gravity_center,
                                               _orientation_effect)])
                            else:
                                for i_a,a in enumerate(A):
                                     V_arcs_dict[a] = np.log(v_real*(
                                         .1 if int( 
                                        centers_[i_a][0]/Delta )%2 == int( 
                                            centers_[i_a][1]/Delta )%2 
                                            else .9)) 
                            # _edge_color_instance_2=[_edge_color___[a] for a in A]
                            # get_plot(G=G_2,edge_color_=_edge_color_instance_2,pp=pp,plt_=plt,Title="Instance")
                            instance__= np.array( [V_arcs_dict[a] 
                                                   for a in A ] )
                    
                        get_plot(G=G_2,edge_color_=-instance__,
                                 pp=pp,plt_=plt,
                                 Title="Instance")    
                        # T_real=.35/v_real #hours
                        
                        # uncorrelated_truth=( np.random.random(len(A)) + 0.1 
                        #                     ) * np.log(v_real)
                        
                        #edge_truth= stats.norm(v_real* np.array(
                        #    [(.1 if (centers_[i][0]<=min_x_y[0] or centers_[i][1]>=max_x_y[1])
                        #                           else .9) for i in range(len(A))])).rvs()
                        
                        
                        
                        
                        sigma = np.ones(len(A)) * np.log(1.2) 
                        v_prior = v_real*.5
                        
                        prior_NormalIG = {}
                        prior_NormalIG["mu"] = np.ones(len(A)
                                                       ) * np.log(v_prior)
                        prior_NormalIG["kappa"] = np.ones(len(A)) 
                        prior_NormalIG["alpha"] = np.ones(len(A)
                                                          ) *_alpha__
                        prior_NormalIG["beta"] = np.ones(len(A)
                                                         ) * _beta__
                        
                        
                        prior_Spatial = {}
                        prior_Spatial["mu"] = np.ones(len(A)
                                                      ) * np.log(v_prior)
                        prior_Spatial["kappa"] = np.ones(len(A)) 
                        prior_Spatial["alpha"] = np.ones(len(A)
                                                         ) *_alpha__
                        prior_Spatial["beta"] = np.ones(len(A)
                                                        )  * _beta__
                        prior_Spatial["theta"] = general_multiplication(
                            __theta_,dist__)
                        prior_Spatial["rho"] = {}
                        prior_Spatial["phi"] = np.zeros(
                            (len(A), len(A)))#
                        # prior_Spatial["phi"]=dok_matrix((len(A), len(A)), dtype=np.float64)
                        prior_Spatial["name"] = kernel_name+"_"+str(
                            prior_Spatial["theta"])
                        prior_Spatial["reversed"] = _reversed__
                        # prior_Spatial["influential_arcs"]={}
                        # kernel = exponential
                        
                        
                        if not any([one_border,borders_faster,
                                    obstacle_slower,snake_opt_path]):
                            for i in range(len(A)):
                                for j in range(len(A)):
                                     if (int( 
                                        centers_[i][0]/Delta )%2 == int( 
                                            centers_[i][1]/Delta )%2 == int( 
                                        centers_[j][0]/Delta )%2 == int( 
                                            centers_[j][1]/Delta )%2) or (int( 
                                        centers_[i][0]/Delta )%2 != int( 
                                            centers_[i][1]/Delta )%2 and int( 
                                        centers_[j][0]/Delta )%2 != int( 
                                            centers_[j][1]/Delta )%2):
                                        #     print(centers_[i],centers_[j])
                                        rho_ij = np.linalg.norm(
                                            centers_[i]*dist__
                                            -centers_[j]*dist__,norm_)
                                        prior_Spatial["phi"][i,j] = (
                                            kernel(__theta_,rho_ij)
                                                    )
                        else:
                            for i in range(len(A)):
                                # prior["influential_arcs"][i]=[]
                                for j in range(len(A)):
                                    # if i!=j: continue
                            #        if not (A[i][1][0]==A[i][0][0]==A[j][1][0]==A[j][0][0] or A[i][1][1]==A[i][0][1]==A[j][1][1]==A[j][0][1]) and i!=j: continue
                            #        i==j or not (centers_[i][1]-centers_[j][1]==0 and ): continue
                                    rho_ij=np.linalg.norm(centers_[i]*
                                                          dist__-centers_[j]*
                                                          dist__,norm_)
                                    if type(prior_Spatial["theta"]) is str:
                                        if not (
                                                np.abs(A[i][1][0]-A[i][0][0])
                                                ==np.abs(A[j][1][0]-A[j][0][0]) 
                                                or np.abs(A[i][1][1]-
                                                          A[i][0][1])==np.abs(
                                                              A[j][1][1]-
                                                              A[j][0][1])
                                                ):
                                            prior_Spatial["phi"][i,j]= kernel(
                                                float(
                                                prior_Spatial["theta"
                                                        ].split("_")[1]),
                                                rho_ij)
                                        else:
                                            prior_Spatial["phi"][i,j] = kernel(
                                                float(
                                                prior_Spatial["theta"
                                                        ].split("_")[0]),
                                                rho_ij)
                                    else:
                                    # if kernel(prior_Spatial["theta"] , rho_ij ) >= np.exp(-21):
                                        prior_Spatial["phi"][i,j]= kernel(
                                            prior_Spatial["theta"],rho_ij)
                                        # prior["influential_arcs"][i].append(j)
                        # distrib_=
                        # min_eigen_v=np.linalg.eigvals( (prior_Spatial["phi"]*np.eye(len(A))) ).min()
                        # if -10**(-6)<min_eigen_v<0:
                        #     for i in range(len(A)):
                        #         prior_Spatial["phi"][i,i]-=1.01*min_eigen_v
                            
                        
                        non_corrated_dist = NormalIG(prior_NormalIG)
                        # spatial_distribution = Spatial_Metric(prior_Spatial)
                        spatial_distribution = Spatial_Metric_2(prior_Spatial)
                        # spatial_distribution = Spatial_Metric3(prior_Spatial)
                        plot_histogram_sigma(non_corrated_dist,
                            pp,plt,"Histogram_sigma2_before_training")
                        plot_histogram_sigma(spatial_distribution,
                            pp,plt,"Histogram_sigma2_before_training")
                        plot_distribution(non_corrated_dist,
                            pp,plt,"mu_density_before_training")
                        plot_distribution(spatial_distribution,
                            pp,plt,"mu_density_before_training")
                        (regret_,Paths_,Observed_values_,
                         Paths_expert,mu_updates,Time_
                         )=TS_Stationary(T_iter,G_2,instance__,sigma,
                                         d_all_arcs,
                                         [non_corrated_dist,
                                          spatial_distribution],
                                         borders_faster,obstacle_slower)#
                        
                        
                        
                        regret_NormalIG = regret_[0]
                        regret_Spatial_ = regret_[1]
                        Paths_NormalIG = Paths_[0]
                        Paths_Spatial_ = Paths_[1]
                        Time_NormalIG = Time_[0]
                        Time_Spatial_ = Time_[1]
                        
                        for d_ in  range(len(regret_)):
                            # index_reg = 0
                            for name__,_pths in Paths_[d_]:
                                kernel_iteration,regret_printted = (
                                    name__.split(", "))
                                dist_name,t_ = kernel_iteration.split("=")
                                dist_name = dist_name.replace("_Iter","")
                                t_ = int(t_.split("/")[0])
                                registered_results["distribution"].append(
                                    dist_name)
                                registered_results["benchmark"].append(
                                    kernel_name)
                                registered_results["theta"].append(__theta_)
                                registered_results["norm"].append(norm_)
                                registered_results["n_sim"].append(it)
                                registered_results["T"].append(T_iter)
                                registered_results["t"].append(t_)
                                registered_results["regret"].append(
                                    regret_[d_][t_-1])
                                # for i in range(len(A)):
                                #     registered_results[f"mu_{A[i]}"].append(mu_updates[d_][index_reg][i])
                                # for i_a____,_a____ in enumerate(_pths):
                                    # registered_paths["distribution"].append(dist_name)
                                    # registered_paths["theta"].append(__theta_)
                                    # registered_paths["norm"].append(norm_)
                                    # registered_paths["n_sim"].append(it)
                                    # registered_paths["T"].append(T_iter)
                                    # registered_paths["t"].append(t_)
                                    # registered_paths["index_arc_in_path"].append(i_a____)
                                    # registered_paths["arc"].append(_a____)
                                    # registered_paths["observed_value"].append(Observed_values_[d_][index_reg][i_a____])
                                # index_reg += 1
    
                        # =TS_Stationary(T_iter,G_2,non_corrated_dist,instance__,sigma,d_all_arcs)
                        
                        plot_histogram_sigma(non_corrated_dist,
                                pp,plt,"Histogram_sigma2_after_training")
                        plot_histogram_sigma(spatial_distribution,
                                pp,plt,"Histogram_sigma2_after_training")
                        plot_distribution(non_corrated_dist,
                                pp,plt,"mu_density_after_training")
                        plot_distribution(spatial_distribution,
                                pp,plt,"mu_density_after_training")
                        
                        
                        if Plot_Instance_1:
                            for j,obs__ in enumerate(Paths_NormalIG):
                                edge_color_instance_2=[ 
                                    ((-mu_updates[0][j][i]#
                                      #   'r'  if int( centers_[i][0]/Delta )%2==
                                      # int( centers_[i][1]/Delta )%2 else 'b' 
                                      ) if A[i] not in obs__[1] else -N/2#'k'
                                     ) for i in range(len(A))]
                                get_plot(G=G_2,
                                         edge_color_=edge_color_instance_2,
                                         pp=pp,plt_=plt,Title=f"{obs__[0]}",
                                         _paths_=obs__[1])
                        
                        
                        if Plot_Instance_2:
                            for j,obs__ in enumerate(Paths_NormalIG):
                                edge_color_instance_2=[ 
                                    (-mu_updates[0][j][i]#'k'
                                     ) for i in range(len(A))]
                                get_plot(G=G_2,
                                         edge_color_=edge_color_instance_2,
                                         pp=pp,plt_=plt,Title=f"{obs__[0]}",
                                         _paths_=obs__[1])
                        
                        
                        if Plot_Instance_1:
                            for j,obs__ in enumerate(Paths_Spatial_):
                                edge_color_instance_2=[ 
                                    ((-mu_updates[1][j][i]#
                                        # 'r'  if int( centers_[i][0]/Delta )%2==
                                        #                     int( centers_[i][1]/Delta )%2 else 'b' 
                                        ) if A[i] not in obs__[1] else -N/2#'k'
                                     ) for i in range(len(A))]
                                get_plot(G=G_2,
                                         edge_color_=edge_color_instance_2,
                                         pp=pp,plt_=plt,Title=f"{obs__[0]}",
                                         _paths_=obs__[1])
                            edge_color_instance_2=[ ((-instance__[i]
                                                      #'r'  if int( centers_[i][0]/Delta )%2==
                                                            # int( centers_[i][1]/Delta )%2 else 'b' 
                                                ) if A[i] not in 
                                                Paths_expert[-1][1] else -N/2#'k'
                                                ) for i in range(len(A))]
                            get_plot(G=G_2,
                                     edge_color_=edge_color_instance_2,
                                     pp=pp,
                                     plt_=plt,Title=f"{Paths_expert[-1][0]}",
                                     _paths_=Paths_expert[-1][1]
                                     )
                        if Plot_Instance_2:
                            for j,obs__ in enumerate(Paths_Spatial_):
                                edge_color_instance_2=[ (-mu_updates[1][j][i]#_edge_color_instance_2[i] 
                                    # if A[i] not in obs__[1] else -N/2#'k'
                                    ) for i in range(len(A))]
                                get_plot(G=G_2,
                                         edge_color_=edge_color_instance_2,
                                         pp=pp,plt_=plt,Title=f"{obs__[0]}",
                                         _paths_=obs__[1])
                            edge_color_instance_2=[ (-instance__[i]#_edge_color_instance_2[i] 
                                    if A[i] not in Paths_expert[-1][1] 
                                    else -instance__[i]#'k'
                                    ) for i in range(len(A))]
                            get_plot(G=G_2,edge_color_=edge_color_instance_2,
                                     pp=pp,plt_=plt,
                                     Title=f"{Paths_expert[-1][0]}",
                                     _paths_=Paths_expert[-1][1])
                        plot_regret_t(prior_Spatial["name"],
                                      regret_NormalIG,
                                      regret_Spatial_,pp,plt)
                        if pp:
                            pp.close()
                if output_data_file:
                    (pd.DataFrame.from_dict(registered_results)).to_csv(
                        os.path.join(output_dir,
                        f"SIMULATION_Regret_{kernel_name}"+
                        f"_norm{norm_}"+
                        f"{'_reversed' if _reversed__ else ''}"+
                        ("_two_border" if borders_faster else (
                            "_one_border" if one_border else (
                            "_with_obstacle" if obstacle_slower else (
                            "_chessboard" if not snake_opt_path else 
                            "_snake"))) 
                            )+".csv"),index=False)
                if not any([one_border,borders_faster,
                            obstacle_slower,snake_opt_path]): break
            if not any([one_border,borders_faster,
                        obstacle_slower,snake_opt_path]): break
                    # (pd.DataFrame.from_dict(registered_paths)).to_csv(os.path.join(output_dir,f"SIMULATION_Observed_Paths_{kernel_name}_norm{norm_}{'_reversed' if _reversed__ else ''}"+
                        # ("_two_border" if borders_faster else ("_one_border" if not obstacle_slower else "_with_obstacle") )+".csv"),index=False)
        if not any([one_border,borders_faster,
                    obstacle_slower,snake_opt_path]): break
    # if not any([one_border,borders_faster,obstacle_slower,snake_opt_path]): break
# print("")
# #
#def long_time(n):
#	for i in range(n):
#		for j in range(100000):
#			i*j
#
#
#long_time(5)




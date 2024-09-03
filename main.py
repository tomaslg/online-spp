#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 16:02:14 2021

@author: tomas
"""

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats
from scipy.sparse import dok_matrix

from time import time
from _meta_util import timer_func
from distrib import (NormalIG,Spatial_Metric,Spatial_Metric_2,Spatial_Metric3,
                     naive_approach,np)
from TS import get_nodes_and_edges,TS_Stationary,nx,truth_sampler
    
from plot_util import (show_plot,plot_distribution,plot_histogram_sigma,get_plot,
                       get_graph_instance_1,get_graph_instance_2,plot_regret_t)


def hyperbolic(theta,rho):
    return 1/(1 + rho/theta)

def gaussian(theta,rho):
    return np.exp(- (rho/theta)**2)

def exponential(theta,rho):
    return np.exp(-rho/theta)

def general_multiplication(__theta_,dist__):
    if type(__theta_) is str:
        __theta_ = [float(i) for i in __theta_.split("_")]
    if type(__theta_) is list:
        st_=""
        for i in range(len(__theta_)):
            __theta_[i] = __theta_[i]*dist__
            st_ += ("_" if i>0 else "") + str(__theta_[i]) 
        return st_
    else:
        return __theta_*dist__

def sample_truth_from_mean(truth,sigma,A):
    return {a : np.exp(truth[j_] + sigma[j_]**2 / 2)  
                                        for j_,a in enumerate(A)}

def artificial_instance_sampler(*args):
    truth , sigma , borders_faster , obstacle_slower = args
    def __function(V):
        source = int( (len(V)**.5 -1 )*len(V)**.5) if obstacle_slower else (
                int( len(V)/3 ) if borders_faster else int(len(V)-1))#target=np.random.randint(len(V))
        target = int( len(V)*2/3 ) if borders_faster else 0
        source,target=V[source],V[target]
        return source,target,np.log(np.random.lognormal(truth , sigma))
    return __function


def get_next_type_arg(next_index, _type=int):
    try:
        if _type is int: return int(sys.argv[next_index])
    except IndexError:
        return 0

pp_bol = False

T_iter_set_size = get_next_type_arg(1)
if T_iter_set_size==0:
    T_iter_set_size = (50, 51)
else:
    T_iter_set_size = (T_iter_set_size, get_next_type_arg(2))
if T_iter_set_size[1]==0:
    T_iter_set_size = (T_iter_set_size[0], T_iter_set_size[0] + 1)

index_instance = get_next_type_arg(3)
index_theta = get_next_type_arg(4)
index_norm_range = get_next_type_arg(5)

T_iter_set = range(T_iter_set_size[0], T_iter_set_size[1])
N = 20
nruns = get_next_type_arg(6)
v_input = get_next_type_arg(7)


_alpha__ = 1
_beta__ = 3

np.random.seed(10)
nruns = nruns if nruns > 0 else 1
output_dir=os.path.join(os.path.dirname(os.getcwd()),
    f"output_beta{_beta__}_alpha{_alpha__}_{nruns}_{T_iter_set_size[0]}"#"output_beta3_alpha1"
    )
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
output_data_file = True
all_instances_ = [#reversed,one-border,start-get-out-of-the-slow-zone-and-come-back,
                  #avoid-an-obstacle,snake##,chessboard-non-smooth-no-theta
    # [True,True,False,False,False],#37(1)
    # [True,False,False,False,True],#37(1)
    # [True,False,False,True,False],#37(2)
    # [True,False,True,False,False],#37(2) 
    # [True,False,False,False,False],#37(1)
    [False,True,False,False,False],#11(1) one border
    [False,False,False,False,True],#11(2) snake
    [False,False,False,True,False],#11(3) obstacle
    [False,False,True,False,False],#11(4) start slow zone leave and enter
    # [False,False,False,False,False]#11(1)
    ]
___theta____ = [2/3, 
                4/3,
                2,
                8/3,
                10/3, 
                4
    ]
if index_instance>0:
    all_instances_ = [all_instances_[index_instance - 1]]
if index_theta>0:
    ___theta____ = [___theta____[index_theta - 1]]

norm_range_ = range(1,3)
Kernels={#"hyperbolic" : hyperbolic , 
    "gaussian" : gaussian , 
    "exponential" : exponential
     }
if 1<=index_norm_range<=2:#1-norm only exponential#2-norm exponential and Gaussian
    norm_range_ = [norm_range_[index_norm_range - 1]]
elif index_norm_range==3:#2-norm Gaussian only
    Kernels = {"gaussian" : gaussian}
elif index_norm_range==4:#2-norm exponential only
    norm_range_ = range(2,3)
    Kernels = {"exponential" : exponential}

for (_reversed__,one_border,borders_faster,obstacle_slower,snake_opt_path
     ) in all_instances_:
    G_2 = get_graph_instance_2(N)
    V,A,centers_ = get_nodes_and_edges(G_2,True)
    
    multipliers_theta = (___theta____
        if any([one_border,borders_faster,obstacle_slower,snake_opt_path]
                 ) else [1000.])
    
    for norm_ in norm_range_:
        Theta__ = {kernel_name : multipliers_theta 
                   for kernel_name in Kernels.keys()}          
        for kernel_name,kernel in Kernels.items():
            if any([one_border,borders_faster,obstacle_slower,snake_opt_path]
                   ) and kernel_name=="hyperbolic": continue
            elif norm_==1 and kernel_name=="gaussian": continue
            elif norm_!=1 and kernel_name=="hyperbolic": continue
            elif not any([one_border,borders_faster,obstacle_slower,
                          snake_opt_path]) and kernel_name!="hyperbolic": continue
            registered_results={}
            registered_results["benchmark"] = []
            registered_results["distribution"] = []
            registered_results["theta"] = []
            registered_results["norm"] = []
            registered_results["n_sim"] = []
            registered_results["T"] = []
            registered_results["t"] = []
            registered_results["regret"] = []
    
            for __theta_p in Theta__[kernel_name]:
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
                        pp = pp_bol
                        M=3#3 for regular chessboard
                        Delta= N / (2**M)
                        
                        _kappa__= (_beta__ / (
                            _alpha__ - 1) / (4 if one_border else (
                            25 if obstacle_slower else (
                                16 if snake_opt_path else 1)))
                                ) if _alpha__>1 else 1
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
                        for a in A:
                            d_all_arcs[a[1],a[0]] = dist__
                        
                        if pp:pp=PdfPages(os.path.join(output_dir,
                            f"test_N{N}_ntest{it}_T{T_iter}_norm{norm_}"+
                            f"{'_reversed' if _reversed__ else ''}"+
                            f"_{kernel_name}"+
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
                            chessboard_truth= np.random.random(
                                len(A))*.1 + np.array([ np.log(v_real*(
                                .1 if int( centers_[i][0]/Delta )%2 == int( 
                                    centers_[i][1]/Delta )%2 
                                else .9
                                )) for (i,a) in enumerate(A) ])
                            instance__= chessboard_truth
                        
                        if Plot_Instance_2:
                            V_arcs_dict={}
                            color_keys = list(mcolors.CSS4_COLORS.keys())
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
                            elif one_border:
                                for i in range(2*N - 1):
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
                            elif snake_opt_path:
                                _slow_gravity_center = [((N-1)/3. , (N-1) /3 ),
                                                       (2*(N-1)/3. , 2*(N-1)/3
                                                        ) ]
                                _orientation_effect = [(10. , 3.5),
                                                      (10. , 3.5)]
                                for i_a,a in enumerate(A):
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
                            instance__= np.array( [V_arcs_dict[a] 
                                                   for a in A ] )
                    
                        get_plot(G=G_2,edge_color_=instance__,
                                 pp=pp,plt_=plt,
                                 Title="")
                        
                        
                        
                        
                        sigma = np.ones(len(A)) * np.log(1.2) 
                        v_prior = v_input if v_input>0 else np.exp(
                            2 if one_border else (
                                3 if obstacle_slower else (
                                    -2 if snake_opt_path else (
                                        4 if borders_faster else
                                        np.log(v_real)))))
                        
                        prior_NormalIG = {}
                        prior_NormalIG["mu"] = np.ones(len(A)
                                                       ) * np.log(v_prior)
                        prior_NormalIG["kappa"] = _kappa__ * np.ones(len(A)) 
                        prior_NormalIG["alpha"] = np.ones(len(A)
                                                          ) *_alpha__
                        prior_NormalIG["beta"] = np.ones(len(A)
                                                         ) * _beta__
                        
                        prior_Naive = {}
                        prior_Naive["mu"] = np.ones(len(A)) * np.log(v_prior)
                        prior_Naive["kappa"] = np.ones(len(A)) 
                        prior_Naive["alpha"] = np.ones(len(A)) *_alpha__
                        prior_Naive["beta"] = np.ones(len(A)) * _beta__
                        prior_Naive["name"] = "Naive"
                        prior_Naive["epsilon"] = .1
                        
                        PB_prior_Spatial = {}
                        PB_prior_Spatial["mu"] = np.ones(len(A)
                                                      ) * np.log(v_prior)
                        PB_prior_Spatial["kappa"] = np.ones(len(A)) 
                        PB_prior_Spatial["alpha"] = np.ones(len(A)
                                                         ) *_alpha__
                        PB_prior_Spatial["beta"] = np.ones(len(A)
                                                        )  * _beta__
                        PB_prior_Spatial["theta"] = general_multiplication(
                            __theta_,dist__)
                        PB_prior_Spatial["rho"] = {}
                        PB_prior_Spatial["phi"] = np.zeros(
                            (len(A), len(A)))
                        PB_prior_Spatial["name"] = "PB_"+kernel_name+"_"+str(
                            np.round(PB_prior_Spatial["theta"]/dist__,1))
                        PB_prior_Spatial["reversed"] = _reversed__
                        PB_prior_Spatial["generate_cov_matrix"] = False
                        
                        PA_prior_Spatial = {}
                        PA_prior_Spatial["mu"] = np.ones(len(A)
                                                      ) * np.log(v_prior)
                        PA_prior_Spatial["kappa"] = np.ones(len(A)) 
                        PA_prior_Spatial["alpha"] = np.ones(len(A)
                                                         ) *_alpha__
                        PA_prior_Spatial["beta"] = np.ones(len(A)
                                                        )  * _beta__
                        PA_prior_Spatial["theta"] = general_multiplication(
                            __theta_,dist__)
                        PA_prior_Spatial["rho"] = {}
                        PA_prior_Spatial["phi"] = np.zeros(
                            (len(A), len(A)))
                        PA_prior_Spatial["name"] = "PA_"+kernel_name+"_"+str(
                            np.round(PA_prior_Spatial["theta"]/dist__,1))
                        PA_prior_Spatial["reversed"] = _reversed__
                        PA_prior_Spatial["generate_cov_matrix"] = False
                        
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
                                        rho_ij = np.linalg.norm(
                                            centers_[i]*dist__
                                            -centers_[j]*dist__,norm_)
                                        PB_prior_Spatial["phi"][i,j] = (
                                            kernel(__theta_,rho_ij)
                                                    )
                                        PA_prior_Spatial["phi"][i,j] = (
                                            kernel(__theta_,rho_ij)
                                                    )
                        else:
                            for i in range(len(A)):
                                for j in range(len(A)):
                                    rho_ij=np.linalg.norm(centers_[i]*
                                                          dist__-centers_[j]*
                                                          dist__,norm_)
                                    if type(PB_prior_Spatial["theta"]) is str:
                                        if not (
                                                np.abs(A[i][1][0]-A[i][0][0])
                                                ==np.abs(A[j][1][0]-A[j][0][0]) 
                                                or np.abs(A[i][1][1]-
                                                          A[i][0][1])==np.abs(
                                                              A[j][1][1]-
                                                              A[j][0][1])
                                                ):
                                            PB_prior_Spatial["phi"][i,j]= kernel(
                                                float(
                                                PB_prior_Spatial["theta"
                                                        ].split("_")[1]),
                                                rho_ij)
                                        else:
                                            PB_prior_Spatial["phi"][i,j] = kernel(
                                                float(
                                                PB_prior_Spatial["theta"
                                                        ].split("_")[0]),
                                                rho_ij)
                                    else:
                                        PB_prior_Spatial["phi"][i,j]= kernel(
                                            PB_prior_Spatial["theta"],rho_ij)
                                    if type(PA_prior_Spatial["theta"]) is str:
                                        if not (
                                                np.abs(A[i][1][0]-A[i][0][0])
                                                ==np.abs(A[j][1][0]-A[j][0][0]) 
                                                or np.abs(A[i][1][1]-
                                                          A[i][0][1])==np.abs(
                                                              A[j][1][1]-
                                                              A[j][0][1])
                                                ):
                                            PA_prior_Spatial["phi"][i,j]= kernel(
                                                float(
                                                PA_prior_Spatial["theta"
                                                        ].split("_")[1]),
                                                rho_ij)
                                        else:
                                            PA_prior_Spatial["phi"][i,j] = kernel(
                                                float(
                                                PA_prior_Spatial["theta"
                                                        ].split("_")[0]),
                                                rho_ij)
                                    else:
                                        PA_prior_Spatial["phi"][i,j]= kernel(
                                            PA_prior_Spatial["theta"],rho_ij)
                                        
                        
                        non_corrated_dist = NormalIG(prior_NormalIG)
                        naive_dist = naive_approach(prior_Naive)
                        PA_spatial_distribution = Spatial_Metric_2(PA_prior_Spatial)
                        PB_spatial_distribution = Spatial_Metric3(PB_prior_Spatial)
                        plot_histogram_sigma(non_corrated_dist,
                            pp,plt,"Histogram_sigma2_before_training")
                        plot_histogram_sigma(PA_spatial_distribution,
                            pp,plt,"Histogram_sigma2_before_training")
                        plot_histogram_sigma(PB_spatial_distribution,
                            pp,plt,"Histogram_sigma2_before_training")
                        plot_distribution(non_corrated_dist,
                            pp,plt,"mu_density_before_training")
                        plot_distribution(PA_spatial_distribution,
                            pp,plt,"mu_density_before_training")
                        plot_distribution(PB_spatial_distribution,
                            pp,plt,"mu_density_before_training")
                        (regret_,Paths_,Observed_values_,
                         Paths_expert,mu_updates,Time_,Exp_obj,delta_
                         ) = TS_Stationary(T_iter,G_2,
                             truth_sampler(artificial_instance_sampler(
                                 instance__ , sigma ,
                                 borders_faster , obstacle_slower),
                             V,sample_truth_from_mean(
                                 instance__,sigma,A)),
                             d_all_arcs,[non_corrated_dist,
                                          PB_spatial_distribution,
                                          PA_spatial_distribution,
                                          naive_dist])
                        
                        
                        
                        regret_NormalIG = regret_[0]
                        PB_regret_Spatial_ = regret_[1]
                        PA_regret_Spatial_ = regret_[2]
                        regret_naive = regret_[3]
                        Paths_NormalIG = Paths_[0]
                        PB_Paths_Spatial_ = Paths_[1]
                        PA_Paths_Spatial_ = Paths_[2]
                        Paths_naive = Paths_[3]
                        Time_NormalIG = Time_[0]
                        PB_Time_Spatial_ = Time_[1]
                        PA_Time_Spatial_ = Time_[2]
                        Time_naive = Time_[3]
                        
                        for d_ in  range(len(regret_)):
                            for name__,_pths in Paths_[d_]:
                                dist_name,t_ = name__.split("$")[
                                    0], int(name__.split("=")[-3].split(
                                        "$")[0])
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
                        
                        plot_histogram_sigma(non_corrated_dist,
                                pp,plt,"Histogram_sigma2_after_training")
                        plot_histogram_sigma(PB_spatial_distribution,
                                pp,plt,"Histogram_sigma2_after_training")
                        plot_histogram_sigma(PA_spatial_distribution,
                                pp,plt,"Histogram_sigma2_after_training")
                        plot_distribution(non_corrated_dist,
                                pp,plt,"mu_density_after_training")
                        plot_distribution(PB_spatial_distribution,
                                pp,plt,"mu_density_after_training")
                        plot_distribution(PA_spatial_distribution,
                                pp,plt,"mu_density_after_training")
                        
                        
                        if Plot_Instance_1:
                            for j___,Paths_Spatial_ in enumerate([Paths_NormalIG,
                                    PB_Paths_Spatial_, PA_Paths_Spatial_,
                                    Paths_naive]):
                                for j,obs__ in enumerate(Paths_Spatial_):
                                    get_plot(G=G_2,
                                             edge_color_=[ 
                                        (np.log(mu_updates[j___][j][i])
                                         ) for i in range(len(A))],
                                             pp=pp,plt_=plt,Title=f"{obs__[0]}",
                                             _paths_=obs__[1])
                            edge_color_instance_2=[ ((instance__[i])
                                                ) for i in range(len(A))]
                            get_plot(G=G_2,
                                     edge_color_=edge_color_instance_2,
                                     pp=pp,
                                     plt_=plt,Title=f"{Paths_expert[-1][0]}",
                                     _paths_=Paths_expert[-1][1]
                                     )
                        if Plot_Instance_2:
                            for j___,Paths_Spatial_ in enumerate([Paths_NormalIG,
                                    PB_Paths_Spatial_, PA_Paths_Spatial_,
                                    Paths_naive]):
                                for j,obs__ in enumerate(Paths_Spatial_):
                                    get_plot(G=G_2,
                                             edge_color_=[ (np.log(
                                        mu_updates[j___][j][i])
                                        ) for i in range(len(A))],
                                             pp=pp,plt_=plt,Title=f"{obs__[0]}",
                                             _paths_=obs__[1])
                            edge_color_instance_2=[ (instance__[i]
                                    ) for i in range(len(A))]
                            get_plot(G=G_2,edge_color_=edge_color_instance_2,
                                     pp=pp,plt_=plt,
                                     Title=f"{Paths_expert[-1][0]}",
                                     _paths_=Paths_expert[-1][1])
                        plot_regret_t(["Naive",
                                        PB_prior_Spatial["name"],
                                        PA_prior_Spatial["name"]
                                       ],
                                      regret_NormalIG,
                                      [regret_naive,
                                       PB_regret_Spatial_,
                                       PA_regret_Spatial_
                                       ],pp,plt)
                        if pp_bol:
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
        if not any([one_border,borders_faster,
                    obstacle_slower,snake_opt_path]): break

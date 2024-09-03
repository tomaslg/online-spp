#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 22:57:31 2021

@author: tomas
"""

import networkx as nx
import numpy as np
from scipy import stats
from time import time
from _meta_util import timer_func
from Dijkstra import heap_dijkstra_


def get_nodes_and_edges(G,get_centers=True):
    if isinstance(G,type(nx.Graph())):
        V,A=list(G.nodes()),list(G.edges())
    elif isinstance(G,type(nx.DiGraph())):
        V,A=list(G.nodes()),list(G.in_edges())+list(G.out_edges())
    if not get_centers:
        return V,A
    else:
        centers_=[]
        for i in range(len(A)):
            centers_.append( (np.array(A[i][0]) + np.array(A[i][1]) )/2)
        return V,A,centers_
    
# (,)

class truth_sampler():
    def __init__(self,sampler_function,args,mean):
        self.sampler_function = sampler_function
        self.args = args
        self.mean = mean

def scientific_notation(x):
    return '{:.2e}'.format(x)

def round_format(val):
    return (scientific_notation(val) if 0<abs(val)<.01 else np.round(
        val,2)) if 0<abs(val)<1. else (np.round(val,2) if 0<abs(val)<10. else (
            scientific_notation(val) if abs(val)>100 else (
                np.round(val,2) if val>0 else val)))

def format_distrib_name(distrib_name):
    st__ = ""
    distrib_name_list = distrib_name.split("_")
    if len(distrib_name_list)==1:
        st__ += f"{distrib_name}"
    elif len(distrib_name_list)==2:
        st__ += (distrib_name.split("_")[0] + r" $(\epsilon = " +
                 distrib_name.split("_")[1] + ")$")
    elif 3 <= len(distrib_name_list) <= 4:
        st__ += f"{distrib_name_list[-3]} {distrib_name_list[-2]} ".replace(
            "gaussian","Gaussian")
        st__ += r"$(\varphi = " + f"{distrib_name_list[-1]})$"
    return st__


@timer_func
def TS_Stationary(T_iter,G,truth_sampler,d_all_arcs,Distrib_):
    nsample = 1
    V,A = get_nodes_and_edges(G,False)
    map_index_A = {a : j_ for j_,a in enumerate(A)}
    if not G.is_directed():
        for j_,a in enumerate(A):
            map_index_A[a[1],a[0]] = j_
    regret=[]
    instant_regret=[]
    Paths=[]
    Observed_values = []
    Paths_expert=[]
    mu_updates=[]
    Time_=[]
    Exp_obj = []
    for d_,distrib_ in enumerate(Distrib_):
        regret.append([])
        instant_regret.append([])
        Paths.append([])
        Observed_values.append([])
        Time_.append([])
        mu_updates.append([])
    delta___0=[0. for d_ in range(len(Distrib_))]
    for t in range(T_iter):
        bol_ = True
        last_source = last_target = last_dist_source_all = last_pred_ = None
        while bol_:
            (source , target , scen_log_speed
             ) = truth_sampler.sampler_function(V)
            
            if (last_source, last_target) == (source , target):
                dist_source_all = last_dist_source_all
                pred_ = last_pred_
            else:
                (dist_source_all,pred_
                 ) = heap_dijkstra_(G,{akey : d_all_arcs[akey]/val_ 
                     for akey,val_ in  truth_sampler.mean.items() },
                     source,True)
                last_source = source
                last_target = target
                last_dist_source_all = dist_source_all
                last_pred_ = pred_
            if dist_source_all[target]==np.float("inf"): continue
            exp_obj = dist_source_all[target]
            print("exp_obj =",exp_obj)
            Exp_obj.append(exp_obj)
            expert_sol = []
            aux = target
            while aux != source:
                if G.is_directed():
                    expert_sol.append((pred_[aux],aux))
                else:
                    expert_sol.append((pred_[aux],aux) if 
                        (pred_[aux],aux) in truth_sampler.mean.keys() 
                        else (aux,pred_[aux]) )
                aux = pred_[aux]
            bol_ = False
        for d_,distrib_ in enumerate(Distrib_):
            t1_=time()
            posterior_sample_log_speed = distrib_.sample_poserior() 
            (dist_source_all,pred_
             ) = heap_dijkstra_(G,{a : d_all_arcs[a]/np.exp(
                posterior_sample_log_speed[map_index_A[a]])
                                 for a in A},source,True)
            obj = dist_source_all[target]
            sol = []
            aux = target
            while aux != source:
                if G.is_directed():
                    sol.append((pred_[aux],aux))
                else:
                    sol.append((pred_[aux],aux) if 
                        (pred_[aux],aux) in truth_sampler.mean.keys() 
                        else (aux,pred_[aux]) )
                aux = pred_[aux]
            if "Naive" in distrib_.name and all([distrib_.N[map_index_A[
                    a]]>0 for a in sol]) and distrib_.epsilon>=np.random.random():
                n_tabu = 15
                try:
                    tabu_arc = sol[np.random.randint(min(n_tabu,len(sol),
                            len(sol)-n_tabu), max(
                                len(sol)-n_tabu,n_tabu) )]
                except IndexError:
                    n_tabu = 3
                    tabu_arc = sol[np.random.randint(min(n_tabu,len(sol),
                            len(sol)-n_tabu), max(
                                len(sol)-n_tabu,n_tabu) )]
                posterior_sample_log_speed[map_index_A[tabu_arc]] = - abs(
                    min(posterior_sample_log_speed)) * 6
                (dist_source_all,pred_
                 ) = heap_dijkstra_(G,{a : d_all_arcs[a]/np.exp(
                    posterior_sample_log_speed[map_index_A[a]])
                                     for a in A},source,True)
                obj = dist_source_all[target]
                sol = []
                aux = target
                while aux != source:
                    if G.is_directed():
                        sol.append((pred_[aux],aux))
                    else:
                        sol.append((pred_[aux],aux) if 
                            (pred_[aux],aux) in truth_sampler.mean.keys() 
                            else (aux,pred_[aux]) )
                    aux = pred_[aux]
                
            cnt_ = 1
            posterior_sample_log_speed = np.exp(
                posterior_sample_log_speed) / nsample
            while cnt_<nsample:
                posterior_sample_log_speed += (
                    np.exp(distrib_.sample_poserior()) / nsample)
                cnt_ += 1
            
            distrib_.update_posterior_one_observation_stationary(
                { map_index_A[a] : scen_log_speed[map_index_A[a]] 
                 for a in sol })
            Time_[d_].append(time()-t1_)
            delta___ = (sum([ d_all_arcs[a]/truth_sampler.mean[a]
                for a in sol]) - sum(
                    [ d_all_arcs[a]/truth_sampler.mean[a]
                        for a in expert_sol]))
            regret[d_].append(( delta___ + 
                          (regret[d_][-1]*t if t>0 else 0)   )/ (t+1) )
            instant_regret[d_].append( delta___  )
            if  True:
                st = (format_distrib_name(distrib_.name) + f" $t = {t+1}$ " + 
                       r"$\frac{\mathcal{R}_t}{t z^*} = " +
                       f"{round_format(regret[d_][-1]/exp_obj)}$ " +
                       r"$\frac{\Delta_t}{ z^*} = " + f"{round_format(delta___/exp_obj)}$")
                print(st)
                Paths[d_].append((st,sol))
                Observed_values[d_].append([ scen_log_speed[map_index_A[a]] for a in sol])
                mu_updates[d_].append(np.copy(posterior_sample_log_speed))
            delta___0[d_]=delta___
        Paths_expert.append(("",expert_sol))
    return regret,Paths,Observed_values,Paths_expert,mu_updates,Time_,Exp_obj,instant_regret




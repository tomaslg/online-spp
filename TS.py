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
from solve_sp import Solve_SP_LP#,gp,GRB,


# def empty_function(*args):
#     pass

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

@timer_func
def TS_Stationary(T_iter,G,truth,sigma,d_all_arcs,Distrib_,
                  borders_faster,obstacle_slower):
    nsample = 1
    V,A=get_nodes_and_edges(G,False)
    map_index_A = {a : j_ for j_,a in enumerate(A)}
    regret=[]
    Paths=[]
    Observed_values = []
    Paths_expert=[]
    mu_updates=[]
    Time_=[]
    # score_iter_=[0 for tm_snp in curve_time_snapshots]
    for d_,distrib_ in enumerate(Distrib_):
        regret.append([])
        Paths.append([])
        Observed_values.append([])
        Paths_expert.append([])
        Time_.append([])
        mu_updates.append([])
    delta___0=[0. for d_ in range(len(Distrib_))]
    for t in range(T_iter):
        source = int( (len(V)**.5 -1 )*len(V)**.5) if obstacle_slower else (
            int( len(V)/3 ) if borders_faster else int(len(V)-1))#target=np.random.randint(len(V))
        # while source==target:
        #     target=np.random.randint(len(V))
        target = int( len(V)*2/3 ) if borders_faster else 0
        source,target=V[source],V[target]
        exp_obj,expert_sol=Solve_SP_LP(G,V,A,
                                       {a : np.exp(truth[j_] + sigma[j_]**2 / 2)  for j_,a in enumerate(A)},d_all_arcs,
                                       source,target)
        
        scen_log_speed = np.log(np.random.lognormal(truth,sigma))
        for d_,distrib_ in enumerate(Distrib_):
            t1_=time()
            posterior_sample_log_speed = distrib_.sample_poserior() 
            if any([np.exp(posterior_sample_log_speed[j_]) <.1**4 
                    for j_,a in enumerate(A)]):
                   print("here")
            try: 
                obj,sol=Solve_SP_LP(G,V,A,
                        {a : np.exp(posterior_sample_log_speed[j_])#,dtype="float128")
                                     for j_,a in enumerate(A)},
                        d_all_arcs,
                        source,target)#,timelimit=3600.0,debug=1,debug_sp=True)
            except TypeError:
                print("here")
                
            cnt_ = 1
            posterior_sample_log_speed = posterior_sample_log_speed / nsample
            while cnt_<nsample:
                posterior_sample_log_speed += distrib_.sample_poserior() / nsample
                cnt_ += 1
            
            # stats.norm(truth,sigma).rvs()
            distrib_.update_posterior_one_observation_stationary({ map_index_A[a] : scen_log_speed[map_index_A[a]] for a in sol })
            Time_[d_].append(time()-t1_)
            # regret[d_].append( sum([ d_all_arcs[a]/scen_log_speed[map_index_A[a]] for a in sol]) - 
            #               sum([ d_all_arcs[a]/scen_log_speed[map_index_A[a]] for a in expert_sol])  )
            delta___=sum([ d_all_arcs[a]/np.exp(truth[map_index_A[a]]) for a in sol]) - sum([ d_all_arcs[a]/np.exp(truth[map_index_A[a]]) for a in expert_sol])
            regret[d_].append(( delta___ + 
                          (regret[d_][-1]*t if t>0 else 0)   )/ (t+1) )
            if (sum([ d_all_arcs[a]/np.exp(truth[map_index_A[a]]) 
                    for a in sol]) - sum([ d_all_arcs[a]/np.exp(truth[map_index_A[a]]) 
                                          for a in expert_sol])<0):
                print("here")
            assert(delta___ >=0)
            if  True:#delta___ > delta___0[d_]*2. or delta___<=.1**4 or t+1==T_iter:#t>0 and (( (t+1)/T_iter) * 4  ) % 1 < 1/T_iter:
                st=f"{distrib_.name}_Iter={t+1}/{T_iter}, (regret[d_])={(regret[d_][-1])}"
                print(st)
                Paths[d_].append((st,sol))
                Observed_values[d_].append([ scen_log_speed[map_index_A[a]] for a in sol])
                mu_updates[d_].append(np.copy(posterior_sample_log_speed#distrib_.mu#
                                      ))
            delta___0[d_]=delta___
        if  t+1==T_iter:
            Paths_expert.append(("Expert",expert_sol))
        # if (t+1) in curve_time_snapshots:
        #     for d_,distrib_ in enumerate(Distrib_):
        #         if regret[d_][-1]==np.min([regret[d____][-1] for d____ in range(len(Distrib_))]):
        #             score_iter_[d_]+=1
    return regret,Paths,Observed_values,Paths_expert,mu_updates,Time_#,score_iter_




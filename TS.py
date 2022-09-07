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
# from solve_sp import Solve_SP_LP#,gp,GRB,
from Dijkstra import heap_dijkstra_

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

class truth_sampler():
    def __init__(self,sampler_function,args,mean):
        self.sampler_function = sampler_function
        self.args = args
        self.mean = mean
    # def sample(self):
    #     return self.sampler_function(self.args)




@timer_func
def TS_Stationary(T_iter,G,truth_sampler,d_all_arcs,Distrib_#,borders_faster,obstacle_slower
                  ):
    nsample = 20
    # limit_in_seconds = 1000 # 0 to shuft off
    V,A = get_nodes_and_edges(G,False)
    map_index_A = {a : j_ for j_,a in enumerate(A)}
    if not G.is_directed():
        for j_,a in enumerate(A):
            map_index_A[a[1],a[0]] = j_
    regret=[]
    Paths=[]
    Observed_values = []
    Paths_expert=[]
    mu_updates=[]
    Time_=[]
    Exp_obj = []
    # score_iter_=[0 for tm_snp in curve_time_snapshots]
    for d_,distrib_ in enumerate(Distrib_):
        regret.append([])
        Paths.append([])
        Observed_values.append([])
        # Paths_expert.append([])
        Time_.append([])
        mu_updates.append([])
    delta___0=[0. for d_ in range(len(Distrib_))]
    for t in range(T_iter):
        # source = int( (len(V)**.5 -1 )*len(V)**.5) if obstacle_slower else (
        #     int( len(V)/3 ) if borders_faster else int(len(V)-1))#target=np.random.randint(len(V))
        # target = int( len(V)*2/3 ) if borders_faster else 0
        # source = target = np.random.randint(len(V))
        # while source==target:
        #     target = np.random.randint(len(V))
        # source,target = V[source],V[target]
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
            # if limit_in_seconds and exp_obj<=limit_in_seconds: continue
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
                # Solve_SP_LP(G,V,A,
                #         {a : np.exp(posterior_sample_log_speed[j_])#,dtype="float128")
                #                      for j_,a in enumerate(A)},
                #         d_all_arcs,
                #         source,target)#,timelimit=3600.0,debug=1,debug_sp=True)
            if distrib_.name=="Naive" and all([distrib_.N[map_index_A[
                    a]]>0 for a in sol]) and distrib_.epsilon>=np.random.random():
                n_tabu = 15
                tabu_arc = sol[np.random.randint(min(n_tabu,len(sol)), max(
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
            # (np.array([d_all_arcs[a]
            #     for a in A]) / )
            while cnt_<nsample:
                posterior_sample_log_speed += (
                    # np.array([d_all_arcs[a] for a in A]) / 
                    np.exp(distrib_.sample_poserior()) / nsample)
                cnt_ += 1
            
            # stats.norm(truth,sigma).rvs()
            distrib_.update_posterior_one_observation_stationary(
                { map_index_A[a] : scen_log_speed[map_index_A[a]] 
                 for a in sol })
            Time_[d_].append(time()-t1_)
            # regret[d_].append( sum([ d_all_arcs[a]/scen_log_speed[map_index_A[a]] for a in sol]) - 
            #               sum([ d_all_arcs[a]/scen_log_speed[map_index_A[a]] for a in expert_sol])  )
            delta___ = sum([ d_all_arcs[a]/#np.exp(
                truth_sampler.mean[#map_index_A[
                    a]#])
                for a in sol]) - sum(
                    [ d_all_arcs[a]/#np.exp(
                        truth_sampler.mean[#map_index_A[
                            a]#]) 
                        for a in expert_sol])
            regret[d_].append(( delta___ + 
                          (regret[d_][-1]*t if t>0 else 0)   )/ (t+1) )
            assert(delta___ >-.1**3 )
            if  True:#delta___ > delta___0[d_]*2. or delta___<=.1**4 or t+1==T_iter:#t>0 and (( (t+1)/T_iter) * 4  ) % 1 < 1/T_iter:
                st=(f"{distrib_.name}_Iter={t+1}/{T_iter}," +
                    f" Acc. Pseudo-Regret={np.round(regret[d_][-1],2)}")
                print(st + f", Marginal Pseudo-Regret={np.round(delta___,2)}")
                Paths[d_].append((st,sol))
                Observed_values[d_].append([ scen_log_speed[map_index_A[a]] for a in sol])
                mu_updates[d_].append(np.copy(posterior_sample_log_speed#distrib_.mu#
                                      ))
            delta___0[d_]=delta___
        # if  t+1==T_iter:
        Paths_expert.append(("Expert",expert_sol))
        # if (t+1) in curve_time_snapshots:
        #     for d_,distrib_ in enumerate(Distrib_):
        #         if regret[d_][-1]==np.min([regret[d____][-1] for d____ in range(len(Distrib_))]):
        #             score_iter_[d_]+=1
    return regret,Paths,Observed_values,Paths_expert,mu_updates,Time_,Exp_obj#,score_iter_




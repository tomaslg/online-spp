# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 23:40:04 2022

@author: tol28
"""

import os
import sys
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.sparse import dok_matrix
import pandas as pd

from distrib import (NormalIG,Spatial_Metric,Spatial_Metric_2,Spatial_Metric3,
                     naive_approach,np)
from TS import TS_Stationary,nx,truth_sampler

from Dijkstra import heap_dijkstra_


def haversine(lon1, lat1, lon2, lat2):
    #Convert decimal degrees to Radians:
    lon1 = np.radians(lon1)
    lat1 = np.radians(lat1)
    lon2 = np.radians(lon2)
    lat2 = np.radians(lat2)

    #Implementing Haversine Formula: 
    dlon = np.subtract(lon2, lon1)
    dlat = np.subtract(lat2, lat1)

    a = np.add(np.power(np.sin(np.divide(dlat, 2)), 2), np.multiply(
        np.cos(lat1), np.multiply(np.cos(lat2), np.power(
            np.sin(np.divide(dlon, 2)), 2))))
    c = np.multiply(2, np.arcsin(np.sqrt(a)))
    r = 6371

    return c*r

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
def sorted_tup(arc_):
    return arc_ if arc_[0]<=arc_[1] else (arc_[1],arc_[0])

def read_from_arcs_nodes(_paths_dir , arcs_paths_fl , nodes_fl,__theta_,
                         run_spatial):
    fl = open(os.path.join(_paths_dir , nodes_fl) , "r")
    lines_array_ = fl.readlines()
    i_ = 0
    header = lines_array_[i_].split(";")
    col_index = {}
    for hd in header:
        col_index[hd.replace("\n","")] = i_
        i_ += 1
    i_ = 1
    V = []
    nodes = {}
    while i_ < len(lines_array_)-1:
        line_atr = lines_array_[i_].replace("\n","").split(";")
        i_ += 1
        latlon___ = np.array([float(line_atr[col_index["latitude"]]),
            float(line_atr[col_index["longitude"]])])
        V.append(int(line_atr[col_index["id"]]))
        nodes[V[-1]] = latlon___
    
    fl = open(os.path.join(_paths_dir , arcs_paths_fl) , "r")
    lines_array_ = fl.readlines()
    i_ = 0
    header = lines_array_[i_].split(";")
    col_index = {}
    for hd in header:
        col_index[hd.replace("\n","")] = i_
        i_ += 1
    i_ = 1
    A = []
    arcs = {}
    speeds_ = {}
    d_all_arcs = {}
    centers_arcs = {}
    max_speed = {'primary': 0.013952060298131714,
     'service': 0.016185045276776066,
     'secondary': 0.013142365972990917,
     'tertiary': 0.011991765896341483,
     'unclassified': 0.016586339453121653,
     'residential': 0.012556273419106957,
     'trunk': 0.028051301428762634,
     'secondary_link': 0.012621462349693325,
     'trunk_link': 0.021720378528622842,
     'primary_link': 0.016184201440056726,
     'motorway': 0.037081843831752594,
     'construction': 0.023001513879863184,
     'tertiary_link': 0.012394071883019542, 
     'cycleway': 0.011092236991396128, 
     'motorway_link': 0.032905218269646444, 
     'living_street': 0.014658973938660752, 
     'track': 0.028663446979922465,
     'proposed': 0.013689777028947817,
     'platform': 0.006422628455614414}
    for key_arc_type in max_speed:
        if max_speed[key_arc_type] > .025:
            max_speed[key_arc_type] = .025
        elif max_speed[key_arc_type]>.02:
            max_speed[key_arc_type] = .02
        elif max_speed[key_arc_type]>.015:
            max_speed[key_arc_type] = .015
        else:
            max_speed[key_arc_type] = .011
    while i_ < len(lines_array_):
        line_atr = lines_array_[i_].replace("\n","").split(";")
        i_ += 1
        arc_ = (int(line_atr[col_index["node1"]])
            , int(line_atr[col_index["node2"]]))
            
        if int(line_atr[col_index["id"]]) not in A:
            try:
                lat1, lon1 = nodes[arc_ [0]]
                lat2, lon2 = nodes[arc_ [1]]
            except KeyError:
                continue
            dist__ = haversine(lon1, lat1, lon2, lat2)
            values_ = line_atr[col_index["times"]].replace(
                "{","").replace("}","").split(",")
            values_ = [dist__/float(val) for val in values_]
            A.append(int(line_atr[col_index["id"]]))
            arcs[A[-1]] = arc_
            
            centers_arcs[arc_] = np.array([lat1 + lat2 , lon1 + lon2])/2
            d_all_arcs[arc_] = dist__
            speeds_[arc_] = []
        for val in values_:
            speeds_[arcs[int(
                    line_atr[col_index["id"]])]
                    ].append(val)

    
    G_ = nx.DiGraph()
    G_.add_nodes_from(nodes)
    G_.add_edges_from(arcs.values())
    mapping_arc_id = {id_arc : arc for id_arc,arc in enumerate(G_.edges())}
    
    dist_for_kernel_cmp = {}
    
    if not run_spatial:
        return G_,V,nodes,A,arcs,d_all_arcs,dist_for_kernel_cmp,speeds_,mapping_arc_id
    
    
    
    
    
    for a in G_.edges():
        inf_zone_a = []
        for id_node,location_lat_lon in nodes.items():
            if haversine(centers_arcs[a][1], centers_arcs[a][0], 
                         location_lat_lon[1], location_lat_lon[0])>3.5*__theta_:
                continue
            inf_zone_a.append(id_node)
        A_a = [a_ for a_ in G_.edges() if a_!=a and
               a_[0] in inf_zone_a or a_[1] in inf_zone_a ]
        G_a = G_.edge_subgraph(A_a).copy()
        G_a.add_nodes_from({a : centers_arcs[a]})
        G_a.add_edges_from([(a[0], a),(a, a[1])])
        c_ = {}
        for a_ in G_a.edges():
            c_[a_] = (d_all_arcs[a_] if a_ in d_all_arcs else
                      (haversine(centers_arcs[a_[0]][1],
                                 centers_arcs[a_[0]][0],
                                 nodes[a_[1]][1],
                                 nodes[a_[1]][0])
                           if type(a_[0]) is tuple else (
                          haversine(centers_arcs[a_[1]][1],
                                 centers_arcs[a_[1]][0],
                                 nodes[a_[0]][1],
                                 nodes[a_[0]][0]))))
        dist_a_i = heap_dijkstra_(G_a,c_,a)
        for a_ in A_a:
            val__ = (
                dist_a_i[a_[0]]+d_all_arcs[a_]/2)
            if val__ == float("inf"):
                continue
            if (a,a_) in dist_for_kernel_cmp:
                dist_for_kernel_cmp[a,a_] = dist_for_kernel_cmp[a_,a] = min(
                    val__,dist_for_kernel_cmp[a,a_])                                                                      
            else:
                dist_for_kernel_cmp[a,a_] = dist_for_kernel_cmp[a_,a] = val__
    

    if False:
        nx.draw(G_)
        plt.draw() 
        plt.show()
    return G_,V,nodes,A,arcs,d_all_arcs,dist_for_kernel_cmp,speeds_,mapping_arc_id

def sample_truth_from_mean(d_all_arcs , speeds_, mapping_arc_id):
    return {mapping_arc_id[index] : np.mean(speeds_[mapping_arc_id[index]])
        for index in mapping_arc_id.keys()}

def real_instance_sampler(*args):
    d_all_arcs , speeds_, mapping_arc_id, nodes, target, source, t_s_index = args
    target_, source_ = target, source
    set_of_pairs = [(29212, 41103),
        (62067, 18216),
        (95313, 59982),
        (63033, 72819),
        (100118, 57423),
        (87250, 33750),
        (13789, 99505),
        (56202, 26804),
        (90382, 17249),
        (6880, 15937)]
    if target_==source_:
        target_, source_ = set_of_pairs[t_s_index]
    def __function(V):
        nodes_names = [key for key in nodes.keys()]
        bol = target_==source_ and t_s_index==-1
        if not bol:
            target, source = target_, source_
        increasing_lon = np.argsort([val[1] for val in list(nodes.values())])
        increasing_lat = np.argsort([val[0] for val in list(nodes.values())])
        while bol:
            target, source = np.random.randint(85796), np.random.randint(85796)
            target, source = nodes_names[
                increasing_lon[target]
                ], nodes_names[increasing_lat[source]]
            bol = haversine(nodes[source][1], nodes[source][0], nodes[target][1],
                            nodes[target][0]) <= 23.
        print(target, source, haversine(nodes[source][1], nodes[source][0], nodes[target][1],
                              nodes[target][0]))
        return source,target,np.log([np.random.choice(
            speeds_[mapping_arc_id[index]])
                for index in mapping_arc_id.keys()])
    return __function

def set_function_as_kernel_on_G(dist_for_kernel_cmp,map_index_A,kernel_func,
                                __theta_):
    def kernel(index_set_1,index_set_2):
        phi___ = np.zeros((len(index_set_1),len(index_set_2)))
        map_index_1 = {i_ : i for i,i_ in enumerate(index_set_1)}
        map_index_2 = {i_ : i for i,i_ in enumerate(index_set_2)}
        for i in index_set_1:
            if i in index_set_2:
                phi___[map_index_1[i],map_index_2[i]] = 1.
        for (a1,a2),dist_a1_a2 in dist_for_kernel_cmp.items():
            i__ = map_index_A[a1]
            j__ = map_index_A[a2]
            if i__ in index_set_1 and j__ in index_set_2:
                phi___[map_index_1[i__],map_index_2[j__]
                       ] =    kernel_func(__theta_,dist_a1_a2)
                if i__ in index_set_2 and j__ in index_set_1:
                    phi___[map_index_1[j__],map_index_2[i__]
                           ] =    kernel_func(__theta_,dist_a1_a2)
        return phi___
    return kernel
if __name__ == "__main__":
    v_real = 40#km/h
    v_real = v_real / 3600
    _paths_dir = os.path.join(os.path.join(
        os.path.dirname(os.getcwd()), "data"), 
        "new_large")
    run_spatial = False
    run_naive = not run_spatial
    try:
        T_iter = int(sys.argv[1])
    except IndexError:
        T_iter = 2
    try:
        target, source = int(sys.argv[2]), int(sys.argv[3])
    except IndexError:
        target = source = 0
    try:
        t_s_index = int(sys.argv[2]) % 10
        id_run = int(sys.argv[2]) + 10
    except IndexError:
        t_s_index = -1
        id_run = np.random.randint(2**31)
    pp = T_iter==2
    np.random.seed(id_run + 10)
    only_show_observed_arcs = not run_spatial
    if pp: pp = PdfPages(os.path.join(_paths_dir,
                        f"real_instance_run_{id_run}"+
                        ".pdf"))
    arcs_paths_fl = "arcs.csv"
    nodes_fl = "nodes.csv"
    period_mod_plot = int(T_iter/4) if T_iter>=4 else 1
    v_prior = v_real
    _alpha__ = 1.
    _beta__ = 3
    __theta_ = 0.25
    kernel_name = "exponential"
    kernel = Kernels[kernel_name]
    (G_,V,nodes,A_,arcs,d_all_arcs,dist_for_kernel_cmp,speeds_,mapping_arc_id
         ) = read_from_arcs_nodes(_paths_dir , 
            arcs_paths_fl , nodes_fl, __theta_, run_spatial)
    map_index_A = {a : j_ for j_,a in enumerate(G_.edges())}
    prior_NormalIG = {}
    prior_NormalIG["mu"] = np.ones(len(A_)) * np.log(v_prior)
    prior_NormalIG["kappa"] = np.ones(len(A_)) 
    prior_NormalIG["alpha"] = np.ones(len(A_)) *_alpha__
    prior_NormalIG["beta"] = np.ones(len(A_)) * _beta__
    non_corrated_dist = NormalIG(prior_NormalIG)
    
    if run_naive:
        prior_Naive = {}
        naive_dist = {}
        init_naive = 0
        end_naive = 5
        for j__ in range(init_naive,end_naive):
            prior_Naive[j__] = {}
            prior_Naive[j__]["mu"] = np.ones(len(A_)) * np.log(v_prior)
            prior_Naive[j__]["kappa"] = np.ones(len(A_)) 
            prior_Naive[j__]["alpha"] = np.ones(len(A_)) *_alpha__
            prior_Naive[j__]["beta"] = np.ones(len(A_)) * _beta__
            prior_Naive[j__]["name"] = f"Naive_{np.round(.1 + (j__ * .2),2)}"
            prior_Naive[j__]["epsilon"] = .1 + (j__ * .2)
            naive_dist[j__] = naive_approach(prior_Naive[j__])
    
    if run_spatial:
        PB_prior_Spatial = {}
        PB_prior_Spatial["mu"] = np.ones(len(A_)) * np.log(v_prior)
        PB_prior_Spatial["kappa"] = np.ones(len(A_)) 
        PB_prior_Spatial["alpha"] = np.ones(len(A_)) *_alpha__
        PB_prior_Spatial["beta"] = np.ones(len(A_))  * _beta__
        PB_prior_Spatial["theta"] = __theta_
        PB_prior_Spatial["rho"] = {}
        PB_prior_Spatial["phi"]=dok_matrix((len(A_), len(A_)), 
                                           dtype=np.float64)
        for i__ in range(len(A_)):
            PB_prior_Spatial["phi"][i__,i__] = 1.
        PB_prior_Spatial["name"] = "PB_"+kernel_name+"_"+str(
            PB_prior_Spatial["theta"])
        PB_prior_Spatial["reversed"] = False
        PB_prior_Spatial["generate_cov_matrix"] = True
        PB_prior_Spatial["kernel"] = set_function_as_kernel_on_G(
            dist_for_kernel_cmp,map_index_A,kernel,__theta_)
        
        PA_prior_Spatial = {}
        PA_prior_Spatial["mu"] = np.ones(len(A_)) * np.log(v_prior)
        PA_prior_Spatial["kappa"] = np.ones(len(A_)) 
        PA_prior_Spatial["alpha"] = np.ones(len(A_)) *_alpha__
        PA_prior_Spatial["beta"] = np.ones(len(A_))  * _beta__
        PA_prior_Spatial["theta"] = __theta_
        PA_prior_Spatial["rho"] = {}
        PA_prior_Spatial["phi"] = dok_matrix((len(A_), len(A_)), 
                                             dtype=np.float64)
        for i__ in range(len(A_)):
            PA_prior_Spatial["phi"][i__,i__] = 1.
        PA_prior_Spatial["name"] = "PA_"+kernel_name+"_"+str(
            PA_prior_Spatial["theta"])
        PA_prior_Spatial["reversed"] = False
        PA_prior_Spatial["generate_cov_matrix"] = True
        PA_prior_Spatial["kernel"] = set_function_as_kernel_on_G(
            dist_for_kernel_cmp,map_index_A,kernel,__theta_)
        
        PA_spatial_distribution = Spatial_Metric_2(PA_prior_Spatial)
        PB_spatial_distribution = Spatial_Metric3(PB_prior_Spatial)
    
    
    
    if pp:
        from plot_util import (plot_distribution,plot_histogram_sigma,
                               plot_regret_t)
    
    distributions_ = [non_corrated_dist]
    if run_spatial:
        distributions_.append(PA_spatial_distribution)
    if run_naive:
        for j__,naive_dist_ in naive_dist.items():
            distributions_.append(naive_dist_)
    (regret_,Paths_,Observed_values_,Paths_expert,mu_updates,Time_,Exp_obj,
     delta_) = TS_Stationary(
         T_iter,G_,truth_sampler(real_instance_sampler(
             d_all_arcs , speeds_, {id_arc : arc for id_arc,arc in enumerate(G_.edges())},
             nodes, target, source, t_s_index),V,sample_truth_from_mean(
                 d_all_arcs , speeds_, {id_arc : arc for id_arc,arc in enumerate(G_.edges())})),
         d_all_arcs,distributions_)
    
    regret_NormalIG = regret_[0]
    Paths_NormalIG = Paths_[0]
    Time_NormalIG = Time_[0]
    
    if run_naive:
        regret_Naive = {}
        Paths_Naive = {}
        Time_Naive = {}
        for j__ in range(1,len(Time_)-1):
            regret_Naive[j__] = regret_[j__]
            Paths_Naive[j__] = Paths_[j__]
            Time_Naive[j__] = Time_[j__]
    
    
    if run_spatial:
        PA_regret_Spatial_ = regret_[-1]
        PA_Paths_Spatial_ = Paths_[-1]
        PA_Time_Spatial_ = Time_[-1]
    
    if pp:
        mean_real_instance = sample_truth_from_mean(d_all_arcs , speeds_,
                        {id_arc : arc for id_arc,arc in enumerate(G_.edges())})
        visited_arcs = [[] for j__ in range(len(Time_))]
        for j,path in enumerate(Paths_expert):
            if not run_spatial and j%period_mod_plot!=0 and j!=len(
                Paths_expert)-1 and (j>0 and regret_NormalIG[j] -
                        regret_NormalIG[j-1]<1.): continue
            if True:
                values = list()
                for i,a in enumerate(G_.edges()):
                    value = min(
                        mean_real_instance[a]*1000
                                 ,20
                                )
                    values.append([0,nodes[a[0]][0],nodes[a[0]][1],
                            nodes[a[1]][0],nodes[a[1]][1],
                            (value 
                              if not a in path[1] and
                              not (a[1],a[0]) in path[1]
                              else "NaN") ])
                values = list()
                for i,a in enumerate(G_.edges()):
                    value = min(mean_real_instance[a]*1000,20
                             )
                    values.append([0,nodes[a[0]][0],nodes[a[0]][1],
                            nodes[a[1]][0],nodes[a[1]][1],
                            value ])
        
            for j__ in range(len(Time_)):
                if not run_spatial and j%period_mod_plot!=0 and j!=len(
                    mu_updates[j__])-1 and (j>0 and regret_NormalIG[j] -
                        regret_NormalIG[j-1]<1.): continue
                for path_prime in Paths_[j__][:(j+1)]:
                    for a in path_prime[1]:
                        if a not in visited_arcs[j__]:
                            visited_arcs[j__].append(a)
                values = list()
                for i in range(len(G_.edges())):
                    a = mapping_arc_id[i]
                    if only_show_observed_arcs and a not in visited_arcs[j__]:
                        continue
                    value = min(mu_updates[j__][j][i]*1000,20
                                ) 
                    values.append([0,nodes[a[0]][0],nodes[a[0]][1],
                                nodes[a[1]][0],nodes[a[1]][1],
                                (value
                                  if not a in Paths_[j__][j][1] and
                                  not (a[1],a[0]) in Paths_[j__][j][1]
                                  else "NaN") ])
                values = list()
                for i in range(len(G_.edges())):
                    a = mapping_arc_id[i]
                    if only_show_observed_arcs and a not in visited_arcs[j__]:
                        continue
                    value = min(mu_updates[j__][j][i]*1000,20
                                ) 
                    values.append([0,nodes[a[0]][0],nodes[a[0]][1],
                                nodes[a[1]][0],nodes[a[1]][1],value ])
            
        plot_regret_t([("Spatial_" + PA_prior_Spatial["name"]
                        ) if run_spatial else prior_Naive[j__]["name"]
                       for j__ in range(init_naive, end_naive)],
                      regret_NormalIG, [PA_regret_Spatial_ if run_spatial else
                      regret_Naive[j__] for j__ in range(1,len(Time_)-1)],pp,plt)
        pp.close()
    (pd.DataFrame.from_dict(
        {**{("Independent" if j__==0 else ("Spatial_PA" if run_spatial else 
            f"Naive_{np.round(.1 + ((init_naive + j__ - 1) * .2),2)}")
          ) : a_r  for j__,a_r in enumerate(regret_) }, **{
              ("Independent" if j__==0 else ("Spatial_PA" if run_spatial else 
            f"Naive_{np.round(.1 + ((init_naive + j__ - 1) * .2),2)}")
          ) + "_delta" : a_r  for j__,a_r in enumerate(delta_)}}
        )).to_csv(
                    os.path.join(_paths_dir,
                    f"SIMULATION_Regret_{T_iter}_"+
                    ("Spatial_PA" if run_spatial else "Naive") + f"{id_run}"+
                    ".csv"),index=False)
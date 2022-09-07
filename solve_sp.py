#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 23:47:11 2021

@author: tomas
"""

import gurobipy as gp
from gurobipy import GRB


import networkx as nx

from time import time

from _meta_util import timer_func

def form_feasible_flow(is_digraph,m,G,V,A,b,x,C):
    for i in V:
        C["null_divergence",i]=m.addConstr(
                gp.quicksum(x[(pred_,i)] for pred_ in (G.predecessors(i) if is_digraph
                            else [(e[0] if e[0]!=i else e[1]) for e in  G.edges(i)] )) - 
                gp.quicksum(x[(i,succ_)] for succ_ in (G.successors(i) if is_digraph
                            else [(e[0] if e[0]!=i else e[1]) for e in  G.edges(i)] )) == (
                                b[i] #-1 if i==source else (1 if i==target else 0) 
                        ) , name=f"flow_balance_{i}");
    m.update();


def Solve_SP_LP(G,V,A,speed,dist,source,target,timelimit=3600.0,debug=0,debug_sp=False):
    m = gp.Model("sp");
    m.setParam( 'OutputFlag', debug_sp);
    m.setParam( 'Threads', 1);
#    m.setParam( 'MIPGap', 0. );
    m.setParam(GRB.Param.TimeLimit, timelimit)
    
    is_digraph=isinstance(G,type(nx.DiGraph()))
    
    x={};
    C={};    
    for a in A:
        # try:
        x[a] = m.addVar(vtype=GRB.CONTINUOUS, lb=0, ub=1,obj=(
            dist[a]/speed[a] if float("inf")!=speed[a] else 0.),name=f"x_{a}");
        if not is_digraph:
            x[(a[1],a[0])] = m.addVar(vtype=GRB.CONTINUOUS, lb=0, ub=1,obj=(
                dist[a]/speed[a] if float("inf")!=speed[a] else 0.),name=f"x_{(a[1],a[0])}");
        # except gp.GurobiError:
        #     print("here")
#        m.addVar(vtype=GRB.BINARY, lb=0);#,obj=1.);
    form_feasible_flow(is_digraph,m,G,V,A,{i : -1 if i==source else (
        1 if i==target else 0) for i in V},x,C)
    
    m.modelSense = GRB.MINIMIZE
#    m.update();
    m.optimize();
    # m.write("modelsp.lp");
    if m.status == GRB.Status.OPTIMAL:
        if debug: print('Optimal objective: %g' % m.objVal)
        return m.objVal ,[a for a in A if (x[a].x>=0.99 or (not is_digraph and  x[(a[1],a[0])].x>=0.99) )];
    elif m.status == GRB.Status.INF_OR_UNBD:
        if debug: print('Model is infeasible or unbounded')
        return None
    elif m.status == GRB.Status.INFEASIBLE:
        if debug: 
            print('Model is infeasible')
            m.computeIIS()
            m.write("sp.ilp")
            print('\nThe following constraint(s) cannot be satisfied:')
            for c in m.getConstrs():
                if c.IISConstr:
                    print('%s' % c.constrName)
                    print('')
        return None
    elif m.status == GRB.Status.UNBOUNDED:
        if debug: print('Model is unbounded')
        return None
    else:
        if debug: print('Optimization ended with status %d' % m.status)
        return None


# def Solve__all_pairs_SP_LP(G,V,A,speed,dist,source,target,timelimit=3600.0,
#                            debug=0,debug_sp=False):
#     m = gp.Model("ap_sp");
#     m.setParam( 'OutputFlag', debug_sp);
#     m.setParam( 'Threads', 1);
# #    m.setParam( 'MIPGap', 0. );
#     m.setParam(GRB.Param.TimeLimit, timelimit)
    
#     is_digraph=isinstance(G,type(nx.DiGraph()))
    
#     x={};
#     C={};    
#     for a in A:
#         # try:
#         x[a] = m.addVar(vtype=GRB.CONTINUOUS, lb=0, ub=1,obj=(
#             dist[a]/speed[a] if float("inf")!=speed[a] else 0.),name=f"x_{a}");
#         if not is_digraph:
#             x[(a[1],a[0])] = m.addVar(vtype=GRB.CONTINUOUS, lb=0, ub=1,obj=(
#                 dist[a]/speed[a] if float("inf")!=speed[a] else 0.),name=f"x_{(a[1],a[0])}");
#         # except gp.GurobiError:
#         #     print("here")
# #        m.addVar(vtype=GRB.BINARY, lb=0);#,obj=1.);
#     # for i in V:
#     #     C["null_divergence",i]=m.addConstr(
#     #             gp.quicksum(x[(pred_,i)] for pred_ in (G.predecessors(i) if is_digraph
#     #                         else [(e[0] if e[0]!=i else e[1]) for e in  G.edges(i)] )) - 
#     #             gp.quicksum(x[(i,succ_)] for succ_ in (G.successors(i) if is_digraph
#     #                         else [(e[0] if e[0]!=i else e[1]) for e in  G.edges(i)] )) == (
#     #                     -1 if i==source else (1 if i==target else 0) ) , name=f"flow_balance_{i}");
            
#     form_feasible_flow(is_digraph,m,G,V,A,[-len(target)+1 if i==source else (
#         1 if i in target else 0) for i in V],x,C)
#     # m.update();
#     m.modelSense = GRB.MINIMIZE
# #    m.update();
#     m.optimize();
#     # m.write("modelsp.lp");
#     if m.status == GRB.Status.OPTIMAL:
#         if debug: print('Optimal objective: %g' % m.objVal)
#         return m.objVal ,[a for a in A if (x[a].x>=0.99 or (not is_digraph and  x[(a[1],a[0])].x>=0.99) )];
#     elif m.status == GRB.Status.INF_OR_UNBD:
#         if debug: print('Model is infeasible or unbounded')
#         return None
#     elif m.status == GRB.Status.INFEASIBLE:
#         if debug: 
#             print('Model is infeasible')
#             m.computeIIS()
#             m.write("sp.ilp")
#             print('\nThe following constraint(s) cannot be satisfied:')
#             for c in m.getConstrs():
#                 if c.IISConstr:
#                     print('%s' % c.constrName)
#                     print('')
#         return None
#     elif m.status == GRB.Status.UNBOUNDED:
#         if debug: print('Model is unbounded')
#         return None
#     else:
#         if debug: print('Optimization ended with status %d' % m.status)
#         return None

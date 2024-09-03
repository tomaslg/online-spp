#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 01:43:33 2021

@author: tomas
"""

from TS import nx,format_distrib_name
from distrib import stats,np


import pandas as pd
import subprocess
from matplotlib.font_manager import FontProperties

show_plot = False

font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
font = FontProperties(size=16)
font.set_family('serif')
font_dict = {'fontsize': 'x-large', 'fontproperties': font}

def replace_infty_by_max(np, ar):
    max_no_infty = max( [a for a in ar if np.abs(a)!=np.inf] )
    min_no_infty = min( [a for a in ar if np.abs(a)!=np.inf] )
    quant95 = np.quantile(ar, .8)
    quant05 = np.quantile(ar, .2)
    return np.array( [(a if quant05 <= a <= quant95 else (
        quant05 if a<quant05 else quant95)
                       )  if np.abs(a)!=np.inf else (
        max_no_infty if np.sign(a)==1 else min_no_infty) for a in ar])
def replace_by_bins(np, ar, number_of_bins = 76):
    sorted_indices = sorted(range(len(ar)), key=lambda i: ar[i])
    edge_color_ = {}
    bin_size = int(np.ceil(len(ar) / number_of_bins))
    bin_average = None
    value_for_min_bins = float("Inf")
    cut_off = .1
    for i in range(number_of_bins):
        bin_elements = [sorted_indices[j] for j in range(
            i * bin_size, min((i + 1) * bin_size, len(ar)))]
        value_for_min_bins = min(np.mean(ar[bin_elements]), value_for_min_bins)
        bin_average = np.mean(ar[bin_elements]) if cut_off <= (
            i / number_of_bins) <= (1 - cut_off) else bin_average
        for j in bin_elements:
            edge_color_[j] = bin_average
    return np.array([edge_color_[i] if edge_color_[i]!=None else 
        value_for_min_bins for i in range(len(ar))])

def plot_distribution(distrib_,pp,plt_,Title=""):
    if pp==False and not show_plot: return None
    _lambda_ = distrib_.sample_lambda_posterior()
    dist_mu_=stats.norm(distrib_.mu,1/(distrib_.kappa*_lambda_)**.5)
    _x_ = np.linspace(dist_mu_.ppf(0.001),dist_mu_.ppf(0.999), 100)
    fig, ax_ = plt_.subplots(1, 1, sharey=True, tight_layout=True)
    ax_.plot(_x_, dist_mu_.pdf(_x_),'r-', lw=5, alpha=0.6, label='mu pdf')
    ax_.set_title(f"{Title}_{distrib_.name}")
    plt_.xlabel('mu', fontdict=font_dict)
    plt_.ylabel('density', fontdict=font_dict)
    if pp==False:
        plt_.show()
    else:
        plt_.savefig(pp, format='pdf')

    plt_.close()
def plot_histogram_sigma(distrib_,pp,plt_,Title=""):
    if pp==False and not show_plot: return None
    Yx=distrib_.sample_lambda_posterior()
    fig, ax_ = plt_.subplots(1, 1, sharey=True, tight_layout=True)
    ax_.hist(Yx, bins=1000)
    ax_.set_title(f"{Title}_{distrib_.name}")
    plt_.xlabel('lambda', fontdict=font_dict)
    plt_.ylabel('frec.', fontdict=font_dict)
    if pp==False:
        plt_.show()
    else:
        plt_.savefig(pp, format='pdf')

    plt_.close()

def get_plot(G,edge_color_,pp,plt_,Title="",_paths_=[]):
    if pp==False and not show_plot: return None
    pos={}
    for (i,j) in G.nodes():
        pos[(i,j)]=( i , j )
    plt_.rcParams['text.usetex'] = True
    ax = plt_.gca()
    ax.set_title(Title, fontdict=font_dict)
    if all([type(ec) is float or isinstance(ec, np.floating)
            for ec in edge_color_]):
        edge_color_ = replace_infty_by_max(np, edge_color_)
    nx.draw_networkx(G, with_labels = False,pos=pos,
                     node_size=1/(len(G.nodes())**.5),
                     arrowsize=60/(len(G.nodes())**.5),
                     width=60/(len(G.nodes())**.5),
                     edge_cmap=plt_.get_cmap("RdYlGn"),
                     edge_color=edge_color_,
                     )
    if len(_paths_)>0:
        path_pos = {(i,j) : pos[(i,j)] for (i,j) in pos.keys()
                    if (i,j) in (
            [a__[0] for a__ in _paths_] + [a__[1] for a__ in _paths_] ) }
        nx.draw_networkx_nodes(G,pos,nodelist=list(path_pos.keys()),
                               node_size=1/(len(G.nodes())**.5),node_color='k')
        nx.draw_networkx_edges(G,pos,edgelist=_paths_,edge_color='k',
                               width=60/(len(G.nodes())**.5))
    
    if pp==False:
        plt_.show()
    else:
        plt_.savefig(pp, format='pdf')

    plt_.close()


def plot_regret_t(prior_Spatial_name,regret_NormalIG,regret_Spatial_,pp,plt_):
    if pp==False and not show_plot: return None
    plt_.plot((list(range(1,len(regret_NormalIG)+1))), 
             regret_NormalIG, 
             label="NormalIG = " +
             f'{np.round(regret_NormalIG[-1],2)}', 
             color='k')
    colors_ = ['b','r', "g", "m", "c"]
    for _color, regret_Spatial, name in zip(colors_, regret_Spatial_, prior_Spatial_name): 
        plt_.plot((list(range(1,len(regret_Spatial)+1))), 
                 regret_Spatial, 
                 label= format_distrib_name(name) +
                  f' = {np.round(regret_Spatial[-1],2)}', 
                 color=_color)
    plt_.legend()
    plt_.title('Travel Time Regret', fontdict=font_dict)
    plt_.xlabel('periods', fontdict=font_dict)
    plt_.ylabel('[hrs]', fontdict=font_dict)
    
    # Tweak spacing to prevent clipping of ylabel
    plt_.subplots_adjust(left=0.15)
    if pp:
        plt_.savefig(pp, format='pdf')
    else:
        plt_.show()
    plt_.close()

def get_graph_instance_1(N=20):
    return nx.grid_2d_graph(N,N,periodic=False,create_using=nx.DiGraph)
    # return G

def get_graph_instance_2(N=20):
    return nx.grid_2d_graph(N,N,periodic=False,create_using=nx.Graph)




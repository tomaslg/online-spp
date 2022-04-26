#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 01:43:33 2021

@author: tomas
"""

from TS import nx
from distrib import stats,np

show_plot = False

def plot_distribution(distrib_,pp,plt_,Title=""):
    # Yx=distrib_.sample_poserior()
    if pp==False and not show_plot: return None
    _lambda_ = distrib_.sample_lambda_posterior()#(stats.invgamma(self.alpha, loc=0, scale=self.beta)).rvs()
    dist_mu_=stats.norm(distrib_.mu,1/(distrib_.kappa*_lambda_)**.5)
    _x_ = np.linspace(dist_mu_.ppf(0.001),dist_mu_.ppf(0.999), 100)
    fig, ax_ = plt_.subplots(1, 1, sharey=True, tight_layout=True)
    ax_.plot(_x_, dist_mu_.pdf(_x_),'r-', lw=5, alpha=0.6, label='mu pdf')
    ax_.set_title(f"{Title}_{distrib_.name}")
    font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
    plt_.xlabel('mu', fontdict=font)
    plt_.ylabel('density', fontdict=font)
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
    font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
    plt_.xlabel('lambda', fontdict=font)
    plt_.ylabel('frec.', fontdict=font)
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
    # nx.set_node_attributes(G, pos, 'coord')
    # nx.draw(G)  # networkx draw()
    # plt.draw()  # pyplot draw()
    # pp = PdfPages(f'graph_speeds_N{N}_M{M}.pdf')
    ax = plt_.gca()
    ax.set_title(Title)
    nx.draw_networkx(G, with_labels = False,pos=pos,
                     node_size=100/(len(G.nodes())**.5),
                     arrowsize=60/(len(G.nodes())**.5),
                     edge_cmap=plt_.get_cmap("RdYlGn_r"),#edge_cmap=plt_.get_cmap("RdYlGn"),#plt.cm.Blues,cmap="RdYlGn",
                     edge_color=edge_color_)
    if len(_paths_)>0:
        path_pos = {(i,j) : pos[(i,j)] for (i,j) in pos.keys()
                    if (i,j) in (
            [a__[0] for a__ in _paths_] + [a__[1] for a__ in _paths_] ) }
        # len_path_ = len(path_pos)
        # for i in range(len_path_):
        #     if path_pos[len_path_-i-1] in path_pos[:(len_path_-i-1)]:
        #         path_pos.pop(len_path_-i-1)
        nx.draw_networkx_nodes(G,pos,nodelist=list(path_pos.keys()),
                               node_size=100/(len(G.nodes())**.5),node_color='k')
        nx.draw_networkx_edges(G,pos,edgelist=_paths_,edge_color='k',
                               width=60/(len(G.nodes())**.5))
        # nx.draw_networkx_edges(
        #     G.edge_subgraph(_paths_),
        #     edgelist=_paths_,
        #                 pos=path_pos,node_size=100/(len(G.nodes())**.5),
        #              arrowsize=60/(len(G.nodes())**.5),
        #              # edge_cmap=plt_.get_cmap("RdYlGn_r"),#edge_cmap=plt_.get_cmap("RdYlGn"),#plt.cm.Blues,cmap="RdYlGn",
        #              edge_color="k")
    
    if pp==False:
        plt_.show()
    else:
        plt_.savefig(pp, format='pdf')

    plt_.close()

def plot_regret_t(prior_Spatial_name,regret_NormalIG,regret_Spatial_,pp,plt_):
    if pp==False and not show_plot: return None
    plt_.plot((list(range(1,len(regret_NormalIG)+1))), 
             regret_NormalIG, 
             label='sum(regret_NormalIG)/T='+
             f'{(regret_NormalIG[-1])}', 
             color='k')
    plt_.plot((list(range(1,len(regret_Spatial_[0])+1))), 
             regret_Spatial_[0], 
             label='sum(regret_Spatial_PB'+
             f'{prior_Spatial_name[0]})/T={regret_Spatial_[0][-1]}', 
             color='b')
    plt_.plot((list(range(1,len(regret_Spatial_[1])+1))), 
             regret_Spatial_[1], 
             label='sum(regret_Spatial_PA'+
             f'{prior_Spatial_name[1][0]})/T={regret_Spatial_[1][-1]}', 
             color='r')
    plt_.legend()
    font = {'family': 'serif',
            'color':  'darkred',
            'weight': 'normal',
            'size': 16,
            }
    plt_.title('Travel Time Regret', fontdict=font)
    plt_.xlabel('periods', fontdict=font)
    plt_.ylabel('sum(regret_t)/t [hrs]', fontdict=font)
    
    # Tweak spacing to prevent clipping of ylabel
    plt_.subplots_adjust(left=0.15)
    # 



    if pp:
        plt_.savefig(pp, format='pdf')
    else:
        plt_.show()
    plt_.close()
# nx.grid_2d_graph(m, n[, periodic, create_using])

def get_graph_instance_1(N=20):
    return nx.grid_2d_graph(N,N,periodic=False,create_using=nx.DiGraph)
    # return G

def get_graph_instance_2(N=20):
    return nx.grid_2d_graph(N,N,periodic=False,create_using=nx.Graph)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 01:43:33 2021

@author: tomas
"""

from TS import nx,format_distrib_name
from distrib import stats,np


# import matplotlib.pyplot as plt
import pandas as pd
import subprocess
# subprocess.check_output('conda activate geo', shell=True)
import geopandas as gpd
# import psycopg2
from shapely.geometry import LineString
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
    return np.array( [a  if np.abs(a)!=np.inf else (
        max_no_infty if np.sign(a)==1 else min_no_infty) for a in ar])

def plot_distribution(distrib_,pp,plt_,Title=""):
    # Yx=distrib_.sample_poserior()
    if pp==False and not show_plot: return None
    _lambda_ = distrib_.sample_lambda_posterior()#(stats.invgamma(self.alpha, loc=0, scale=self.beta)).rvs()
    dist_mu_=stats.norm(distrib_.mu,1/(distrib_.kappa*_lambda_)**.5)
    _x_ = np.linspace(dist_mu_.ppf(0.001),dist_mu_.ppf(0.999), 100)
    fig, ax_ = plt_.subplots(1, 1, sharey=True, tight_layout=True)
    ax_.plot(_x_, dist_mu_.pdf(_x_),'r-', lw=5, alpha=0.6, label='mu pdf')
    ax_.set_title(f"{Title}_{distrib_.name}")
    # font = {'family': 'serif',
    #     'color':  'darkred',
    #     'weight': 'normal',
    #     'size': 16,
    #     }
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
    # font = {'family': 'serif',
    #     'color':  'darkred',
    #     'weight': 'normal',
    #     'size': 16,
    #     }
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
    # nx.set_node_attributes(G, pos, 'coord')
    # nx.draw(G)  # networkx draw()
    # plt.draw()  # pyplot draw()
    # pp = PdfPages(f'graph_speeds_N{N}_M{M}.pdf')
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
                     edge_cmap=plt_.get_cmap("RdYlGn"),#edge_cmap=plt_.get_cmap("RdYlGn"),#plt.cm.Blues,cmap="RdYlGn",
                     edge_color=edge_color_,
                     # missing_kwds= dict(color = "k",)
                     )
    if len(_paths_)>0:
        path_pos = {(i,j) : pos[(i,j)] for (i,j) in pos.keys()
                    if (i,j) in (
            [a__[0] for a__ in _paths_] + [a__[1] for a__ in _paths_] ) }
        # len_path_ = len(path_pos)
        # for i in range(len_path_):
        #     if path_pos[len_path_-i-1] in path_pos[:(len_path_-i-1)]:
        #         path_pos.pop(len_path_-i-1)
        nx.draw_networkx_nodes(G,pos,nodelist=list(path_pos.keys()),
                               node_size=1/(len(G.nodes())**.5),node_color='k')
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
             label="NormalIG = " +#'sum(regret_NormalIG)/T='+
             f'{np.round(regret_NormalIG[-1],2)}', 
             color='k')
    colors_ = ['b','r', "g", "m", "c"]
    for _color, regret_Spatial, name in zip(colors_, regret_Spatial_, prior_Spatial_name): 
        plt_.plot((list(range(1,len(regret_Spatial)+1))), 
                 regret_Spatial, 
                 label= format_distrib_name(name) + # 'sum(regret_'+
                  f' = {np.round(regret_Spatial[-1],2)}', 
                 color=_color)
    plt_.legend()
    # font = {'family': 'serif',
    #         'color':  'darkred',
    #         'weight': 'normal',
    #         'size': 16,
    #         }
    plt_.title('Travel Time Regret', fontdict=font_dict)
    plt_.xlabel('periods', fontdict=font_dict)
    plt_.ylabel('[hrs]', fontdict=font_dict)
    
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




    
def plot_real_life_instance(values,title,pp,plt_,color_scale_='RdYlGn'):
    # data 
    # sql = """select w.id as arc, n1.latitude as lat1, n1.longitude as lon1, 
    #      n2.latitude as lat2, n2.longitude as lon2, s.speed as speed
    #      from way w join node as n1 on w.node1 = n1.id 
    #      join node as n2 on w.node2 = n2.id join speed as s
    #      on w.id = s.id where s.hour = 8 and (valores is not NULL and array_length(s.valores, 1) > 1)"""

    data = dict()
    data['LineString_obj'] = list()
    for i in range(len(values)):
        data['LineString_obj'].append(LineString([[float(values[i][2]), 
                                                    float(values[i][1])], 
                                                  [float(values[i][4]), 
                                                    float(values[i][3])]]))
    data['speed'] = [(1 if i==0 else (-1 if i==1 else 0
        )) if not np.isnan(float(values[i][
        5])) and color_scale_=="gist_yarg" else float(
        values[i][5]) for i in range(len(values))]
    print(title)
    print(f"min()={min(data['speed'])}",
          f"max()={max(data['speed'])}",
          f"mean()={np.mean(data['speed'])}")
    df = pd.DataFrame(data)
    geo_df = gpd.GeoDataFrame(df, geometry = 'LineString_obj')
    
    # create figure and axes, assign to subplot
    fig, ax = plt_.subplots(figsize=(7,6))
    geo_df.plot(column='speed', ax=ax, legend=color_scale_!="gist_yarg",
                cmap=color_scale_, legend_kwds={'shrink': 0.45, 'label': 'meters/second'},
                missing_kwds= dict(color = "k", linewidth=2)) #lightgrey brg
    # font = FontProperties()
    # font.set_family('serif')
    plt_.title(title#"Beijing"
                , fontdict={'fontsize': 'large', 'fontproperties': font})
    plt_.ylabel("Latitude", fontdict={'fontsize': 'large', 'fontproperties': font})
    plt_.xlabel("Longitude", fontdict={'fontsize': 'large', 'fontproperties': font})
    # fig.savefig("../../real_instance.pdf", bbox_inches='tight')
    if pp:
        plt_.savefig(pp, format='pdf')
    else:
        plt_.show()
    plt_.close()



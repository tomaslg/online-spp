#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 28 23:23:48 2023

@author: tomas
"""


import geopandas as gpd
from shapely.geometry import LineString
from distrib import stats,np
import pandas as pd

font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
    
def plot_real_life_instance(values,title,pp,plt_,color_scale_='RdYlGn'):

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
    
    fig, ax = plt_.subplots(figsize=(7,6))
    geo_df.plot(column='speed', ax=ax, legend=color_scale_!="gist_yarg",
                cmap=color_scale_, legend_kwds={'shrink': 0.45, 'label': 'meters/second'},
                missing_kwds= dict(color = "k", linewidth=2))
    
    plt_.title(title
                , fontdict={'fontsize': 'large', 'fontproperties': font})
    plt_.ylabel("Latitude", fontdict={'fontsize': 'large', 'fontproperties': font})
    plt_.xlabel("Longitude", fontdict={'fontsize': 'large', 'fontproperties': font})
    if pp:
        plt_.savefig(pp, format='pdf')
    else:
        plt_.show()
    plt_.close()



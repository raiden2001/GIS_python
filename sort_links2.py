# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 15:57:33 2021

@author: LESNIAM
"""

import geopandas as gpd
import mapping_module6 as MM
"""
imports the IPython Library for display data function
Trois-Rivieres2.shp
"""
from IPython.display import display
#from arcpy.display import display

gdf_corr = gpd.read_file('C:/Users/Jason/Downloads/New folder (35)/Trois_N.shp')
corridor_name = 'Trois_N'

status,message,gdf_corridor_sorted = MM.Order_the_Links(corridor_name,gdf_corr)
print(message)


#gdf_node = MM.pop_row(gdf_corridor_sorted,gdf_corr)
#if gdf_corridor_sorted is not None:
if status == 1:
    display(gdf_corridor_sorted)
    gdf_corridor_sorted = gpd.GeoDataFrame(gdf_corridor_sorted,geometry='geometry') # print the sorted data into the GEODATAFRAME on to ARCMAP
    gdf_corridor_sorted.crs ={'init': 'epsg:4326'}  # prints out the projection plane
    gdf_corridor_sorted.to_file('C:/Users/Jason/Downloads/New folder (35)/Trois_N_done.shp',driver = 'ESRI Shapefile',index=True)

    print("successful sorting")

else:
    #if the status is not 1, the sorting was not successful, so we just print the appropriate message and do not proceed further
    print("not the direction shapefile needs, it will not be sorted")
    

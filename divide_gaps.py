import geopandas as gpd
import pandas as pd
import os


shapefile_dir = 'C:/Users/Jason/Downloads/New folder (35)/Trois_update.shp'


shapefile = gpd.read_file(shapefile_dir)

#splits the data into two dataframes, one for 'E' and other for 'W'
shapefile_n = shapefile[shapefile['DIR'] == 'N']
shapefile_s = shapefile[shapefile['DIR'] == 'S']

#Write the two dataframes to separate shapeifiles

output_file_n = os.path.join(os.path.dirname(shapefile_dir),'Trois_N.shp')
shapefile_n.to_file(output_file_n, driver = 'ESRI Shapefile')




output_file_s = os.path.join(os.path.dirname(shapefile_dir),'Trois_S.shp')
shapefile_s.to_file(output_file_s, driver = 'ESRI Shapefile')
print('All splited successfully')

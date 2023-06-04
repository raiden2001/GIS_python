import geopandas as gpd

new_gdf = gpd.read_file('C:/Users/Jason/Downloads/sorted/SortedCorridors/Toronto3_N_done.shp')
#new_gdf = gpd.read_file('G:/Jason Ma/SortedCorridors/Edmonton1_N_done.shp')


#Create a list of position values
positions = list(range(1,len(new_gdf)+1))


#Add the 'POS_ALONG' column to the GEODATAFRAME and lable the rows
new_gdf['POS_ALONG_'] = positions

#save the update file
new_gdf.to_file('C:/Users/Jason/Downloads/sorted/SortedCorridors/Toronto3_N_done.shp',driver='ESRI Shapefile')
print("ouput success")

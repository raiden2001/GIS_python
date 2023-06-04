import geopandas as gpd

gdf = gpd.read_file('C:/Users/Jason/Downloads/sorted/SortedCorridors/2022Q1_Corridors7.shp')


#drop the columns not required
gdf = gdf.drop(columns=['FID','Shape'])


gdf.to_file('C:/Users/Jason/Downloads/sorted/SortedCorridors/2022Q1_Corridors7.shp',driver='ESRI Shapefile')
print('columns dropped')

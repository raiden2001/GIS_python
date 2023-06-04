import geopandas as gpd
import os
import numpy as np

#Specifies the path to the sahpefile 
shapefile_dir = r'C:/Users/Jason/Downloads/sorted/SortedCorridors/Edmonton1_N_done.shp'

shapefile = gpd.read_file(shapefile_dir)




#Extract the city and corridor name from the shapefile's name
file_name = shapefile_dir.split('/')[-1]
city = file_name.split('_')[0]
corridor = '_'.join(file_name.split('_')[:2])

#loop through each shapefiles from the list


'''
for shapefile in shapefile_dir:
    if 'CORR' not in shapefile.columns or 'PROV' not in shapefile.columns:
        print(f'skipping{shapefile} as it does not contain required columns')
        continue
'''

#Set the 'CITY' , 'PROV' and 'CORR' columns to their respective values
shapefile['CITY'] = 'Edmonton'
shapefile['PROV'] = 'AB'
shapefile['CORR'] = 'Edmonton1'


#Creates new Field
shapefile['LINK_DIR'] = shapefile['LINK_ID'].astype(str) + shapefile['DIR_TRAVEL']
shapefile['CORR_ID']= shapefile['CORR'] + '_' +shapefile['DIR']
#shapefile['Length_km'] = shapefile[' Length_Km'].fillna(0.0)
#shapefile['AADT']= 0
#shapefile['AADT_Year']= 0


#Save the modified shapefile
shapefile.to_file(shapefile_dir, driver='ESRI Shapefile')
print("adding columns sucessfully")

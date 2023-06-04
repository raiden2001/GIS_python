import geopandas as gpd
import os
import pandas as pd
import shapely


#specify the location of the folder containing the corridors of  shapefiles 
shapefiles_path = r'C:/Users/Jason/Downloads/sorted/SortedCorridors/'

#Get a list of all the shapefiles in the folder
shapefiles_list = [gpd.read_file(os.path.join(shapefiles_path,shapefile),encoding='utf-8' )
                   for shapefile in os.listdir(shapefiles_path)
                   if shapefile.endswith('_E_done.shp') or
                      shapefile.endswith('_W_done.shp') or
                      shapefile.endswith('_N_done.shp') or
                      shapefile.endswith('_S_done.shp')]

#creates an empty list of the shapefiles
master_df_list = []

#loop through each of the shapefiles form the list
for shapefile,file_name in zip(shapefiles_list,os.listdir(shapefiles_path)):
    #Check if the LINK_ID is in the columns of the shapefiles 
    if 'LINK_ID' not in shapefile.columns:  # prevents the Keyerror , but looking for specific column field name 
        continue

        #keep only the specified columnsin the master shapefiles
        shapefile = ['LINK_ID','REF_IN_ID','NREF_IN_ID','DIR_TRAVEL','FUNC_CLASS','CONTRACC','POS_ALONG_','LINK_DIR','PROV','CITY','CORR','DIR','CORR_ID']

        # Drop all columns except the ones to keep
        master_df.drop(columns=[col for col in master_df.columns if col not in shapefile], inplace=True)
    
    #Creates the new columns in the shapefile  one Dataframe to another
    #shapefile = shapefile.copy()
    #shapefile['geometry']  = None
   # shapefile['LINK_DIR'] = ''
    #shapefile['PROV'] = ''
    #shapefile['CITY'] =''
    #shapefile['CORR'] =''
    #shapefile['DIR'] = ''
    #shapefile['CORR_ID'] = ''
    shapefile['Length_Km'] = 0.0
    shapefile['AADT'] = 0
    shapefile['AADT_Year'] = 0
   #shapefile['POS_ALONG'] = 0
        
         
    
#set the coordinate reference system (CRS) of the GeoDataFrame to that of the first shapefile in the list 
#master_shapefile.crs = shapefiles_list[0].crs



#calculate the length of each link in the shapefile 
    shapefile['length_km'] = shapefile['geometry'].length / 1000.00

#Append the shapefile to the master list
    master_df_list.append(shapefile)

#Concatenate all the individual shapefiles into the one large shapefile
master_df = pd.concat(master_df_list,ignore_index=True,sort=False)

#create a list of position of values
#positions = list(range(1,len(master_df)+ 1))

#Add the 'POS_ALONG' column to the master GEODataFrame and label the rows
# master_df['POS_ALONG'] = positions

#Creates new fields
#master_df['LINK_DIR'] = master_df['LINK_ID'].astype(str)+ master_df['DIR_TRAVEL']
#master_df['CORR_ID'] = master_df['CORR'] + '_' + master_df['DIR']
#master_df['PROV'] = ''
#master_df['CITY'] = ''
#master_df['DIR'] = ''
#master_df['Length_Km'] = master_df['Length_Km'].fillna(0.0)
#master_df['AADT'] = 0
#master_df['AADT_Year'] = 0
#master_df['POS_ALONG'] = master_df['POS_ALONG'].fillna(0)

#Set the coordinate referencec system (CRS) of the GeoDataFrame to that of the first shapefile in the list
master_gdf = gpd.GeoDataFrame(master_df,crs=shapefiles_list[0].crs)
master_gdf.set_index('POS_ALONG',inplace=True)

#Save the GeoDataFrame to a shapefile
master_gdf.to_file('C:/Users/Jason/Downloads/sorted/SortedCorridors/2022Q1_Corridors8.shp',driver='ESRI Shapefile')
print("finish merging the files") 




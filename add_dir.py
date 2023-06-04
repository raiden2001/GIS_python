import geopandas as gpd
import os
import numpy as np


#Load the shapefiles
gatineau1 = gpd.read_file('C:/Users/Jason/Downloads/trois/Trois-Rivieres2.shp')
gatineau1_final = gpd.read_file('C:/Users/Jason/Downloads/New folder (35)/new_roads.shp')


#montreal5_final = montreal5_final.assign(DIR=np.nan)
#montreal5_final['DIR'] = np.nan
#Copy the 'DIR' field frm the montreal5 to 'montreal5_final' with default value of np.nan
#montreal5_final = montreal5_final.assign(DIR=np.nan)

#Creates the column in the montreal5_final with the default values from the montreal3
gatineau1_final['DIR'] = np.nan

#loop through the rows of montreal5
for index,row in gatineau1.iterrows():
    #Get coresponding row in montreal5_final
    final_row = gatineau1_final.loc[gatineau1_final['LINK_ID'] == row['LINK_ID']]
    print(final_row)
    #Assign the value of 'DIR' form montreal5 to montreal5_final
    gatineau1_final.at[final_row.index,'DIR'] = row['DIR']
    print(gatineau1_final)
#montreal5_final = montreal5_final.merge(montreal5[['LINK_ID','DIR']],on='LINK_ID',how='left')



#save to update
#output_file = os.path.join(os.path.dirname(shapefile_dir),'

gatineau1_final.to_file('C:/Users/Jason/Downloads/New folder (35)/Trois_update.shp',driver='ESRI Shapefile')
print("copt success")

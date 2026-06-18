# check spatial references 

import arcpy
import random

arcpy.env.overwriteOutput = True



#############Configuration ############Get list 
premise_layer = "Premise_History"   # Make sure this matches exactly
zone_layer = "JointMunicipal1"
output_layer = r"C:\Users\raide\OneDrive\Documents\ArcGIS\Projects\NRPS_CAD_GIS\NRPS_CAD_GIS.gdb\Zone_Pop_crime_fire6"

print("step 1: Performing SPatial JOin to assign zones...")

###################3 Spatial Join: Assign Zone Name based on location ====


for lyr in [premise_layer, zone_layer, output_layer]:
    desc = arcpy.Describe(lyr)
    sr = desc.spatialReference
    print(f"{lyr}")
    print(f"  Spatial Reference: {sr.name}")
    print(f"  Factory Code: {sr.factoryCode}")
    print(f"{output_layer}")
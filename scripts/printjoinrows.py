import arcpy
import random

arcpy.env.overwriteOutput = True

arcpy.env.overwriteOutput = True

print(f"Populating sample data...")
print(f"Populating 25 premises with each zone information")

#############Configuration ############Get list 
premise_layer = "Premise_History"   # Make sure this matches exactly
zone_layer = "JointMunicipal1"
output_layer = r"C:\Users\raide\OneDrive\Documents\ArcGIS\Projects\NRPS_CAD_GIS\NRPS_CAD_GIS.gdb\Zone_Pop_crime_fire6"

print("step 1: Performing SPatial JOin to assign zones...")

###################3 Spatial Join: Assign Zone Name based on location ====

try:
    
    arcpy.analysis.SpatialJoin(
        target_features=premise_layer,
        join_features=zone_layer,
        out_feature_class=output_layer,
        join_operation="JOIN_ONE_TO_ONE",
        join_type="KEEP_ALL",
        match_option="INTERSECT" # change within to intersect 

    )

    print("Spatial JOin completed! Zones assigned and completed successfully")

#Step 2: 
except Exception as e:
    print(f" Spatial join fail: {e}")
    print(f"Tip: Make sure both layers are added to the map and have valid geometry")
    exit()


#-------- step 2: Populate realistic sample data ===
print("step : fixing Name and Mayer fields +  Populating sample data...")

#Get the list of fields in output_layer
fields = [f.name for f in arcpy.ListFields(output_layer)]

#Determine the correct joined field names (ArcGIS adds _1)
name_field = "Name_1" if "Name_1" in fields else ("Name" if "Name" in fields else None)
mayor_field = "Mayor_1" if "Mayor_1" in fields else ("Mayor" if "Mayor" in fields else None)


sample_fields = ["Join_Count", name_field, mayor_field]

print("\nJoined examples only:")

printed = 0

with arcpy.da.SearchCursor(output_layer, sample_fields) as cursor:
    for row in cursor:
        if row[0] and row[0] > 0:
            print(row)
            printed += 1

        if printed >= 10:
            break
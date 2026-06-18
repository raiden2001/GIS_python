import arcpy
import random

arcpy.env.overwriteOutput = True

print(f"Populating sample data...")
print(f"Populating 25 premises with each zone information")

#############Configuration ############Get list 
premise_layer = "Premise_History"   # Make sure this matches exactly
zone_layer = "JointMunicipal1"
output_layer = r"C:\Users\raide\OneDrive\Documents\ArcGIS\Projects\NRPS_CAD_GIS\NRPS_CAD_GIS.gdb\Zone_Pop_crime_fire7"

print("step 1: Performing SPatial JOin to assign zones...")

###################3 Spatial Join: Assign Zone Name based on location ====

try:
    
    arcpy.analysis.SpatialJoin(
        target_features=premise_layer,
        join_features=zone_layer,
        out_feature_class=output_layer,
        join_operation="JOIN_ONE_TO_ONE",
        join_type="KEEP_ALL",
        match_option="INTERSECT", # change within to intersect 
        search_radius = "100 meters"

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

if not name_field:
    raise RuntimeError("Could not find source zone field: Name or Name_1")
if not mayor_field:
    raise RuntimeError("Could not find source zone field: Mayor or Mayor_1")

#Use of different final field names so they do not conflict with source Name/Mayor 
final_name_field = "Zone_Name"
final_mayor_field = "Zone_Mayor"

#Add fields if they dont exists
for fld in [final_mayor_field,final_name_field]:
    if fld not in fields:
        arcpy.management.AddField(output_layer,fld,"Text",field_length=100)

cursor_fields = [
"Join_Count",
"Risk_Level",
"HAS_ALARM",
"WARNING_FLAG",
"Premise_Type",
"ADDRESS",
final_name_field,
final_mayor_field,
name_field,
mayor_field


]      

idx = {f: i for i, f in enumerate(cursor_fields)} # puts the rwos of columns in each columns of index in each field 



 #final_mayor_field,final_name_field,name_field,mayor_field      ["Risk_Level", "HAS_ALARM", "WARNING_FLAG", "Premise_Type","ADDRESS","Name"] 

with arcpy.da.UpdateCursor(output_layer,cursor_fields) as cursor:
    count=0
    for row in cursor:
        
        if not row[idx["Join_Count"]] or row[idx["Join_Count"]] == 0 : # this means telling line does not skip 
            continue
        
        if count >= 137:
           break
        
        has_alarm = random.choice([0, 1])
        premise_type = random.choice(["Residential", "Commercial", "Industrial","school"])
        #Get zone anme fomr spatial join (adjust index if needed)
        #zone_name = row[0] if row[0] else "Unknown Zone" # this will comes from the SPatial JOin(the 'Name' field)
        
        if has_alarm == 1:
            risk = "High"
            warning = 1 
        else:
            risk = random.choice(["Low","Medium","High"])
            #warning = 1 if risk == "Low"  else 0
            warning = 1 if risk in ('Medium', 'High') else 0
            
        street_number = random.randint(100,9999)
        fake_address = f"{street_number}{random.choice(['Main st','Lake Rd','Queen st','Niagara Blvd'])}"
        
        #Copy Name and Mayor from joined zone ===
        #zone_name = row[7] if row[7] else "Unknown" #index name_field
        #mayor_name = row[8] if row[8] else "unknown" # index of mayor_field
        zone_name = row[idx[name_field]] if row[idx[name_field]] else "Unknown"
        zone_mayor = row[idx[mayor_field]] if row[idx[mayor_field]] else "Unknown"
        
        
        #final_name_field = "zone_name"
        #final_mayor_field = "zone_mayor"
        #Assignes a zone within the premise within 100 meters for the polygon of a zone 
        
        row[idx["Risk_Level"]] = risk
        row[idx["HAS_ALARM"]] = has_alarm
        row[idx["WARNING_FLAG"]] = warning
        row[idx["Premise_Type"]] = premise_type
        row[idx["ADDRESS"]] = fake_address
        row[idx[final_name_field]] = zone_name
        row[idx[final_mayor_field]] = zone_mayor
        
        #row[idx[[0]] = risk
        #row[1] = has_alarm
        #row[2] = warning
        #row[3] = premise_type
        #row[4] = fake_address 
        #row[5] = zone_name  # < -- now actually write
        #row[6] = mayor_name 
        #row[5] is comment out due to the "name" field exisit
        
        cursor.updateRow(row)
        count += 1
        
print(f"Successfully procesed {count} premises with zone information")

        
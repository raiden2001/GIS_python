import arcpy
import random



print("\n=== Spatial JOin Dianostics ===")

#############Configuration ############Get list 
premise_layer = "Premise_History"   # Make sure this matches exactly
zone_layer = "JointMunicipal1"
output_layer = r"C:\Users\raide\OneDrive\Documents\ArcGIS\Projects\NRPS_CAD_GIS\NRPS_CAD_GIS.gdb\Zone_Pop_crime_fire6"

print("step 1: Performing SPatial join zones...")

total_rows  = int(arcpy.management.GetCount(output_layer)[0])
print(f"Total output rows: {total_rows}")

join_count_positive = 0 
join_count_zero = 0 

with arcpy.da.SearchCursor(output_layer,["Join_Count"]) as cursor:
    for row in cursor:
        if row[0] and row[0] > 0:
            join_count_positive += 1
        else:
            join_count_zero += 1
            
print(f"Rows with join_Count > 0: {join_count_positive}")
print(f"Rows with Join_Count = 0 or NULL:  {join_count_zero}")
 
fields = [f.name for f in arcpy.ListFields(output_layer)]
 
print(f"\n Fields found:")
for check_field in ["Name","Name_1","Mayor","Mayor_1","Join_Count"]:
    print(f" {check_field}: {check_field in fields}")
     
name_field = "Name_1" if "Name_1" in fields else ("Name" if "Name" in fields else None)
mayor_field = "Mayor_1" if "Mayor_1" in fields else ("Mayor" if "Mayor" in fields else None)

print("\n Using source fields:")
print(f" name_field:{name_field}")
print(f" mayor_field: {mayor_field}")

if not name_field:
    raise RuntimeError("Could not find source zone field: Name or Name_1")
    
if not mayor_field:
    raise RuntimeError("COuld not find source mayor field: Mayor or Mayor_1")
    
    

if name_field and mayor_field:
    print("\nFirst 10 sample joined rows")
    sample_fields = ["Join_Count",name_field,mayor_field]
    
    with arcpy.da.SearchCursor(output_layer,sample_fields) as cursor:
        for i,row, in enumerate(cursor):
            print(row)
            if i >=9:
                break
import arcpy

arcpy.env.overwriteOutput = True

# === Define your paths (adjust if needed) ===
target_features = "Roads"     # The roads layer
join_features = "Municipal_Boundaries" #The polygon layer with district names
output_fc = r"C:\Users\raide\OneDrive\Documents\ArcGIS\Projects\NRPS_CAD_GIS\NRPS_CAD_GIS.gdb\jointMunicipal1"

#Create the filed field map and field mapping objects 
field_mappings = arcpy.FieldMappings()

#Add_mapping all fields from the Target (roads)
field_mappings.addTable(target_features)

#Add_mappings.addTable(join_feature)
field_mappings.addTable(join_features)

#keep the fields : "Name" 
fields_to_keep = ["Name","CSD_CODE","Mayor","LandArea"]

#Loop through field maps(not fields)
for i in reversed(range(field_mappings.fieldCount)):
    field_map = field_mappings.getFieldMap(i)
    out_field = field_map.outputField
    
    if out_field.name not in fields_to_keep and out_field.name not in ["OBJECTID","Shape"]:
        field_mappings.removeFieldMap(i)
    
#Loop through
#for field in field_mappings.fields:
 #   if field.name not in fields_to_keep and field.origin == join_features:
  #      field_mappings.removeFieldMap(field_mappings.findFieldMapIndex(field.name))
        
print("Field mappings created. Running Spatial Join")  

#Run the spatail Join
arcpy.analysis.SpatialJoin(
    target_features = target_features,
    join_features=join_features,
    out_feature_class=output_fc,
    join_operation = "JOIN_ONE_TO_ONE",
    join_type = "KEEP_ALL",
    field_mapping=field_mappings, # key part 
    match_option = "WITHIN"    

) 

print("Spatial JOin is complete successful")
print(f"Output saved to:{output_fc}")

#print("Fields in output:")
#for f in arcpy.ListFields(output_fc):
 #   print(f.name)
    
#open the attribute table automatically (optional)
#arcpy.management.AddJoin(output_fc,"OBJECTID",output_fc,"OBJECTID") #justto watch to refresh
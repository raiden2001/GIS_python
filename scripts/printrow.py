import zone_crime_fire

import arcpy

output_layer = r"C:\Users\raide\OneDrive\Documents\ArcGIS\Projects\NRPS_CAD_GIS\NRPS_CAD_GIS.gdb\Zone_Pop_crime_fire7"

fields = [f.name for f in arcpy.ListFields(output_layer)]

name_field = "Name_1" if "Name_1" in fields else ("Name" if "Name" in fields else None)
mayor_field = "Mayor_1" if "Mayor_1" in fields else ("Mayor" if "Mayor" in fields else None)

if not name_field:
    raise RuntimeError("Could not find Name or Name_1")

if not mayor_field:
    raise RuntimeError("Could not find Mayor or Mayor_1")

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

print(f"\nPrinted {printed} joined rows.")
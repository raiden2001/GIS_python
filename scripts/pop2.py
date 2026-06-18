import arcpy
import random

layer = "Premise_History"   # Make sure this matches exactly

print("Populating sample data...")

with arcpy.da.UpdateCursor(layer, ["Risk_Level", "HAS_ALARM", "WARNING_FLAG", "Premise_Type"]) as cursor:
    for row in cursor:
        has_alarm = random.choice([0, 1])
        premise_type = random.choice(["Residential", "Commercial", "Industrial"])
        
        if not has_alarm == 1:
            risk = "High"
            warning = 1
        else:
            risk = random.choice(["Low", "Medium", "High"])
            warning = 1 if risk == "High" else 0
        
        row[0] = risk
        row[1] = has_alarm
        row[2] = warning
        row[3] = premise_type
        cursor.updateRow(row)

print(f"✅ Sample data populated successfully in{layer}!")
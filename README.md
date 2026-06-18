# GIS Python

ArcGIS Pro / ArcPy scripts and sample outputs for the NRPS CAD GIS workflow.

## Contents

- `scripts/` - Python scripts copied from the ArcGIS geodatabase workspace.
- `outputs/printed_10_joined_preview.png` - Preview image for the first 10 joined output features.
- `outputs/printed_10_joined_shapefile/` - Shapefile components for the first 10 joined output features.

## ArcGIS Pro Python

These scripts require ArcGIS Pro's Python environment with `arcpy` available.

Example local interpreter used during testing:

```powershell
C:\Users\raide\AppData\Local\Programs\ArcGIS\Pro\bin\Python\envs\arcgis_clean\python.exe
```

Example run:

```powershell
& "C:\Users\raide\AppData\Local\Programs\ArcGIS\Pro\bin\Python\envs\arcgis_clean\python.exe" scripts\zone_crime_fire.py
```

The full file geodatabase is not committed here. Use Git LFS or cloud storage for large `.gdb` data if needed.

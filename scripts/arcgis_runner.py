"""Simple runner for loading ArcGIS-exported DEM terrain into Genesis.

This script mirrors the terrain preview flow in ``raster.py`` so the
repository has a stable ``arcgis_runner.py`` entrypoint.
"""

from pathlib import Path

import genesis as gs
import rasterio


def main() -> None:
    gs.init()

    terrain_tif = Path("assets/terrain.tif")
    if not terrain_tif.exists():
        raise FileNotFoundError(
            "Missing DEM file at assets/terrain.tif. "
            "Export terrain from ArcGIS Pro and place it there."
        )

    with rasterio.open(terrain_tif) as src:
        height_field = src.read(1)
        downsample_factor = 16
        height_field = height_field[::downsample_factor, ::downsample_factor]
        horizontal_scale = src.res[0] * downsample_factor

    scene = gs.Scene(show_viewer=True)
    scene.add_entity(
        gs.morphs.Terrain(
            height_field=height_field,
            horizontal_scale=horizontal_scale,
            vertical_scale=1.0,
        )
    )
    scene.add_entity(gs.morphs.Sphere(pos=(10, 15, 15), radius=1))
    scene.build(n_envs=1000)

    for _ in range(1000):
        scene.step()


if __name__ == "__main__":
    main()

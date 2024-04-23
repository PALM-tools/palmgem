follow [preprocesor documentations](../../docs/run_prerocessor.md)
1) from OpenStreetMap download Prague. Use user [QGIS procedure](../../docs/user_preproces.md) to clip unnecessary polygons.
2) from UrbanAtlas download Prague. Use user [QGIS procedure](../../docs/user_preproces.md) to clip unnecessary polygons.
3) from EU-DEM website download zip file. Use user [QGIS procedure](../../docs/user_preproces.md) to clip unnecessary raster tiles.
4) from UrbanAtlas building height download Prague. Use user [QGIS procedure](../../docs/user_preproces.md) to clip unnecessary raster tiles.
5) Configure importing shell script: [import_pgsql_prague.sh](import_pgsql_prague.sh). Use correct path for shapefiles and raster files. 
6) Configure user credential and password, database name and connection in [default_share.yaml](../../config/default_share.yaml) or directly in preprocessor config [prague.yaml](prague.yaml). Finish [prague.yaml](prague.yaml) configuration, change if necessary input domain extent. 
7) Run preprocessing script: "python3 main_urban_atlas_preprocess.py -c examples/prague/prague.yaml". See log in log folder.
8) Using QGIS check imported dataset, follow [visualization](../../docs/visuallization.md)
9) Check [default configuration](../../config/default_config.yaml) and example [Prague config](prague_palm.yaml). Run "python3 main_palm_static.py -c examples/prague/prague_palm.yaml". See log in logs. 
10) Using QGIS check created grid, follow [visualization](../../docs/visuallization.md). Check created png files in visual_check/prague folder. See created static driver in output folder.
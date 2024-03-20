# Run preprocessor
## Data Downloading
### OpenStreetMaps (OSM) download
OSM datasets can be downloaded from OSM geofabrik website https://download.geofabrik.de/europe.html. Download desired layer as .shp.zip file. Extract afterwards and copy to your working dir. For simplicity and algorithm speed up, there is an option to use some GIS software to crop dateset before importing into PostgreSQL database. Please see [crop shapefile in QGIS](user_preproces.md)

### Urban Atlas (UA) Landcover
Create a free account in UA website https://land.copernicus.eu/. From UA website https://land.copernicus.eu/local/urban-atlas/urban-atlas-2018?tab=download find Prague (Praha) and download data. Extract the files after download. Extract zip file inside. In folder data, process gpkg file with ogr2ogr command in terminal: ogr2ogr -f "ESRI Shapefile" . ${file_name}.gpkg. Copy files to your working dir. For simplicity and algorithm speed up, there is an option to use some GIS software to crop dateset before importing into PostgreSQL database. Please see [crop shapefile in QGIS](user_preproces.md)

### EU-DEM 
FROM EU-DEM website download Open Topo Data. Follow instruction at the website. Extract downloaded dataset, copy .tif file you require for your domain, please use N***E*** specification, to you working dir. For simplicity and algorithm speed up, there is an option to use some GIS software to clip dateset before importing into PostgreSQL database. Please see [crop raster in QGIS](user_preproces.md)

### Urban Atlas Building height
From UA website https://land.copernicus.eu/local/urban-atlas/building-height-2012?tab=download (with previous created free account) download buildings heights for your specified city. Extract after download, extract zip file inside, copy files from Dataset to your working dir.For simplicity and algorithm speed up, there is an option to use some GIS software to clip dateset before importing into PostgreSQL database. Please see [crop raster in QGIS](user_preproces.md)

### Data import to PostgreSQL database
Use prepared shell script in import_to_sql folder import_pgsql.sh to import downloaded data into postgresql database (please check names and paths of the importing files). After successful import you can check imported geodata using QGIS, detailed instruction can be found [here](../../docs/visuallization.md).

## Configuration and running preprocessor
Check default configurations in config folder. Create you configuration in config folder, name_of_case.yaml. Fill in your PostgreSQL database connection, you connection username and password, name of input schema. Specify if you want to crop to user defined extend (see example_preproc.yaml). Then run preprocessing python script **python3 main_urban_atlas_preprocess.py -c brno.yaml &*. See log in logs, the name of log is name+scenario.log. You can visualize processed tables using QGIS, detailed instruction can be found [here](../../docs/visuallization.md).
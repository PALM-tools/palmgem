# Brno testcase
TODO: crop the text only do names of websites and name of datasets.
## Data Downloading
### OpenStreetMaps (OSM) download
From OSM geofabrik website https://download.geofabrik.de/europe.html download layer of Czech Republic as .shp.zip file. Extract the file after download. Copy gis_osm_buildings_a_free_1.* files to /examples/brno/downloads. 

### Urban Atlas (UA) Landcover
Create a free account in UA website https://land.copernicus.eu/. From UA website https://land.copernicus.eu/local/urban-atlas/urban-atlas-2018?tab=download find Brno and download data. Extract the files after download. Extract zip file inside. In folder data, process gpkg file with ogr2ogr command in terminal: ogr2ogr -f "ESRI Shapefile" . ${file_name}.gpkg. Copy CZ002L2_BRNO_UA2018.* files to /examples/brno/downloads. 

### Urban Atlas Digital Elevation Model (DEM)
From UA website https://land.copernicus.eu/imagery-in-situ/eu-dem/eu-dem-v1.1 (with previous creating of free account) inspect in interactive window area of interest. E40N20 for Brno location and download the data. Extract after download, extract zip file inside, copy files inside to /examples/brno/downloads.

### Urban Atlas Building height
From UA website https://land.copernicus.eu/local/urban-atlas/building-height-2012?tab=download (with previous creating of free account) download buildings heights for Brno. Extract after download, extract zip file inside, copy files from Dataset to /examples/brno/downloads.

### Data import to PostgreSQL database
Use prepared shell script in import_to_sql folder import_pgsql_brno.sh to import downloaded data into postgresql database (please check names and paths of the importing files). After successful import you can check imported geodata using QGIS, detailed instruction can be found [here](../../docs/visuallization.md).

## Configuration and running preprocessor
Check configuration for Brno testcase [default_config](../../config/default_config_preproc.yaml) and [brno_config](../../config/brno.yaml). Run preprocessing python script python3 main_urban_atlas_preprocess.py -c brno.yaml &. See log in logs. You can visualize processed tables using QGIS, detailed instruction can be found [here](../../docs/visuallization.md).

## Configuration and running PALM preprocessor
TODO: note about loading all psl function to server, or make it automatic?
Check [configuration](../../config/brno_palm.yaml) for PALM static driver preprocessor for selected Brno domain. Run PALM static driver with configuration python3 main_palm_static.py -c brno_palm.yaml. In logs see progress and in visual_check folder see output png figures. Or use QGIS to visualized PALM grid with properties, detailed instruction can be found [here](../../docs/visuallization.md).

## Example configuration of PALM simulation? Some comparison figures from manuscript?
Show we do that? 
# Tutorials
Please find an examples how PALM-GEM can process various Level-of-Details (LODs) using differently detailed inputs. \
In our project we define several layers of LOD. \
In each case, fill in your credentials and copy .yaml config into root config file. Import into dtb, run the scenario, see created .nc file.
## LOD0
Most simple, only landcover is provided (with/without buidling information). Terrain is defined as flat (must be specified in config). \
Building height is taken as default height based on building type.
## LOD1
Most often used case. Landcover in combination with terrain height (and building height) is processed into PALM input netcdf file. \
This case is based on using open source data from Copernicus and OSM.
## LOD2 
Detailed parametrization of landcover, building roof, walls, (trees). This option utilizes almost all PALM configurations to obtain as detailed as possible simulation.
## OTHER OPTIONAL / ADDITIONAL
### LAI / CANOPY HEIGHT
In lod1_lai, we combine downloaded Leaf Area Index and Canopy Height into Leaf Area Density - Palm Input variable.
### Global Building Atlas
GBA dataset with building height from Prague. Only a snipped, not working example. Replace building height in lod1 and run the case. \
TOBE added.
### Impervious
An additional raster variable that can be combined with coarse resolution - open source copernicus data to narrow down and differentiate pavement/building/vegetation grids. \
A code snipped, add impervious impervious.tif into lod1 import script and run the case.
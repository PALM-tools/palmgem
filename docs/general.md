# Logging system
PALM-GEM provides standard logging tool. There are 6 different log levels separately for python and SQL outputs. 
* 1 - EXTRA VERBOSE - prints some sql queries and detailed outputs. For development purposes.
* 2 - VERBOSE - Detailed log, log inside loops, etc.
* 3 - DEBUG - Print lower detailed logs.
* 4 - PROGRESS - Print which part of process is currently being processed.
* 5 - WARNING - Print Warning, if something non critical is not missing.
* 6 - ERROR - Print Error and broke the run.

All standard and error outputs are redirected into log file in log folder. The name of log file is **name**+**scenario**.**log**

# Configuration of PostgreSQL and PostGIS, loading special SQL function
Before the first run of python scripts, you have to configure the PostgreSQL and PostGIS server and  upload PALM-GEM SQL functions into your PostgreSQL database. The procedure is described in [instalation guide](install.md).

# Creation of final landcover - separation into smaller segments
Due to optimization of spatial operation, we developed a segmentation of larger polygons into smaller using fishnet. Configuration of fishned size can be found in [general configuration](configuration_docs.md).

# Types
To standard PALM surface types, we use integer value to identify PALM type. For particular surface categories the following ranges of values are used: vegetation_type type=(100, 199), pavement_type=(200, 299), water_type=(300, 399), building_type=(900,999).   

# mt table
UA = UrbanAtlas \
OSM = OpenStreetMap

| **UrbanAtlas type** | **PALM type, if not OSM** | **PALM type, if OSM** |
|:------------------|:---------------|:---------------|
| 11100 | 203 | 901 | 
| 11210 | 203 | 902 |
| 11220 | 203 | 903 |
| 11230 | 203 | 902 |
| 11240 | 203 | 903 |
| 11300 | 203 | 903 |
| 12100 | 203 | 906 |
| 12210 | 201 | 906 |
| 12220 | 202 | 906 |
| 12230 | 209 | 906 |
| 12300 | 203 | 906 |
| 12400 | 103 | 906 |
| 13100 | 101 | 906 |
| 13300 | 101 | 906 |
| 13400 | 108 | 906 |
| 14100 | 118 | 906 |
| 14200 | 103 | 906 |
| 21000 | 102 | 906 |
| 22000 | 101 | 906 |
| 23000 | 103 | 906 |
| 24000 | 102 | 906 |
| 25000 | 102 | 906 |
| 31000 | 117 | 906 |
| 32000 | 110 | 906 |
| 33000 | 112 | 906 |
| 40000 | 110 | 906 |
| 50000 | 301 | 906 |

# TODO work in progress
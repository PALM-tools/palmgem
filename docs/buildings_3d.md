# Buildings 3d description
In Palm-Gem, 3d buildings including bridges, overhanging structures and passages is included. Placing of 3D buildings is done according to Palm official documentation. In special case with bridges, bottom of bridge grid cells are placed above terrain with height defined by extras raster. Top of bridge is placed above bottom with configurable height build_3d.bridge_width. In case of overhanging structures and passages, the buildings 3d are firstly created as 2d, subsequently in places where extras_shp are defined, the height of building's bottom is calculated by joining grid with raster extras. \
With this approach, only structures between ground and building's bottom are available. In project it is prepared to include multiple empty zones in buildings in future (e.g., include cascade of balconies).\
This module requires extras raster and extras_shp polygons tables, all parameters are include in tables bellow. 

| Attribute | Type  | Values | Desription |
|:----------|:------|:-------|:----------------------|
| gid       | int   | >1     | unique identified for each polygon |
| type      | int   | > 900  | type of building |
| typeu     | int   | > 1    | type of polygon under bridge |
| typed     | int   | > 1    | type of polygon above bridge | 
| class3d   | char  | char   | type of structure, available 'bridge', 'passage', 'overhang' |

Raster table: 

| Attribute | Type  | Values | Desription |
|:----------|:------|:-------|:----------------------|
| gid       | int   | >1     | unique identified for each raster |
| rast      | real  | > 1.0  | height between ground and first building grid |


During processing of 3d structure, checking routines are used. E.g., gaps with size of 1 grid are omitted. \
In current version, construction of bridges are done to create as close as possible structures. But main problems are mainly in data availability and resolution. Out suggestion is to use 3D strucure only in < 3 m resolutions. 
Note: Cut-Cell-Topograghy is not compatible with buildings 3D, resulting to turn off buildings 3D feature in case do_cct is True.
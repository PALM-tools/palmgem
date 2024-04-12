#Configuration 
## Shared configuration
| Configuration Item | Default value | Describtion |
|:------------------|:---------------|:---------------|
| **pg_host** |  localhost | name of database host |
| **pg_port** |  5432 | port number |
| **pg_user** |  unknown | name user in postgresql database |
| **pg_password** |  unknown | password for pg_user in postgresql database |
| **database** |  palm_static | name of database schema |
| **pg_owner** |  palm | name of general group in postgresql database |

| **log.level** | 1 | Level of logging. See [docs](general.md) |

| **srid** | 32633 | SRID reference system for final schema |
| **srid_input** | 32633 | SRID reference system for final schema |
| **srid_palm** | 32633 | SRID reference system for final schema |
| **dem_srid** | 3035 | SRID reference system for imported DEM |
| **srid_utm** | 32633 | SRID reference system for UTM |
| **srid_wgs84** | 4326 | SRID reference system for WGS84 |



## Preprocessor tool
| Configuration Item | Default value | Describtion |
|:------------------|:---------------|:---------------|
| **domain.name** | input_general | name of case |
| **domain.scenario** | "" | name of case scenario, **domain.name+domain.scerio** create name of postgresql schema |
| **domain.xl** | -1e10 | x lower coordinate of cropping domain |
| **domain.xh** | +1e10 | x top coordinate of cropping domain |
| **domain.yl** | -1e10 | y lower coordinate of cropping domain |
| **domain.yh** | +1e10 | y top coordinate of cropping domain |
| **domain.crop_domain** | False | Option that will crop originally uploaded datasets |
| **domain.fill_boundary** | True | Filling missing places, e.g., sea, ocean |
| **fishnet.nx** | 200 | number of grid cell in x-direction in fishnet |
| **fishnet.ny** | 200 | number of grid cell in y-direction in fishnet, resolution is (higher corner - lower corner) / fishnet.ny |
| **tables.im_landcover_or** | imported_landcover | table name in case schema |
| **tables.im_landcover** | landcover | table name in case schema |
| **tables.streetmaps_or** | imported_streetmap | table name in case schema |
| **tables.streetmaps** | streetmap | table name in case schema |
| **tables.buildings_or** | imported_buildings | table name in case schema |
| **tables.buildings** | buildings | table name in case schema |
| **tables.dem_or** | imported_dem | table name in case schema |
| **tables.dem** | dem | table name in case schema |
| **tables.fishnet** | fishnet |  table name in case schema |
| **rast_tr.tilewidth** | 32 | Configuration of raster transformation |
| **rast_tr.tileheight** |  32|  |
| **rast_tr.scalex** | 25 |  |
| **rast_tr.scaley** | -25 |  |
| **rast_tr.skewx** | 0 |  |
| **rast_tr.skewy** | 0 |  |
| **max_stl_area** | 50.0 | area of urban object that are not considered in intersection |
| **max_fishnet_split_area** | Maximal area of polygons that are not intersected by fishnet |  |
| **mt_default** | 103 | Default type that fills places in case of **domain.fill_boundary**. See Types in [docs](general.md)  |
| **clean_up** | False | Delete unused table when they are no longer needed. Use carefully.|
| **mt** | | Dictionary for transformation from UrbanAtlas and OpenStreetMaps into PALM types. Table can be found in [mt table](general.md) |
| **** | | |



## Static driver preprocessor
| Configuration Item                      | Default value                      | Describtion |
|:----------------------------------------|:-----------------------------------|:--------------------|
| **input_schema**                        | inputs_prague                      | Name of PostgreSQL input schema |
| **domain.name**                         | case_name                          | Name of case |
| **domain.scenario**                     | scenario_00                        | Scenario name, schema name is name + scenario |
| **domain.parent_domain_schema**         | ""                                 | Name of parent domain schema, in case of nesting setup |
| **domain.dx**                           | 10.0                               | Resolution in x-direction |
| **domain.dy**                           | 10.0                               | Resolution in y-direction |
| **domain.dz**                           | 10.0                               | Resolution in z-direction |
| **domain.nx**                           | 128                                | Number of grid cell in x-direction |
| **domain.ny**                           | 128                                | Number of grid cell in x-direction  |
| **domain.cent_x**                       | 459180                             | Center of domain in palm_srid, x-coordinate |
| **domain.cent_y**                       | 5547800                            | Center of domain in palm_srid, y-coordinate |
| **default_height.0**                    | 6.0                                | Default height of building type 0 (According to PALM PIDS), in case of buildings height are missing at building location. |
| **default_height.1**                    | 6.0                                | Default height of building type 1 (According to PALM PIDS), in case of buildings height are missing at building location. |
| **default_height.2**                    | 6.0                                | Default height of building type 2 (According to PALM PIDS), in case of buildings height are missing at building location. |
| **default_height.3**                    | 6.0                                | Default height of building type 3 (According to PALM PIDS), in case of buildings height are missing at building location. |
| **default_height.4**                    | 12.0                               | Default height of building type 4 (According to PALM PIDS), in case of buildings height are missing at building location. |
| **default_height.5**                    | 12.0                               | Default height of building type 5 (According to PALM PIDS), in case of buildings height are missing at building location. |
| **default_height.6**                    | 12.0                               | Default height of building type 6 (According to PALM PIDS), in case of buildings height are missing at building location. |
| **visual_check.enabled**                | True                               | Create .png picture of final static driver file |
| **visual_check.show_plots**             | False                              | Interactivelly show plots during generation |
| **tables.grid**                         | grid                               | Name of grid table |
| **tables.dem**                          | dem                                | Name of DEM table |
| **tables.buildings_height**             | buildings                          | Name of buildings height raster |
| **tables.buildings_grid**               | buildings_grid                     | Name of buildings grid table |
| **tables.landcover**                    | landcover                          | Name of landcover table |
| **tables.buildings_offset**             | buildings_offset                   | Name of buildings offset, for filtration routines |
| **tables.aspect**                       | buildings_aspect_deg               | Name of Buldings's aspect table |
| **tables.slope**                        | buildings_slope_deg                | Name of Buildings's slope table |
| **tables.slanted_wall**                 | slanted_wall                       | Name of slanted wall table |
| **tables.slanted_roof**                 | slanted_roof                       | Name of slanted roof table |
| **tables.slanted_terr_wall**            | slanted_terr_wall                  | Name of slanted terrain wall table |
| **tables.slanted_terrain**              | slanted_terrain                    | Name of slanted terrain table |
| **tables.slanted_wall_gridded_temp**    | slanted_wall_gridded_temp          | Name of slanted temporal gridded wall table |
| **tables.slanted_wall_gridded**         | slanted_wall_gridded               | Name of slanted gridded wall table |
| **tables.slanted_roof_gridded_temp**    | slanted_roof_gridded_temp          | Name of slanted temporal gridded roof table |
| **tables.slanted_roof_gridded**         | slanted_roof_gridded               | Name of slanted gridded roof table |
| **tables.slanted_terrain_gridded_temp** | slanted_terrain_gridded_temp       | Name of slanted temporal gridded terrain table |
| **tables.slanted_terrain_gridded**      | slanted_terrain_gridded            | Name of slanted gridded terrain table |
| **tables.walls_outer**                  | walls_outer                        | Name of outer walls table |
| **tables.slanted_faces**                | slanted_faces                      | Name of slanted faces table |
| **tables.vertices**                     | vertices                           | Name of vertices table |
| **tables.land_build**                   | landcover_buildings                | Name of landcover only buildings terrain |
| **tables.build_new**                    | buildings_new                      | Name of new buildings table |
| **tables.height_corrected**             | building_height_correct            | Name of buildings height correction's table |
| **tables.height_terr_corrected**        | terrain_height_correct             | Name of terrain's height correction's table |
| **tables.slanted_wall_points**          | slanted_wall_points                | Name of slanted points in wall table |
| **idx.landcover**                       | lid                                | Landcover unique index |
| **idx.walls**                           | wid                                | Walls unique index |
| **idx.roofs**                           | rid                                | Roofs unique index |
| **origin_time**                         | 1970-01-01 00:00:00                | Origin time of simulation |
| **maxbuildingdisance**                  | 5.0                                | Distance to search in neighbors in case of missing height in buildings grid cell |
| **cortyard_fill.apply**                 | False                              | Option that fill cortyard |
| **cortyard_fill.count**                 | 30                                 | Cortyards smaller that count grid cell are filled |
| **topo_fill_v2.apply**                  | False | PALM replicated filter for topo filtering |
| **topo_fill_v2.apply**                  | 9 | Number of adjacent gridcell tobe filtered |
| **force_lsm_only**                          | False                              | Option to force static driver into land surface only |
| **force_cyclic**                            | False                              | Option to force static driver to cyclic boundary condition |
| **force_cyclic_nc**                         | 4                                  | Number of boundary cell to apply cyclic boundary condition |
| **force_building_boundary**                 | True                               | Option to delete buildings that are close to domain boarder |
| **force_building_boundary_dist**            | 3                                  | Distance in grid points for **force_buildings_boundary**, real distance =* dx |
| **zlib**                                    | False                              | Option that compresses static driver netcdf file |
| **complevel**                               | 4                                  | Compression level |
| **fill_values.b**                           | -127                               | Fill value for netCDF file, byte |
| **fill_values.i**                           | -9999                              | Fill value for netCDF file, integer |
| **fill_values.i4**                          | -9999                              | Fill value for netCDF file, integer |
| **fill_values.f4**                          | -9999.0                            | Fill value for netCDF file, float |
| **fill_values.f8**                          | -9999.0                            | Fill value for netCDF file, float8 |
| **ncprops.acronym**                         | ICS                                | Accronim of Institution |
| **ncprops.author**                          | Jaroslav Resler \ resler@cs.cas.cz | Name of static driver file author |
| **ncprops.campaign**                        | ""                                 | Name of campaign |
| **ncprops.comment**                         | ""                                 | Comment |
| **ncprops.contact_person**                  | resler@cs.cas.cz                   | Contact author |
| **ncprops.data_content**                    | ""                                 |  |
| **ncprops.dependencies**                    | ""                                 |  |
| **ncprops.institution**                     | ""                                 |  |
| **ncprops.keywords**                        | ""                                 |  |
| **ncprops.location**                        | ""                                 |  |
| **ncprops.palm_version**                    | ""                                 |  |
| **ncprops.references**                      | ""                                 |  |
| **ncprops.rotation_angle**                  | ""                                 |  |
| **ncprops.site**                            | ""                                 |  |
| **ncprops.source**                          | UrbanAtlas, OpenStreetMaps         |  |
| **ncprops.version**                         | ""                                 |  |
| **ndims.nwater_pars**                       | 7                                  | Number of elements in water_pars |
| **ground.soil_type_default**                | 3                                  | Default soil type medium-fine in terms of porosity |
| **type_range.vegetation_min**               | 100                                | See [Types](general.md) |
| **type_range.vegetation_max**               | 199                                |  |
| **type_range.pavement_min**                 | 200                                |  |
| **type_range.pavement_max**                 | 299                                |  |
| **type_range.water_min**                    | 300                                |  |
| **type_range.water_max**                    | 399                                |  |
| **type_range.building_min**                 | 900                                |  |
| **type_range.building_max**                 | 999                                |  |
| **water_pars_temp.1**                       | 283.15                             | Water body temperature for type 1 |
| **water_pars_temp.2**                       | 283.15                             | Water body temperature for type 2 |
| **water_pars_temp.3**                       | 283.15                             | Water body temperature for type 3 |
| **water_pars_temp.4**                       | 283.15                             | Water body temperature for type 4 |
| **water_pars_temp.5**                       | 283.15                             | Water body temperature for type 5 |
| **water_pars_temp.6**                       | 283.15                             | Water body temperature for type 6 |
| **** |  |  |
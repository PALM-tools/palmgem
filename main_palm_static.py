#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright 2018-2024 Institute of Computer Science of the Czech Academy of
# Sciences, Prague, Czech Republic. Authors: Martin Bures, Jaroslav Resler.
#
# This file is part of PALM-GeM.
#
# PALM-GeM is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# PALM-GeM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# PALM-GeM. If not, see <https://www.gnu.org/licenses/>.

import psycopg2
from utils.palm_static_pg_lib import *
from utils.palm_static_pg_lib_cct import *
from utils.consistency_checks import *
from utils.tree_lib import *
from config.config import load_config, cfg
from argparse import ArgumentParser
import getpass
from config.logger import *

progress('Reading configuration')
argp = ArgumentParser(description=__doc__)
argp.add_argument('-c', '--config', help='configuration file')

# load config file
argv = argp.parse_args()
load_config(argv)
logging_level(cfg)

debug('Config file: {}', argv.config)
debug('Domain: {}', cfg.domain.name)
debug('Scenario: {}', cfg.domain.scenario)
debug('Input schema: {}', cfg.input_schema)
debug('Origin time: {}', cfg.origin_time)

# calculate dz in case dz is not supplied (dz<=0)
if cfg.domain.dz <= 0.0:
    debug("dz = 0.0: set dz = dx")
    cfg.domain._settings['dz'] = cfg.domain.dx
debug('Domain resolution (dx,dy,dz):({},{},{})', cfg.domain.dx, cfg.domain.dy, cfg.domain.dz)

# calculate origin of the domain (not multiplied)
origin_x = cfg.domain.cent_x - cfg.domain.nx * cfg.domain.dx / 2.0
origin_y = cfg.domain.cent_y - cfg.domain.ny * cfg.domain.dy / 2.0
cfg.domain._settings['origin_x'] = origin_x
cfg.domain._settings['origin_y'] = origin_y
debug('Domain centre (x,y):({},{}) and origin (x,y):({},{})', cfg.domain.cent_x, cfg.domain.cent_y, origin_x, origin_y)

# create connection to the postgresql server
if cfg.pg_password is None or cfg.pg_password == '':
    pg_password = getpass.getpass()
connection = psycopg2.connect(database=cfg.database, host=cfg.pg_host,
                              port=cfg.pg_port, user=cfg.pg_user, password=cfg.pg_password)
connection.set_client_encoding('UTF8')
cur = connection.cursor()
debug('Connection status: {}', connection.status)
# check and create case schema
progress('CREATE NEW CASE schema: {}', cfg.domain.case_schema)
sqltext = 'drop schema if exists "{}" cascade'.format(cfg.domain.case_schema)
cur.execute(sqltext)
sql_debug(connection)
connection.commit()

sqltext = 'CREATE schema if not exists "{}"'.format(cfg.domain.case_schema)
debug('Drop old schema, if existing')
cur.execute(sqltext)
sql_debug(connection)
connection.commit()

debug('Create new schema')
sqltext = 'ALTER SCHEMA "{}" OWNER TO {}'.format(cfg.domain.case_schema, cfg.pg_owner)
cur.execute(sqltext)
sql_debug(connection)
connection.commit()

progress('Creating grid')
create_grid(cfg, connection, cur)

progress('Calculate extent of the grid')
grid_ext = calculate_grid_extend(cfg, connection, cur)

progress('Coping and transforming data from inputs, vector data')
vtabs = copy_vectors_from_input(grid_ext, cfg, connection, cur)

progress('Coping and transforming data from inputs, raster data')
rtabs = copy_rasters_from_input(grid_ext, cfg, connection, cur)

check_surface_params(cfg, connection, cur)

check_buildings(cfg, connection, cur, rtabs, vtabs, grid_ext)

check_landcover(cfg)

if cfg.do_cct and cfg.tables.buildings_height in rtabs:
    debug('Cut cell topo, Modifying landcover, simplified buildings')
    preprocess_building_landcover(cfg, connection, cur)

progress('Calculate terrain height in grid')
calculate_terrain_height(cfg, connection, cur)

progress('Calculate origin_z and ori_min')
calculate_origin_z_oro_min(cfg, connection, cur)

progress('Connect landcover to grid')
connect_landcover_grid(cfg, connection, cur)

if cfg.cortyard_fill.apply:
    progress('Fill cortyards smaller than {} grid cells', cfg.cortyard_fill.count)
    fill_cortyard(cfg, connection, cur)

if cfg.fill_missing_holes:
    progress('Filling missing building holes')
    fill_missing_holes_in_grid(cfg, connection, cur)

progress('Processing building elevation model BEM')
connect_buildings_height(cfg, connection, cur)

if cfg.lod2:
    progress('Connection building roofs')
    connect_roofs(cfg, connection, cur)
    connect_walls(cfg, connection, cur, vtabs)

if cfg.force_cyclic:
    update_force_cyclic(cfg, connection, cur)

if cfg.has_trees:
    progress('Process trees')
    process_trees(cfg, connection, cur)

if cfg.canopy.using_lai:
    progress('Processing LAI and canopy height into LAD')
    process_lai(cfg, connection, cur)

progress('Done with preparation of geo inputs')
progress('Process data into netCDF4 static driver according to PALM Input Data Standard')

progress('Initialize netCDF4 static driver file')
ncfile = nc_create_file(cfg.domain.static_driver_file)

progress('Prepare domain origins (x,y,lat,lon)')
prepare_domain_extends(cfg, connection, cur)

progress('Initialize static driver global attributes')
nc_write_global_attributes(ncfile, cfg)

progress('Write Coordinate Reference System to netcdf4 file')
nc_write_crs(ncfile, cfg, connection, cur)

progress('Prepare netcdf4 dimensions')
create_dim_xy(ncfile, cfg, connection, cur)

progress('Writing terrain type')
write_terrain(ncfile, cfg, connection, cur)

if cfg.landcover.surface_fractions:
    progress('Write fractions')
    write_surface_fractions(ncfile, cfg, connection, cur)

progress('Writing pavement type')
write_pavements(ncfile, cfg, connection, cur)

progress('Writing water type')
write_water(ncfile, cfg, connection, cur)

progress('Writing vegetation type')
write_vegetation(ncfile, cfg, connection, cur)

progress('Writing soil type')
write_soil(ncfile, cfg, connection, cur)

progress('Writing buildings type')
if not cfg.slurb:
    write_buildings(ncfile, cfg, connection, cur)
    if cfg.lod2:
        test_building_insulation(cfg, connection, cur)
        # write_building_pars(ncfile, cfg, connection, cur)
        write_building_pars_depricated(ncfile, cfg, connection, cur)
        write_building_surface_pars(ncfile, cfg, connection, cur, vtabs)
        write_albedo_pars(ncfile, cfg, connection, cur)
else:
    write_mask_usm(ncfile, cfg, connection, cur)

if cfg.has_trees:
    progress('Writing lad, bad')
    write_trees_grid(ncfile, cfg, connection, cur)

if cfg.canopy.using_lai:
    progress('Writing lad, bad')
    write_lad_grid(ncfile, cfg, connection, cur)

if cfg.do_cct:
    slanted_surface_init(cfg, connection, cur)
    slanted_write_nc(ncfile, cfg, connection, cur)
    check_cct_consistency(ncfile, cfg, connection, cur)
    cct_continuity_check(ncfile, cfg)

progress('Check consistency of static driver file')
check_consistency(ncfile, cfg)

ncfile.close()
debug('File {} was closed', cfg.domain.static_driver_file)

if cfg.slurb:
    progress('Initialize netCDF4 slurb static driver')
    ncfile_slurb = nc_create_file(cfg.domain.slurb_driver_file)

    progress('Initilize SLURb driver')
    nc_write_global_attributes_slurb(ncfile_slurb, cfg)

    progress('Write crs for slurb file')
    nc_write_crs(ncfile_slurb, cfg, connection, cur)

    progress('Write SLURb dimensions')
    create_slurb_dims(ncfile_slurb, cfg, connection, cur)

    progress('Now create driver for slurb')
    create_slurb_vars(ncfile_slurb, cfg, connection, cur)

    progress('SLURB driver done')
    ncfile_slurb.close()



progress('PALM GEM finished OK')
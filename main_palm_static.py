#!/usr/bin/python3B308
# -*- coding: utf-8 -*-
import psycopg2
from utils.palm_static_pg_lib import *
from config.config import load_config, cfg
from argparse import ArgumentParser
import getpass
from config.logger import *

########################################
#
# Authors: Martin Bures + Jaroslav Resler
#
# Institute of Computer Science, Prague
#
########################################

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

progress('Calculate terrain height in grid')
calculate_terrain_height(cfg, connection, cur)

progress('Calculate origin_z and ori_min')
calculate_origin_z_oro_min(cfg, connection, cur)

progress('Connect landcover to grid')
connect_landcover_grid(cfg, connection, cur)

progress('Filling missing building holes')
fill_missing_holes_in_grid(cfg, connection, cur)

progress('Processing building elevation model BEM')
connect_buildings_height(cfg, connection, cur)

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

progress('Writing pavement type')
write_pavements(ncfile, cfg, connection, cur)

progress('Writing water type')
write_water(ncfile, cfg, connection, cur)

progress('Writing vegetation type')
write_vegetation(ncfile, cfg, connection, cur)

progress('Writing soil type')
write_soil(ncfile, cfg, connection, cur)

progress('Writing terrain type')
write_buildings(ncfile, cfg, connection, cur)

progress('Check consistency of static driver file')
check_consistency(ncfile, cfg)

ncfile.close()
debug('File {} was closed', cfg.domain.static_driver_file)
progress('Generating of static driver was successfully completed')
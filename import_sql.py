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
from pandas.plotting import table

from config.config import load_config, cfg
from argparse import ArgumentParser
import getpass
import os
import geopandas as gpd

from config.logger import *
from utils.import_utils import *

geom_name = 'geometry'

progress('Reading configuration')
argp = ArgumentParser(description=__doc__)
argp.add_argument('-c', '--config', help='configuration file')

# load config file
argv = argp.parse_args()

# TODO: solve default config
load_config(argv, v=2, cfg_default_path='config/default_import.yaml')
logging_level(cfg)

progress('Importing file to postgreSQL database, VERSION: python script')

debug('Config file: {}', argv.config)

# create connection to the postgresql server
if cfg.pg_password is None or cfg.pg_password == '':
    pg_password = getpass.getpass()
connection = psycopg2.connect(database=cfg.database, host=cfg.pg_host,
                              port=cfg.pg_port, user=cfg.pg_user, password=cfg.pg_password)
connection.set_client_encoding('UTF8')
cur = connection.cursor()
sql_debug(connection)

# TODO: move this to .yaml config
vector_files = {
    'imported_landcover': {
        'table_name': 'imported_landcover',
        'file': cfg.upload_files.urban_atlas,
        'columns2copy': ['code_2018']
    },
    'imported_streetmap': {
        'table_name': 'imported_streetmaps_py',
        'file': cfg.upload_files.street_maps,
        'columns2copy': ['osm_id']
    },
    'landcover': {
        'table_name': 'landcover',
        'file': cfg.upload_files.landcover,
        'columns2copy': ['lid', 'type']
    },
    'roofs': {
        'table_name': 'roofs',
        'file': cfg.upload_files.roofs,
        'columns2copy': ['']
    },
    'walls': {
        'table_name': 'walls',
        'file': cfg.upload_files.walls,
        'columns2copy': ['']
    },
    'trees': {
        'table_name': 'trees',
        'file': cfg.upload_files.trees,
        'columns2copy': ['']
    },
}
raster_files = {
    'imported_dem': {
        'table_name': 'imported_dem',
        'file': cfg.upload_files.imported_dem,
    },
    'imported_buildings': {
        'table_name': 'imported_buildings',
        'file': cfg.upload_files.imported_buildings,
    },
    'dem': {
        'table_name': 'dem',
        'file': cfg.upload_files.dem,
    },
    'buildings': {
        'table_name': 'buildings',
        'file': cfg.upload_files.buildings,
    },
    'extras': {
        'table_name': 'extras',
        'file': cfg.upload_files.extras,
    },
}

# CREATE SCHEMA
create_schema(cfg, connection, cur)

# PROCESS SHAPEFILE
for vector_file in vector_files.keys():
    progress('Processing {}', vector_file)

    shp_file = r"{}".format(vector_files[vector_file]['file'])
    if not os.path.isfile(shp_file):
        debug('\tFile does not exists, skip. \n{}', shp_file)
        continue
    # read shape file
    gdf = gpd.read_file(shp_file)

    table_name = vector_files[vector_file]['table_name']

    columns2copy = check_column(vector_files[vector_file]['columns2copy'], gdf)
    columns2copy += [geom_name]

    debug('shp_file: {} \n table_name: {} \n colums2copy: {}', shp_file, table_name, columns2copy)

    # create empty table based on column types
    create_empty_table(gdf, cfg, connection, cur, table_name, columns2copy, geom_name)

    # upload data from shapefile into postgreSQL
    upload_data(gdf, cfg, connection, cur, table_name, columns2copy, geom_name)

# PROCESS RASTERFILE
for raster_file in raster_files.keys():
    progress('Processing {}', raster_file)
    tiff_file = r"{}".format(raster_files[raster_file]['file'])
    if not os.path.isfile(tiff_file):
        debug('\tFile does not exists, skip. \n{}', tiff_file)
        continue
    table_name = raster_files[raster_file]['table_name']
    upload_raster(cfg, connection, cur, tiff_file, table_name)
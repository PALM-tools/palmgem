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

from config.logger import *
from shapely import wkb
import platform
import os

def create_schema(cfg, connection, cur):
    """
    CREATE SCHEMA
    """
    progress('Create a schema: {} in database', cfg.input_schema)
    if cfg.scratch_import:
        debug('Drop previous schema')
        sqltext = 'DROP SCHEMA IF EXISTS "{}" CASCADE'.format(cfg.input_schema)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

    debug('Creating new schema')
    sqltext = 'CREATE SCHEMA IF NOT EXISTS "{}"'.format(cfg.input_schema)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    debug('Changing owner to {}', cfg.pg_owner)
    sqltext = 'ALTER SCHEMA "{}" OWNER TO {}'.format(cfg.input_schema, cfg.pg_owner)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

def check_column(column2copy, gdf):
    """ A simple function to check if column exists in import file """
    gdf_column = gdf.columns
    for column in column2copy:
        if not column in gdf_column:
            debug('Column: {} does not exists in this imported file', column)
            column2copy.remove(column)

    return column2copy


def create_empty_table(gdf, cfg, connection, cur, table_name, columns2copy, geom_name):
    """
        Create an empty table with defined columns
    """
    debug('Preparing table: {}', table_name)
    sqltext = 'DROP TABLE IF EXISTS "{0}"."{1}" CASCADE'.format(cfg.input_schema, table_name)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    srid = gdf.crs.to_epsg()

    sqltext = 'CREATE TABLE "{0}"."{1}" (gid SERIAL PRIMARY KEY '

    for column in columns2copy:
        column_type = str(gdf.dtypes[column])
        if column == geom_name:
            sqltext += f', geom GEOMETRY(GEOMETRY, {srid})'
        elif column_type.find('int') != -1:
            sqltext += ', ' + column + ' INTEGER'
        elif column_type.find('float') != -1:
            sqltext += ', ' + column + ' NUMERIC'
        else:
            sqltext += ', ' + column + ' TEXT'

    sqltext += ')'

    cur.execute(sqltext.format(cfg.input_schema, table_name))
    sql_debug(connection)
    connection.commit()

def upload_data(gdf, cfg, connection, cur, table_name, columns2copy, geom_name):
    """
    Function to upload geometry data into postgreSQL database
    """
    debug('Uploading values into table {}', table_name)
    srid = gdf.crs.to_epsg()

    # convert polygon string from shapefile to WKB - POSTGIS standard
    # geom = to_wkb(
    #     set_srid(gdf[geom_name].values, srid=srid), hex=True, include_srid=True
    # )
    geom = gdf[geom_name].to_wkb(hex=True, include_srid=True)
    geom = wkb.dumps(gdf[geom_name], hex=True, srid=srid)
    gdf = gdf.drop(geom_name, axis=1)
    gdf[geom_name] = geom  # .values

    tuples = [tuple(x) for x in gdf.to_numpy()]
    columns2copy = [col.replace('geometry', 'geom') for col in columns2copy]
    cols = ','.join(columns2copy)
    sqltext = ('INSERT INTO "{0}"."{1}" ({2}) VALUES ({3})'
               .format(cfg.input_schema, table_name, cols, ','.join(['%s'] * (len(columns2copy))), ))
    cur.executemany(sqltext, tuples)
    sql_debug(connection)
    connection.commit()

def upload_raster(cfg, connection, cur, raster_file, table_name):
    """
    A function to upload raster tiff file to postgresql database
    """
    # check the platform
    if platform.system() == "Windows":
        exe_str = '.exe'
    elif platform.system() == "Linux":
        exe_str = ''
    else:
        error('Unknown system')
        exit(1)

    sqltext = 'DROP TABLE IF EXISTS "{0}"."{1}" CASCADE'.format(cfg.input_schema, table_name)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    pw = cfg.pg_password
    file_path = raster_file
    schema = cfg.input_schema
    file_name = table_name
    user = cfg.pg_user
    dbname = cfg.database
    host = cfg.pg_host
    port = cfg.pg_port
    bin_path = r"{}".format(cfg.bin_path)
    cmd = f"""
    set PGPASSWORD={pw}&&SET PATH=%PATH%;{bin_path}&& "raster2pgsql{exe_str}" -I -C -M -t auto "{file_path}" -q "{schema}"."{file_name}" | psql.exe -U {user} -d {dbname} -h {host} -p {port}
            """
    debug(cmd)
    os.system(cmd)
    debug('\tUpload DONE')

    sqltext = 'ALTER TABLE "{}"."{}" OWNER TO {}'.format(cfg.input_schema, table_name, cfg.pg_owner)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()
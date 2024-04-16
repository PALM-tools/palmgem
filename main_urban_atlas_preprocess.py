#!/usr/bin/python3
# -*- coding: utf-8 -*-
import psycopg2
from argparse import ArgumentParser
import getpass
import os
import sys
from argparse import ArgumentParser
import getopt
from config.logger import *
from config.config import load_config, cfg
from utils.palm_static_pg_lib import *

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
load_config(argv, cfg_default_path='config/default_config_preproc.yaml')
logging_level(cfg)

# print configuration into log file
progress('Configuration: {}', argv.config)
progress('pg User: {}', cfg.pg_user )
progress('pg Password: ', 'Weak Password! Thanks!' )
progress('pg Host: {}', cfg.pg_host)
progress('pg port: {}', cfg.pg_port )
progress('database: {}', cfg.database)
progress('pg owner: {}', cfg.pg_owner)
progress('case schema: {}', cfg.domain.case_schema)
progress('srid: {}', cfg.srid)
progress('srid dem: {}', cfg.dem_srid)
progress('clean up: {}', cfg.clean_up)
progress('log path: {}', cfg.logs.path)
progress('level of logging: {}', cfg.logs.level)

# create connection to the postgresql server
if cfg.pg_password is None or cfg.pg_password == '':
    pg_password = getpass.getpass()
connection = psycopg2.connect(database=cfg.database, host=cfg.pg_host,
                              port=cfg.pg_port, user=cfg.pg_user, password=cfg.pg_password)
connection.set_client_encoding('UTF8')
cur = connection.cursor()
sql_debug(connection)

if cfg.domain.crop_domain:
    progress('Creating user defined envelope')
    sqltext = 'select ST_MakeEnvelope(%s, %s, %s, %s, %s)'
    cur.execute(sqltext, (cfg.domain.xl, cfg.domain.yl, cfg.domain.xh, cfg.domain.yh, cfg.srid,))
    envelope = cur.fetchone()[0]
    connection.commit()
else:
    progress('Creating envelope of imported landcover')
    sqltext = 'SELECT ST_Union(geom) AS geom FROM "{0}"."{1}"' \
        .format(cfg.domain.case_schema, cfg.tables.im_landcover_or)
    cur.execute(sqltext)
    union = cur.fetchone()[0]
    sql_debug(connection)
    connection.commit()

    verbose('Creating rectangular envelop of Union')
    sqltext = 'SELECT ST_Transform(ST_Envelope(%s::geometry), %s)'
    cur.execute(sqltext, (union, cfg.srid))
    envelope = cur.fetchone()[0]
    sql_debug(connection)
    connection.commit()

    verbose('Obtaining envelope corners')
    sqltext = 'SELECT ST_XMin(%s::geometry), ST_YMin(%s::geometry), ' \
              '       ST_XMax(%s::geometry), ST_YMax(%s::geometry)'
    cur.execute(sqltext, (envelope, envelope, envelope, envelope,))
    envelope_bounds = cur.fetchall()
    sql_debug(connection)
    connection.commit()
    xl, yl, xh, yh = envelope_bounds[0]
    cfg.domain._settings['xl'] = xl
    cfg.domain._settings['yl'] = yl
    cfg.domain._settings['xh'] = xh
    cfg.domain._settings['yh'] = yh

""" DEM """
progress('Processing DEM')
debug('Updating DEM raster with SRID')
sqltext = 'SELECT UpdateRasterSRID(%s, %s, %s, %s)'
cur.execute(sqltext, (cfg.domain.case_schema, cfg.tables.dem_or, 'rast', cfg.dem_srid, ))
sql_debug(connection)
connection.commit()

debug('Delete all entries from DEM imported table that are outside envelope')
sqltext = 'DELETE FROM "{0}"."{1}" ' \
          'WHERE NOT ST_Intersects(ST_Transform(%s::geometry, %s), rast)'.format(cfg.domain.case_schema, cfg.tables.dem_or)
cur.execute(sqltext, (envelope, cfg.dem_srid,))
sql_debug(connection)
connection.commit()

debug('Create new table with transformed coordinates')
sqltext = 'DROP TABLE IF EXISTS "{0}"."{1}"; ' \
          'CREATE TABLE "{0}"."{1}" AS ' \
          'SELECT ROW_NUMBER() OVER(ORDER BY t.rast::geometry) AS rid, ' \
          '       ST_Union(ST_Clip( ST_Transform( r.rast, t.rast), t.rast::geometry ), {3}) AS rast ' \
          'FROM (SELECT ST_Transform(ST_SetSRID(ST_Extent(rast::geometry), %s), %s) AS geom ' \
          '      FROM "{0}"."{2}") AS g, ' \
          'ST_MakeEmptyCoverage(tilewidth => (SELECT (ST_MetaData(rast)).width FROM "{0}"."{2}" LIMIT 1), ' \
          '                     tileheight => (SELECT (ST_MetaData(rast)).height FROM "{0}"."{2}" LIMIT 1), ' \
          '                     width => (ST_XMax(g.geom) - ST_XMin(g.geom))::integer,' \
          '                     height => (ST_YMax(g.geom) - ST_YMin(g.geom))::integer,' \
          '                     upperleftx => ST_XMin(g.geom), ' \
          '                     upperlefty => ST_YMax(g.geom), ' \
          '                     scalex =>  (SELECT (ST_MetaData(rast)).scalex FROM "{0}"."{2}" LIMIT 1),' \
          '                     scaley => (SELECT (ST_MetaData(rast)).scaley FROM "{0}"."{2}" LIMIT 1),' \
          '                     skewx => (SELECT (ST_MetaData(rast)).skewx FROM "{0}"."{2}" LIMIT 1), ' \
          '                     skewy => (SELECT (ST_MetaData(rast)).skewy FROM "{0}"."{2}" LIMIT 1),' \
          '                     srid => %s) AS t(rast) ' \
          'INNER JOIN "{0}"."{2}" AS r ON ST_Transform(t.rast::geometry, %s) && r.rast ' \
          'GROUP BY t.rast;'.format(cfg.domain.case_schema, cfg.tables.dem, cfg.tables.dem_or, "'MAX'")
cur.execute(sqltext, (cfg.dem_srid, cfg.srid, cfg.srid, cfg.dem_srid))
sql_debug(connection)
connection.commit()

verbose('Changing owner of DEM table')
sqltext = 'ALTER TABLE "{}"."{}" OWNER TO {}'.format(cfg.domain.case_schema, cfg.tables.dem, cfg.pg_owner)
cur.execute(sqltext)
sql_debug(connection)
connection.commit()

if cfg.clean_up:
    debug('Deleting original imported DEM table')
    sqltext = 'DROP TABLE "{0}"."{1}" CASCADE;'.format(cfg.domain.case_schema, cfg.tables.dem_or)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

""" DEM BUILDINGS """
progress('Processing buildings DEM')
debug('Create new table with transformed coordinates')
sqltext = 'DROP TABLE IF EXISTS "{0}"."{1}" CASCADE; ' \
          'CREATE TABLE "{0}"."{1}" AS ' \
          'SELECT ROW_NUMBER() OVER(ORDER BY t.rast::geometry) AS rid, ' \
          '       ST_Union(ST_Clip( ST_Transform( r.rast, t.rast), t.rast::geometry ), {3}) AS rast ' \
          'FROM (SELECT ST_Transform(ST_SetSRID(ST_Extent(rast::geometry), %s), %s) AS geom ' \
          '      FROM "{0}"."{2}") AS g, ' \
          'ST_MakeEmptyCoverage(tilewidth => (SELECT (ST_MetaData(rast)).width FROM "{0}"."{2}" LIMIT 1), ' \
          '                     tileheight => (SELECT (ST_MetaData(rast)).height FROM "{0}"."{2}" LIMIT 1), ' \
          '                     width => (ST_XMax(g.geom) - ST_XMin(g.geom))::integer,' \
          '                     height => (ST_YMax(g.geom) - ST_YMin(g.geom))::integer,' \
          '                     upperleftx => ST_XMin(g.geom), ' \
          '                     upperlefty => ST_YMax(g.geom), ' \
          '                     scalex =>  (SELECT (ST_MetaData(rast)).scalex FROM "{0}"."{2}" LIMIT 1),' \
          '                     scaley => (SELECT (ST_MetaData(rast)).scaley FROM "{0}"."{2}" LIMIT 1),' \
          '                     skewx => (SELECT (ST_MetaData(rast)).skewx FROM "{0}"."{2}" LIMIT 1), ' \
          '                     skewy => (SELECT (ST_MetaData(rast)).skewy FROM "{0}"."{2}" LIMIT 1),' \
          '                     srid => %s) AS t(rast) ' \
          'INNER JOIN "{0}"."{2}" AS r ON ST_Transform(t.rast::geometry, %s) && r.rast ' \
          'GROUP BY t.rast;'.format(cfg.domain.case_schema, cfg.tables.buildings, cfg.tables.buildings_or, "'MAX'")
cur.execute(sqltext, (cfg.dem_srid, cfg.srid, cfg.srid, cfg.dem_srid))
sql_debug(connection)
connection.commit()

verbose('Changing owner of the table')
sqltext = 'ALTER TABLE "{}"."{}" OWNER TO {}'.format(cfg.domain.case_schema, cfg.tables.buildings, cfg.pg_owner)
cur.execute(sqltext)
sql_debug(connection)

debug('Adding serial rid index')
sqltext = 'ALTER TABLE "{0}"."{1}" DROP COLUMN IF EXISTS rid; ' \
          'ALTER TABLE "{0}"."{1}" ADD COLUMN IF NOT EXISTS rid SERIAL' \
    .format(cfg.domain.case_schema, cfg.tables.buildings)
cur.execute(sqltext)
sql_debug(connection)
connection.commit()

if cfg.clean_up:
    debug('Cleaning up imported buildings DEM')
    sqltext = 'DROP TABLE "{0}"."{1}" CASCADE;'.format(cfg.domain.case_schema, cfg.tables.buildings_or)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

""" ORIGINAL LANDCOVER """
progress('Processing Landcover')
# verbose('Adding geometry index to landcover table')
# sqltext = 'create index if not exists {1}_geom_idx on "{0}"."{1}" using gist(geom)' \
#     .format(cfg.domain.case_schema, cfg.tables.im_landcover_or)
# cur.execute(sqltext)
# sql_debug(connection)
# connection.commit()

debug('Creating, transforming and clipping landcover from originally imported one')
tranform_original_landcover(cfg, connection, cur, envelope)

# create background, mainly for cases with sea
# obtain union
debug('Processing background to landcover')
if cfg.domain.fill_boundary:
    verbose('Obtaining Union of landcover')
    sqltext = 'SELECT ST_Union(geom) AS geom FROM "{0}"."{1}"'\
              .format(cfg.domain.case_schema, cfg.tables.im_landcover)
    cur.execute(sqltext)
    union = cur.fetchone()[0]
    sql_debug(connection)
    connection.commit()

    verbose('Creating rectangular envelop of Union')
    sqltext = 'SELECT ST_Envelope(%s::geometry)'
    cur.execute(sqltext, (union,))
    envelope_background = cur.fetchone()[0]
    sql_debug(connection)
    connection.commit()

    verbose('Creating background as difference of envelope and landcover')
    sqltext = 'SELECT ST_Difference(%s::geometry, %s::geometry)'
    cur.execute(sqltext, (envelope_background, union))
    background = cur.fetchone()[0]
    sql_debug(connection)
    connection.commit()

    verbose('Inserting background into landcover')
    sqltext = 'INSERT INTO "{0}"."{1}" (geom) SELECT %s::geometry'.format(cfg.domain.case_schema, cfg.tables.im_landcover)
    cur.execute(sqltext, (background,))
    sql_debug(connection)
    connection.commit()

""" PROCESS STREETMAPS"""
progress('Processing Streetmaps layer')
# create from import
debug('Creating, clipping and transforming originally streetmap imported table')
sqltext = 'DROP TABLE IF EXISTS "{0}"."{1}";' \
          'CREATE TABLE "{0}"."{1}" AS ' \
          'SELECT osm_id, ST_Transform((ST_Dump(geom)).geom, %s) AS geom ' \
          'FROM "{0}"."{2}" ' \
          'WHERE (ST_Area(ST_Transform(geom, %s)) > %s) AND' \
          '       ST_Intersects(ST_Transform(geom, %s), %s::geometry)'.format(cfg.domain.case_schema, cfg.tables.streetmaps, cfg.tables.streetmaps_or)
cur.execute(sqltext, (cfg.srid, cfg.srid, cfg.max_stl_area, cfg.srid, envelope, ))
sql_debug(connection)
connection.commit()

verbose('Changing owner of the new streetmap table')
sqltext = 'ALTER TABLE "{}"."{}" OWNER TO {}'.format(cfg.domain.case_schema, cfg.tables.streetmaps, cfg.pg_owner)
cur.execute(sqltext)
sql_debug(connection)

# create gist spatial index
verbose('Adding geometry index on streetmap table')
sqltext = 'create index if not exists {1}_geom_idx on "{0}"."{1}" using gist(geom)' \
    .format(cfg.domain.case_schema, cfg.tables.streetmaps)
cur.execute(sqltext)
sql_debug(connection)
connection.commit()

verbose('Add code_2018 attribute into streepmaps table, for further use')
sqltext = 'ALTER TABLE "{0}"."{1}" ' \
          'ADD COLUMN IF NOT EXISTS code_2018 integer'\
          .format(cfg.domain.case_schema, cfg.tables.streetmaps)
cur.execute(sqltext)
sql_debug(connection)
connection.commit()

verbose('Add center of streetmaps buildings to optimize joining')
sqltext = 'ALTER TABLE "{0}"."{1}" ' \
          'ADD COLUMN IF NOT EXISTS cent geometry("POINT", %s)'\
          .format(cfg.domain.case_schema, cfg.tables.streetmaps)
cur.execute(sqltext, (cfg.srid,))
sql_debug(connection)
connection.commit()

sqltext = 'UPDATE "{0}"."{1}" ' \
          'SET cent = ST_Centroid(geom)'\
          .format(cfg.domain.case_schema, cfg.tables.streetmaps)
cur.execute(sqltext)
sql_debug(connection)
connection.commit()

verbose('Join code_2018 in streetmaps with landcover')
sqltext = 'UPDATE "{0}"."{1}" AS s SET ' \
          'code_2018 = (SELECT CAST(code_2018 AS integer) ' \
          '             FROM  "{0}"."{2}" AS l ' \
          '             WHERE ST_Intersects(l.geom, s.cent) LIMIT 1)'\
          .format(cfg.domain.case_schema, cfg.tables.streetmaps, cfg.tables.im_landcover)
cur.execute(sqltext)
sql_debug(connection)
connection.commit()

if cfg.clean_up:
    debug('Deleting original import streetmap table')
    sqltext = 'DROP TABLE "{0}"."{1}" CASCADE;'.format(cfg.domain.case_schema, cfg.tables.streetmaps_or)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

""" FINAL LANDCOVER """
progress('Finalizing Landcover table')
joined_streetmaps = 'joined_streetmaps'
debug('Preparing streetmap union')
sqltext = 'DROP TABLE IF EXISTS "{0}"."{1}";' \
          'CREATE TABLE "{0}"."{1}" AS ' \
          'SELECT ST_Union(geom) AS geom FROM "{0}"."{2}"'.format(cfg.domain.case_schema, joined_streetmaps, cfg.tables.streetmaps)
cur.execute(sqltext,)
sql_debug(connection)
connection.commit()

verbose('Changing table owner')
sqltext = 'ALTER TABLE "{}"."{}" OWNER TO {}'.format(cfg.domain.case_schema, joined_streetmaps, cfg.pg_owner)
cur.execute(sqltext)
sql_debug(connection)
connection.commit()

debug('Adding geom index')
# create gist spatial index
sqltext = 'create index if not exists {1}_geom_idx on "{0}"."{1}" using gist(geom)' \
          .format(cfg.domain.case_schema, joined_streetmaps)
cur.execute(sqltext)
sql_debug(connection)
connection.commit()

progress('Creating difference between landcover and unioned streetmaps, NOTE: it might take few hours')
sqltext = 'WITH s AS (SELECT geom FROM "{0}"."{1}") ' \
          'UPDATE "{0}"."{2}" AS l SET ' \
          'geom = ST_Difference(l.geom, s.geom) ' \
          'FROM s ' \
          'WHERE ST_Intersects(l.geom, s.geom) '\
          .format(cfg.domain.case_schema, joined_streetmaps, cfg.tables.im_landcover)
cur.execute(sqltext)
sql_debug(connection)
connection.commit()

verbose('Adding additional information into landcover table')
sqltext = 'ALTER TABLE "{0}"."{1}" ' \
          'RENAME COLUMN code_2018 TO code_2018_char'\
          .format(cfg.domain.case_schema, cfg.tables.im_landcover)
cur.execute(sqltext)
sql_debug(connection)
connection.commit()

verbose('Add extra column into landcover table')
sqltext = 'ALTER TABLE "{0}"."{1}" ' \
          'ADD COLUMN IF NOT EXISTS osm_id integer, ' \
          'ADD COLUMN IF NOT EXISTS type integer, ' \
          'ADD COLUMN IF NOT EXISTS code_2018 integer '.format(cfg.domain.case_schema, cfg.tables.im_landcover)
cur.execute(sqltext)
sql_debug(connection)
connection.commit()

verbose('Update code_2018 from char into integer')
sqltext = 'UPDATE "{0}"."{1}" SET ' \
          'code_2018 = CAST(code_2018_char AS integer) '\
          .format(cfg.domain.case_schema, cfg.tables.im_landcover)
cur.execute(sqltext)
sql_debug(connection)
connection.commit()

debug('Inserting streetmaps into landcover')
sqltext = 'INSERT INTO "{0}"."{1}" ' \
          '(type, osm_id, code_2018, geom) ' \
          'SELECT NULL, CAST(osm_id AS integer), code_2018, geom ' \
          'FROM "{0}"."{2}"'.format(cfg.domain.case_schema, cfg.tables.im_landcover, cfg.tables.streetmaps)
cur.execute(sqltext, )
sql_debug(connection)
connection.commit()

debug('Adding serial lid index to landcover table')
sqltext = 'ALTER TABLE "{0}"."{1}" ADD COLUMN IF NOT EXISTS lid SERIAL' \
    .format(cfg.domain.case_schema, cfg.tables.im_landcover)
cur.execute(sqltext)
sql_debug(connection)
connection.commit()

verbose('Make lid unique key')
sqltext = 'ALTER TABLE "{0}"."{1}" ADD PRIMARY KEY (lid)'.format(cfg.domain.case_schema, cfg.tables.im_landcover)
cur.execute(sqltext)
sql_debug(connection)
connection.commit()

if cfg.clean_up:
    sqltext = 'DROP TABLE "{0}"."{1}" CASCADE;'.format(cfg.domain.case_schema, joined_streetmaps)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    debug('Deleting original imported landcover table')
    sqltext = 'DROP TABLE "{0}"."{1}" CASCADE;'.format(cfg.domain.case_schema, cfg.tables.im_landcover_or)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

progress('Joining user defined class and ladcover classes into PALM types')
debug('Updating PALM type')
sqltext = 'UPDATE "{0}"."{1}" SET type = CASE '\
          .format(cfg.domain.case_schema, cfg.tables.im_landcover)
for m in cfg.mt:
    sqltext += 'WHEN code_2018 = {0} THEN CASE WHEN osm_id IS NULL THEN {1} ELSE {2} END '.format(m[0], m[1], m[2])
sqltext += 'ELSE {0} END '.format(cfg.mt_default)
cur.execute(sqltext)
sql_debug(connection)
connection.commit()
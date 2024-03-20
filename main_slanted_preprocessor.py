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

exit(1)

if cfg.domain.crop_domain:
    progress('Creating user defined envelope')
    sqltext = 'select ST_MakeEnvelope(%s, %s, %s, %s, %s)'
    cur.execute(sqltext, (cfg.domain.xl, cfg.domain.yl, cfg.domain.xh, cfg.domain.yh, cfg.srid,))
    envelope = cur.fetchone()[0]
    connection.commit()
else:
    progress('Creating envelope of imported landcover')
    sqltext = 'SELECT ST_Union(geom) AS geom FROM "{0}"."{1}"' \
        .format(cfg.domain.case_schema, cfg.tables.landcover)
    cur.execute(sqltext)
    envelope = cur.fetchone()[0]
    sql_debug(connection)
    connection.commit()

# create slanted mask, without buildings, if necessary, include option in configuration
progress('Created slanted mask')
sqltext = 'DROP TABLE IF EXISTS "{0}"."{2}"; ' \
          'CREATE TABLE "{0}"."{2}" AS ' \
          'WITH lb AS (SELECT ST_Buffer(ST_Union(ST_MakeValid(geom)), 10, {3}) AS geom FROM "{0}"."{1}" WHERE type >= 900)  ' \
          'SELECT ST_Difference(%s::geometry, lb.geom) AS geom ' \
          'FROM lb, ll'.format(cfg.domain.case_schema, cfg.tables.landcover, cfg.tables.slanted_mask, 'quad_segs=2, endcap=square')
cur.execute(sqltext, (envelope, ))
sql_debug(connection)
connection.commit()

progress('Add unique index to slanted mask table')
sqltext = 'UPDATE TABLE "{0}"."{1}" ' \
          'ADD COLUMN IF NOT EXISTS gid SERIAL'.format(cfg.domain.case_schema, cfg.tables.slanted_mask)
cur.execute(sqltext)
sql_debug(connection)
connection.commit()

# create roof table for slanted faces
progress('Creating roof table for slanted faces usage')
sqltext = 'DROP TABLE IF EXISTS "{0}"."{2}"; ' \
          'CREATE TABLE "{0}"."{2}" AS ' \
          'SELECT geom  ' \
          'FROM "{0}"."{1}"  ' \
          'WHERE type >= 900'.format(cfg.domain.case_schema, cfg.tables.landcover, cfg.tables.roofs)
cur.execute(sqltext)
sql_debug(connection)
connection.commit()

progress('Add unique index to roofs table')
sqltext = 'UPDATE TABLE "{0}"."{1}" ' \
          'ADD COLUMN IF NOT EXISTS rid SERIAL'.format(cfg.domain.case_schema, cfg.tables.roofs)
cur.execute(sqltext)
sql_debug(connection)
connection.commit()

# create walls table for slanted faces
progress('Creating walls table for slanted faces usage')
sqltext = 'DROP TABLE IF EXISTS "{0}"."{2}"; ' \
          'CREATE TABLE "{0}"."{2}" AS ' \
          'WITH lb AS (SELECT geom FROM "{0}"."{1}" WHERE type >= 900) ' \
          'SELECT ST_Boundary(lb.geom) AS geom ' \
          'FROM lb'.format(cfg.domain.case_schema, cfg.tables.landcover, cfg.tables.walls)
cur.execute(sqltext)
sql_debug(connection)
connection.commit()

progress('Add unique index to walls table')
sqltext = 'UPDATE TABLE "{0}"."{1}" ' \
          'ADD COLUMN IF NOT EXISTS wid SERIAL'.format(cfg.domain.case_schema, cfg.tables.walls)
cur.execute(sqltext)
sql_debug(connection)
connection.commit()
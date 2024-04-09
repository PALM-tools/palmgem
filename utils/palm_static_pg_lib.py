import sys
from datetime import datetime
from pathlib import Path
from math import ceil
import numpy as np
from netCDF4 import Dataset
import pandas as pd
import matplotlib.pyplot  as plt
import os
from config.logger import *

def create_fishnet(cfg, connection, cur):
    """ CREATE fishnet
        Regular gridnet in defined extend with defined number of grid fishnet.nx fishnet.ny
    """
    progress('Creating fishnet')
    # TODO: do case when nx is not set, do some optimum dx, dy
    dx = cfg.fishnet.dx
    dy = cfg.fishnet.dy
    xl = cfg.domain.xl
    xh = cfg.domain.xh
    yl = cfg.domain.yl
    yh = cfg.domain.yh
    nx = int(np.ceil((xh - xl) / dx))
    ny = int(np.ceil((yh - yl) / dy))
    xh = xl + nx * dx
    yh = yl + ny * dy

    xs = (xl + xh) / 2.0
    ys = (yl + yh) / 2.0

    debug('Creating Fishnet using palm_create_grid psql function')
    cur.callproc('palm_create_grid',
                 [cfg.domain.case_schema, cfg.tables.fishnet,
                  nx, ny, dx, dy,
                  xs, ys,
                  cfg.srid, cfg.srid_wgs84, cfg.srid_utm, cfg.pg_owner,
                  cfg.logs.level])
    connection.commit()
    debug('Fishnet created')


def tranform_original_landcover(cfg, connection, cur, envelope):
    """
    Function that crop imported landcover by defined fishnet
    The reason is to split larger polygons into smaller and subsequently speed up spatial joins.
    """
    debug('Creating fishnet')
    create_fishnet(cfg, connection, cur)

    progress('Clipping landcover by fishnet')
    debug('Creating new cliped landcover table')
    sqltext = 'DROP TABLE IF EXISTS "{0}"."{3}" CASCADE;' \
              'CREATE TABLE "{0}"."{3}" AS  ' \
              'WITH f AS (SELECT geom AS geom FROM "{0}"."{1}"), ' \
              '     l AS (SELECT code_2018, ST_Transform((ST_Dump(geom)).geom, %s) AS geom ' \
              '            FROM "{0}"."{2}" ' \
              '            WHERE ST_Intersects(ST_Transform(geom, %s), %s::geometry)' \
              '            )' \
              'SELECT code_2018, (ST_Dump(ST_Intersection(l.geom, f.geom))).geom AS geom ' \
              'FROM l, f ' \
              'WHERE ST_Intersects(l.geom, f.geom) AND ST_Area(l.geom) > {4} ' \
              'UNION ALL ' \
              'SELECT code_2018, l.geom ' \
              'FROM l, f ' \
              'WHERE ST_Intersects(l.geom, f.geom) AND ST_Area(l.geom) <= {4}' \
              .format(cfg.domain.case_schema, cfg.tables.fishnet, cfg.tables.im_landcover_or,
                      cfg.tables.im_landcover, cfg.max_fishnet_split_area)
    cur.execute(sqltext, (cfg.srid, cfg.srid, envelope, ))
    sql_debug(connection)
    connection.commit()

    if cfg.clean_up:
        debug('Drop original landcover table')
        sqltext = 'DROP TABLE "{0}"."{1}" CASCADE;'.format(cfg.domain.case_schema, cfg.tables.im_landcover_or)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

    verbose('Changing owner of new landcover table into {}', cfg.pg_owner)
    sqltext = 'ALTER TABLE "{}"."{}" OWNER TO {}'.format(cfg.domain.case_schema, cfg.tables.im_landcover, cfg.pg_owner)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    verbose('Adding geometry index to landcover table')
    # create gist spatial index
    sqltext = 'create index if not exists {1}_geom_idx on "{0}"."{1}" using gist(geom)' \
        .format(cfg.domain.case_schema, cfg.tables.im_landcover)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

def create_grid(cfg, connection, cur):
    """ Create grid from configuration, using psql function utils/palm_create_grid.sql

    :param cfg: configuration object that hold all user specified and default configurations
    :param connection: psql object that connects to postgreSQL database
    :param cur: cursor to postgreSQL database
    """
    debug('Create regular case grid: {}', cfg.tables.grid)
    cur.callproc('palm_create_grid',
                 [cfg.domain.case_schema, cfg.tables.grid,
                  cfg.domain.nx, cfg.domain.ny, cfg.domain.dx, cfg.domain.dy,
                  cfg.domain.cent_x, cfg.domain.cent_y,
                  cfg.srid_palm, cfg.srid_wgs84, cfg.srid_utm, cfg.pg_owner,
                  cfg.logs.level])
    sql_debug(connection)
    connection.commit()
    debug('Grid created')

    # add nz and height column to grid
    debug('Adding height and nz field to grid')
    sqltext = 'alter table "{}"."{}" ' \
              'add if not exists height double precision, ' \
              'add if not exists nz integer, ' \
              'add if not exists lid integer' \
        .format(cfg.domain.case_schema, cfg.tables.grid)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

def calculate_grid_extend(cfg, connection, cur):
    """ Calculate grid envelope for further use in inputs clipping.
        Envelope is created as a rectangle around grid.
    """
    sqltext = 'select min(xmi), max(xma), min(ymi), max(yma) from "{}"."{}"' \
        .format(cfg.domain.case_schema, cfg.tables.grid)
    cur.execute(sqltext)
    gxmin, gxmax, gymin, gymax = cur.fetchone()[0:4]
    sqltext = 'select ST_MakeEnvelope(%s, %s, %s, %s, %s)'
    cur.execute(sqltext, (gxmin, gymin, gxmax, gymax, cfg.srid_palm,))
    connection.commit()
    grid_ext = cur.fetchone()[0]

    return grid_ext

def copy_vectors_from_input(grid_ext, cfg, connection, cur):
    """ Copy inputs vector layers from input schema.
        Transform and clip them to grid extend and grid coordinate system.

        params: grid_ext: rectangle polygon around grid, created in calculate_grid_extend function
    """
    vtables = [cfg.tables.landcover]
    vidx = [cfg.idx.landcover]
    vtabs = []
    for rel, idx in zip(vtables, vidx):
        # check if table exists in input source schema
        verbose('Transforming {} table', rel)
        sqltext = 'select exists(select * from information_schema.tables where table_schema=%s and table_name=%s)'
        cur.execute(sqltext, (cfg.input_schema, rel,))
        sql_debug(connection)
        has_rel = cur.fetchone()[0]
        if has_rel:
            vtabs.append(rel)
        else:
            error('Table {} does not exist in input schema', rel)
            exit(1)

        debug('Transform srid of input {} and limit it to grid', rel)
        # drop table if exists
        sqltext = 'drop table if exists "{}"."{}"'.format(cfg.domain.case_schema, rel)
        cur.execute(sqltext)
        sqltext = 'create table "{}"."{}" (like "{}"."{}" including all)' \
            .format(cfg.domain.case_schema, rel, cfg.input_schema, rel)
        cur.execute(sqltext)
        # check if table is empty or not
        cur.execute('SELECT COUNT(*) FROM "{}"."{}"'.format(cfg.input_schema, rel))
        count = cur.fetchone()[0]
        if count == 0:
            warning('There is 0 items in the input table: {}, pop, delete this table from case', rel)
            vtabs.remove(rel)
            continue
        # test srid of the input table
        sqltext = 'select ST_SRID(geom) from "{}"."{}" limit 1'.format(cfg.input_schema, rel)
        cur.execute(sqltext)
        srid_rel = cur.fetchone()[0]
        if srid_rel == 0:
            srid_rel = cfg.srid_input
            # update relation SRID
            sqltext = 'select SELECT UpdateGeometrySRID(%s, %s, %s)'
            cur.execute(sqltext, (cfg.input_schema, rel, 'geom', cfg.srid_input,))
        sqltext = 'insert into "{}"."{}" select * from "{}"."{}" where ST_Intersects(ST_Transform(geom, %s), %s)' \
            .format(cfg.domain.case_schema, rel, cfg.input_schema, rel)
        cur.execute(sqltext, (cfg.srid_palm, grid_ext,))
        sqltext = 'update "{}"."{}" set geom = ST_Transform(geom, %s)'.format(cfg.domain.case_schema, rel)
        cur.execute(sqltext, (cfg.srid_palm,))

        # check if table is empty or not
        cur.execute('SELECT COUNT(*) FROM "{}"."{}"'.format(cfg.domain.case_schema, rel))
        count = cur.fetchone()[0]
        if count == 0:
            warning('There is 0 items in the table: {}, pop, delete this table from case', rel)
            vtabs.remove(rel)
            cur.execute('DROP TABLE "{}"."{}"'.format(cfg.domain.case_schema, rel))
            sql_debug(connection)
            connection.commit()
            continue
        # create gist spatial index
        sqltext = 'create index if not exists {1}_geom_idx on "{0}"."{1}" using gist(geom)' \
            .format(cfg.domain.case_schema, rel)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()
        sqltext = 'ALTER TABLE "{}"."{}" OWNER TO {}'.format(cfg.domain.case_schema, rel, cfg.pg_owner)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()
        # check if unique index is unique at all
        sqltext = 'SELECT COUNT(*) FROM (SELECT COUNT(*) FROM "{0}"."{1}" ' \
                  ' GROUP BY {2} HAVING COUNT(*) > 1) AS a'.format(cfg.domain.case_schema, rel, idx)

        cur.execute(sqltext)
        sql_debug(connection)
        idx_uniques = cur.fetchone()[0]
        connection.commit()
        if idx_uniques > 0:
            warning('Index {} in table {} is not unique, moving original index to old column and creating new one',
                    idx, rel)
            sqltext = 'ALTER TABLE "{0}"."{1}" RENAME COLUMN {2} TO {3}'.format(cfg.domain.case_schema, rel, idx,
                                                                                idx + '_old')
            cur.execute(sqltext)
            sql_debug(connection)
            connection.commit()

            sqltext = 'ALTER TABLE "{0}"."{1}" ADD COLUMN {2} SERIAL' \
                      ''.format(cfg.domain.case_schema, rel, idx)
            cur.execute(sqltext)
            sql_debug(connection)
            connection.commit()

        # Assign unique index to idx
        # drop previous unique index
        sqltext = "SELECT indexname FROM pg_indexes" \
                  " WHERE schemaname = '{0}' and tablename = '{1}'  and indexname not like '%geom%' " \
            .format(cfg.domain.case_schema, rel)
        cur.execute(sqltext)
        sql_debug(connection)
        prev_ui = cur.fetchone()[0]
        connection.commit()
        sqltext = 'ALTER TABLE "{0}"."{1}" DROP CONSTRAINT {2}'.format(cfg.domain.case_schema, rel, prev_ui)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()
        sqltext = 'ALTER TABLE "{0}"."{1}" ADD PRIMARY KEY ({2})'.format(cfg.domain.case_schema, rel, idx)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        # check if table is empty or not
        cur.execute('SELECT COUNT(*) FROM "{}"."{}"'.format(cfg.domain.case_schema, rel))
        count = cur.fetchone()[0]
        if count == 0:
            error('There is 0 items in the table: {}, pop, delete this table from case', rel)
            exit(1)
    return vtabs

def copy_rasters_from_input(grid_ext, cfg, connection, cur):
    """ Copy raster tables from input schema.
        Transform and clip them to grid extend and grid coordinate system

        params: grid_ext: rectangle polygon around grid, created in calculate_grid_extend function
    """
    rtables = [cfg.tables.dem, cfg.tables.buildings_height, ]
    rtabs = []
    for rel in rtables:
        # check if table exists in input source schema
        verbose('Transforming {} table', rel)
        sqltext = 'select exists(select * from information_schema.tables where table_schema=%s and table_name=%s)'
        cur.execute(sqltext, (cfg.input_schema, rel))
        sql_debug(connection)
        has_rel = cur.fetchone()[0]
        if has_rel:
            rtabs.append(rel)
        else:
            error('Table {} does not exist in input schema', rel)
            exit(1)

        debug('Transform srid of input {} and limit it to grid', rel)
        # drop table if exists
        sqltext = 'drop table if exists "{}"."{}"'.format(cfg.domain.case_schema, rel)
        cur.execute(sqltext)
        sqltext = 'create table "{}"."{}" (like "{}"."{}" including all excluding constraints)' \
            .format(cfg.domain.case_schema, rel, cfg.input_schema, rel)
        cur.execute(sqltext)
        sql_debug(connection)
        sqltext = 'ALTER TABLE "{}"."{}" OWNER TO {}'.format(cfg.domain.case_schema, rel, cfg.pg_owner)
        cur.execute(sqltext)
        sql_debug(connection)
        # test srid of the input table
        sqltext = 'select ST_SRID(rast) from "{}"."{}" limit 1'.format(cfg.input_schema, rel)
        # alternative query:
        # select srid from "raster_columns" where r_table_catalog = 'palm_static' and r_table_schema = 'inputs_tunnel'
        # and r_table_name = 'buildings' and r_raster_column = 'rast'
        cur.execute(sqltext)
        srid_rel = cur.fetchone()[0]
        if srid_rel == 0:
            srid_rel = cfg.srid_input
            # update relation SRID
            sqltext = 'select UpdateRasterSRID(%s, %s, %s, %s)'
            cur.execute(sqltext, (cfg.input_schema, rel, 'rast', cfg.srid_input,))
        # transform and clip the raster layer
        if not srid_rel == cfg.srid_input:
            sqltext = 'insert into "{0}"."{2}" (rid, rast) ' \
                      'select rid, ST_Transform(rast, %s) ' \
                      ' from "{1}"."{2}" ' \
                      ' where ST_Intersects(ST_Transform(rast, %s), %s::geometry)' \
                .format(cfg.domain.case_schema, cfg.input_schema, rel)
            cur.execute(sqltext, (cfg.srid_palm, cfg.srid_palm, grid_ext,))
        else:
            sqltext = 'insert into "{0}"."{2}" (rid, rast) ' \
                      'select rid, rast ' \
                      ' from "{1}"."{2}" ' \
                      ' where ST_Intersects(rast, %s::geometry)' \
                .format(cfg.domain.case_schema, cfg.input_schema, rel)
            cur.execute(sqltext, (grid_ext,))

        sql_debug(connection)
        connection.commit()

        # check if table is empty or not
        cur.execute('SELECT COUNT(*) FROM "{}"."{}"'.format(cfg.domain.case_schema, rel))
        count = cur.fetchone()[0]
        if count == 0:
            error('There is 0 items in the table: {}, pop, delete this table from case', rel)
            exit(1)
    return rtabs

def calculate_terrain_height(cfg, connection, cur):
    """ Calculate terrain height for each grid cell in PALM domain.
        Height is calculated from raster table DEM in grid center.
    """
    sqltext = 'CREATE TEMP TABLE "temp" AS ' \
              'SELECT g.i, g.j, r.height, cast(0 AS INTEGER) AS nz, g.geom ' \
              'FROM (SELECT i, j, xcen, ycen, geom FROM "{0}"."{1}") AS g ' \
              'JOIN LATERAL ( ' \
              'SELECT ST_NearestValue(rast, ST_SetSRID(ST_Point(g.xcen,g.ycen), %s)) AS height ' \
              'FROM "{0}"."{2}" WHERE ST_Intersects(rast, ST_SetSRID(ST_Point(g.xcen,g.ycen), %s))' \
              'LIMIT 1) r on true ' \
              'WHERE r.height IS NOT NULL'
    sqltext = sqltext.format(cfg.domain.case_schema, cfg.tables.grid, cfg.tables.dem)
    cur.execute(sqltext, (cfg.srid_palm, cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()
    sqltext = 'ALTER TABLE "temp" ADD primary key (i,j)'
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()
    sqltext = 'UPDATE "{0}"."{1}" g ' \
              'SET height = (SELECT d.height FROM "temp" d WHERE d.i = g.i AND d.j = g.j)'
    sqltext = sqltext.format(cfg.domain.case_schema, cfg.tables.grid)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()
    # fill remaining
    sqltext = 'SELECT COUNT(*) FROM "{0}"."{1}"' \
              ' WHERE height IS NULL'.format(cfg.domain.case_schema, cfg.tables.grid)
    cur.execute(sqltext)
    missing_grid_heights = cur.fetchone()[0]
    sql_debug(connection)
    connection.commit()
    if missing_grid_heights > 0:
        warning('There is {} missing height in the grid, filling from the nearest raster', missing_grid_heights)
        # sqltext = 'UPDATE "{0}"."{1}" g ' \
        #           'SET height = (SELECT gg.height FROM "{0}"."{1}" AS gg' \
        #           '  WHERE gg.height IS NOT NULL ' \
        #           '  ORDER BY ST_Distance(ST_SetSrid(ST_Point(g.xcen,  g.ycen), %s), ' \
        #           '                       ST_SetSrid(ST_Point(gg.xcen,gg.ycen), %s))) ' \
        #           ' WHERE g.height IS NULL'.format(cfg.domain.case_schema, cfg.tables.grid)
        # cur.execute(sqltext, (cfg.srid_palm, cfg.srid_palm,))
        # sql_debug(connection)
        # connection.commit()

        sqltext = 'UPDATE "{0}"."{1}" g ' \
                  'SET height = (SELECT ST_Value(rast, ST_SetSRID(ST_Point(g.xcen,g.ycen), %s)) ' \
                  ' FROM "{0}"."{2}" LIMIT 1) ' \
                  'WHERE g.height IS NULL'.format(cfg.domain.case_schema, cfg.tables.grid, cfg.tables.dem)
        cur.execute(sqltext, (cfg.srid_palm,))
        sql_debug(connection)
        connection.commit()

        sqltext = 'select min(height) from "{}".{}'. \
            format(cfg.domain.case_schema, cfg.tables.grid)
        cur.execute(sqltext)
        min_fill_height = cur.fetchone()[0]

        sqltext = 'UPDATE "{0}"."{1}" g ' \
                  'SET height = {2} ' \
                  'WHERE height IS NULL'.format(cfg.domain.case_schema, cfg.tables.grid, min_fill_height)
        cur.execute(sqltext, (cfg.srid_palm,))
        sql_debug(connection)
        connection.commit()

def calculate_origin_z_oro_min(cfg, connection, cur):
    """ Calculate domain origin_z and oro_min.
        Origin_z is essential for domain location and in nested relation.
        Oro_min stores information about minimum terrain height in current domain.
        In case of single domain configuration origin_z == oro_min.
        In case of shared origin_z between nesting domains origin_z == oro_min.
        Calculation of gridded height in table grid.
    """
    progress('Calculation of origin_z and oro_min')
    sqltext = 'select min(height) from "{}".{}'. \
        format(cfg.domain.case_schema, cfg.tables.grid)
    cur.execute(sqltext)
    cfg.domain._settings['oro_min'] = cur.fetchone()[0]
    debug('Oro_min is in height level: {}', cfg.domain.oro_min)

    # calculate grid nz height
    sqltext = 'UPDATE "{0}"."{1}" SET nz = CAST((height - {2})/{3} AS INTEGER)'.format(
        cfg.domain.case_schema, cfg.tables.grid, cfg.domain.oro_min, cfg.domain.dz)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # calculate origin_z as bottom of the current or parent domain
    sqltext = 'select min(height) from "{}".{}'. \
        format(cfg.domain.case_schema if cfg.domain.parent_domain_schema == ''
               else cfg.domain.parent_domain_schema, cfg.tables.grid)
    cur.execute(sqltext)
    cfg.domain._settings['origin_z'] = cur.fetchone()[0]
    debug('Origin_z is in height level: {}', cfg.domain.origin_z)

    sqltext = 'select min(height) from "{}".{}'. \
        format(cfg.domain.case_schema, cfg.tables.grid)
    cur.execute(sqltext)
    cfg.domain._settings['oro_min'] = cur.fetchone()[0]
    debug('Oro_min is in height level: {}', cfg.domain.oro_min)

def connect_landcover_grid(cfg, connection, cur):
    """ Join landcover with grid """
    change_log_level(cfg.logs.level_landcover)
    sqltext = 'UPDATE "{0}"."{1}" ' \
              'SET lid = (SELECT l.lid FROM "{0}"."{2}" AS l ' \
              'WHERE ST_Within( ST_SetSRID(ST_Point(xcen,ycen),%s), l.geom) LIMIT 1)'
    sqltext = sqltext.format(cfg.domain.case_schema, cfg.tables.grid, cfg.tables.landcover)
    cur.execute(sqltext, (cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()
    # Remark: lid can be None in case that parts of the domain are not covered by landcover shapes!
    sqltext = 'update "{0}"."{1}" g ' \
              ' set lid = (select l.lid from "{0}"."{2}" l ' \
              ' order by ST_Distance(l.geom, ST_SetSRID(ST_Point(g.xcen,g.ycen),%s)) limit 1) ' \
              ' where g.lid is null'
    sqltext = sqltext.format(cfg.domain.case_schema, cfg.tables.grid, cfg.tables.landcover)
    cur.execute(sqltext, (cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()
    # test number of missing lid (should by none)
    sqltext = 'SELECT COUNT(*) FROM "{0}"."{1}" WHERE lid IS NULL'.format(cfg.domain.case_schema, cfg.tables.grid)
    cur.execute(sqltext)
    if cur.fetchone()[0] > 0:
        warning('There are some grid that dont have LID specified')
        sqltext = 'SELECT id, i, j, xcen, ycen FROM "{0}"."{1}" WHERE lid IS NULL' \
            .format(cfg.domain.case_schema, cfg.tables.grid)
        cur.execute(sqltext)
        missing_lids = cur.fetchall()
        for ms in missing_lids:
            warning('Missing grid [id,i,j,xcen,ycen] = [{},{},{},{},{}]', *ms)
    sql_debug(connection)
    connection.commit()

def fill_missing_holes_in_grid(cfg, connection, cur):
    """ Fill empty 3d grid cell that has more than 3 neighbors, replicating process in PALM program to omit creating grid cell with default configuration """
    cur.callproc('palm_fill_building_holes', [cfg.domain.case_schema, cfg.tables.grid, cfg.tables.landcover,
                                              cfg.type_range.building_min, cfg.type_range.building_max,
                                              cfg.domain.nx, cfg.domain.ny, cfg.logs.level])
    sql_debug(connection)
    connection.commit()
    restore_log_level(cfg)

def filling_grid(cfg, connection, cur):
    """ filling grid """

    # try adding height+max then round or add each rounds
    debug('Create nz_temp')
    sqltext = 'CREATE TEMP TABLE "nz_temp" AS ' \
	            'SELECT g.i AS i, g.j AS j,  ' \
	            'CASE WHEN b.height IS NOT NULL THEN b.nz+bo.max_int' \
		        '     ELSE g.nz END AS nz, ' \
                'CASE WHEN b.height IS NOT NULL THEN true ' \
                '     ELSE false END AS is_building, ' \
                'CASE WHEN b.height IS NOT NULL THEN b.height+bo.max_int + {6} ' \
                '     ELSE g.height END AS height ' \
                'FROM "{2}"."{3}" AS g ' \
	            'LEFT OUTER JOIN "{2}"."{4}" AS b ON g.id=b.id  ' \
	            'LEFT OUTER JOIN "{2}"."{5}" AS bo ON b.lid = bo.lid'.format(
                cfg.domain.oro_min, cfg.domain.dz, cfg.domain.case_schema, cfg.tables.grid,
                cfg.tables.buildings_grid, cfg.tables.buildings_offset, cfg.domain.oro_min)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    debug('Add primary keys to nz_temp')
    sqltext = 'ALTER TABLE "nz_temp" ADD primary key (i,j)'
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    sqltext = 'SELECT nt.i,nt.j FROM nz_temp AS nt ' \
                'WHERE (SELECT COUNT(*) FROM nz_temp AS ntt  ' \
	                'WHERE ((ntt.i = nt.i+1 AND ntt.j   = nt.j) OR ' \
                          ' (ntt.i = nt.i-1 AND ntt.j   = nt.j) OR ' \
                          ' (ntt.i = nt.i   AND ntt.j   = nt.j+1) OR ' \
                          ' (ntt.i = nt.i   AND ntt.j   = nt.j-1)) AND ' \
				          ' (ntt.nz > nt.nz)) > 2'
    debug('Fetching all places that need to be filled')
    sqlupdate = 'UPDATE "nz_temp" AS nt SET (nz, height) = ' \
                '   (SELECT ntt.nz, ntt.height FROM "nz_temp" AS ntt' \
                '   WHERE ((ntt.i = nt.i+1 AND ntt.j   = nt.j) OR ' \
                         ' (ntt.i = nt.i-1 AND ntt.j   = nt.j) OR ' \
                         ' (ntt.i = nt.i   AND ntt.j   = nt.j+1) OR ' \
                         ' (ntt.i = nt.i   AND ntt.j   = nt.j-1)) AND ' \
				         ' (ntt.nz > nt.nz) ' \
                          'ORDER BY ntt.nz LIMIT 1) ' \
                'WHERE nt.i = %s AND nt.j = %s'

    cur.execute(sqltext)
    missings = cur.fetchall()
    verbose('All missing place: {}', missings)
    empty_grids = len(missings)
    debug('Empty grids tobe filled: ', empty_grids)
    filled_grids = []
    if empty_grids > 0:
        empty_grids_old = 10*empty_grids
        debug('There is {} grid to fill at the beginning', empty_grids)
        iter = 0
        while empty_grids > 0:
            iter += 1
            empty_grids_old = empty_grids
            for i, j in missings:
                extra_verbose('In filling grid; filling [j, i] = [{},{}] ', j, i)
                cur.execute(sqlupdate, (i,j,))
                filled_grids.append([i,j])
            cur.execute(sqltext)
            missings = cur.fetchall()
            empty_grids = len(missings)
            # if empty_grids == 0:
            #     debug('All empty grid were filled')
            #     break
            debug('There is still {} grid to fill in iteration: {}', empty_grids, iter)
            if empty_grids == 0:
                debug('Nothing to fill, ending')
                break

    debug('Done with filling')
    sql_debug(connection)
    connection.commit()

    # update nz, height in grid
    debug('Updating nz, height in grid')
    sqltext = 'UPDATE "{0}"."{1}" AS g SET (height, nz) = (nt.height, nt.nz)' \
              '   FROM "nz_temp" AS nt ' \
              '   WHERE g.i=nt.i AND g.j = nt.j AND NOT is_building'.format(
                cfg.domain.case_schema, cfg.tables.grid)
    cur.execute(sqltext)

    # update nz, height in buildings_grid
    debug('Updating nz, height in buildings grid')
    sqltext = 'UPDATE "{0}"."{1}" AS b SET (height, nz) = (SELECT nt.height - bo.max_int - {2}, nt.nz - bo.max_int' \
              '   FROM "nz_temp" AS nt ' \
              '   LEFT OUTER JOIN "{0}"."{3}"  AS bo ON b.lid = bo.lid' \
              '   WHERE nt.is_building AND b.i=nt.i AND b.j = nt.j)'.format(
                cfg.domain.case_schema, cfg.tables.buildings_grid, cfg.domain.oro_min, cfg.tables.buildings_offset)
    cur.execute(sqltext)

    sql_debug(connection)
    connection.commit()

    debug('Dropping temp table')
    cur.execute('drop table "nz_temp"')

    filled_grids.sort()
    return filled_grids

def connect_buildings_height(cfg, connection, cur):
    """ Connection of raster building heights with grid and creating special buildings_grid"""

    change_log_level(cfg.logs.level_buildings)

    progress('Calculate building heights')
    debug('Creating table of buildings grid')
    sqltext = 'CREATE TABLE "{0}"."{1}" AS SELECT gg.id, gg.i, gg.j, gg.xcen, gg.ycen,' \
              '      {4} AS azimuth, 0.0 AS zenith, gg.geom, gg.lid, gg.type' \
              '      FROM (select g.id, g.i, g.j, g.xcen, g.ycen, g.geom, g.lid, l.type' \
              '      FROM "{0}"."{2}" as g left outer join "{0}"."{3}" as l on l.lid = g.lid' \
              '      WHERE type BETWEEN {5} AND {6}) AS gg'
    sqltext = sqltext.format(cfg.domain.case_schema, cfg.tables.buildings_grid,
                             cfg.tables.grid, cfg.tables.landcover,
                             cfg.fill_values.f8, cfg.type_range.building_min, cfg.type_range.building_max)
    cur.execute(sqltext)
    sql_debug(connection)
    verbose('Altering of table owner')
    sqltext = 'ALTER TABLE "{0}"."{1}" OWNER TO {2}'.format(
              cfg.domain.case_schema, cfg.tables.buildings_grid, cfg.pg_owner)
    cur.execute(sqltext)
    sql_debug(connection)
    verbose('Adding primary key')
    sqltext = 'ALTER TABLE "{0}"."{1}" ADD PRIMARY KEY (i,j,lid)'.format(
              cfg.domain.case_schema, cfg.tables.buildings_grid)
    cur.execute(sqltext)
    sql_debug(connection)
    verbose('Adding geom index')
    sqltext = 'CREATE INDEX buildings_geom_idx ON "{0}"."{1}" USING gist(geom)'.format(
              cfg.domain.case_schema, cfg.tables.buildings_grid)
    cur.execute(sqltext)
    verbose('Adding new column (height, nz), tobe filled')
    sqltext = 'alter table "{0}"."{1}" ' \
              ' add if not exists height double precision, ' \
              ' add if not exists nz integer'.format(
               cfg.domain.case_schema, cfg.tables.buildings_grid)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # connect raster height with building_heights grid
    # alternatively "where ST_Intersects(bp.geom, bg.geom)" instead of "where ST_Within(bp.geom, bg.geom)"
    debug('Updating building_grids heights from buildings raster')
    sqltext = 'update "{0}"."{1}" as bg ' \
              ' set height = ( select val from ( ' \
              ' select (ST_PixelAsPoints(b.rast)).geom as geom, (ST_PixelAsPoints(b.rast)).val as val ' \
              ' from "{0}"."{2}" b where ST_intersects(b.rast,  bg.geom) ) bp ' \
              ' where ST_Intersects(bp.geom, bg.geom) ' \
              ' order by ST_Distance(bp.geom, ST_SetSRID(ST_Point(bg.xcen,bg.ycen), %s)) ' \
              ' limit 1 ) where bg.height is null or bg.height = 0'\
              .format(cfg.domain.case_schema, cfg.tables.buildings_grid, cfg.tables.buildings_height)
    cur.execute(sqltext, (cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

    # fill missing heights by the value from the nearest filled gridpoint of the same building
    debug('Filling missing heights')
    sqltext = 'update "{0}"."{1}" as bg ' \
              ' set height = (' \
              '  select bn.height from "{0}"."{2}" as bn ' \
              '   where bn.lid = bg.lid and (bn.height != 0 or (bn.height is not null and bn.height != 0))' \
              '   order by ST_Distance(ST_SetSrid(ST_Point(bn.xcen,bn.ycen), %s), ' \
              '                        ST_SetSrid(ST_Point(bg.xcen,bg.ycen), %s)) ' \
              '   limit 1 ) ' \
              ' where bg.height is null or bg.height = 0'\
              .format(cfg.domain.case_schema, cfg.tables.buildings_grid, cfg.tables.buildings_grid)
    cur.execute(sqltext, (cfg.srid_palm, cfg.srid_palm, ))
    sql_debug(connection)
    connection.commit()

    # fill remaining missing heights by the value from the nearest filled gridpoint of any building to some distance
    debug('Fill remaining heights')
    sqltext2 = 'update "{0}"."{1}" as bg ' \
              ' set height = (' \
              '  select bn.height from "{0}"."{2}" as bn ' \
              '   where (bn.height != 0 or (bn.height is not null and bn.height != 0)) and ' \
              '         ST_Distance(ST_SetSrid(ST_Point(bn.xcen,bn.ycen), %s), ' \
              '                     ST_SetSrid(ST_Point(bg.xcen,bg.ycen), %s)) <= %s ' \
              '   order by ST_Distance(ST_SetSrid(ST_Point(bn.xcen,bn.ycen), %s), ' \
              '                        ST_SetSrid(ST_Point(bg.xcen,bg.ycen), %s)) ' \
              '   limit 1 ) ' \
              ' where bg.height is null or bg.height = 0 '\
              .format(cfg.domain.case_schema, cfg.tables.buildings_grid, cfg.tables.buildings_grid)
    cur.execute(sqltext2, (cfg.srid_palm, cfg.srid_palm, cfg.maxbuildingdisance, cfg.srid_palm, cfg.srid_palm, ))
    sql_debug(connection)
    connection.commit()

    # repeat once more filling of missing heights by the value
    # from the nearest filled gridpoint of the same building (to complete the buildings filled in previous step)
    verbose('Do it one more to fill grid points')
    cur.execute(sqltext, (cfg.srid_palm, cfg.srid_palm, ))
    sql_debug(connection)
    connection.commit()

    # # fill remain with default height
    # sqltext = 'UPDATE "{0}"."{1}" SET height = {2} WHERE height IS NULL OR height = 0.0'.format(
    #             cfg.domain.case_schema, cfg.tables.buildings_grid, cfg.default_height)
    # cur.execute(sqltext)
    # sql_debug(connection)
    # connection.commit()

    # fill remain missing with default according to its type
    debug('Fill all remaining missing heights using user configuration')
    sqltext = 'SELECT COUNT(*) FROM "{0}"."{1}" WHERE height IS NULL OR height = 0.0'.format(
                            cfg.domain.case_schema, cfg.tables.buildings_grid)
    cur.execute(sqltext)
    empty_buildings = cur.fetchone()[0]
    sql_debug(connection)
    connection.commit()
    if empty_buildings > 0:
        textdh = 'CASE '
        for t, dh in cfg.default_height._settings.items():
            textdh += 'WHEN type = {0} THEN {1} '.format(t+cfg.type_range.building_min, dh)
        textdh += ' ELSE {} END '.format(cfg.domain.dz)
        sqltext = 'UPDATE "{0}"."{1}" SET height = {2} WHERE height IS NULL OR height = 0'.format(
            cfg.domain.case_schema, cfg.tables.buildings_grid, textdh)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()


    # calculate nz
    debug('Calculation of nz')
    sqltext = 'UPDATE "{0}"."{1}" SET nz = cast(round(height/{2}+0.001) as integer)'.format(
                     cfg.domain.case_schema, cfg.tables.buildings_grid, cfg.domain.dz)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    verbose('add new column nz_min')
    sqltext = 'ALTER TABLE "{0}"."{1}" ' \
              ' ADD IF NOT EXISTS nz_min integer'\
              .format(cfg.domain.case_schema, cfg.tables.buildings_grid)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    verbose('update nz_min')
    sqltext = 'UPDATE "{0}"."{1}" SET nz_min = cast(0 as integer)'\
              .format(cfg.domain.case_schema, cfg.tables.buildings_grid)

    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()
    progress('2d buildings done')

    # terrain improvement
    progress('Terrain improvement')
    sqltext = 'drop table if exists "{0}"."{1}"'.format(cfg.domain.case_schema, cfg.tables.buildings_offset)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    sqltext = 'CREATE TABLE "{0}"."{1}" AS SELECT b.lid, ' \
              ' MAX(b.height-({3})) AS max, ' \
              ' CAST(ROUND(MAX(b.height-({3}))/{4}+0.001) AS INTEGER) AS max_int,' \
              ' MIN(b.height-({3})) AS min, ' \
              ' CAST(ROUND(MIN(b.height-({3}))/{4}+0.001) AS INTEGER) AS min_int, ' \
              ' MAX(b.height) - MIN(b.height) AS difference, ' \
              ' CAST(ROUND((MAX(b.height) - MIN(b.height))/{4}+0.001) AS INTEGER) AS difference_int, ' \
              '	COUNT(*) AS area, ' \
              ' SUM(b.nz) AS terr_sum, ' \
              ' ROUND(SUM(b.nz)*1.0/COUNT(*)*1.0) AS avg_terrain,' \
              ' ROUND(SUM(b.nz)*1.0/COUNT(*)*1.0)-CAST(ROUND(MIN(b.height-({3}))/{4} + 0.001) AS INTEGER) AS art_elev' \
              ' FROM "{0}"."{2}" AS b ' \
              ' GROUP BY b.lid ORDER BY b.lid' \
        .format(cfg.domain.case_schema, cfg.tables.buildings_offset, cfg.tables.grid, cfg.domain.oro_min, cfg.domain.dz)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    sqltext = 'ALTER TABLE "{}"."{}" OWNER TO {}' \
        .format(cfg.domain.case_schema, cfg.tables.buildings_offset, cfg.pg_owner)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # add art_elev to build height
    sqltext = 'UPDATE "{0}"."{1}" AS b SET nz = nz + bo.art_elev, height = height + bo.art_elev*{3} ' \
              ' FROM "{0}"."{2}" AS bo ' \
              ' WHERE b.lid = bo.lid AND bo.art_elev > 0'.format(cfg.domain.case_schema,
              cfg.tables.buildings_grid, cfg.tables.buildings_offset, cfg.domain.dz)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # Adjust terrain height according to building_offset
    sqltext = 'UPDATE "{0}"."{1}" AS g SET nz = bo.min_int, height = bo.min + {4} ' \
              '  FROM "{0}"."{2}"  AS bo ' \
              '  LEFT JOIN "{0}"."{3}"  AS l ON l.lid = bo.lid ' \
              '  WHERE l.type BETWEEN %s AND %s AND g.lid = bo.lid'.format(
              cfg.domain.case_schema, cfg.tables.grid, cfg.tables.buildings_offset, cfg.tables.landcover, cfg.domain.oro_min )
    cur.execute(sqltext, (cfg.type_range.building_min, cfg.type_range.building_max, ))
    sql_debug(connection)
    connection.commit()

    # fill holes in grid
    debug('Filling holes in grid')
    filled_grids = [0,0]
    all_filled = []
    while len(filled_grids) > 0:
        filled_grids = filling_grid(cfg, connection, cur)
        all_filled.append(filled_grids)

    restore_log_level(cfg)

# Procedure creating static driver
def nc_create_file(filename):
    fp = Path(filename)
    if fp.exists():
        if fp.is_dir():
            error("Error: {} is existing directory!", filename)
            sys.exit(1)
        else:
            debug("Delete existing file: {}", filename)
            fp.unlink()
    # create new file
    try:
        ncfile = Dataset(filename, "w", format="NETCDF4")
        debug("Created: {}.", filename)
    except FileNotFoundError:
        error("Error. Could not create file: {}!", filename)
        ncfile = None
    return ncfile

def prepare_domain_extends(cfg, connection, cur):
    """ Calculate origins, transform origins to latitude and longitude """
    origin_x = cfg.domain.cent_x - cfg.domain.nx * cfg.domain.dx / 2.0
    cfg.domain._settings['origin_x'] = origin_x
    origin_y = cfg.domain.cent_y - cfg.domain.ny * cfg.domain.dy / 2.0
    cfg.domain._settings['origin_y'] = origin_y

    # calculate origin lat and lon
    sqltext = 'select ST_X(ST_Transform(ST_SetSRID(ST_Point(%s,%s),%s),%s)), ' \
              'ST_Y(ST_Transform(ST_SetSRID(ST_Point(%s,%s),%s),%s)) '
    cur.execute(sqltext, (origin_x, origin_y, cfg.srid_palm, cfg.srid_wgs84,
                          origin_x, origin_y, cfg.srid_palm, cfg.srid_wgs84,))
    origin_lon, origin_lat = list(cur.fetchone())
    cfg.domain._settings['origin_lon'] = origin_lon
    cfg.domain._settings['origin_lat'] = origin_lat
    debug('Domain origin x,y: {}, {}', origin_x, origin_y)
    debug('Domain origin lon,lat: {}, {}', origin_lon, origin_lat)

def nc_write_global_attributes(ncfile, cfg):
    debug("Writing global attributes to file...")
    ncfile.setncattr('Conventions', "CF-1.7")
    ncfile.setncattr("origin_x", cfg.domain.origin_x)
    ncfile.setncattr("origin_y", cfg.domain.origin_y)
    ncfile.setncattr("origin_z", cfg.domain.origin_z)
    ncfile.setncattr("origin_time", cfg.origin_time)
    ncfile.setncattr("origin_lat", cfg.domain.origin_lat)
    ncfile.setncattr("origin_lon", cfg.domain.origin_lon)
    ncfile.setncattr("acronym", cfg.ncprops.acronym)
    ncfile.setncattr("author", cfg.ncprops.author)
    ncfile.setncattr("campaign", cfg.ncprops.campaign)
    ncfile.setncattr("contact_person", cfg.ncprops.contact_person)
    ncfile.setncattr("creation_time", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    ncfile.setncattr("comment", cfg.ncprops.comment)
    ncfile.setncattr("data_content", cfg.ncprops.data_content)
    ncfile.setncattr("dependencies", cfg.ncprops.dependencies)
    ncfile.setncattr("institution", cfg.ncprops.institution)
    ncfile.setncattr("keywords", cfg.ncprops.keywords)
    ncfile.setncattr("location", cfg.ncprops.location)
    ncfile.setncattr("palm_version", cfg.ncprops.palm_version)
    ncfile.setncattr("references", cfg.ncprops.references)
    ncfile.setncattr("rotation_angle", cfg.ncprops.rotation_angle)
    ncfile.setncattr("site", cfg.ncprops.site)
    ncfile.setncattr("source", cfg.ncprops.source)
    ncfile.setncattr("version", cfg.ncprops.version)

def nc_write_crs(ncfile, cfg, connection, cur):
    debug("Writing crs to file...")
    sqltext = 'select srtext from "spatial_ref_sys" where srid = Find_SRID(%s, %s, %s)'
    cur.execute(sqltext,(cfg.domain.case_schema, cfg.tables.grid, 'geom'))
    srtext = cur.fetchone()[0]
    temp = ncfile.createVariable("crs", "i")

    # default info
    temp.long_name = "coordinate reference system"
    temp.grid_mapping_name = "transverse_mercator"
    temp.semi_major_axis = 6378137.0
    temp.inverse_flattening = 298.257222101
    temp.longitude_of_prime_meridian = 0.0
    temp.longitude_of_central_meridian = 15.0
    temp.latitude_of_projection_origin = 0.0
    temp.scale_factor_at_central_meridian = 0.9996
    temp.false_easting = 500000.0
    temp.false_northing = 0.0
    temp.units = 'm'
    temp.epsg_code = 'EPSG: 25833'

    # overwrite default info
    temp.long_name = "coordinate reference system"
    # PROJCS
    proj = srtext.split(sep="PROJECTION[")[1]
    temp.grid_mapping_name = proj.split('],')[0].strip('"')
    # PARAMETERs
    params = proj.split('PARAMETER[')
    for i in range(1, len(params)):
        param = params[i].split(']')[0].split(',')
        if param[0].strip('"') == 'latitude_of_origin':
            temp.latitude_of_projection_origin = float(param[1].strip('"'))
        elif param[0].strip('"') == 'scale_factor':
            temp.scale_factor_at_central_meridian = float(param[1].strip('"'))
        else:
            temp.setncattr(param[0].strip('"'), float(param[1].strip('"')))
    # UNIT
    temp.units = proj.split('UNIT[')[1].split('],')[0].split(',')[0].strip('"')
    # AUTHORITY
    auths = proj.split('AUTHORITY[')
    auth = auths[len(auths) - 1].split(']')[0].split(',')
    authstr = auth[0].replace('"', '') + ':' + auth[1].replace('"', '')
    temp.epsg_code = authstr


def nc_create_dimension(ncfile, dimname, dimlen):
    try:
        debug("Creating dimension {}", dimname)
        ncfile.createDimension(dimname, dimlen)
        return 0
    except:
        return 0

def nc_create_variable(ncfile, var_name, precision, dims, fill_value = None ):
    try:
        ncfile.createVariable(var_name, precision, dims, fill_value = fill_value )
    except:
        pass
    return ncfile.variables[var_name]

def nc_write_attribute(ncfile, variable, attribute, value):
    if not hasattr(ncfile[variable], attribute):
        var = ncfile.variables[variable]
        var.setncattr(attribute, value)
    return 0


def create_dim_xy(ncfile, cfg, connection, cur):
    """ """
    sqltext = 'select distinct xcen from "{0}"."{1}" order by xcen'.format(cfg.domain.case_schema, cfg.tables.grid)
    cur.execute(sqltext)
    x1d = [x[0] - cfg.domain.origin_x for x in cur.fetchall()]
    sqltext = 'select distinct ycen from "{0}"."{1}" order by ycen'.format(cfg.domain.case_schema, cfg.tables.grid)
    cur.execute(sqltext)
    sql_debug(connection)
    y1d = [y[0] - cfg.domain.origin_y for y in cur.fetchall()]
    nxm, nym = len(x1d), len(y1d)

    debug("Writing 2D variables x, y to file...")
    nc_create_dimension(ncfile, 'x', nxm)
    nc_create_dimension(ncfile, 'y', nym)
    vt = 'f8'
    temp_x = ncfile.createVariable('x', vt, 'x')
    temp_y = ncfile.createVariable('y', vt, 'y')
    temp_x[:] = x1d[:]
    temp_y[:] = y1d[:]
    del x1d, y1d
    nc_write_attribute(ncfile, 'x', 'long_name', 'x')
    nc_write_attribute(ncfile, 'x', 'standard_name', 'projection_x_coordinate')
    nc_write_attribute(ncfile, 'x', 'units', 'm')
    nc_write_attribute(ncfile, 'y', 'long_name', 'y')
    nc_write_attribute(ncfile, 'y', 'standard_name', 'projection_y_coordinate')
    nc_write_attribute(ncfile, 'y', 'units', 'm')

    ################################
    # transform xcen, ycen coordinates to lat,lon and E_UTM, N_UTM coordinates
    # and store into coresponding netcdfs variables
    debug("Writing 2D variables lon, lat, E_UTM, N_UTM to file...")
    sqltext = 'select lon, lat, "E_UTM", "N_UTM" from "{0}"."{1}" order by j,i'.format(cfg.domain.case_schema, cfg.tables.grid)
    cur.execute(sqltext)
    res = cur.fetchall()
    vt = 'f4'
    nc_create_dimension(ncfile, 'lon', nxm)
    var = ncfile.createVariable('lon', vt, ('y', 'x'), fill_value=cfg.fill_values[vt])
    var[:, :] = np.reshape(np.asarray([x[0] for x in res], dtype=vt), (nym, nxm))
    nc_create_dimension(ncfile, 'lat', nym)
    var = ncfile.createVariable('lat', vt, ('y', 'x'), fill_value=cfg.fill_values[vt])
    var[:, :] = np.reshape(np.asarray([x[1] for x in res], dtype=vt), (nym, nxm))
    vt = 'f8'
    nc_create_dimension(ncfile, 'E_UTM', nxm)
    var = ncfile.createVariable('E_UTM', vt, ('y', 'x'), fill_value=cfg.fill_values[vt])
    var[:, :] = np.reshape(np.asarray([x[2] for x in res], dtype=vt), (nym, nxm))
    nc_create_dimension(ncfile, 'N_UTM', nym)
    var = ncfile.createVariable('N_UTM', vt, ('y', 'x'), fill_value=cfg.fill_values[vt])
    var[:, :] = np.reshape(np.asarray([x[3] for x in res], dtype=vt), (nym, nxm))
    del res

    # write needed attributes
    nc_write_attribute(ncfile, 'lat', 'long_name', 'latitude')
    nc_write_attribute(ncfile, 'lat', 'standard_name', 'latitude')
    nc_write_attribute(ncfile, 'lat', 'units', 'degrees_north')
    nc_write_attribute(ncfile, 'lon', 'long_name', 'longitude')
    nc_write_attribute(ncfile, 'lon', 'standard_name', 'longitude')
    nc_write_attribute(ncfile, 'lon', 'units', 'degrees_east')
    nc_write_attribute(ncfile, 'E_UTM', 'long_name', 'easting')
    nc_write_attribute(ncfile, 'E_UTM', 'standard_name', 'projection_x_coorindate')
    nc_write_attribute(ncfile, 'E_UTM', 'units', 'm')
    nc_write_attribute(ncfile, 'N_UTM', 'long_name', 'northing')
    nc_write_attribute(ncfile, 'N_UTM', 'standard_name', 'projection_y_coorindate')
    nc_write_attribute(ncfile, 'N_UTM', 'units', 'm')

    for name in cfg.ndims._settings.keys():
        d = cfg.ndims[name]
        nc_create_dimension(ncfile, name, d)
        temp = ncfile.createVariable(name, 'i', name)
        temp[:] = np.asarray(list(range(d)))

def write_terrain(ncfile, cfg, connection, cur):
    # write topography
    nx = cfg.domain.nx
    ny = cfg.domain.ny
    # artificial adjustment, according to fact that child domain and parent domain share the same origin
    if cfg.domain.origin_z > cfg.domain.oro_min:
        error('Origin z [{} m] is higher that oro min [{} m]', cfg.domain.origin_z, cfg.domain.oro_min)
    nesting_adjust = ceil((cfg.domain.oro_min - cfg.domain.origin_z) / cfg.domain.dz)

    sqltext = 'select (nz+{3})*{2} from "{0}"."{1}" order by j,i'.format(
              cfg.domain.case_schema, cfg.tables.grid, cfg.domain.dz, nesting_adjust)
    cur.execute(sqltext)
    res = cur.fetchall()
    sql_debug(connection)
    vn = 'zt'
    vt = 'f4'
    var = ncfile.createVariable(vn, vt, ('y', 'x'), fill_value=cfg.fill_values[vt])
    var[:, :] = np.reshape(np.asarray([x[0] for x in res], dtype=vt), (ny, nx))
    del res
    if cfg.visual_check.enabled:
        variable_visualization(var=var,
                               x=np.asarray(ncfile.variables['x']), y=np.asarray(ncfile.variables['y']),
                               var_name=vn, par_id='', text_id='terrain_height', path=cfg.visual_check.path,
                                   show_plots=cfg.visual_check.show_plots)


def write_type_variable(ncfile, cfg, vn, vln, vt, fill, lod, connection, cur):
    res = cur.fetchall()
    sql_debug(connection)
    ncfile.createVariable(vn, vt, ('y', 'x'), fill_value=fill[vt])
    res = [fill[vt] if x[0] is None else x[0] for x in res]
    var = np.reshape(np.asarray([x for x in res], dtype=vt), (cfg.domain.ny, cfg.domain.nx))
    del res
    var = np.nan_to_num(var, copy=False, nan=fill[vt], posinf=fill[vt], neginf=fill[vt])
    ncfile[vn][...] = var
    nc_write_attribute(ncfile, vn, 'long_name', vln)
    nc_write_attribute(ncfile, vn, 'units', '')
    nc_write_attribute(ncfile, vn, 'res_orig', cfg.domain.dz)
    if lod is not None and lod != '' and lod != 0:
        nc_write_attribute(ncfile, vn, 'lod', lod)
    nc_write_attribute(ncfile, vn, 'coordinates', 'E_UTM N_UTM lon lat')
    nc_write_attribute(ncfile, vn, 'grid_mapping', 'E_UTM N_UTM lon lat')
    debug('Variable {} ({}) has been written.', vn, vln)

def write_pavements(ncfile, cfg, connection, cur):
    # write landcover
    vn = 'pavement_type'
    vt = 'b'
    sqltext = 'select l.type-%s from "{0}"."{1}" g ' \
              'left outer join "{0}"."{2}" l on l.lid = g.lid and l.type >= %s and l.type < %s ' \
              'order by g.j, g.i' \
        .format(cfg.domain.case_schema, cfg.tables.grid, cfg.tables.landcover)
    cur.execute(sqltext, (cfg.type_range.pavement_min, cfg.type_range.pavement_min, cfg.type_range.pavement_max,))
    sql_debug(connection)
    write_type_variable(ncfile, cfg, vn, 'pavement type', vt, cfg.fill_values, 0, connection, cur)

    if cfg.visual_check.enabled:
        variable_visualization(var=ncfile[vn],
                               x=np.asarray(ncfile.variables['x']), y=np.asarray(ncfile.variables['y']),
                               var_name=vn, par_id='', text_id='pavement_type', path=cfg.visual_check.path,
                               show_plots=cfg.visual_check.show_plots)

def write_water(ncfile, cfg, connection, cur):
    # write water surfaces
    vn = "water_type"
    vt = 'b'
    sqltext = 'select l.type-%s from "{0}"."{1}" g left outer join "{0}"."{2}" l ' \
              ' on l.lid = g.lid and l.type >= %s and l.type < %s order by g.j, g.i' \
        .format(cfg.domain.case_schema, cfg.tables.grid, cfg.tables.landcover)
    cur.execute(sqltext, (cfg.type_range.water_min, cfg.type_range.water_min,
                          cfg.type_range.water_max,))
    sql_debug(connection)
    write_type_variable(ncfile, cfg, vn, 'water type', vt, cfg.fill_values, 0, connection, cur)
    if cfg.visual_check.enabled:
        variable_visualization(var=ncfile[vn],
                               x=np.asarray(ncfile.variables['x']), y=np.asarray(ncfile.variables['y']),
                               var_name=vn, par_id='', text_id='water_type', path=cfg.visual_check.path,
                               show_plots = cfg.visual_check.show_plots)

    # process water pars - water body temperature
    vn = "water_pars"
    vt = 'f4'
    ncfile.createVariable(vn, vt, ('nwater_pars', 'y', 'x'), fill_value=cfg.fill_values[vt])
    nc_write_attribute(ncfile, vn, 'long_name', 'water parameters')
    nc_write_attribute(ncfile, vn, 'units', '1')
    nc_write_attribute(ncfile, vn, 'source', '')
    nc_write_attribute(ncfile, vn, 'res_orig', cfg.domain.dz)
    nc_write_attribute(ncfile, vn, 'coordinates', 'E_UTM N_UTM lon lat')
    nc_write_attribute(ncfile, vn, 'grid_mapping', 'E_UTM N_UTM lon lat')

    # specific treatment of water body
    sqltext = 'CASE  '
    for wtype, wtemp in cfg.water_pars_temp._settings.items():
        sqltext += 'WHEN type = {0} THEN {1} '.format(wtype + cfg.type_range.water_min, wtemp)
    sqltext += 'ELSE NULL END'

    sqltext = 'SELECT {0} ' \
              'FROM "{1}"."{2}" AS g ' \
              'LEFT OUTER JOIN "{1}"."{3}" AS l ON l.lid = g.lid AND l.type >= %s AND l.type < %s ' \
              'ORDER BY g.j, g.i' \
        .format(sqltext, cfg.domain.case_schema, cfg.tables.grid,
                cfg.tables.landcover, )
    cur.execute(sqltext, (cfg.type_range.water_min, cfg.type_range.water_max,))
    res = cur.fetchall()
    sql_debug(connection)
    res = [cfg.fill_values[vt] if x[0] is None else x[0] for x in res]
    var = np.reshape(np.asarray(res, dtype=vt), (cfg.domain.ny, cfg.domain.nx))
    del res
    var = np.nan_to_num(var, copy=False, nan=cfg.fill_values[vt],
                        posinf=cfg.fill_values[vt], neginf=cfg.fill_values[vt])
    ncfile.variables[vn][0, ...] = var
    debug('Variable {}, parameter {} written.', vn, 0)

    if cfg.visual_check.enabled:
        variable_visualization(var=ncfile.variables[vn][0, ...],
                               x=np.asarray(ncfile.variables['x']), y=np.asarray(ncfile.variables['y']),
                               var_name=vn, par_id=0, text_id='water_pars', path=cfg.visual_check.path,
                               show_plots = cfg.visual_check.show_plots)

def write_vegetation(ncfile, cfg, connection, cur):
    # write vegetation surfaces
    vn = "vegetation_type"
    vt = 'b'
    sqltext = 'select l.type-%s from "{0}"."{1}" g left outer join "{0}"."{2}" l ' \
              ' on l.lid = g.lid and l.type >= %s and l.type < %s order by g.j, g.i' \
        .format(cfg.domain.case_schema, cfg.tables.grid, cfg.tables.landcover)
    cur.execute(sqltext, (cfg.type_range.vegetation_min, cfg.type_range.vegetation_min, cfg.type_range.vegetation_max,))
    sql_debug(connection)
    write_type_variable(ncfile, cfg, vn, 'vegetation type', vt, cfg.fill_values, 0, connection, cur)
    if cfg.visual_check.enabled:
        variable_visualization(var=ncfile[vn],
                               x=np.asarray(ncfile.variables['x']), y=np.asarray(ncfile.variables['y']),
                               var_name=vn, par_id='', text_id='vegetation_type', path=cfg.visual_check.path,
                               show_plots = cfg.visual_check.show_plots)

def write_soil(ncfile, cfg, connection, cur):
    # write soil type for vegetation surfaces
    # TODO: used one default type of the soil so far - get specific soil type for landcover database !!!
    vn = "soil_type"
    vt = 'b'
    sqltext = 'select case when l.type is not null then %s else null end ' \
              'from "{0}"."{1}" g left outer join "{0}"."{2}" l ' \
              ' on l.lid = g.lid and ((l.type >= %s and l.type < %s) or (l.type >= %s and l.type < %s)) ' \
              'order by g.j, g.i' \
        .format(cfg.domain.case_schema, cfg.tables.grid, cfg.tables.landcover)
    cur.execute(sqltext, (cfg.ground.soil_type_default, cfg.type_range.vegetation_min, cfg.type_range.vegetation_max,
                          cfg.type_range.pavement_min, cfg.type_range.pavement_max,))
    sql_debug(connection)
    write_type_variable(ncfile, cfg, vn, 'soil type', vt, cfg.fill_values, 1, connection, cur)
    if cfg.visual_check.enabled:
        variable_visualization(var=ncfile[vn],
                               x=np.asarray(ncfile.variables['x']), y=np.asarray(ncfile.variables['y']),
                               var_name=vn, par_id='', text_id='soil_type', path=cfg.visual_check.path,
                               show_plots = cfg.visual_check.show_plots)

def write_buildings(ncfile, cfg, connection, cur):
    # write building_height (buildings_2d), building_id and building_type into netcdf file
    # building id
    sqltext = 'select b.lid from "{0}"."{1}" g ' \
              'left outer join "{0}"."{2}" b on b.id = g.id order by g.j, g.i' \
        .format(cfg.domain.case_schema, cfg.tables.grid, cfg.tables.buildings_grid)
    cur.execute(sqltext)
    sql_debug(connection)
    vn = 'building_id'
    vt = 'i'
    write_type_variable(ncfile, cfg, vn, vn, vt, cfg.fill_values, 0, connection, cur)
    if cfg.visual_check.enabled:
        variable_visualization(var=ncfile[vn],
                               x=np.asarray(ncfile.variables['x']), y=np.asarray(ncfile.variables['y']),
                               var_name=vn, par_id='', text_id='building_id', path=cfg.visual_check.path,
                               show_plots = cfg.visual_check.show_plots)

    # building height
    sqltext = 'select b.nz*{3} from "{0}"."{1}" g ' \
              'left outer join "{0}"."{2}" b on b.id = g.id order by g.j, g.i' \
        .format(cfg.domain.case_schema, cfg.tables.grid, cfg.tables.buildings_grid, cfg.domain.dz)
    cur.execute(sqltext)
    sql_debug(connection)
    vn = 'buildings_2d'
    vt = 'f4'
    write_type_variable(ncfile, cfg, vn, vn, vt, cfg.fill_values, 1, connection, cur)
    if cfg.visual_check.enabled:
        variable_visualization(var=ncfile[vn],
                               x=np.asarray(ncfile.variables['x']), y=np.asarray(ncfile.variables['y']),
                               var_name=vn, par_id='', text_id='building_2d', path=cfg.visual_check.path,
                               show_plots = cfg.visual_check.show_plots)

    # building_type
    sqltext = 'select case when b.id is not null then ' \
              ' case when l.type >= %s and l.type < %s then l.type-%s ' \
              ' else 1 end else null end as type ' \
              'from "{0}"."{1}" g ' \
              'left outer join "{0}"."{2}" b on b.id = g.id ' \
              'left outer join "{0}"."{3}" l on l.lid = g.lid ' \
              'order by g.j, g.i ' \
        .format(cfg.domain.case_schema, cfg.tables.grid, cfg.tables.buildings_grid, cfg.tables.landcover)
    cur.execute(sqltext, (cfg.type_range.building_min, cfg.type_range.building_max, cfg.type_range.building_min,))
    sql_debug(connection)
    vn = 'building_type'
    vt = 'b'
    write_type_variable(ncfile, cfg, vn, vn, vt, cfg.fill_values, 0, connection, cur)
    if cfg.visual_check.enabled:
        variable_visualization(var=ncfile[vn],
                               x=np.asarray(ncfile.variables['x']), y=np.asarray(ncfile.variables['y']),
                               var_name=vn, par_id='', text_id='building_type', path=cfg.visual_check.path,
                               show_plots = cfg.visual_check.show_plots)

def variable_visualization(var, x, y, var_name, par_id, text_id, path, scale=15, color_name='gist_rainbow', show_plots=False): #hsv
    """ create 2d plot from var 2d matrix"""
    # TODO: add posibility to directly add X, Y = lat, lon
    if not show_plots:
        plt.ioff()
    ny, nx = var.shape
    fig = plt.figure()
    ax = plt.subplot(111)
    X, Y = np.meshgrid(x, y)

    ax_sl = ax.imshow(var, aspect='equal', cmap=plt.get_cmap(color_name), origin='lower',
                      )
    fig.colorbar(ax_sl, extend='max', orientation = 'horizontal')

    plt.title('{}: {}\n'
              '[{}], '
              'min: {:.2f}, max: {:.2f}, avg: {:.2f}, stdev: {:.2f}'.format(
                var_name, par_id, text_id,
                np.min(var), np.max(var), np.average(var), np.std(var)))

    fig.savefig(os.path.join(path,'{}_{}.png'.format(var_name, par_id)),
                dpi=400)
    debug('{}: {} was successfully plotted', var_name, par_id)

def check_consistency(ncfile, cfg):
    """ Check whether there are no grid points without any type """
    change_log_level(cfg.logs.level_check_consistency)
    progress('Checking consistency ...')
    pavement_type_default = 2
    mask = (ncfile.variables['vegetation_type'][:,:].mask & ncfile.variables['pavement_type'][:,:].mask & \
            ncfile.variables['building_type'][:,:].mask & ncfile.variables['water_type'][:,:].mask )
    missing_values = np.sum(mask)
    if missing_values > 0:
        warning('There are {} missing values that were filled with default pavement type {}',
                missing_values, pavement_type_default)
        extra_verbose('Missing grid: {}', np.where(mask))

    pt = ncfile.variables['pavement_type'][:, :]
    pt = np.where(mask, pavement_type_default, pt)
    ncfile.variables['pavement_type'][:, :] = pt

    # soil type
    soil_type_default = 3
    mask = np.logical_and(np.logical_or(~ncfile.variables['vegetation_type'][:,:].mask, ~ncfile.variables['pavement_type'][:,:].mask),
                                         ncfile.variables['soil_type'][:,:].mask)
    missing_soils = np.sum(mask)
    if missing_soils > 0:
        warning('There are {} missing soil values that were filled with default soil type {}',
             missing_soils, soil_type_default)
        extra_verbose('Missing grids: {}', np.where(mask))
    pt = ncfile.variables['soil_type'][:, :]
    pt[mask] = soil_type_default
    ncfile.variables['soil_type'][:, :] = pt
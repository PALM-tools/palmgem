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

import numpy as np
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
    sqltext = 'SELECT ST_MakeEnvelope(%s, %s, %s, %s, %s)'
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

def check_buildings(cfg, connection, cur, rtabs, grid_ext):
    """ Check if USM are present in domain """
    verbose('Checking if buildings raster is present in inputs')
    cfg._settings['has_buildings'] = False
    if cfg.tables.buildings_height in rtabs:
        cfg._settings['has_buildings'] = True
    debug('Modification of buildings that intersect with domain boundary or are adjacet to domain extent')
    if cfg.force_building_boundary:
        # TODO: join with nearest adjacent landcover
        sqltext = 'UPDATE "{0}"."{1}" AS l SET type = 202 ' \
                  'WHERE type >= {2} AND ST_Distance(l.geom, ST_Boundary(%s::geometry)) < {3}'\
                  .format(cfg.domain.case_schema, cfg.tables.landcover, cfg.type_range.building_min,
                          cfg.force_building_boundary_dist * cfg.domain.dx )
        cur.execute(sqltext, (grid_ext,))
        sql_debug(connection)
        connection.commit()

    debug('Check if buildings grid will be defined in more that 3 grid cells')
    # TODO: join with nearest adjacent landcover
    sqltext = 'UPDATE "{0}"."{1}" AS l SET type = 202 ' \
              'WHERE type >= {2} AND ST_Area(geom) < {3}'\
              .format(cfg.domain.case_schema, cfg.tables.landcover, cfg.type_range.building_min,
                      3.0 ** 2 * cfg.domain.dx ** 2 )
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    verbose('Checking if buildings (type >= 900) are present in landover')
    sqltext = 'SELECT COUNT(*) FROM "{0}"."{1}" WHERE type BETWEEN {2} AND {3}'\
              .format(cfg.domain.case_schema, cfg.tables.landcover,
                      cfg.type_range.building_min, cfg.type_range.building_max,)
    cur.execute(sqltext)
    buildings_counts = cur.fetchone()[0]
    sql_debug(connection)
    connection.commit()
    if buildings_counts > 0:
        cfg._settings['has_buildings'] = True

    if cfg.force_lsm_only:
        if cfg.tables.buildings_height in rtabs:
            tab = cfg.tables.buildings_height
            warning('Because of force_lsm_only, delete table {} from case', tab)
            rtabs.remove(tab)
            cur.execute('DROP TABLE "{}"."{}"'.format(cfg.domain.case_schema, tab))
            sql_debug(connection)
            connection.commit()

        debug('Modification of all building types to LSM type')
        sqltext = 'UPDATE "{0}"."{1}" AS l SET type = 202 ' \
                  'WHERE type >= {2}'.format(cfg.domain.case_schema, cfg.tables.landcover, cfg.type_range.building_min)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

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
    # calculate origin_z as bottom of the current or parent domain
    sqltext = 'select min(height) from "{}".{}'. \
        format(cfg.domain.case_schema if cfg.domain.parent_domain_schema == ''
               else cfg.domain.parent_domain_schema, cfg.tables.grid)
    cur.execute(sqltext)
    cfg.domain._settings['origin_z'] = cur.fetchone()[0]
    cfg.domain._settings['oro_min'] = cfg.domain.origin_z
    debug('Origin_z is in height level: {}', cfg.domain.origin_z)

    # calculate grid nz height
    sqltext = 'UPDATE "{0}"."{1}" SET nz = CAST((height - {2})/{3} AS INTEGER)'.format(
        cfg.domain.case_schema, cfg.tables.grid, cfg.domain.oro_min, cfg.domain.dz)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

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

def fill_cortyard(cfg, connection, cur):
    """ fill the cortyard that are surrounded by buildings and their size in grid is lower that user defined value """
    progress('Filling cortyard')
    debug('Finding all suspicious polygons')
    sqltext = 'WITH ll AS ( ' \
              ' SELECT lid  ' \
              ' FROM "{0}"."{1}" AS l' \
              ' WHERE (SELECT COUNT(*) FROM "{0}"."{1}" AS nb WHERE type BETWEEN 0 AND 899 AND ST_Touches(nb.geom, l.geom) AND nb.lid != l.lid) = 0  ' \
              '        ' \
              '        AND l.type BETWEEN 0 AND 899' \
              ' ) ' \
              'SELECT g.lid ' \
              'FROM "{0}"."{2}" AS g, ll ' \
              'WHERE g.lid IN (ll.lid) ' \
              'GROUP BY g.lid ' \
              'HAVING COUNT(*) < {3}'.format(cfg.domain.case_schema, cfg.tables.landcover, cfg.tables.grid, cfg.cortyard_fill.count)
    cur.execute(sqltext)
    lids = cur.fetchall()
    sql_debug(connection)
    connection.commit()

    lids_list = [x[0] for x in lids]

    if len(lids_list) == 0:
        return

    debug('Modification of those cortyard into nearest building type')
    verbose('In landcover')
    sqltext = 'UPDATE "{0}"."{1}" AS l SET (type, katland, albedo, emisivita) = ' \
              ' (SELECT type, katland, albedo, emisivita FROM "{0}"."{1}" AS b ' \
              '   WHERE type BETWEEN 900 AND 999 ' \
              '   ORDER BY ST_Distance(l.geom, b.geom) LIMIT 1) ' \
              'WHERE lid IN {2}'.format(cfg.domain.case_schema, cfg.tables.landcover, tuple(lids_list))
    cur.execute(sqltext,)
    sql_debug(connection)
    connection.commit()

def fill_missing_holes_in_grid(cfg, connection, cur):
    """ Fill empty 3d grid cell that has more than 3 neighbors, replicating process in PALM program to omit creating grid cell with default configuration """
    if cfg.force_lsm_only:
        return
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

def fill_topo_v2(cfg, connection, cur):
    """ Fill all topologies that satisfies condition:
        Connected ${topo_fill_v2_count} air grid cell in one k-th layer
        Algorithm:
        While something was adjusted loop:
            Loop k = 0 .. nz_max+1
                Create temp table with joined terrain and buildings
                Filter and Union all adjacent grids with nz = k -> table grid_k
                Filter and Union all adjacent grids with nz > k -> table grid_higher
                Select all polygons with ST_Area() <= dx*dy*{topo_fill_v2_count} that are surrounded by only grid_higher
                    Do intersection with temp table and increase height in those {j,i} grid cells.
                Save info about which ones {j,i} were filled

        Move temp table back do grid and buildings grid
        Delete temp table and grid_k, grid_higher
    """
    progress('Filling topology using fill_topo_v2')
    debug('Create nz_temp')
    sqltext = 'CREATE TABLE "{2}"."nz_temp" AS ' \
              'SELECT g.id, g.i AS i, g.j AS j, g.geom AS geom,  ' \
              'CASE WHEN b.height IS NOT NULL THEN b.nz+bo.max_int' \
              '     ELSE g.nz END AS nz, ' \
              'CASE WHEN b.height IS NOT NULL THEN true ' \
              '     ELSE false END AS is_building ' \
              'FROM "{2}"."{3}" AS g ' \
              'LEFT OUTER JOIN "{2}"."{4}" AS b ON g.id=b.id  ' \
              'LEFT OUTER JOIN "{2}"."{5}" AS bo ON b.lid = bo.lid'.format(
        cfg.domain.oro_min, cfg.domain.dz, cfg.domain.case_schema, cfg.tables.grid,
        cfg.tables.buildings_grid, cfg.tables.buildings_offset, cfg.domain.oro_min)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    debug('Add primary keys to nz_temp')
    sqltext = 'ALTER TABLE "{0}"."nz_temp" ADD primary key (id)'.format(cfg.domain.case_schema)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    debug('Create geom index on nz_temp')
    sqltext = 'CREATE INDEX IF NOT EXISTS nz_temp_geom_idx on "{0}"."nz_temp" using gist(geom)' \
        .format(cfg.domain.case_schema)
    cur.execute(sqltext)
    sql_debug(connection)

    debug('Start filling grids')
    verbose('Find max nz')
    sqltext = 'SELECT MAX(nz) FROM "{0}"."nz_temp"'.format(cfg.domain.case_schema)
    cur.execute(sqltext)
    nz_max = cur.fetchone()[0]
    sql_debug(connection)
    connection.commit()

    sqltext_k =      'DROP TABLE IF EXISTS "{0}".k; ' \
                     'CREATE TABLE "{0}".k AS  ' \
                     'SELECT (ST_Dump(ST_Union(geom))).geom ' \
                     'FROM "{0}"."nz_temp" ' \
                     'WHERE nz={1}'
    sqltext_kplus =  'DROP TABLE IF EXISTS "{0}".kplus; ' \
                     'CREATE TABLE "{0}".kplus AS  ' \
                     'SELECT (ST_Dump(ST_Union(geom))).geom ' \
                     'FROM "{0}"."nz_temp" ' \
                     'WHERE nz>{1}'
    sqltext_kminus = 'DROP TABLE IF EXISTS "{0}".kminus; ' \
                     'CREATE TABLE "{0}".kminus AS  ' \
                     'SELECT (ST_Dump(ST_Union(geom))).geom ' \
                     'FROM "{0}"."nz_temp" ' \
                     'WHERE nz<{1}'
    sqltext_gist   = 'CREATE INDEX IF NOT EXISTS {1}_geom_idx on "{0}"."{1}" using gist(geom)'
    sqltext_find   = 'WITH gg AS (SELECT geom FROM "{0}".k AS gkk ' \
                     '            WHERE ST_Area(gkk.geom) < {1} AND ' \
                     '                  (SELECT COUNT(*) AS count ' \
                     '                   FROM "{0}".kminus AS gkm ' \
                     '                   WHERE ST_Touches(gkk.geom, gkm.geom)) = 0 ) ' \
                     'SELECT id, i, j, nz FROM "{0}"."nz_temp" AS g, gg ' \
                     'WHERE ST_Intersects(ST_Centroid(g.geom), gg.geom)'\
                     .format(cfg.domain.case_schema, cfg.domain.dx * cfg.domain.dy * cfg.topo_fill_v2.count+0.1)

    sqltext_update = 'UPDATE "{1}"."nz_temp" SET ' \
                     'nz = nz + 1 ' \
                     'WHERE id IN {0}'
    sqltext_update_1 = 'UPDATE "{1}"."nz_temp" SET ' \
                       'nz = nz + 1 ' \
                       'WHERE id = {0}'
    sqltext_count  = 'SELECT COUNT(*) FROM "{0}"."{1}"'

    filled_grids = []
    fillings = 0
    for k in range(nz_max+2):
        verbose('\tk: {}', k)
        extra_verbose('\tCreate k table, with an geom index')
        cur.execute(sqltext_k.format(cfg.domain.case_schema, k))
        cur.execute(sqltext_gist.format(cfg.domain.case_schema, 'k'))
        cur.execute(sqltext_count.format(cfg.domain.case_schema, 'k'))
        k_count = cur.fetchone()[0]
        if k_count == 0:
            continue
        extra_verbose('\tCreate kminus table')
        cur.execute(sqltext_kminus.format(cfg.domain.case_schema, k))
        cur.execute(sqltext_gist.format(cfg.domain.case_schema, 'kminus'))
        # extra_verbose('\tCreate kplus table')
        # cur.execute(sqltext_kplus.format(cfg.domain.case_schema, k))
        # cur.execute(sqltext_gist.format(cfg.domain.case_schema, 'kplus'))
        sql_debug(connection)
        connection.commit()
        extra_verbose('\tFind grids to update')
        cur.execute(sqltext_find, (cfg.srid_palm,))
        missings = cur.fetchall()
        if len(missings) == 0:
            continue
        filled_grids.append([missings])
        fillings += len(missings)
        id, i, j, nz = [x[0] for x in missings], [x[1] for x in missings], [x[2] for x in missings], [x[3] for x in missings]
        for ii, jj in zip(i, j):
            extra_verbose('In filling grid; filling [j, i] = [{},{}] ', jj, ii)
        extra_verbose('\tCorrection of missing grid height in {} grid points', len(missings))
        if len(i) == 1:
            cur.execute(sqltext_update_1.format(id[0], cfg.domain.case_schema))
        else:
            cur.execute(sqltext_update.format(tuple(id), cfg.domain.case_schema))
        sql_debug(connection)
        connection.commit()

    verbose('{} number of grid has been filled', fillings)

    debug('Done with filling')
    sql_debug(connection)
    connection.commit()

    # debug('Add primary keys to nz_temp')
    # sqltext = 'ALTER TABLE "{0}"."nz_temp" ADD primary key (i,j)'.format(cfg.domain.case_schema)
    # cur.execute(sqltext)
    # sql_debug(connection)
    # connection.commit()

    # sqltext = 'ALTER TABLE "{0}"."nz_temp" DROP CONSTRAINT id'.format(cfg.domain.case_schema, rel, prev_ui)

    # update nz, height in grid
    debug('Updating nz, height in grid')
    sqltext = 'UPDATE "{0}"."{1}" AS g SET nz = nt.nz' \
              '   FROM "{0}"."nz_temp" AS nt ' \
              '   WHERE g.id=nt.id AND NOT is_building'.format(
                cfg.domain.case_schema, cfg.tables.grid)
    cur.execute(sqltext)

    # update nz, height in buildings_grid
    debug('Updating nz, height in buildings grid')
    sqltext = 'UPDATE "{0}"."{1}" AS b SET nz = (SELECT nt.nz - bo.max_int' \
              '   FROM "{0}"."nz_temp" AS nt ' \
              '   LEFT OUTER JOIN "{0}"."{3}"  AS bo ON b.lid = bo.lid' \
              '   WHERE nt.is_building AND b.id = nt.id) ' \
              .format(cfg.domain.case_schema, cfg.tables.buildings_grid, cfg.domain.oro_min, cfg.tables.buildings_offset)
    cur.execute(sqltext)

    sql_debug(connection)
    connection.commit()

    debug('Dropping temp table')
    cur.execute('DROP TABLE "{0}"."nz_temp"'.format(cfg.domain.case_schema))
    cur.execute('DROP TABLE "{0}".k'.format(cfg.domain.case_schema))
    cur.execute('DROP TABLE "{0}".kminus'.format(cfg.domain.case_schema))
    # cur.execute('DROP TABLE "{0}".kplus'.format(cfg.domain.case_schema))
    sql_debug(connection)
    connection.commit()

    filled_grids.sort()
    return filled_grids

def connect_buildings_height(cfg, connection, cur):
    """ Connection of raster building heights with grid and creating special buildings_grid"""
    if cfg.force_lsm_only:
        return
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
    progress('Filling holes in grid')
    filled_grids = [0, 0]
    all_filled = []

    if cfg.topo_fill_v2.apply:
        filled_grids = fill_topo_v2(cfg, connection, cur)
        all_filled.append(filled_grids)

    while len(filled_grids) > 0:
        filled_grids = filling_grid(cfg, connection, cur)
        all_filled.append(filled_grids)

    if cfg.do_cct:
        sqltext = 'UPDATE "{0}"."{1}" ' \
                  'SET nz = 3, height = 3 * {2}' \
                  'WHERE nz < 3 OR height < 3 * {2}'\
                  .format(cfg.domain.case_schema,
                          cfg.tables.buildings_grid, cfg.domain.dz)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

    restore_log_level(cfg)

def update_force_cyclic(cfg, connection, cur):
    """ In case of force cyclic option, force terrain changes to satisfy cyclic boundary condition """
    progress('Updating nz, height in grid table in order to fulfill cyclic boundary condition')
    verbose('Update front - back cycling bc')
    j_to_floor = cfg.force_cyclic_nc
    i_to_floor = cfg.force_cyclic_nc
    j_to_modifies = [j for j in range(j_to_floor)] + [cfg.domain.ny - 1 - j for j in range(j_to_floor + 1)]
    i_to_modifies = [i for i in range(i_to_floor)] + [cfg.domain.nx - 1 - i for i in range(i_to_floor + 1)]
    for j_to_modify in j_to_modifies:
        verbose('Updating j level: {}', j_to_modify)
        sqltext = 'WITH gf AS (SELECT i, height, nz FROM "{0}"."{1}" WHERE j = {3}) ' \
                  'UPDATE "{0}"."{1}" AS g SET (nz, height) =  ' \
                  '   (gf.nz, gf.height)  ' \
                  'FROM gf  ' \
                  'WHERE gf.i = g.i AND j = {2}' \
            .format(cfg.domain.case_schema, cfg.tables.grid, j_to_modify, j_to_floor)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

    verbose('Update left - right cycling bc')
    for i_to_modify in i_to_modifies:
        verbose('Updating i level: {}', i_to_modify)
        sqltext = 'WITH gl AS (SELECT j, height, nz FROM "{0}"."{1}" WHERE i = {3}) ' \
                  'UPDATE "{0}"."{1}" AS g SET (nz, height) =  ' \
                  '   (gl.nz, gl.height)  ' \
                  'FROM gl  ' \
                  'WHERE gl.j = g.j AND i = {2}'\
                  .format(cfg.domain.case_schema, cfg.tables.grid, i_to_modify, i_to_floor)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()


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

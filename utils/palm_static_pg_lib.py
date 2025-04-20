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
import sys
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
    vtables = [cfg.tables.landcover, cfg.tables.roofs, cfg.tables.walls, cfg.tables.trees, cfg.tables.extras_shp,
               cfg.tables.centerline, cfg.tables.building_area]
    vidx = [cfg.idx.landcover, cfg.idx.roofs, cfg.idx.walls, cfg.idx.trees, cfg.idx.extras_shp,
            cfg.idx.centerline, cfg.idx.building_area]
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
            warning('Table {} does not exist in input schema', rel)
            continue

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
        if rel == 'landcover' and not cfg.multipolygon:
            sqltext = 'SELECT column_name FROM information_schema.columns WHERE table_schema = %s AND table_name = %s'
            cur.execute(sqltext, (cfg.input_schema, rel,))
            colums = cur.fetchall()
            sql_debug(connection)
            connection.commit()

            columns_list = [x[0] for x in colums]
            columns_list.remove('geom')
            try:
                columns_list.remove('gid')
            except:
                pass
            try:
                columns_list.remove('lid')
            except:
                pass

            # change column type from multipolygon to polygon
            sqltext = 'ALTER TABLE "{0}"."{1}" ALTER COLUMN geom type geometry(Polygon, {2})' \
                .format(cfg.domain.case_schema, rel, srid_rel)
            cur.execute(sqltext, (cfg.srid_palm, grid_ext,))

            sqltext = 'insert into "{0}"."{1}" ({4}, geom) select {4}, (ST_Dump(geom)).geom::geometry(Polygon,{3})' \
                      ' from "{2}"."{1}" where ST_Intersects(ST_Transform(geom, %s), %s)' \
                .format(cfg.domain.case_schema, rel, cfg.input_schema, srid_rel, ','.join(columns_list))
            cur.execute(sqltext, (cfg.srid_palm, grid_ext,))

        else:
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
            sys.exit(1)

    if cfg.tables.trees in vtabs:
        cfg._settings['has_trees'] = True
    else:
        cfg._settings['has_trees'] = False

    if cfg.slurb:
        if not cfg.tables.centerline in vtabs or not cfg.tables.building_area in vtabs:
            warning('Centerline or building_area table is missing in input schema, supress option for slurb driver')
            cfg._settings['slurb'] = False
    return vtabs

def copy_rasters_from_input(grid_ext, cfg, connection, cur):
    """ Copy raster tables from input schema.
        Transform and clip them to grid extend and grid coordinate system

        params: grid_ext: rectangle polygon around grid, created in calculate_grid_extend function
    """
    rtables = [cfg.tables.dem, cfg.tables.buildings_height, cfg.tables.extras,
               cfg.tables.lai, cfg.tables.canopy_height]
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
            continue

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
            sys.exit(1)
    return rtabs

def check_surface_params(cfg, connection, cur):
    """ Check if surface params are present in Input Schema """
    debug('Checking surface parameters existence in Input schema')
    cur.execute('SELECT EXISTS(SELECT * FROM information_schema.tables '
                'WHERE table_schema=%s and table_name=%s)',
                (cfg.input_schema, cfg.tables.surface_params,))
    rel_exists = cur.fetchone()[0]

    cur.execute('SELECT EXISTS('
                ' SELECT column_name '
                ' FROM information_schema.columns '
                ' WHERE table_schema=%s AND table_name=%s AND column_name=%s)',
                (cfg.input_schema, cfg.tables.landcover, cfg.landcover_params_var,))
    rel_exists_katland = cur.fetchone()[0]
    if rel_exists and rel_exists_katland:
        progress('Surface params and katland parameter in landcover detected in Input Schema, apply LOD2 routines')
        cfg._settings['has_surface_params'] = True
        sqltext = 'CREATE TABLE "{}"."{}" (LIKE "{}"."{}" INCLUDING ALL)' \
            .format(cfg.domain.case_schema, cfg.tables.surface_params, cfg.input_schema, cfg.tables.surface_params)
        cur.execute(sqltext)
        sqltext = 'INSERT INTO "{}"."{}" SELECT * FROM "{}"."{}"' \
            .format(cfg.domain.case_schema, cfg.tables.surface_params, cfg.input_schema, cfg.tables.surface_params)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()
        sqltext = 'ALTER TABLE "{}"."{}" OWNER TO {}'\
                  .format(cfg.domain.case_schema, cfg.tables.surface_params, cfg.pg_owner)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()
    else:
        debug('Surface params was not detected in Input Schema')
        cfg._settings['has_surface_params'] = False

def check_buildings(cfg, connection, cur, rtabs, vtabs, grid_ext):
    """ Check if USM are present in domain """
    verbose('Checking if buildings raster is present in inputs')
    cfg._settings['has_buildings'] = False
    if cfg.tables.buildings_height in rtabs:
        cfg._settings['has_buildings'] = True

    verbose('Checking if LOD 2 will be applied')
    if cfg.has_surface_params and \
            cfg.tables.roofs in vtabs and \
            cfg.tables.walls in vtabs:
        cfg._settings['lod2'] = True
    else:
        cfg._settings['lod2'] = False

    if cfg.cortyard_fill.apply_polygon:
        debug('Filling courtyards in polygon format ')
        fill_cortyard_polygon(cfg, connection, cur)

    debug('Modification of buildings that intersect with domain boundary or are adjacet to domain extent')
    if cfg.force_building_boundary:
        if cfg.lod2:
            # modify Katland also in case of lod2
            sqltext = f"""
                UPDATE "{cfg.domain.case_schema}"."{cfg.tables.landcover}" AS l 
                SET type = 202, 
                    katland = 32,
                    albedo = 0.1,
                    emisivita = 0.93
                WHERE type >= {cfg.type_range.building_min} 
                    AND ST_Distance(l.geom, ST_Boundary(%s::geometry)) < {cfg.force_building_boundary_dist * cfg.domain.dx}
                """
            cur.execute(sqltext, (grid_ext,))
            sql_debug(connection)
            connection.commit()
        else:
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
    if cfg.crop_small_buildings:
        sqltext = 'UPDATE "{0}"."{1}" AS l SET type = 202 ' \
                  'WHERE type >= {2} AND ST_Area(geom) < {3}'\
                  .format(cfg.domain.case_schema, cfg.tables.landcover, cfg.type_range.building_min,
                          cfg.small_buildings_area * cfg.domain.dx ** 2 )
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

    if cfg.tables.extras_shp in vtabs and cfg.tables.extras in rtabs:
        debug('Buildings 3d are present in inputs, will be processed')
        cfg._settings['has_3d_buildings'] = True
    else:
        cfg._settings['has_3d_buildings'] = False
        verbose('Buildings 3d are not present in inputs')
        cur.execute('DROP TABLE IF EXISTS "{0}"."{1}"'.format(cfg.domain.case_schema, cfg.tables.extras_shp))
        cur.execute('DROP TABLE IF EXISTS "{0}"."{1}"'.format(cfg.domain.case_schema, cfg.tables.extras))
        sql_debug(connection)
        connection.commit()

    debug('Checking canopy lai')
    if cfg.canopy.using_lai:
        if not cfg.tables.lai in rtabs:
            warning('LAI is not in inputs, suppress canopy.using_lai option')
            cfg.canopy._settings['using_lai'] = False

        if not cfg.tables.canopy_height in rtabs:
            warning('Canopy_height is not in inputs, suppress canopy.using_lai option')
            cfg.canopy._settings['using_lai'] = False

def check_landcover(cfg):
    """ Placeholder for landcover related checks """
    if cfg.landcover.surface_fractions and cfg.lod2:
        debug('Lod2 is not yet implemented with surface fractions')
        cfg._settings['lod2'] = False

def calculate_terrain_height(cfg, connection, cur):
    """ Calculate terrain height for each grid cell in PALM domain.
        Height is calculated from raster table DEM in grid center.
    """
    if cfg.flat_terrain.force:
        debug('Forcing flat terrain')
        sqltext = 'UPDATE "{0}"."{1}" g ' \
                  'SET height = {2}'
        sqltext = sqltext.format(cfg.domain.case_schema, cfg.tables.grid, cfg.flat_terrain.height)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()
    else:
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
    if cfg.slurb:
        sqltext = 'UPDATE "{0}"."{1}" ' \
                  'SET type = 101 ' \
                  'where type < 900'
        sqltext = sqltext.format(cfg.domain.case_schema, cfg.tables.landcover)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

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

    if cfg.slurb:
        # Process building only with coverage > default
        debug('Calculate building fraction')
        sqltext = f"""
            create temp table building2drop as 
            select g.id
            from "{cfg.domain.case_schema}"."{cfg.tables.grid}" g
                join "{cfg.domain.case_schema}"."{cfg.tables.landcover}" l on l.lid = g.lid
            where l.type between {cfg.type_range.building_min} and {cfg.type_range.building_max}
            group by g.id
            having SUM(ST_Area(ST_Intersection(g.geom, l.geom))) / {cfg.domain.dx * cfg.domain.dy} < {cfg.min_plan_area}
            ;
        """
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        debug('Replace grids not sufficiently covered by buildings')
        sqltext = f"""
        update  "{cfg.domain.case_schema}"."{cfg.tables.grid}" g 
        set lid = (select l.lid from "{cfg.domain.case_schema}"."{cfg.tables.landcover}" l 
                     where not l.type between {cfg.type_range.building_min} and {cfg.type_range.building_max}
                     order by ST_Distance(l.geom, ST_SetSRID(ST_Point(g.xcen,g.ycen), %s)) limit 1)
        where id in (select id from building2drop)
        """
        cur.execute(sqltext, (cfg.srid_palm,))
        sql_debug(connection)
        connection.commit()

    # surface fractions
    if cfg.landcover.surface_fractions:
        debug('Calculation of surface fraction')
        # fractions for vegetation, pavement, water
        verbose('Add surface fraction columns to grid table')
        sqltext = ('ALTER TABLE "{0}"."{1}" '
                   'ADD COLUMN IF NOT EXISTS veg_fraction double precision DEFAULT 0.0, '
                   'ADD COLUMN IF NOT EXISTS wat_fraction double precision DEFAULT 0.0, '
                   'ADD COLUMN IF NOT EXISTS pav_fraction double precision DEFAULT 0.0, '
                   'ADD COLUMN IF NOT EXISTS veg_fract_type integer, '
                   'ADD COLUMN IF NOT EXISTS wat_fract_type integer, '
                   'ADD COLUMN IF NOT EXISTS pav_fract_type integer, '
                   'ADD COLUMN IF NOT EXISTS build_fraction boolean DEFAULT FALSE')
        cur.execute(sqltext.format(cfg.domain.case_schema, cfg.tables.grid))
        sql_debug(connection)
        connection.commit()

        verbose('Mark building grid')
        sqltext = (('UPDATE "{0}"."{1}" g SET '
                    'build_fraction = TRUE '
                    'FROM "{0}"."{2}" l '
                    'WHERE g.lid = l.lid AND l.type BETWEEN 900 AND 999 ')
                   .format(cfg.domain.case_schema, cfg.tables.grid, cfg.tables.landcover))
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        verbose('For vegetation')
        sqltext = (('UPDATE "{0}"."{1}" g SET '
                    'veg_fraction = s.sum_area '
                    'FROM ( '
                    '       SELECT g.id AS gid, SUM(ST_Area(ST_Intersection(g.geom, l.geom))) / {5} AS sum_area '
                    '       FROM "{0}"."{1}" g'
                    '       JOIN "{0}"."{2}" l ON ST_Intersects(l.geom, g.geom) '
                    '       WHERE l.type BETWEEN {3} AND {4} '
                    '       GROUP BY g.id '
                    ') AS s '
                    'WHERE g.id = s.gid '
                    '      AND NOT g.build_fraction;')
                   .format(cfg.domain.case_schema, cfg.tables.grid, cfg.tables.landcover,
                           cfg.type_range.vegetation_min, cfg.type_range.vegetation_max,
                           cfg.domain.dx * cfg.domain.dy))
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        sqltext = (('UPDATE "{0}"."{1}" g SET '
                    'veg_fract_type = ( '
                    ' SELECT l.type '
                    ' FROM "{0}"."{2}" l '
                    ' WHERE ST_Intersects(l.geom, g.geom) '
                    '       AND l.type BETWEEN {3} AND {4} '
                    ' ORDER BY ST_Area(ST_Intersection(g.geom, l.geom)) DESC '
                    ' LIMIT 1'
                    ') '
                    'WHERE NOT g.build_fraction ')
                    .format(cfg.domain.case_schema, cfg.tables.grid, cfg.tables.landcover,
                            cfg.type_range.vegetation_min, cfg.type_range.vegetation_max))
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        verbose('For water')
        sqltext = (('UPDATE "{0}"."{1}" g SET '
                    'wat_fraction = s.sum_area '
                    'FROM ( '
                    '       SELECT g.id AS gid, SUM(ST_Area(ST_Intersection(g.geom, l.geom))) / {5} AS sum_area '
                    '       FROM "{0}"."{1}" g'
                    '       JOIN "{0}"."{2}" l ON ST_Intersects(l.geom, g.geom) '
                    '       WHERE l.type BETWEEN {3} AND {4} '
                    '       GROUP BY g.id '
                    ') AS s '
                    'WHERE g.id = s.gid '
                    '      AND NOT g.build_fraction;')
                   .format(cfg.domain.case_schema, cfg.tables.grid, cfg.tables.landcover,
                           cfg.type_range.water_min, cfg.type_range.water_max,
                           cfg.domain.dx * cfg.domain.dy))
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        sqltext = (('UPDATE "{0}"."{1}" g SET '
                    'wat_fract_type = ( '
                    ' SELECT l.type '
                    ' FROM "{0}"."{2}" l '
                    ' WHERE ST_Intersects(l.geom, g.geom) '
                    '       AND l.type BETWEEN {3} AND {4} '
                    ' ORDER BY ST_Area(ST_Intersection(g.geom, l.geom)) DESC '
                    ' LIMIT 1'
                    ') '
                    'WHERE NOT g.build_fraction ')
                    .format(cfg.domain.case_schema, cfg.tables.grid, cfg.tables.landcover,
                            cfg.type_range.water_min, cfg.type_range.water_max))
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        verbose('For pavement')
        sqltext = (('UPDATE "{0}"."{1}" g SET '
                    'pav_fraction = s.sum_area '
                    'FROM ( '
                    '       SELECT g.id AS gid, SUM(ST_Area(ST_Intersection(g.geom, l.geom))) / {5} AS sum_area '
                    '       FROM "{0}"."{1}" g'
                    '       JOIN "{0}"."{2}" l ON ST_Intersects(l.geom, g.geom) '
                    '       WHERE l.type BETWEEN {3} AND {4} '
                    '       GROUP BY g.id '
                    ') AS s '
                    'WHERE g.id = s.gid '
                    '      AND NOT g.build_fraction;')
                   .format(cfg.domain.case_schema, cfg.tables.grid, cfg.tables.landcover,
                           cfg.type_range.pavement_min, cfg.type_range.pavement_max,
                           cfg.domain.dx * cfg.domain.dy))
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        sqltext = (('UPDATE "{0}"."{1}" g SET '
                    'pav_fract_type = ( '
                    ' SELECT l.type '
                    ' FROM "{0}"."{2}" l '
                    ' WHERE ST_Intersects(l.geom, g.geom) '
                    '       AND l.type BETWEEN {3} AND {4} '
                    ' ORDER BY ST_Area(ST_Intersection(g.geom, l.geom)) DESC '
                    ' LIMIT 1'
                    ') '
                    'WHERE NOT g.build_fraction ')
                   .format(cfg.domain.case_schema, cfg.tables.grid, cfg.tables.landcover,
                           cfg.type_range.pavement_min, cfg.type_range.pavement_max))
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        debug('Check minumum fraction {}', cfg.landcover.min_fraction)
        sqltext = (('UPDATE "{0}"."{1}" set '
                    'veg_fraction = 0.0 '
                    'WHERE veg_fraction <= {2}'
                    '      AND NOT build_fraction')
                   .format(cfg.domain.case_schema, cfg.tables.grid,
                           cfg.landcover.min_fraction))
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        sqltext = (('UPDATE "{0}"."{1}" set '
                    'wat_fraction = 0.0 '
                    'WHERE wat_fraction <= {2}'
                    '      AND NOT build_fraction')
                   .format(cfg.domain.case_schema, cfg.tables.grid,
                           cfg.landcover.min_fraction))
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        sqltext = (('UPDATE "{0}"."{1}" set '
                    'pav_fraction = 0.0 '
                    'WHERE pav_fraction <= {2}'
                    '      AND NOT build_fraction')
                   .format(cfg.domain.case_schema, cfg.tables.grid,
                           cfg.landcover.min_fraction))
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        debug('Checking sum of fractions to 1.0')
        sqltext = (('UPDATE "{0}"."{1}" set '
                    '(veg_fraction, wat_fraction, pav_fraction) = '
                    '(veg_fraction / (veg_fraction + wat_fraction + pav_fraction), '
                    ' wat_fraction / (veg_fraction + wat_fraction + pav_fraction), '
                    ' pav_fraction / (veg_fraction + wat_fraction + pav_fraction)) '
                    'WHERE NOT build_fraction')
                   .format(cfg.domain.case_schema, cfg.tables.grid))
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()
        # Crop minimum fractions and check fraction to 100%. Also invalidate grids where building is located.
        # TODO: include building detection in fraction and type checking

        # Last thing is to download it to netcdf


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
    if cfg.lod2:
        sqltext = 'UPDATE "{0}"."{1}" AS l SET (type, katland, albedo, emisivita) = ' \
                  ' (SELECT type, katland, albedo, emisivita FROM "{0}"."{1}" AS b ' \
                  '   WHERE type BETWEEN 900 AND 999 ' \
                  '   ORDER BY ST_Distance(l.geom, b.geom) LIMIT 1) ' \
                  'WHERE lid IN {2}'.format(cfg.domain.case_schema, cfg.tables.landcover, tuple(lids_list))
    else:
        sqltext = 'UPDATE "{0}"."{1}" AS l SET (type) = ' \
                  ' (SELECT type FROM "{0}"."{1}" AS b ' \
                  '   WHERE type BETWEEN 900 AND 999 ' \
                  '   ORDER BY ST_Distance(l.geom, b.geom) LIMIT 1) ' \
                  'WHERE lid IN {2}'.format(cfg.domain.case_schema, cfg.tables.landcover, tuple(lids_list))
    cur.execute(sqltext,)
    sql_debug(connection)
    connection.commit()

def fill_cortyard_polygon(cfg, connection, cur):
    """Function to fill Courtyards based on polygon search """
    progress('Filling polygon courtyards')
    while True:
        sqltext = f"""
            DROP TABLE IF EXISTS temp_c_lids;
            CREATE TEMP TABLE temp_c_lids as 
            select 
                l.lid as llid,
                count(*) filter(where ln.type between 900 and 999) as building_count,
                count(*) as count,
                min(ln.type) as new_type
            from "{cfg.domain.case_schema}"."{cfg.tables.landcover}" l
                join "{cfg.domain.case_schema}"."{cfg.tables.landcover}" ln ON  st_intersects(ST_Boundary(st_buffer(l.geom, 0.5)), ln.geom) and not l.lid = ln.lid
            where l.type < {cfg.type_range.building_min}
                and ST_Area(l.geom) < {cfg.cortyard_fill.polygon_area}
                --and ln.type >= {cfg.type_range.building_min}
            group by 1
            having count(*) filter(where ln.type between 900 and 999) = count(*)
        """
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        debug('Fetch number of modified polygons')
        sqltext = """
            select coalesce(count(*), 0) from temp_c_lids;
        """
        cur.execute(sqltext)
        modified = cur.fetchone()[0]

        if modified == 0:
            debug('{} polygon will be modified', modified)
            break

        debug('Replacing [{}] Courtyards to nearest building', modified)
        sqltext = f"""
            update "{cfg.domain.case_schema}"."{cfg.tables.landcover}" l
            set type = tcl.new_type
            from temp_c_lids tcl
            where tcl.llid = l.lid 
        """
        cur.execute(sqltext)
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
	            'CASE WHEN b.height IS NOT NULL AND NOT b.is_bridge THEN b.nz+bo.max_int' \
		        '     ELSE g.nz END AS nz, ' \
                'CASE WHEN b.height IS NOT NULL AND NOT b.is_bridge THEN true ' \
                '     ELSE false END AS is_building, ' \
                'CASE WHEN b.height IS NOT NULL AND NOT b.is_bridge THEN b.height+bo.max_int + {6} ' \
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
              '   WHERE nt.is_building AND b.i=nt.i AND b.j = nt.j) ' \
              ' WHERE NOT b.is_bridge'.format(
                cfg.domain.case_schema, cfg.tables.buildings_grid, cfg.domain.oro_min, cfg.tables.buildings_offset)
    cur.execute(sqltext)

    sql_debug(connection)
    connection.commit()

    debug('Dropping temp table')
    cur.execute('drop table "nz_temp"')

    filled_grids.sort()
    return filled_grids

def fill_near_boundary(cfg, connection, cur):
    """ Fill grid cell near boundary, that creates canyon structures. """
    if cfg.boundary_fill_1:
        progress('Correcting grid cell near boundary in first row around.')
        sqltext = f"""
        with potential_ji as (
            select 
                g1.id as id,   
                g2.nz as new_nz,
                g2.height as new_height
            from "{cfg.domain.case_schema}"."{cfg.tables.grid}" g1
                join "{cfg.domain.case_schema}"."{cfg.tables.grid}" g2 on g1.j=g2.j and g1.i + 1 = g2.i
            where g1.i = 0
                and g2.nz > g1.nz
        )
        update "{cfg.domain.case_schema}"."{cfg.tables.grid}" g
        set nz = pji.new_nz,
            height = pji.new_height
        from potential_ji pji 
        where pji.id = g.id;
        
        with potential_ji as (
            select 
                g1.id as id,   
                g2.nz as new_nz,
                g2.height as new_height
            from "{cfg.domain.case_schema}"."{cfg.tables.grid}" g1
                join "{cfg.domain.case_schema}"."{cfg.tables.grid}" g2 on g1.j=g2.j and g1.i - 1 = g2.i
            where g1.i = {cfg.domain.nx - 1}
                and g2.nz > g1.nz
        )
        update "{cfg.domain.case_schema}"."{cfg.tables.grid}" g
        set nz = pji.new_nz,
            height = pji.new_height
        from potential_ji pji 
        where pji.id = g.id;
        
        with potential_ji as (
            select 
                g1.id as id,   
                g2.nz as new_nz,
                g2.height as new_height
            from "{cfg.domain.case_schema}"."{cfg.tables.grid}" g1
                join "{cfg.domain.case_schema}"."{cfg.tables.grid}" g2 on g1.j + 1 = g2.j and g1.i = g2.i
            where g1.j = 0
                and g2.nz > g1.nz
        )
        update "{cfg.domain.case_schema}"."{cfg.tables.grid}" g
        set nz = pji.new_nz,
            height = pji.new_height
        from potential_ji pji 
        where pji.id = g.id;
        
        with potential_ji as (
            select 
                g1.id as id,   
                g2.nz as new_nz,
                g2.height as new_height
            from "{cfg.domain.case_schema}"."{cfg.tables.grid}" g1
                join "{cfg.domain.case_schema}"."{cfg.tables.grid}" g2 on g1.j - 1 = g2.j and g1.i = g2.i
            where g1.j = {cfg.domain.ny - 1}
                and g2.nz > g1.nz
        )
        update "{cfg.domain.case_schema}"."{cfg.tables.grid}" g
        set nz = pji.new_nz,
            height = pji.new_height
        from potential_ji pji 
        where pji.id = g.id;
        """
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

    # Extend search for grid cell above one more grid
    if cfg.boundary_fill_2:
        progress('Correcting grid cell near boundary up to second row around the boundary')
        sqltext = f"""
        with potential_ids as (
            select 
                g1.id as id1, g2.id as id2, 
                g3.nz as new_nz, g3.height as new_height
            from "{cfg.domain.case_schema}"."{cfg.tables.grid}" g1
                join "{cfg.domain.case_schema}"."{cfg.tables.grid}" g2 on g1.j = g2.j and g1.i + 1 = g2.i
                join "{cfg.domain.case_schema}"."{cfg.tables.grid}" g3 on g1.j = g3.j and g1.i + 2 = g3.i
            where g1.i = 0
                and g2.nz >= g1.nz and g3.nz > g1.nz
        )
        update "{cfg.domain.case_schema}"."{cfg.tables.grid}" g 
        set nz = pi.new_nz, height = pi.new_height
        from potential_ids pi
        where g.id = pi.id1 or g.id = pi.id2;
        
        with potential_ids as (
            select 
                g1.id as id1, g2.id as id2, 
                g3.nz as new_nz, g3.height as new_height
            from "{cfg.domain.case_schema}"."{cfg.tables.grid}" g1
                join "{cfg.domain.case_schema}"."{cfg.tables.grid}" g2 on g1.j = g2.j and g1.i - 1 = g2.i
                join "{cfg.domain.case_schema}"."{cfg.tables.grid}" g3 on g1.j = g3.j and g1.i - 2 = g3.i
            where g1.i = {cfg.domain.nx - 1}
                and g2.nz >= g1.nz and g3.nz > g1.nz
        )
        update "{cfg.domain.case_schema}"."{cfg.tables.grid}" g 
        set nz = pi.new_nz, height = pi.new_height
        from potential_ids pi
        where g.id = pi.id1 or g.id = pi.id2;
        
        with potential_ids as (
            select 
                g1.id as id1, g2.id as id2, 
                g3.nz as new_nz, g3.height as new_height
            from "{cfg.domain.case_schema}"."{cfg.tables.grid}" g1
                join "{cfg.domain.case_schema}"."{cfg.tables.grid}" g2 on g1.j + 1 = g2.j and g1.i = g2.i
                join "{cfg.domain.case_schema}"."{cfg.tables.grid}" g3 on g1.j + 2 = g3.j and g1.i = g3.i
            where g1.j = 0
                and g2.nz >= g1.nz and g3.nz > g1.nz
        )
        update "{cfg.domain.case_schema}"."{cfg.tables.grid}" g 
        set nz = pi.new_nz, height = pi.new_height
        from potential_ids pi
        where g.id = pi.id1 or g.id = pi.id2;
        
        with potential_ids as (
            select 
                g1.id as id1, g2.id as id2, 
                g3.nz as new_nz, g3.height as new_height
            from "{cfg.domain.case_schema}"."{cfg.tables.grid}" g1
                join "{cfg.domain.case_schema}"."{cfg.tables.grid}" g2 on g1.j - 1 = g2.j and g1.i = g2.i
                join "{cfg.domain.case_schema}"."{cfg.tables.grid}" g3 on g1.j - 2 = g3.j and g1.i = g3.i
            where g1.j = {cfg.domain.ny - 1}
                and g2.nz >= g1.nz and g3.nz > g1.nz
        )
        update "{cfg.domain.case_schema}"."{cfg.tables.grid}" g 
        set nz = pi.new_nz, height = pi.new_height
        from potential_ids pi
        where g.id = pi.id1 or g.id = pi.id2;
        """
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

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
              'CASE WHEN b.height IS NOT NULL AND NOT b.is_bridge THEN b.nz+bo.max_int' \
              '     ELSE g.nz END AS nz, ' \
              'CASE WHEN b.height IS NOT NULL AND NOT b.is_bridge THEN true ' \
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
              ' WHERE NOT b.is_bridge' \
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
    cur.execute('SELECT EXISTS(SELECT * FROM information_schema.tables '
                'WHERE table_schema=%s and table_name=%s)',
                (cfg.domain.case_schema, cfg.tables.buildings_height,))
    rel_exists = cur.fetchone()[0]
    if rel_exists:
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

    verbose('add 3d buildings related columns')
    sqltext = 'alter table "{0}"."{1}" ' \
              ' ADD IF NOT EXISTS has_bottom    BOOLEAN DEFAULT FALSE, ' \
              ' ADD IF NOT EXISTS height_bottom DOUBLE PRECISION, ' \
              ' ADD IF NOT EXISTS lid_extra     INTEGER, ' \
              ' ADD IF NOT EXISTS is_bridge     BOOLEAN DEFAULT FALSE, ' \
              ' ADD IF NOT EXISTS upper         BOOLEAN DEFAULT FALSE, ' \
              ' ADD IF NOT EXISTS under         BOOLEAN DEFAULT FALSE'.format(
        cfg.domain.case_schema, cfg.tables.buildings_grid)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    if cfg.has_3d_buildings:
        progress('Process 3d buildings')
        debug('Join tables with extras_shp')
        sqltext = 'UPDATE "{0}"."{1}" ' \
                  ' SET has_bottom = True, lid_extra = shp.gid ' \
                  ' FROM "{0}"."{2}" AS shp' \
                  ' WHERE ST_Within( ST_SetSRID(ST_Point(xcen,ycen),%s), shp.geom) AND shp.class3d IN (%s, %s)'.format(
            cfg.domain.case_schema, cfg.tables.buildings_grid, cfg.tables.extras_shp)
        cur.execute(sqltext, (cfg.srid_palm, cfg.build_3d.overhanging, cfg.build_3d.passage))
        sql_debug(connection)
        connection.commit()

        debug('Process bottom height of buildings')
        sqltext = 'UPDATE "{0}"."{1}" AS bd ' \
                  ' SET height_bottom = ( SELECT ST_Value(rast, ST_SetSRID( ST_Point(bd.xcen,bd.ycen) , %s)) ' \
                  ' FROM "{0}"."{2}" AS b WHERE ST_intersects(b.rast,  bd.geom) LIMIT 1) ' \
                  ' WHERE bd.has_bottom' \
            .format(cfg.domain.case_schema, cfg.tables.buildings_grid, cfg.tables.extras)
        cur.execute(sqltext, (cfg.srid_palm,))
        sql_debug(connection)
        connection.commit()

        # update missing ones from the nearest neighbors
        sqltext = 'UPDATE "{0}"."{1}" AS bd ' \
                  ' SET height_bottom = ( SELECT bdd.height_bottom FROM "{0}"."{1}" AS bdd ' \
                  '   WHERE (bdd.height_bottom != 0 OR (bdd.height_bottom IS NOT NULL AND bdd.height_bottom != 0)) ' \
                  '   ORDER BY ST_Distance(ST_SetSrid(ST_Point(bdd.xcen,bdd.ycen), %s), ' \
                  '                        ST_SetSrid(ST_Point(bd.xcen, bd.ycen),  %s)) ' \
                  '   LIMIT 1) ' \
                  ' WHERE bd.has_bottom AND bd.height_bottom IS NULL' \
            .format(cfg.domain.case_schema, cfg.tables.buildings_grid)
        cur.execute(sqltext, (cfg.srid_palm, cfg.srid_palm,))
        sql_debug(connection)
        connection.commit()

        # calculate nz_bottom
        sqltext = 'UPDATE "{0}"."{1}" SET nz_min = CAST(ROUND(height_bottom/{2}+0.001) AS INTEGER)' \
                  ' WHERE has_bottom'.format(
            cfg.domain.case_schema, cfg.tables.buildings_grid, cfg.domain.dz)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        debug('3d buildings bottom done')
        progress('Process bridges')

        # remove all height from building_grid, which are bridge
        debug('Update buildings grid where bridges are placed')
        sqltext = 'UPDATE "{0}"."{1}" SET (height, nz) = (Null ,Null)' \
                  'WHERE type = 907 '.format(cfg.domain.case_schema, cfg.tables.buildings_grid)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        # connect lid with bridge shp
        debug('Connect bridge extras_shp')
        sqltext = 'UPDATE "{0}"."{1}" ' \
                  ' SET is_bridge = TRUE, lid_extra = shp.gid FROM "{0}"."{2}" AS shp' \
                  ' WHERE ST_Within( ST_SetSRID(ST_Point(xcen,ycen),%s), shp.geom) AND shp.class3d = %s' \
                  ' AND lid_extra IS NULL'.format(
            cfg.domain.case_schema, cfg.tables.buildings_grid, cfg.tables.extras_shp)
        cur.execute(sqltext, (cfg.srid_palm, cfg.build_3d.bridge))

        debug('Update index is bridge for buildings grid')
        sqltext = 'UPDATE "{0}"."{1}" AS b ' \
                  ' SET (is_bridge, lid_extra) = (SELECT TRUE, bb.lid_extra FROM "{0}"."{1}" AS bb ' \
                  ' WHERE bb.lid_extra IS NOT NULL ' \
                  ' ORDER BY ST_Distance(ST_SetSRID(ST_Point(b.xcen,  b.ycen), %s),' \
                  '                      ST_SetSRID(ST_Point(bb.xcen,bb.ycen), %s)) ' \
                  ' LIMIT 1) ' \
                  ' WHERE b.type = 907 AND lid_extra IS NULL'
        sqltext = sqltext.format(cfg.domain.case_schema, cfg.tables.buildings_grid)
        cur.execute(sqltext, (cfg.srid_palm, cfg.srid_palm,))
        sql_debug(connection)
        connection.commit()

        # fill with heights
        debug('Update height and bottom height of bridges, different approach compared to buildings')
        sqltext = 'UPDATE "{0}"."{1}" AS br ' \
                  ' SET (height, height_bottom) = ( SELECT val, val - {3} FROM ( ' \
                  ' SELECT (ST_PixelAsPoints(b.rast)).geom AS geom, (ST_PixelAsPoints(b.rast)).val AS val ' \
                  ' FROM "{0}"."{2}" AS b WHERE ST_intersects(b.rast,  br.geom) ) AS bp ' \
                  ' WHERE ST_Intersects(bp.geom, br.geom) ' \
                  ' ORDER BY ST_Distance(bp.geom, ST_SetSRID(ST_Point(br.xcen,br.ycen), %s)) ' \
                  ' LIMIT 1 ) ' \
                  ' WHERE is_bridge' \
            .format(cfg.domain.case_schema, cfg.tables.buildings_grid, cfg.tables.extras, cfg.build_3d.bridge_width)
        cur.execute(sqltext, (cfg.srid_palm,))
        sql_debug(connection)
        connection.commit()

        debug('Fill remaining heights from nearest ones or default')
        sqltext = 'UPDATE "{0}"."{1}" AS br ' \
                  ' SET (height, height_bottom) = (' \
                  '  SELECT brn.height, brn.height_bottom FROM "{0}"."{2}" AS brn ' \
                  '   WHERE (brn.height != 0 OR (brn.height IS NOT NULL AND brn.height != 0)) ' \
                  '   ORDER BY ST_Distance(ST_SetSrid(ST_Point(brn.xcen,brn.ycen), %s), ' \
                  '                        ST_SetSrid(ST_Point(br.xcen, br.ycen), %s)) ' \
                  '   LIMIT 1 ) ' \
                  ' WHERE (br.height IS NULL OR br.height = 0) AND br.is_bridge ' \
            .format(cfg.domain.case_schema, cfg.tables.buildings_grid, cfg.tables.buildings_grid)
        cur.execute(sqltext, (cfg.srid_palm, cfg.srid_palm,))
        sql_debug(connection)
        connection.commit()

        debug('Calculate nz, nz_min for bridges')
        sqltext = 'UPDATE "{0}"."{1}" SET (nz_min, nz) = ( CAST(ROUND(height_bottom/{2}+0.001) AS INTEGER), ' \
                  '                                        CAST(ROUND(height/{2}+0.001)        AS INTEGER)) ' \
                  'WHERE is_bridge'.format(
            cfg.domain.case_schema, cfg.tables.buildings_grid, cfg.domain.dz)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        verbose('Correct cases where nz_min == 1, 1 grid empty space, is filled in PALM')
        sqltext = 'UPDATE "{0}"."{1}" ' \
                  ' SET nz_min = 0 WHERE nz_min <= 1 ' \
            .format(cfg.domain.case_schema, cfg.tables.buildings_grid, cfg.domain.dz)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        verbose('Do a correction for buildings with nz < 1, under dz resolution')
        sqltext = 'UPDATE "{0}"."{1}" ' \
                  ' SET (upper, under) = (false, true) ' \
            .format(cfg.domain.case_schema, cfg.tables.buildings_grid, cfg.domain.dz)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        verbose('Do a correction for buildings with nz < 1, under dz resolution')
        sqltext = 'UPDATE "{0}"."{1}" ' \
                  ' SET (upper, under) = (true, false) ' \
                  ' WHERE nz < 1 ' \
            .format(cfg.domain.case_schema, cfg.tables.buildings_grid, cfg.domain.dz)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        verbose('Do a correction for buildings with nz < 1, under dz resolution')
        sqltext = 'UPDATE "{0}"."{1}" ' \
                  ' SET under = false ' \
                  ' WHERE nz < 1 OR nz_min IS NULL OR nz_min = 0' \
            .format(cfg.domain.case_schema, cfg.tables.buildings_grid, cfg.domain.dz)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        debug('Bridges are done')

    # TODO: include is_bridge in all filtering routines
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
              ' WHERE b.lid = bo.lid AND bo.art_elev > 0 AND NOT b.is_bridge'.format(cfg.domain.case_schema,
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

    # Fill grid cells near boundary
    fill_near_boundary(cfg, connection, cur)

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

def connect_roofs(cfg, connection, cur):
    """ Join roofs polygons with grid """
    sqltext = 'ALTER TABLE "{0}"."{1}" ADD IF NOT EXISTS rid INTEGER' \
        .format(cfg.domain.case_schema, cfg.tables.buildings_grid)
    cur.execute(sqltext)

    sqltext = 'CREATE INDEX buildings_rid_idx ON "{0}"."{1}" (rid)'.format(
        cfg.domain.case_schema, cfg.tables.buildings_grid)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    sqltext = 'UPDATE "{0}"."{1}" b SET rid = ' \
              '(SELECT gid FROM "{0}"."{2}" AS r ' \
              ' ORDER BY ST_Distance(r.geom, ST_SetSRID(ST_Point(b.xcen,b.ycen),%s)) LIMIT 1)'
    sqltext = sqltext.format(cfg.domain.case_schema, cfg.tables.buildings_grid, cfg.tables.roofs)
    cur.execute(sqltext, (cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

def connect_walls(cfg, connection, cur, vtabs):
    """ Connect walls lines with created building walls """
    change_log_level(cfg.logs.level_build_walls)
    progress('Processing building walls')
    debug('Prepare table')
    sqltext = 'CREATE TABLE "{0}"."{1}" ' \
              '(id INTEGER not null, direction INTEGER, azimuth DOUBLE PRECISION, zenith DOUBLE PRECISION, ' \
              'xs DOUBLE PRECISION, ys DOUBLE PRECISION, xcen DOUBLE PRECISION, ycen DOUBLE PRECISION, ' \
              'nz_min INTEGER, nz_min_art INTEGER DEFAULT 0,  nz_max INTEGER, isroof BOOLEAN, inner_wall BOOLEAN, ' \
              'wid INTEGER, rid INTEGER, geom geometry(linestring, %s), primary key (id, direction, inner_wall) )'
    sqltext = sqltext.format(cfg.domain.case_schema, cfg.tables.building_walls)
    cur.execute(sqltext, (cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()
    sqltext = 'ALTER TABLE "{}"."{}" OWNER TO {}' \
        .format(cfg.domain.case_schema, cfg.tables.building_walls, cfg.pg_owner)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    debug('Insert outer walls ')
    for d in range(len(cfg.walls.wall_directions)):
        wdx, wdy = cfg.walls.wall_directions[d][0], cfg.walls.wall_directions[d][1]
        wa = cfg.walls.wall_azimuth[d]
        sqltext = 'insert into "{0}"."{1}" ' \
                  '(id, direction, azimuth, zenith, nz_min, nz_min_art, nz_max, isroof, xs, ys, xcen, ycen, wid, rid, inner_wall, geom) ' \
                  'select b.id, %s, %s, %s, b.nz_min, CASE WHEN b.nz_min = 0 THEN ng.nz-g.nz ELSE 0 END, b.nz, false, ' \
                  'g.xcen+%s*%s-%s, g.ycen+%s*%s-%s, g.xcen+%s*%s, g.ycen+%s*%s, null, null, false, ' \
                  'ST_SetSRID(ST_MakeLine(' \
                  'ST_MakePoint(' \
                  'case when %s = 0 then (g.xmi+g.xma+(%s)*%s)/2 else g.xmi end, ' \
                  'case when %s = 0 then (g.ymi+g.yma+(%s)*%s)/2 else g.ymi end), ' \
                  'ST_MakePoint(' \
                  'case when %s = 0 then (g.xmi+g.xma+(%s)*%s)/2 else g.xma end, ' \
                  'case when %s = 0 then (g.ymi+g.yma+(%s)*%s)/2 else g.yma end)), %s) ' \
                  'from "{0}"."{2}" b ' \
                  'join "{0}"."{3}" g on g.id = b.id ' \
                  'join "{0}"."{3}" ng on ng.i = g.i + %s and ng.j = g.j + %s ' \
                  'where b.nz > 0 and not exists (select * from "{0}"."{2}" where i = b.i+(%s) and j = b.j+(%s))'
        sqltext = sqltext.format(cfg.domain.case_schema, cfg.tables.building_walls, cfg.tables.buildings_grid,
                                 cfg.tables.grid)
        cur.execute(sqltext, (d + 1, wa, 90, wdx, cfg.domain.dx / 2.0, cfg.domain.origin_x, wdy, cfg.domain.dy / 2.0, cfg.domain.origin_y,
                              wdx, cfg.domain.dx / 2.0, wdy, cfg.domain.dy / 2.0, wdy, wdx, cfg.domain.dx,
                              wdx, wdy, cfg.domain.dy, wdy, wdx, cfg.domain.dx, wdx, wdy,
                              cfg.domain.dy, cfg.srid_palm, wdx, wdy, wdx, wdy,))
        sql_debug(connection)
        connection.commit()

    debug('Insert inner walls')
    for d in range(len(cfg.walls.wall_directions)):
        wdx, wdy = cfg.walls.wall_directions[d][0], cfg.walls.wall_directions[d][1]
        wa = cfg.walls.wall_azimuth[d]
        sqltext = 'insert into "{0}"."{1}" (id, direction, azimuth, zenith, nz_min, nz_max, isroof, ' \
                  '                         xs, ys, xcen, ycen, wid, rid, inner_wall, geom) ' \
                  'select b.id, %s, %s, %s, CASE WHEN bn.has_bottom THEN b.nz_min ELSE 0 END, bn.nz_min, false, ' \
                  '       g.xcen+%s*%s-%s, g.ycen+%s*%s-%s, g.xcen+%s*%s, g.ycen+%s*%s, null, null, true, ' \
                  'ST_SetSRID(' \
                  'ST_MakeLine(' \
                  'ST_MakePoint(' \
                  'case when %s = 0 then (g.xmi+g.xma+(%s)*%s)/2 else g.xmi end, ' \
                  'case when %s = 0 then (g.ymi+g.yma+(%s)*%s)/2 else g.ymi end), ' \
                  'ST_MakePoint(' \
                  'case when %s = 0 then (g.xmi+g.xma+(%s)*%s)/2 else g.xma end, ' \
                  'case when %s = 0 then (g.ymi+g.yma+(%s)*%s)/2 else g.yma end)), %s) ' \
                  'from "{0}"."{2}" b ' \
                  'join "{0}"."{2}" bn on bn.i = b.i+(%s) and bn.j = b.j+(%s) ' \
                  'join "{0}"."{3}" g on g.id = b.id ' \
                  'WHERE b.nz_min < bn.nz_min OR ' \
                  '((NOT b.has_bottom OR b.has_bottom IS NULL) AND bn.has_bottom)'
        sqltext = sqltext.format(cfg.domain.case_schema, cfg.tables.building_walls,
                                 cfg.tables.buildings_grid, cfg.tables.grid)
        cur.execute(sqltext, (d + 1, wa, 90, wdx, cfg.domain.dx / 2.0, cfg.domain.origin_x, wdy, cfg.domain.dy / 2.0, cfg.domain.origin_y,
                              wdx, cfg.domain.dx / 2.0, wdy, cfg.domain.dy / 2.0,
                              wdy, wdx, cfg.domain.dx, wdx, wdy, cfg.domain.dy, wdy, wdx, cfg.domain.dx,
                              wdx, wdy, cfg.domain.dy, cfg.srid_palm, wdx, wdy,))
        sql_debug(connection)
        connection.commit()

    debug('Insert roof vertical surfaces ')
    for d in range(len(cfg.walls.wall_directions)):
        wdx, wdy = cfg.walls.wall_directions[d][0], cfg.walls.wall_directions[d][1]
        wa = cfg.walls.wall_azimuth[d]
        sqltext = 'insert into "{0}"."{1}" (id, direction, azimuth, zenith, nz_min, nz_max, isroof, ' \
                  '                         xs, ys, xcen, ycen, wid, rid, inner_wall, geom) ' \
                  'select b.id, %s, %s, %s, d.nz-bob.difference_int, b.nz, true, ' \
                  '       g.xcen+%s*%s-%s, g.ycen+%s*%s-%s, g.xcen+%s*%s, g.ycen+%s*%s, null, null, false, ' \
                  'ST_SetSRID(' \
                  'ST_MakeLine(' \
                  'ST_MakePoint(' \
                  'case when %s = 0 then (g.xmi+g.xma+(%s)*%s)/2 else g.xmi end, ' \
                  'case when %s = 0 then (g.ymi+g.yma+(%s)*%s)/2 else g.ymi end), ' \
                  'ST_MakePoint(' \
                  'case when %s = 0 then (g.xmi+g.xma+(%s)*%s)/2 else g.xma end, ' \
                  'case when %s = 0 then (g.ymi+g.yma+(%s)*%s)/2 else g.yma end)), %s) ' \
                  'from "{0}"."{2}" b ' \
                  'join "{0}"."{2}" d on d.i = b.i+(%s) and d.j = b.j+(%s) ' \
                  'join "{0}"."{4}" bob on b.lid = bob.lid ' \
                  'join "{0}"."{4}" bod on d.lid = bod.lid ' \
                  'join "{0}"."{3}" g on g.id = b.id ' \
                  'WHERE b.nz + bob.max_int > d.nz + bod.max_int ' \
                  ' AND NOT b.is_bridge'
        sqltext = sqltext.format(cfg.domain.case_schema, cfg.tables.building_walls,
                                 cfg.tables.buildings_grid, cfg.tables.grid, cfg.tables.buildings_offset)
        cur.execute(sqltext, (d + 1, wa, 90, wdx, cfg.domain.dx / 2.0, cfg.domain.origin_x, wdy, cfg.domain.dy / 2.0, cfg.domain.origin_y,
                              wdx, cfg.domain.dx / 2.0, wdy, cfg.domain.dy / 2.0,
                              wdy, wdx, cfg.domain.dx, wdx, wdy, cfg.domain.dy, wdy, wdx, cfg.domain.dx,
                              wdx, wdy, cfg.domain.dy, cfg.srid_palm, wdx, wdy,))
        sql_debug(connection)
        connection.commit()

    debug('Connect outer walls with their properties')
    sqltext = 'update "{0}"."{1}" b ' \
              'set wid = (select gid from "{0}"."{2}" w ' \
              ' order by ST_Distance(w.geom, ST_SetSRID(ST_Point(b.xcen,b.ycen),%s)) limit 1) ' \
              'where not b.isroof'
    sqltext = sqltext.format(cfg.domain.case_schema, cfg.tables.building_walls, cfg.tables.walls)
    cur.execute(sqltext, (cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

    debug('Connect roof vertical surfaces with their properties')
    sqltext = 'update "{0}"."{1}" b ' \
              'set rid = (select gid from "{0}"."{2}" r ' \
              ' order by ST_Distance(r.geom, ST_SetSRID(ST_Point(b.xcen,b.ycen),%s)) limit 1) ' \
              'where b.isroof'
    sqltext = sqltext.format(cfg.domain.case_schema, cfg.tables.building_walls, cfg.tables.roofs)
    cur.execute(sqltext, (cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

    # create and populate individual surfaces
    progress('Processing individual building surfaces')
    debug('Create an empty table')
    sqltext = 'create table if not exists "{0}"."{1}" ' \
              '(sid integer not null, direction integer, zs double precision, ' \
              ' xs double precision, ys double precision, ' \
              ' azimuth double precision, zenith double precision, ' \
              ' lons double precision, lats double precision, ' \
              ' "Es_UTM" double precision, "Ns_UTM" double precision, ' \
              ' ishorizontal boolean, isroof boolean, gid integer, rid integer, wid integer{2},' \
              ' primary key (sid, direction, zs) )' \
        .format(cfg.domain.case_schema, cfg.tables.surfaces,
                '' if not cfg.tables.extras_shp in vtabs else ', eid integer ')
    cur.execute(sqltext)

    sqltext = 'ALTER TABLE "{0}"."{1}" OWNER TO "{2}"' \
        .format(cfg.domain.case_schema, cfg.tables.surfaces, cfg.pg_owner)
    cur.execute(sqltext)

    debug('Insert horizontal surfaces, upward facing')
    direction = 0  # upward
    sqltext = 'insert into "{0}"."{1}" ' \
              '(sid, xs, ys, zs, direction, azimuth, zenith, ishorizontal, isroof, gid, rid, wid) ' \
              'select row_number() over (order by g.id)-1 as sid, ' \
              'g.xcen-%s as xs, g.ycen-%s as ys, b.nz*%s as zs, %s, b.azimuth as azimuth, b.zenith as zenith, ' \
              'true, true, b.id, b.rid, null from "{0}"."{2}" b ' \
              'left outer join "{0}"."{3}" g on g.id = b.id ' \
              'WHERE b.nz IS NOT null AND b.type != 907' \
        .format(cfg.domain.case_schema, cfg.tables.surfaces, cfg.tables.buildings_grid, cfg.tables.grid)
    cur.execute(sqltext, (cfg.domain.origin_x, cfg.domain.origin_y, cfg.domain.dz, direction,))
    connection.commit()

    if cfg.has_3d_buildings:
        # upward horizontal grid face in bridges (type = 907)
        progress('Insert horizontal bridge surfaces, upward facing')
        direction = 0  # upward
        sqltext = 'select max(sid) from "{0}"."{1}"'.format(cfg.domain.case_schema, cfg.tables.surfaces)
        cur.execute(sqltext)
        max_sid = cur.fetchone()[0]
        sqltext = 'insert into "{0}"."{1}" ' \
                  '(sid, xs, ys, zs, direction, azimuth, zenith, ishorizontal, isroof, gid, rid, wid, eid) ' \
                  'select row_number() over (order by g.id)-1 + %s as sid, ' \
                  'g.xcen-%s as xs, g.ycen-%s as ys, b.nz*%s as zs, %s, b.azimuth as azimuth, b.zenith as zenith, ' \
                  'true, true, b.id, null, null, b.lid_extra from "{0}"."{2}" b ' \
                  'left outer join "{0}"."{3}" g on g.id = b.id ' \
                  'WHERE b.nz IS NOT null AND b.type = 907' \
            .format(cfg.domain.case_schema, cfg.tables.surfaces, cfg.tables.buildings_grid, cfg.tables.grid)
        cur.execute(sqltext, (max_sid, cfg.domain.origin_x, cfg.domain.origin_y, cfg.domain.dz, direction,))
        connection.commit()

        # downward horizontal grid faces
        progress('Insert horizontal surfaces, downward facing ')
        direction = 5  # downward
        sqltext = 'select max(sid) from "{0}"."{1}"'.format(cfg.domain.case_schema, cfg.tables.surfaces)
        cur.execute(sqltext)
        max_sid = cur.fetchone()[0]
        sqltext = 'insert into "{0}"."{1}" ' \
                  '(sid, xs, ys, zs, direction, azimuth, zenith, ishorizontal, isroof, gid, rid, wid, eid) ' \
                  'select row_number() over (order by g.id)-1 + %s as sid, ' \
                  'g.xcen-%s as xs, g.ycen-%s as ys, ' \
                  'b.nz_min*%s as zs, %s, b.azimuth as azimuth, 180 as zenith, ' \
                  'true, false, b.id, null, null, b.lid_extra from "{0}"."{2}" b ' \
                  'left outer join "{0}"."{3}" g on g.id = b.id ' \
                  'where b.has_bottom'.format(
            cfg.domain.case_schema, cfg.tables.surfaces, cfg.tables.buildings_grid,
            cfg.tables.grid)
        cur.execute(sqltext, (max_sid, cfg.domain.origin_x, cfg.domain.origin_y, cfg.domain.dz, direction,))
        sql_debug(connection)
        connection.commit()

    debug('Insert vertical surfaces')
    cur.callproc('palm_vertical_surfaces',
                 [cfg.domain.case_schema, cfg.tables.surfaces,
                  cfg.tables.building_walls, cfg.domain.dz, cfg.logs.level])
    sql_debug(connection)
    connection.commit()

    # calculate lons, lats and Es_UTM and Ns_UTM
    debug('Calculate lons, lats, Es_UTM, Ns_UTM ')
    sqltext = 'update "{0}"."{1}" set ' \
              'lons = ST_X(ST_Transform(ST_SetSRID(ST_Point(xs+%s,ys+%s),%s),%s)), ' \
              'lats = ST_Y(ST_Transform(ST_SetSRID(ST_Point(xs+%s,ys+%s),%s),%s)), ' \
              '"Es_UTM" = ST_X(ST_Transform(ST_SetSRID(ST_Point(xs+%s,ys+%s),%s),%s)), ' \
              '"Ns_UTM" = ST_Y(ST_Transform(ST_SetSRID(ST_Point(xs+%s,ys+%s),%s),%s))'\
              .format(cfg.domain.case_schema, cfg.tables.surfaces)
    cur.execute(sqltext, (cfg.domain.origin_x, cfg.domain.origin_y,
                          cfg.srid_palm, cfg.srid_wgs84,
                          cfg.domain.origin_x, cfg.domain.origin_y,
                          cfg.srid_palm, cfg.srid_wgs84,
                          cfg.domain.origin_x, cfg.domain.origin_y,
                          cfg.srid_palm, cfg.srid_utm,
                          cfg.domain.origin_x, cfg.domain.origin_y,
                          cfg.srid_palm, cfg.srid_utm,))
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
    # calculate origin lat and lon
    sqltext = 'select ST_X(ST_Transform(ST_SetSRID(ST_Point(%s,%s),%s),%s)), ' \
              'ST_Y(ST_Transform(ST_SetSRID(ST_Point(%s,%s),%s),%s)) '
    cur.execute(sqltext, (cfg.domain.origin_x, cfg.domain.origin_y, cfg.srid_palm, cfg.srid_wgs84,
                          cfg.domain.origin_x, cfg.domain.origin_y, cfg.srid_palm, cfg.srid_wgs84,))
    origin_lon, origin_lat = list(cur.fetchone())
    cfg.domain._settings['origin_lon'] = origin_lon
    cfg.domain._settings['origin_lat'] = origin_lat
    debug('Domain origin lon,lat: {}, {}', origin_lon, origin_lat)
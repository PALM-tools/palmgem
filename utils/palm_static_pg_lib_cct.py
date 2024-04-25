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

import sys
from datetime import datetime
from pathlib import Path
from math import ceil
import numpy as np
from netCDF4 import Dataset
import pandas as pd
import matplotlib.pyplot  as plt
import os
from pyproj import CRS, Transformer
from config.logger import *
from scipy.linalg import lstsq, norm
from utils.palm_static_pg_lib_netcdf import *
from utils.visualization import create_slanted_vtk

def preprocess_terrain_height(cfg, connection, cur):
    """ correct terrain height """
    progress('Starting preprocessing of terrain height')
    # create grid with, i,j from grid
    sqltext = 'DROP TABLE IF EXISTS "{0}"."{1}" CASCADE ; ' \
              'CREATE TABLE "{0}"."{1}" AS ' \
              'SELECT xmi, ymi, i, j, lid, CAST(NULL AS double precision) AS height, ' \
              '       CAST(NULL AS double precision) AS height_dummy, ' \
              '       FALSE as iswall, FALSE as inside, ' \
              '       ST_SetSRID(ST_MakePoint(xmi, ymi), %s) AS geom ' \
              'FROM "{0}"."{2}" ' \
              .format(cfg.domain.case_schema, cfg.tables.height_terr_corrected, cfg.tables.grid)
    cur.execute(sqltext, (cfg.srid_palm, ))
    sql_debug(connection)
    connection.commit()

    debug('Add integer index on table')
    sqltext = 'CREATE INDEX j_i_index ON "{0}"."{1}" (i asc, j asc)'.format(cfg.domain.case_schema, cfg.tables.height_terr_corrected)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    sqltext = 'INSERT INTO "{0}"."{1}" ' \
              'SELECT xmi, yma, i, {3}+1, lid, NULL, NULL, FALSE, FALSE, ST_SetSRID(ST_MakePoint(xmi, yma), %s) ' \
              'FROM "{0}"."{2}" ' \
              'WHERE j = {3}'.format(cfg.domain.case_schema, cfg.tables.height_terr_corrected, cfg.tables.grid, cfg.domain.ny-1)
    cur.execute(sqltext, (cfg.srid_palm, ))
    sql_debug(connection)
    connection.commit()

    sqltext = 'INSERT INTO "{0}"."{1}" ' \
              'SELECT xma, ymi, i + 1, j, lid, NULL, NULL, FALSE, FALSE, ST_SetSRID(ST_MakePoint(xma, ymi), %s) ' \
              'FROM "{0}"."{2}" ' \
              'WHERE i = {3}'.format(cfg.domain.case_schema, cfg.tables.height_terr_corrected, cfg.tables.grid, cfg.domain.nx-1)
    cur.execute(sqltext, (cfg.srid_palm, ))
    sql_debug(connection)
    connection.commit()

    sqltext = 'INSERT INTO "{0}"."{1}" ' \
              'SELECT xma, yma, i + 1, j + 1, lid, NULL, NULL, FALSE, FALSE, ST_SetSRID(ST_MakePoint(xma, yma), %s) ' \
              'FROM "{0}"."{2}" ' \
              'WHERE i = {3} AND j = {4}'.format(cfg.domain.case_schema, cfg.tables.height_terr_corrected,
                                                 cfg.tables.grid, cfg.domain.nx-1, cfg.domain.ny-1)
    cur.execute(sqltext, (cfg.srid_palm, ))
    sql_debug(connection)
    connection.commit()

    sqltext = 'INSERT INTO "{0}"."{1}" ' \
              'SELECT xmi, yma, i, j + 1, lid, NULL, NULL, FALSE, FALSE, ST_SetSRID(ST_MakePoint(xmi, yma), %s) ' \
              'FROM "{0}"."{2}" AS g ' \
              'WHERE NOT EXISTS (SELECT 1 FROM "{0}"."{1}" AS tc WHERE tc.i = g.i AND tc.j = g.j + 1)'\
              .format(cfg.domain.case_schema, cfg.tables.height_terr_corrected, cfg.tables.grid)
    cur.execute(sqltext, (cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

    sqltext = 'INSERT INTO "{0}"."{1}" ' \
              'SELECT xma, ymi, i + 1, j, lid, NULL, NULL, FALSE, FALSE, ST_SetSRID(ST_MakePoint(xma, ymi), %s) ' \
              'FROM "{0}"."{2}" AS g ' \
              'WHERE NOT EXISTS (SELECT 1 FROM "{0}"."{1}" AS tc WHERE tc.i = g.i + 1 AND tc.j = g.j)'\
              .format(cfg.domain.case_schema, cfg.tables.height_terr_corrected, cfg.tables.grid)
    cur.execute(sqltext, (cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

    sqltext = 'INSERT INTO "{0}"."{1}" ' \
              'SELECT xma, yma, i + 1, j + 1, lid, NULL, NULL, FALSE, FALSE, ST_SetSRID(ST_MakePoint(xma, yma), %s) ' \
              'FROM "{0}"."{2}" AS g ' \
              'WHERE NOT EXISTS (SELECT 1 FROM "{0}"."{1}" AS tc WHERE tc.i = g.i + 1 AND tc.j = g.j + 1)'\
              .format(cfg.domain.case_schema, cfg.tables.height_terr_corrected, cfg.tables.grid)
    cur.execute(sqltext, (cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

    # fill with height
    sqltext = 'UPDATE "{0}"."{1}" AS h ' \
              'SET height = (SELECT AVG(height) - {3} FROM "{0}"."{2}" AS g ' \
              '              WHERE g.i BETWEEN h.i-1 AND h.i+1 AND ' \
              '                    g.j BETWEEN h.j-1 AND h.j+1 )' \
              .format(cfg.domain.case_schema, cfg.tables.height_terr_corrected,
                      cfg.tables.grid, cfg.domain.origin_z)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # spacial index
    debug('Creating spacial index')
    sqltext = 'CREATE INDEX terr_height_geom_index ON "{0}"."{1}" USING gist(geom) '\
              .format(cfg.domain.case_schema, cfg.tables.height_terr_corrected)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # update lid from nearest building polygon grid
    # TODO: ST_within
    if cfg.has_buildings:
        sqltext = 'UPDATE "{0}"."{1}" AS h ' \
                  'SET lid = (SELECT bg.lid FROM "{0}"."{2}" AS bg ' \
                  '           WHERE type BETWEEN 900 AND 999 AND ST_Intersects(h.geom, lb.geom) ' \
                  '           ORDER BY ST_DISTANCE(bg.geom, h.geom) LIMIT 1) ' \
                  'FROM "{0}"."{3}" AS lb  ' \
                  'WHERE ST_Intersects(h.geom, lb.geom)' \
                  .format(cfg.domain.case_schema, cfg.tables.height_terr_corrected,
                          cfg.tables.buildings_grid, cfg.tables.build_new)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        # update information that point is inside
        sqltext = 'UPDATE "{0}"."{1}" AS h ' \
                  'SET inside = TRUE ' \
                  'FROM "{0}"."{3}" AS lb  ' \
                  'WHERE ST_Intersects(h.geom, lb.geom)' \
                  .format(cfg.domain.case_schema, cfg.tables.height_terr_corrected,
                          cfg.tables.buildings_grid, cfg.tables.build_new)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        sqltext = 'WITH bo AS (SELECT * FROM "{0}"."{2}") ' \
                  'UPDATE "{0}"."{1}" AS h ' \
                  'SET height = bo.max ' \
                  'FROM bo, "{0}"."{4}" AS lb ' \
                  'WHERE ST_Intersects(h.geom, lb.geom) AND bo.lid = h.lid ' \
                  .format(cfg.domain.case_schema, cfg.tables.height_terr_corrected,
                          cfg.tables.buildings_offset, cfg.domain.origin_z, cfg.tables.build_new)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

    sqltext = 'UPDATE "{0}"."{1}" AS h ' \
              'SET height_dummy = height'  \
              .format(cfg.domain.case_schema, cfg.tables.height_terr_corrected)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # put an indexes on the table
    sqltext = 'CREATE INDEX terr_height_corrected_j_i_idx ON "{0}"."{1}" (i asc, j asc)'\
               .format(cfg.domain.case_schema, cfg.tables.height_terr_corrected)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    sqltext = 'CREATE INDEX terr_height_corrected_geom_idx ON "{0}"."{1}" USING gist(geom)'\
               .format(cfg.domain.case_schema, cfg.tables.height_terr_corrected)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # put an indexes on the table
    sqltext = 'CREATE INDEX terr_height_corrected_ji_idx ON "{0}"."{1}" (j,i)'\
               .format(cfg.domain.case_schema, cfg.tables.height_terr_corrected)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # now filter special cases
    sqltext = 'WITH hc_diag  AS (SELECT i,j,height FROM "{0}"."{1}"), ' \
              '     hc_top   AS (SELECT i,j,height FROM "{0}"."{1}"), ' \
              '     hc_right AS (SELECT i,j,height FROM "{0}"."{1}") ' \
              '' \
              'SELECT h.i,h.j,h.height AS h_height, hc_diag.height AS h_diag, hc_top.height AS h_top, hc_right.height AS h_right,' \
              '       (h.height < hc_top.height AND h.height < hc_right.height AND hc_diag.height < hc_top.height AND hc_diag.height < hc_right.height),' \
              '       (h.height > hc_top.height AND h.height > hc_right.height AND hc_diag.height > hc_top.height AND hc_diag.height > hc_right.height) ' \
              'FROM "{0}"."{1}" h ' \
              'LEFT JOIN hc_diag  ON hc_diag.i  = h.i + 1 AND hc_diag.j  = h.j + 1 AND hc_diag.height  IS NOT NULL AND NOT iswall ' \
              'LEFT JOIN hc_top   ON hc_top.i   = h.i     AND hc_top.j   = h.j + 1 AND hc_top.height   IS NOT NULL AND NOT iswall ' \
              'LEFT JOIN hc_right ON hc_right.i = h.i + 1 AND hc_right.j = h.j     AND hc_right.height IS NOT NULL AND NOT iswall ' \
              'WHERE' \
              '     (h.height < hc_top.height AND h.height < hc_right.height AND hc_diag.height < hc_top.height AND hc_diag.height < hc_right.height) OR ' \
              '     (h.height > hc_top.height AND h.height > hc_right.height AND hc_diag.height > hc_top.height AND hc_diag.height > hc_right.height) ' \
              'ORDER BY h.i, h.j'.format(cfg.domain.case_schema, cfg.tables.height_terr_corrected)
    cur.execute(sqltext)
    ij2corr = cur.fetchall()
    sql_debug(connection)
    connection.commit()

    sqlcorr = 'UPDATE "{0}"."{1}" SET ' \
              'height = {2} ' \
              'WHERE (i = {3} AND j = {4}) OR (i = {5} AND j = {6}) '
    for i, j, hh, hdiag, htop, hright, sml, hgh in ij2corr:
        verbose('Height correction for [i,j]:[{},{}]', i, j)
        verbose('\th:{}, hdiag:{}, htop:{}, hright:{}, sml: {}, hgh: {}', hh, hdiag, htop, hright, sml, hgh)
        if sml:
            # h and hdiag are both lower, correct height in both of them to min(htop,hright)
            cur.execute(
                sqlcorr.format(cfg.domain.case_schema, cfg.tables.height_terr_corrected, min(htop, hright), i, j, i + 1,
                               j + 1))
            sql_debug(connection)
            connection.commit()
        else:
            # top, right are lower, correct height in both of them to min(hh, hdiag)
            cur.execute(
                sqlcorr.format(cfg.domain.case_schema, cfg.tables.height_terr_corrected, min(hh, hdiag), i + 1, j, i, j + 1))
            sql_debug(connection)
            connection.commit()

def calculate_aspect_slope(cfg, connection, cur):
    """ Create slope and aspect on building roofs"""
    progress('Creating Aspect from buildings')
    # FIXME: rather use corrected building heights
    sqltext = 'DROP TABLE IF EXISTS "{0}"."{1}";' \
              'CREATE TABLE "{0}"."{1}" AS SELECT ST_Aspect(rast) FROM "{0}"."{2}"'\
              .format(cfg.domain.case_schema, cfg.tables.aspect, cfg.tables.buildings_height)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # FIXME: rather use corrected building heights
    progress('Creating Slope from buildings')
    sqltext = 'DROP TABLE IF EXISTS "{0}"."{1}";' \
              'CREATE TABLE "{0}"."{1}" AS SELECT ST_Slope(rast) FROM "{0}"."{2}"'\
              .format(cfg.domain.case_schema, cfg.tables.slope, cfg.tables.buildings_height)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

def preprocess_building_height(cfg, connection, cur):
    """ correct building heights """
    progress('Starting preprocessing building height')

    # create grid with, i,j from grid
    sqltext = 'DROP TABLE IF EXISTS "{0}"."{1}" CASCADE ; ' \
              'CREATE TABLE "{0}"."{1}" AS ' \
              'SELECT xmi, ymi, i, j, ' \
              'CAST(NULL AS double precision) AS height, ' \
              'CAST(NULL AS double precision) AS height_bottom, ' \
              'ST_SetSRID(ST_MakePoint(xmi, ymi), %s) AS geom ' \
              'FROM "{0}"."{2}" ' \
              .format(cfg.domain.case_schema, cfg.tables.height_corrected, cfg.tables.grid)
    cur.execute(sqltext, (cfg.srid_palm, ))
    sql_debug(connection)
    connection.commit()

    sqltext = 'INSERT INTO "{0}"."{1}" ' \
              'SELECT xmi, yma, i, {3}+1, NULL, NULL, ST_SetSRID(ST_MakePoint(xmi, yma), %s) ' \
              'FROM "{0}"."{2}" ' \
              'WHERE j = {3}'.format(cfg.domain.case_schema, cfg.tables.height_corrected, cfg.tables.grid, cfg.domain.ny-1)
    cur.execute(sqltext, (cfg.srid_palm, ))
    sql_debug(connection)
    connection.commit()

    sqltext = 'INSERT INTO "{0}"."{1}" ' \
              'SELECT xma, ymi, i + 1, j, NULL, NULL, ST_SetSRID(ST_MakePoint(xma, ymi), %s) ' \
              'FROM "{0}"."{2}" ' \
              'WHERE i = {3}'.format(cfg.domain.case_schema, cfg.tables.height_corrected, cfg.tables.grid, cfg.domain.nx-1)
    cur.execute(sqltext, (cfg.srid_palm, ))
    sql_debug(connection)
    connection.commit()

    sqltext = 'INSERT INTO "{0}"."{1}" ' \
              'SELECT xma, yma, i + 1, j + 1, NULL, NULL, ST_SetSRID(ST_MakePoint(xma, yma), %s) ' \
              'FROM "{0}"."{2}" ' \
              'WHERE i = {3} AND j = {4}'.format(cfg.domain.case_schema, cfg.tables.height_corrected,
                                                 cfg.tables.grid, cfg.domain.nx-1, cfg.domain.ny-1)
    cur.execute(sqltext, (cfg.srid_palm, ))
    sql_debug(connection)
    connection.commit()


    sqltext = 'INSERT INTO "{0}"."{1}" ' \
              'SELECT xmi, yma, i, j + 1, NULL, NULL, ST_SetSRID(ST_MakePoint(xmi, yma), %s) ' \
              'FROM "{0}"."{2}" AS g ' \
              'WHERE NOT EXISTS (SELECT 1 FROM "{0}"."{1}" AS tc WHERE tc.i = g.i AND tc.j = g.j + 1)'\
              .format(cfg.domain.case_schema, cfg.tables.height_corrected, cfg.tables.grid)
    cur.execute(sqltext, (cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

    sqltext = 'INSERT INTO "{0}"."{1}" ' \
              'SELECT xma, ymi, i + 1, j, NULL, NULL, ST_SetSRID(ST_MakePoint(xma, ymi), %s) ' \
              'FROM "{0}"."{2}" AS g ' \
              'WHERE NOT EXISTS (SELECT 1 FROM "{0}"."{1}" AS tc WHERE tc.i = g.i + 1 AND tc.j = g.j)'\
              .format(cfg.domain.case_schema, cfg.tables.height_corrected, cfg.tables.grid)
    cur.execute(sqltext, (cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

    sqltext = 'INSERT INTO "{0}"."{1}" ' \
              'SELECT xma, yma, i + 1, j + 1, NULL, NULL, ST_SetSRID(ST_MakePoint(xma, yma), %s) ' \
              'FROM "{0}"."{2}" AS g ' \
              'WHERE NOT EXISTS (SELECT 1 FROM "{0}"."{1}" AS tc WHERE tc.i = g.i + 1 AND tc.j = g.j + 1)'\
              .format(cfg.domain.case_schema, cfg.tables.height_corrected, cfg.tables.grid)
    cur.execute(sqltext, (cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

    # fill with height
    debug('Filling table {} with heights', cfg.tables.height_corrected)
    sqltext = 'UPDATE "{0}"."{1}" AS bh ' \
              'SET height = (SELECT AVG(height) ' \
              '              FROM "{0}"."{2}" AS g ' \
              '              WHERE g.i BETWEEN bh.i-1 AND bh.i+1 AND ' \
              '                    g.j BETWEEN bh.j-1 AND bh.j+1 ) + th.height_dummy,' \
              '    height_bottom = th.height_dummy ' \
              'FROM "{0}"."{3}" AS th, "{0}"."{4}" AS lb ' \
              'WHERE th.i = bh.i AND th.j = bh.j AND ST_Intersects(bh.geom, lb.geom)' \
              .format(cfg.domain.case_schema, cfg.tables.height_corrected,
                      cfg.tables.buildings_grid, cfg.tables.height_terr_corrected,
                      cfg.tables.build_new)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # put an indexes on the table
    debug('Creating geom indexes on table {}', cfg.tables.height_corrected)
    sqltext = 'CREATE INDEX height_corrected_geom_idx ON "{0}"."{1}" USING gist(geom)'\
               .format(cfg.domain.case_schema, cfg.tables.height_corrected)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # put an indexes on the table
    sqltext = 'CREATE INDEX height_corrected_ji_idx ON "{0}"."{1}" (j,i)'\
               .format(cfg.domain.case_schema, cfg.tables.height_corrected)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # now filter special cases
    sqltext = 'WITH hc_diag  AS (SELECT i,j,height FROM "{0}"."{1}"), ' \
              '     hc_top   AS (SELECT i,j,height FROM "{0}"."{1}"), ' \
              '     hc_right AS (SELECT i,j,height FROM "{0}"."{1}") ' \
              '' \
              'SELECT h.i,h.j,h.height AS h_height, hc_diag.height AS h_diag, hc_top.height AS h_top, hc_right.height AS h_right,' \
              '       (h.height < hc_top.height AND h.height < hc_right.height AND hc_diag.height < hc_top.height AND hc_diag.height < hc_right.height),' \
              '       (h.height > hc_top.height AND h.height > hc_right.height AND hc_diag.height > hc_top.height AND hc_diag.height > hc_right.height) ' \
              'FROM "{0}"."{1}" h ' \
              'LEFT JOIN hc_diag  ON hc_diag.i  = h.i + 1 AND hc_diag.j  = h.j + 1 AND hc_diag.height  IS NOT NULL ' \
              'LEFT JOIN hc_top   ON hc_top.i   = h.i     AND hc_top.j   = h.j + 1 AND hc_top.height   IS NOT NULL ' \
              'LEFT JOIN hc_right ON hc_right.i = h.i + 1 AND hc_right.j = h.j     AND hc_right.height IS NOT NULL ' \
              'WHERE' \
              '     (h.height < hc_top.height AND h.height < hc_right.height AND hc_diag.height < hc_top.height AND hc_diag.height < hc_right.height) OR ' \
              '     (h.height > hc_top.height AND h.height > hc_right.height AND hc_diag.height > hc_top.height AND hc_diag.height > hc_right.height) ' \
              'ORDER BY h.i, h.j'.format(cfg.domain.case_schema, cfg.tables.height_corrected)
    cur.execute(sqltext)
    ij2corr = cur.fetchall()
    sql_debug(connection)
    connection.commit()

    sqlcorr = 'UPDATE "{0}"."{1}" SET ' \
              'height = {2} ' \
              'WHERE (i = {3} AND j = {4}) OR (i = {5} AND j = {6}) '
    for i, j, hh, hdiag, htop, hright, sml, hgh in ij2corr:
        verbose('Height correction for [i,j]:[{},{}]', i, j)
        verbose('\th:{}, hdiag:{}, htop:{}, hright:{}, sml: {}, hgh: {}', hh, hdiag, htop, hright, sml, hgh)
        if sml:
            # h and hdiag are both lower, correct height in both of them to min(htop,hright)
            cur.execute(sqlcorr.format(cfg.domain.case_schema, cfg.tables.height_corrected, min(htop, hright), i,j,i+1,j+1 ))
            sql_debug(connection)
            connection.commit()
        else:
            # top, right are lower, correct height in both of them to min(hh, hdiag)
            cur.execute(sqlcorr.format(cfg.domain.case_schema, cfg.tables.height_corrected, min(hh, hdiag), i+1,j,i,j+1 ))
            sql_debug(connection)
            connection.commit()

def preprocess_building_landcover(cfg, connection, cur):
    """ join adjacent buildings, using Convex/Concave Hull postgis function """
    progress('Starting preprocessing building landcover')

    # JOIN ON DISTANCE SHORTEN THEN max_dist
    debug('Joining tables polygons when distance is lower that set max_dist')
    sqltext = 'DROP TABLE IF EXISTS "{0}"."{1}"'.format(cfg.domain.case_schema, cfg.tables.land_build)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    sqltext = 'CREATE TABLE "{0}"."{1}" AS ' \
              'SELECT ST_Simplify((ST_Dump(ST_Union(ST_Buffer(geom, 0.0000001)))).geom, {5}) AS geom FROM "{0}"."{2}"' \
              'WHERE type BETWEEN {3} AND {4}'.format(cfg.domain.case_schema, cfg.tables.land_build, cfg.tables.landcover,
                                                      cfg.type_range.building_min, cfg.type_range.building_max,
                                                      cfg.slanted_pars.simplify_dist)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    cur.execute('ALTER TABLE "{0}"."{1}" ADD COLUMN lid SERIAL'
                .format(cfg.domain.case_schema, cfg.tables.land_build))
    sql_debug(connection)
    connection.commit()

    sqltext = 'SELECT lid FROM "{0}"."{1}"'.format(cfg.domain.case_schema, cfg.tables.land_build)
    cur.execute(sqltext)
    blids = cur.fetchall()
    sql_debug(connection)
    connection.commit()

    blids = [blid[0] for blid in blids]

    # find all adjacent buildings
    b_neig = []
    for blid in blids:
        sqltext = 'WITH bll AS (SELECT lid, geom FROM "{0}"."{1}" WHERE lid = {2}) ' \
                  'SELECT bl.lid ' \
                  'FROM "{0}"."{1}" AS bl, bll ' \
                  'WHERE ST_DWithin(bl.geom, bll.geom, {3}) AND bll.lid != bl.lid'.format(cfg.domain.case_schema, cfg.tables.land_build,
                                                                    blid, cfg.slanted_pars.max_dist)
        cur.execute(sqltext)
        b_ads = cur.fetchall()
        sql_debug(connection)
        connection.commit()
        b_ads = [b_ad[0] for b_ad in b_ads]
        # if len(b_ads) == 1:
        #     b_neig.append([blid, b_neig])
        for b_ad in b_ads:
            b_neig.append([blid, b_ad])
            verbose('new pair of close building pars is lids: [{},{}]', blid, b_ad)

    # TODO: remove duplicits {[a,b] ~ [b,a]}

    """ NOW CONVEX HULL """
    debug('Creating convex hull between close pairs of building polygons, using segmentized outer wall')
    build_conv = 'buildings_convex_hull'
    build_seg = 'buildings_wall_segments'

    sqltext = 'DROP TABLE IF EXISTS "{0}"."{1}"'.format(cfg.domain.case_schema, build_seg)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    sqltext = 'CREATE TABLE "{0}"."{1}" AS ' \
              'WITH segments AS ' \
              '    ( SELECT lid, ST_MakeLine(lag((pt).geom, 1, NULL) OVER (PARTITION BY lid ORDER BY lid, (pt).path), (pt).geom) AS geom ' \
              '      FROM ' \
              '          (SELECT lid, ST_DumpPoints(ST_Segmentize(geom, {3})) AS pt ' \
              '           FROM "{0}"."{2}") as dumps ) ' \
              'SELECT * FROM segments WHERE geom IS NOT NULL AND ST_Length(geom) < 2.0*{3};'.format(cfg.domain.case_schema, build_seg,
                                                                                              cfg.tables.land_build, cfg.slanted_pars.min_seg)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    cur.execute('ALTER TABLE "{0}"."{1}" ADD COLUMN llid SERIAL'
                .format(cfg.domain.case_schema, build_seg))
    sql_debug(connection)
    connection.commit()

    sqltext = 'CREATE INDEX segment_geom_idx ON "{0}"."{1}" USING gist(geom)'.format(
        cfg.domain.case_schema, build_seg)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    sqltext = 'DROP TABLE IF EXISTS "{0}"."{1}"'.format(cfg.domain.case_schema, build_conv)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    sqltext = 'CREATE TABLE "{0}"."{1}"' \
              '(lid1 integer, lid2 integer, ser integer, ' \
              ' geom geometry("POLYGON", %s))'.format(cfg.domain.case_schema, build_conv)
    cur.execute(sqltext, (cfg.srid_palm, ))
    sql_debug(connection)
    connection.commit()

    # TODO: Similar algorithm, but find all pair (multiples) and after that look into segments and join them to ConvexHull
    join_id = -1
    for b_n in b_neig:
        join_id += 1
        # do concave
        verbose('Convex Hull joining of pair lid: [{},{}]', b_n[0], b_n[1])
        sqltext = 'INSERT INTO "{0}"."{1}" ' \
                  'SELECT {5}, {3}, {4}, ST_ConvexHull(ST_Collect(ARRAY_AGG(points))) FROM ( ' \
                  '     WITH bl1 AS (SELECT geom, lid FROM "{0}"."{2}" WHERE lid = {3}),' \
                  '          bl2 AS (SELECT geom, lid FROM "{0}"."{2}" WHERE lid = {4})' \
                  '     SELECT (ST_DumpPoints(ST_Collect(bl1.geom, bl2.geom))).geom AS points ' \
                  '     FROM bl1, bl2 ' \
                  '     WHERE ST_DWithin(bl1.geom, bl2.geom, {6})' \
                  ') AS s'.format(cfg.domain.case_schema, build_conv, build_seg, b_n[0], b_n[1], join_id, cfg.slanted_pars.max_dist)
        cur.execute(sqltext)
        sql_debug(connection)
    connection.commit()

    # join This convex structures with buildings
    # if cfg.slanted_pars.clean_up:
    sqltext = 'DROP TABLE IF EXISTS "{0}"."{1}"'.format(cfg.domain.case_schema, cfg.tables.build_new)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    if join_id == -1:
        # there were no connected buildings
        debug('Original buildings landcover created with forced order')
        sqltext = 'CREATE TABLE "{0}"."{1}" AS ' \
                  'SELECT ST_ForceRHR((ST_Dump(ST_UNION(ST_Buffer(lb.geom, 0.001)))).geom) AS geom FROM "{0}"."{2}" AS lb'\
                  .format(cfg.domain.case_schema, cfg.tables.build_new, cfg.tables.land_build, build_conv)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()
    else:
        # there were buildings connected
        debug('Joining created convex hulls and original buildings landcover')
        sqltext = 'CREATE TABLE "{0}"."{1}" AS ' \
                  'SELECT ST_ForceRHR((ST_Dump(ST_UNION(ST_Buffer(geom, 0.001)))).geom) AS geom FROM ' \
                  '(SELECT geom FROM "{0}"."{2}" ' \
                  ' UNION ALL ' \
                  ' SELECT geom FROM "{0}"."{3}") AS l'\
                  .format(cfg.domain.case_schema, cfg.tables.build_new, cfg.tables.land_build, build_conv)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()
    # ST_MakePolygon(ST_ExteriorRing(geom))

    cur.execute('ALTER TABLE "{0}"."{1}" ADD COLUMN lid SERIAL'
                .format(cfg.domain.case_schema, cfg.tables.build_new))
    sql_debug(connection)
    connection.commit()

    progress('Generating outer wall')
    sqltext = 'DROP TABLE IF EXISTS "{0}"."{1}";' \
              'CREATE TABLE "{0}"."{1}" AS SELECT ' \
              'ST_SetSRID(ST_Boundary(ST_Union(geom)), %s) AS geom ' \
              'FROM "{0}"."{2}"'.format(cfg.domain.case_schema, cfg.tables.walls_outer,
                                        cfg.tables.build_new,)
    cur.execute(sqltext, (cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

    sqltext = 'ALTER TABLE "{0}"."{1}" ADD COLUMN {2} SERIAL' \
              ''.format(cfg.domain.case_schema, cfg.tables.walls_outer, cfg.idx.walls)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    sqltext = 'CREATE INDEX wall_outer_geom_index ON "{0}"."{1}" USING gist(geom)'.format(cfg.domain.case_schema, cfg.tables.walls_outer)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # generate outer walls
    progress('Generating building outer walls')
    verbose('Check if wall table is present (lod2 case)')
    sqltext = 'SELECT EXISTS(SELECT * FROM information_schema.tables WHERE table_schema=%s and table_name=%s)'
    cur.execute(sqltext, (cfg.domain.case_schema, cfg.tables.walls))
    sql_debug(connection)
    has_wall = cur.fetchone()[0]
    if has_wall:
        pass
        # TODO if exists join with landcover lid
    else:
        sqltext = 'DROP TABLE IF EXISTS "{0}"."{1}"; ' \
                  'CREATE TABLE "{0}"."{1}" AS  ' \
                  'WITH segments AS ( ' \
                  '     SELECT ST_MakeLine(lag((pt).geom, 1, NULL) OVER (PARTITION BY lid ORDER BY lid, (pt).path), (pt).geom) AS geom ' \
                  '     FROM (SELECT lid, ST_DumpPoints(ST_Boundary(geom)) AS pt FROM "{0}"."{2}") as dumps' \
                  '                 ) ' \
                  'SELECT ST_SetSRID(geom, %s) AS geom ' \
                  'FROM segments WHERE geom IS NOT NULL; ' \
                  'ALTER TABLE "{0}"."{1}" ADD COLUMN wid SERIAL'\
                 .format(cfg.domain.case_schema, cfg.tables.walls,
                          cfg.tables.build_new,)
        cur.execute(sqltext, (cfg.srid_palm,))
        sql_debug(connection)
        connection.commit()

        sqltext = 'CREATE INDEX wall_geom_index ON "{0}"."{1}" USING gist(geom)'\
                  .format(cfg.domain.case_schema, cfg.tables.walls)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

    # TODO: Drop tables convex buildings and segmentized walls
    debug('Droping temporal tables: {},{}', build_conv, build_seg)

    # temp tables:
    landcover_intersect = 'landcover_intersect'
    landcover_dif = 'landcover_dif'

    debug('Joining corrected tables with original landcover')
    sqltext = 'DROP TABLE IF EXISTS "{0}"."{1}";' \
              'CREATE TABLE "{0}"."{1}" AS ' \
              'SELECT ll.lid AS lid, 906 AS type, ST_Union(ST_Buffer(ll.geom, 0.00001)) AS geom ' \
              'FROM ' \
              '   (SELECT lb.lid, ' \
              '       (ST_Dump(ST_Intersection(l.geom, lb.geom))).geom AS geom ' \
              '    FROM "{0}"."{2}" AS l, "{0}"."{3}" AS lb ' \
              '    WHERE ST_Intersects(l.geom, lb.geom)) AS ll ' \
              'GROUP BY ll.lid ' \
              .format(cfg.domain.case_schema, landcover_intersect, cfg.tables.landcover, cfg.tables.build_new, )
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    debug('Deleting non polygon rows')
    sqltext = 'DELETE FROM "{0}"."{1}" ' \
              'WHERE NOT ST_GeometryType(geom) = %s'\
              .format(cfg.domain.case_schema, landcover_intersect)
    cur.execute(sqltext,("ST_Polygon",) )
    sql_debug(connection)
    connection.commit()

    # TODO: there might be problem with multipolygons - use ST_Dump with ST_difference -- will require new indexing, unique indexing
    debug('Creating difference between corrected building and original landcover')
    sqltext = 'DROP TABLE IF EXISTS "{0}"."{1}";' \
              'CREATE TABLE "{0}"."{1}" AS ' \
              'WITH lb AS (SELECT ST_Union(ST_Buffer(geom, 0.00001)) AS geom FROM "{0}"."{2}") ' \
              'SELECT l.lid, l.type, ' \
              '       (ST_Dump(ST_Difference(l.geom, lb.geom))).geom AS geom ' \
              'FROM "{0}"."{3}" AS l, lb'\
              .format(cfg.domain.case_schema, landcover_dif, landcover_intersect, cfg.tables.landcover)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    debug('Updating corrected landcover')
    sqltext = 'UPDATE "{0}"."{1}" AS l ' \
              'SET (type, lid) = ' \
              '(SELECT lo.type, lo.lid ' \
              ' FROM "{0}"."{2}" lo ' \
              ' WHERE NOT lo.type BETWEEN 900 AND 999 ' \
              ' ORDER BY ST_Distance(lo.geom, l.geom) ' \
              ' LIMIT 1) ' \
              'WHERE l.type BETWEEN 900 AND 999 '.format(cfg.domain.case_schema, landcover_dif, cfg.tables.landcover)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    sqltext = 'UPDATE "{0}"."{1}" AS l ' \
              'SET (type, lid) = ' \
              '(SELECT lo.type, lo.lid ' \
              ' FROM "{0}"."{2}" lo ' \
              ' WHERE lo.type BETWEEN 900 AND 999 ' \
              ' ORDER BY ST_Distance(lo.geom, l.geom) ' \
              ' LIMIT 1) ' \
              'WHERE NOT l.type BETWEEN 900 AND 999 '.format(cfg.domain.case_schema, landcover_intersect, cfg.tables.landcover)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    debug('Rename landcover to old_landcover and create new using dif and intersect')
    sqltext = 'DROP TABLE IF EXISTS "{0}"."{2}"; ' \
              'ALTER TABLE "{0}"."{1}" RENAME TO "{2}"'\
               .format(cfg.domain.case_schema, cfg.tables.landcover, cfg.tables.landcover+'_old')
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    sqltext = "SELECT indexname FROM pg_indexes" \
              " WHERE schemaname = '{0}' and tablename = '{1}'  and indexname not like '%geom%' " \
        .format(cfg.domain.case_schema, cfg.tables.landcover+'_old')
    cur.execute(sqltext)
    sql_debug(connection)
    prev_ui = cur.fetchone()[0]

    sqltext = 'ALTER TABLE "{0}"."{1}" DROP CONSTRAINT {2}' \
        .format(cfg.domain.case_schema, cfg.tables.landcover+'_old', prev_ui)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    sqltext = 'DROP TABLE IF EXISTS "{0}"."{1}";' \
              'CREATE TABLE "{0}"."{1}" AS ' \
              'SELECT * FROM "{0}"."{2}" ' \
              'UNION ALL ' \
              'SELECT * FROM "{0}"."{3}" '\
              .format(cfg.domain.case_schema, cfg.tables.landcover, landcover_dif, landcover_intersect)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    sqltext = 'SELECT COUNT(*) FROM (SELECT COUNT(*) FROM "{0}"."{1}" ' \
              ' GROUP BY {2} HAVING COUNT(*) > 1) AS a'\
              .format(cfg.domain.case_schema, cfg.tables.landcover,  cfg.idx.landcover)
    cur.execute(sqltext)
    sql_debug(connection)
    idx_uniques = cur.fetchone()[0]
    connection.commit()
    if idx_uniques > 0:
        warning('Index {} in table {} is not unique, moving original index to old column and creating new one',
                cfg.idx.landcover, cfg.tables.landcover )

        sqltext = 'ALTER TABLE "{0}"."{1}" RENAME COLUMN {2} TO {3}'\
                  .format(cfg.domain.case_schema, cfg.tables.landcover,  cfg.idx.landcover,
                          cfg.idx.landcover + '_old')
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        sqltext = 'ALTER TABLE "{0}"."{1}" ADD COLUMN {2} SERIAL' \
                  ''.format(cfg.domain.case_schema, cfg.tables.landcover,  cfg.idx.landcover)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        sqltext = 'ALTER TABLE "{0}"."{1}" ADD PRIMARY KEY ({2})'\
                  .format(cfg.domain.case_schema, cfg.tables.landcover, cfg.idx.landcover)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

    # generate roof
    progress('Generating building roofs')
    verbose('Check if roofs table is present (lod2 case)')
    sqltext = 'SELECT EXISTS(SELECT * FROM information_schema.tables WHERE table_schema=%s and table_name=%s)'
    cur.execute(sqltext, (cfg.domain.case_schema, cfg.tables.roofs))
    sql_debug(connection)
    has_roof = cur.fetchone()[0]
    if has_roof:
        debug('update lid index on roofs table to connect type')
        sqltext = 'ALTER TABLE "{0}"."{1}" DROP COLUMN IF EXISTS lid INTEGER' \
            .format(cfg.domain.case_schema, cfg.tables.roofs)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        sqltext = 'ALTER TABLE "{0}"."{1}" ADD COLUMN lid INTEGER' \
            .format(cfg.domain.case_schema, cfg.tables.roofs)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        sqltext = 'UPDATE "{0}"."{1}" AS r SET lid = ' \
                  '(SELECT lid FROM "{0}"."{2}" AS l ' \
                  ' WHERE type BETWEEN {3} AND {4} AND ST_Intersect(l.geom, r.geom) ' \
                  ' ORDER BY ST_Distance(r.geom, l.geom) ' \
                  ' LIMIT 1)' \
            .format(cfg.domain.case_schema, cfg.tables.roofs, cfg.tables.landcover,
                    cfg.type_range.building_min, cfg.type_range.building_max)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()
    else:
        sqltext = 'DROP TABLE IF EXISTS "{0}"."{1}";' \
                  'CREATE TABLE "{0}"."{1}" AS SELECT ' \
                  'lid AS lid, geom AS geom ' \
                  'FROM "{0}"."{2}" ' \
                  'WHERE type BETWEEN {3} AND {4}'\
                  .format(cfg.domain.case_schema, cfg.tables.roofs,cfg.tables.landcover,
                          cfg.type_range.building_min, cfg.type_range.building_max)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        sqltext = 'ALTER TABLE "{0}"."{1}" ADD COLUMN {2} SERIAL' \
                  ''.format(cfg.domain.case_schema, cfg.tables.roofs, cfg.idx.roofs)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        sqltext = 'ALTER TABLE "{0}"."{1}" ADD PRIMARY KEY ({2})'\
                  .format(cfg.domain.case_schema, cfg.tables.roofs, cfg.idx.roofs)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        sqltext = 'CREATE INDEX roof_geom_index ON "{0}"."{1}" USING gist(geom)'\
                  .format(cfg.domain.case_schema, cfg.tables.roofs)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

    #
    debug('Create lid index on walls table to connect type')
    sqltext = 'ALTER TABLE "{0}"."{1}" DROP COLUMN IF EXISTS lid' \
        .format(cfg.domain.case_schema, cfg.tables.walls)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    sqltext = 'ALTER TABLE "{0}"."{1}" ADD COLUMN lid INTEGER' \
        .format(cfg.domain.case_schema, cfg.tables.walls)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    sqltext = 'UPDATE "{0}"."{1}" AS w SET lid = ' \
              '(SELECT lid FROM "{0}"."{2}" AS l ' \
              ' WHERE type BETWEEN {3} AND {4} AND ST_DWithin(l.geom, w.geom, 30.0) ' \
              ' ORDER BY ST_Distance(w.geom, l.geom) ' \
              ' LIMIT 1)' \
        .format(cfg.domain.case_schema, cfg.tables.walls, cfg.tables.landcover,
                cfg.type_range.building_min, cfg.type_range.building_max)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()


    # TODO: drop unnecessary tables

def create_slanted_walls_terrain(cfg, connection, cur):
    """ So far dummy function, just idea """
    # SELECT intersecting point from wall
    progress('Processing connection between wall and adjacent terrain')
    sqltext = 'DROP TABLE IF EXISTS "{0}"."{1}"; ' \
              'CREATE TABLE "{0}"."{1}" AS ' \
              'SELECT ROW_NUMBER() OVER () AS id, point, FALSE AS j_line, FALSE AS i_line, ' \
              '       CAST(NULL AS double precision) AS height, ' \
              '       CAST(NULL AS integer) AS i,' \
              '       CAST(NULL AS integer) AS j ' \
              'FROM (SELECT point1 AS point FROM "{0}"."{2}" ' \
              '      UNION ALL ' \
              '      SELECT point2 AS point FROM "{0}"."{2}") AS s ' \
              'GROUP BY point' \
              .format(cfg.domain.case_schema, cfg.tables.slanted_terr_wall, cfg.tables.slanted_wall)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    debug('Finding if point is at i-line or j-line')
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'j_line = CASE WHEN ABS((ST_X(point) - {4}) / {2} - ROUND((ST_X(point) - {4}) / {2})) < 1e-10 THEN TRUE ELSE FALSE END, ' \
              'i_line = CASE WHEN ABS((ST_Y(point) - {5}) / {3} - ROUND((ST_Y(point) - {5}) / {3})) < 1e-10 THEN TRUE ELSE FALSE END '\
              .format(cfg.domain.case_schema, cfg.tables.slanted_terr_wall, cfg.domain.dx, cfg.domain.dy,
                      cfg.domain.origin_x, cfg.domain.origin_y)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    debug('Calculating of i,j')
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'i = FLOOR((ST_X(point) - {4}) / {2}), ' \
              'j = FLOOR((ST_Y(point) - {5}) / {3})'\
              .format(cfg.domain.case_schema, cfg.tables.slanted_terr_wall, cfg.domain.dx, cfg.domain.dy,
                      cfg.domain.origin_x, cfg.domain.origin_y)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    sqltext = 'UPDATE "{0}"."{1}" tw SET ' \
              'height = ' \
              'CASE WHEN i_line THEN (FLOOR((SELECT MAX(height) FROM "{0}"."{2}" AS th ' \
              '                             WHERE (th.i = tw.i OR th.i = tw.i+1) AND th.j = tw.j) / {3}) + 1)*{3} ' \
              '     ELSE             (FLOOR((SELECT MAX(height) FROM "{0}"."{2}" AS th ' \
              '                             WHERE (th.j = tw.j OR th.j = tw.j+1) AND th.i = tw.i) / {3}) + 1)*{3} ' \
              '     END'.format(cfg.domain.case_schema, cfg.tables.slanted_terr_wall,
                                cfg.tables.height_terr_corrected, cfg.domain.dz)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # add those point into terrain correction
    sqltext = 'INSERT INTO "{0}"."{1}" SELECT ' \
              'ST_X(point), ST_Y(point), i, j, NULL, height, height, TRUE, FALSE, point ' \
              'FROM "{0}"."{2}"'.format(cfg.domain.case_schema, cfg.tables.height_terr_corrected,
                                        cfg.tables.slanted_terr_wall)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # Do some filtering in height_terr_corrected table
    debug('Small filtering inside {} table', cfg.tables.height_terr_corrected)
    sqltext_1 = 'WITH tr AS (SELECT i,j,height,iswall FROM "{0}"."{1}" WHERE NOT iswall) ' \
              'SELECT tc.i, tc.j, tr.height ' \
              'FROM "{0}"."{1}" AS tc, tr ' \
              'WHERE ((tc.i + 1 = tr.i AND tc.j     = tr.j)     OR  ' \
              '       (tc.i     = tr.i AND tc.j + 1 = tr.j)     OR ' \
              '       (tc.i     = tr.i AND tc.j     = tr.j)   ) AND ' \
              '	      tc.height < tr.height AND tc.iswall ' \
              'ORDER BY tc.i, tc.j'.format(cfg.domain.case_schema, cfg.tables.height_terr_corrected)

    sqltext_2 = 'WITH hc_diag  AS (SELECT i,j,height FROM "{0}"."{1}" WHERE NOT iswall), ' \
                '     hc_top   AS (SELECT i,j,height FROM "{0}"."{1}" WHERE NOT iswall), ' \
                '     hc_right AS (SELECT i,j,height FROM "{0}"."{1}" WHERE NOT iswall) ' \
                '' \
                'SELECT h.i,h.j,h.height AS h_height, hc_diag.height AS h_diag, hc_top.height AS h_top, hc_right.height AS h_right,' \
                '       (h.height < hc_top.height AND h.height < hc_right.height AND hc_diag.height < hc_top.height AND hc_diag.height < hc_right.height),' \
                '       (h.height > hc_top.height AND h.height > hc_right.height AND hc_diag.height > hc_top.height AND hc_diag.height > hc_right.height) ' \
                'FROM "{0}"."{1}" h ' \
                'LEFT JOIN hc_diag  ON hc_diag.i  = h.i + 1 AND hc_diag.j  = h.j + 1 AND hc_diag.height  IS NOT NULL AND NOT iswall ' \
                'LEFT JOIN hc_top   ON hc_top.i   = h.i     AND hc_top.j   = h.j + 1 AND hc_top.height   IS NOT NULL AND NOT iswall ' \
                'LEFT JOIN hc_right ON hc_right.i = h.i + 1 AND hc_right.j = h.j     AND hc_right.height IS NOT NULL AND NOT iswall ' \
                'WHERE' \
                '     (h.height < hc_top.height AND h.height < hc_right.height AND hc_diag.height < hc_top.height AND hc_diag.height < hc_right.height) OR ' \
                '     (h.height > hc_top.height AND h.height > hc_right.height AND hc_diag.height > hc_top.height AND hc_diag.height > hc_right.height) AND NOT h.iswall ' \
                'ORDER BY h.i, h.j'.format(cfg.domain.case_schema, cfg.tables.height_terr_corrected)


    sqlcorr_1 = 'UPDATE "{0}"."{1}" SET ' \
                'height = {2} ' \
                'WHERE i = {3} AND j = {4} AND iswall'

    sqlcorr_2 = 'UPDATE "{0}"."{1}" SET ' \
                'height = {2} ' \
                'WHERE (i = {3} AND j = {4}) OR (i = {5} AND j = {6}) AND NOT iswall'

    sqltext_inside_update = \
              'UPDATE "{0}"."{1}" AS tc SET ' \
              'height = (SELECT MAX(height)-0.1 FROM "{0}"."{1}" AS tr ' \
              '          WHERE NOT inside AND ST_DWithin(tc.geom, tr.geom, {2}) LIMIT 1    ' \
              '	         ) ' \
              'WHERE inside'\
              .format(cfg.domain.case_schema, cfg.tables.height_terr_corrected, cfg.slanted_pars.dist2edge_filter)


    ij2corr = [1]
    while len(ij2corr) > 0:
        cur.execute(sqltext_inside_update)
        sql_debug(connection)
        connection.commit()

        cur.execute(sqltext_1)
        ij2corr = cur.fetchall()
        sql_debug(connection)
        connection.commit()

        debug('Updating places')
        for i, j, height in ij2corr:
            height = (np.floor(height / cfg.domain.dz) + 1) * cfg.domain.dz
            verbose('Height correction for [i,j]:[{},{}], height = {} ', i, j, height)
            cur.execute(
                sqlcorr_1.format(cfg.domain.case_schema, cfg.tables.height_terr_corrected, height, i, j))
            sql_debug(connection)
            connection.commit()

        # now filter special cases

        cur.execute(sqltext_2)
        ij2corr = cur.fetchall()
        sql_debug(connection)
        connection.commit()

        for i, j, hh, hdiag, htop, hright, sml, hgh in ij2corr:
            verbose('Height correction for [i,j]:[{},{}]', i, j)
            verbose('\th:{}, hdiag:{}, htop:{}, hright:{}, sml: {}, hgh: {}', hh, hdiag, htop, hright, sml, hgh)
            if sml:
                # h and hdiag are both lower, correct height in both of them to min(htop,hright)
                cur.execute(
                    sqlcorr_2.format(cfg.domain.case_schema, cfg.tables.height_terr_corrected, min(htop, hright), i, j,
                                   i + 1,
                                   j + 1))
                sql_debug(connection)
                connection.commit()
            else:
                # top, right are lower, correct height in both of them to min(hh, hdiag)
                cur.execute(
                    sqlcorr_2.format(cfg.domain.case_schema, cfg.tables.height_terr_corrected, min(hh, hdiag), i + 1, j,
                                   i, j + 1))
                sql_debug(connection)
                connection.commit()

    debug('Correcting the ones near the wall (just in case)')
    dist2edge = cfg.slanted_pars.dist2edge
    sqltext = 'UPDATE "{0}"."{1}" AS hc SET height = ' \
              '(SELECT height FROM "{0}"."{1}" AS wo ' \
              ' WHERE iswall AND ST_DWithin(wo.geom, hc.geom, 2.0 * {3}) ' \
              ' ORDER BY ST_Distance(wo.geom, hc.geom) ' \
              ' LIMIT 1)' \
              'FROM "{0}"."{2}" AS wo ' \
              'WHERE NOT iswall AND ST_Distance(hc.geom, wo.geom) < {3}'.format(cfg.domain.case_schema, cfg.tables.height_terr_corrected, cfg.tables.walls_outer, dist2edge)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # now filter special cases
    sqltext = 'WITH hc_diag  AS (SELECT i,j,height FROM "{0}"."{1}" WHERE NOT iswall), ' \
              '     hc_top   AS (SELECT i,j,height FROM "{0}"."{1}" WHERE NOT iswall), ' \
              '     hc_right AS (SELECT i,j,height FROM "{0}"."{1}" WHERE NOT iswall) ' \
              '' \
              'SELECT h.i,h.j,h.height AS h_height, hc_diag.height AS h_diag, hc_top.height AS h_top, hc_right.height AS h_right,' \
              '       (h.height < hc_top.height AND h.height < hc_right.height AND hc_diag.height < hc_top.height AND hc_diag.height < hc_right.height),' \
              '       (h.height > hc_top.height AND h.height > hc_right.height AND hc_diag.height > hc_top.height AND hc_diag.height > hc_right.height) ' \
              'FROM "{0}"."{1}" h ' \
              'LEFT JOIN hc_diag  ON hc_diag.i  = h.i + 1 AND hc_diag.j  = h.j + 1 AND hc_diag.height  IS NOT NULL AND NOT iswall ' \
              'LEFT JOIN hc_top   ON hc_top.i   = h.i     AND hc_top.j   = h.j + 1 AND hc_top.height   IS NOT NULL AND NOT iswall ' \
              'LEFT JOIN hc_right ON hc_right.i = h.i + 1 AND hc_right.j = h.j     AND hc_right.height IS NOT NULL AND NOT iswall ' \
              'WHERE' \
              '     (h.height < hc_top.height AND h.height < hc_right.height AND hc_diag.height < hc_top.height AND hc_diag.height < hc_right.height) OR ' \
              '     (h.height > hc_top.height AND h.height > hc_right.height AND hc_diag.height > hc_top.height AND hc_diag.height > hc_right.height) ' \
              'ORDER BY h.i, h.j'.format(cfg.domain.case_schema, cfg.tables.height_terr_corrected)
    ij2corr = [1]
    while len(ij2corr) > 0:
        cur.execute(sqltext)
        ij2corr = cur.fetchall()
        sql_debug(connection)
        connection.commit()

        sqlcorr = 'UPDATE "{0}"."{1}" SET ' \
                  'height = {2} ' \
                  'WHERE (i = {3} AND j = {4}) OR (i = {5} AND j = {6}) AND NOT iswall '
        for i, j, hh, hdiag, htop, hright, sml, hgh in ij2corr:
            verbose('Height correction for [i,j]:[{},{}]', i, j)
            verbose('\th:{}, hdiag:{}, htop:{}, hright:{}, sml: {}, hgh: {}', hh, hdiag, htop, hright, sml, hgh)
            if sml:
                # h and hdiag are both lower, correct height in both of them to min(htop,hright)
                cur.execute(
                    sqlcorr.format(cfg.domain.case_schema, cfg.tables.height_terr_corrected, min(htop, hright), i, j, i + 1,
                                   j + 1))
                sql_debug(connection)
                connection.commit()
            else:
                # top, right are lower, correct height in both of them to min(hh, hdiag)
                cur.execute(
                    sqlcorr.format(cfg.domain.case_schema, cfg.tables.height_terr_corrected, min(hh, hdiag), i + 1, j, i, j + 1))
                sql_debug(connection)
                connection.commit()

    #### PROCESS CONNECTION BETWEEN WALL AND BUILDING WALL
    progress('Processing connection between wall and adjacent terrain')
    sqltext = 'DROP TABLE IF EXISTS "{0}"."{1}"; ' \
              'CREATE TABLE IF NOT EXISTS "{0}"."{1}" AS ' \
              'SELECT id, i, j, point1 AS w_point1, point2 as w_point2, ' \
              '       split_terr_wall, norm ' \
              'FROM "{0}"."{2}" ' \
        .format(cfg.domain.case_schema, cfg.tables.slanted_terr_wall, cfg.tables.slanted_wall)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # add indexes and unique indexes
    sqltext = 'CREATE INDEX terr_wall_id_idx ON "{0}"."{1}" (id); ' \
              'CREATE INDEX terr_wall_ji_idx ON "{0}"."{1}" (i asc, j asc) ' \
        .format(cfg.domain.case_schema, cfg.tables.slanted_terr_wall)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # Add extra columns to table
    debug('Adding extra column in terr_wall table')
    sqltext = 'ALTER TABLE "{0}"."{1}" ' \
              'ADD COLUMN IF NOT EXISTS cent   geometry("POINT", %s), ' \
              'ADD COLUMN IF NOT EXISTS norm   geometry("POINT", %s), ' \
              'ADD COLUMN IF NOT EXISTS point1 geometry("POINT", %s), ' \
              'ADD COLUMN IF NOT EXISTS point2 geometry("POINT", %s), ' \
              'ADD COLUMN IF NOT EXISTS point3 geometry("POINT", %s), ' \
              'ADD COLUMN IF NOT EXISTS point4 geometry("POINT", %s), ' \
              'ADD COLUMN IF NOT EXISTS point5 geometry("POINT", %s), ' \
              'ADD COLUMN IF NOT EXISTS point6 geometry("POINT", %s), ' \
              'ADD COLUMN IF NOT EXISTS z1 double precision, ' \
              'ADD COLUMN IF NOT EXISTS z2 double precision, ' \
              'ADD COLUMN IF NOT EXISTS z3 double precision, ' \
              'ADD COLUMN IF NOT EXISTS z4 double precision, ' \
              'ADD COLUMN IF NOT EXISTS z5 double precision, ' \
              'ADD COLUMN IF NOT EXISTS z6 double precision, ' \
              'ADD COLUMN IF NOT EXISTS is_wall1 boolean, ' \
              'ADD COLUMN IF NOT EXISTS is_wall2 boolean, ' \
              'ADD COLUMN IF NOT EXISTS is_wall3 boolean, ' \
              'ADD COLUMN IF NOT EXISTS is_wall4 boolean, ' \
              'ADD COLUMN IF NOT EXISTS is_wall5 boolean, ' \
              'ADD COLUMN IF NOT EXISTS is_wall6 boolean, ' \
              'ADD COLUMN IF NOT EXISTS n_edges integer,' \
              'ADD COLUMN IF NOT EXISTS geom3d geometry("POLYGONZ", %s) '. \
        format(cfg.domain.case_schema, cfg.tables.slanted_terr_wall)
    cur.execute(sqltext, (
        cfg.srid_palm, cfg.srid_palm, cfg.srid_palm, cfg.srid_palm, cfg.srid_palm, cfg.srid_palm, cfg.srid_palm,
        cfg.srid_palm,cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

    # Assing points to polygon corners
    debug('Assigning building corners')
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'n_edges = ST_NPoints(split_terr_wall)-1, ' \
              'point1 = CASE WHEN ST_NPoints(split_terr_wall) > 1 THEN ST_PointN(ST_Boundary(split_terr_wall),1) ELSE NULL END, ' \
              'point2 = CASE WHEN ST_NPoints(split_terr_wall) > 2 THEN ST_PointN(ST_Boundary(split_terr_wall),2) ELSE NULL END, ' \
              'point3 = CASE WHEN ST_NPoints(split_terr_wall) > 3 THEN ST_PointN(ST_Boundary(split_terr_wall),3) ELSE NULL END, ' \
              'point4 = CASE WHEN ST_NPoints(split_terr_wall) > 4 THEN ST_PointN(ST_Boundary(split_terr_wall),4) ELSE NULL END, ' \
              'point5 = CASE WHEN ST_NPoints(split_terr_wall) > 5 THEN ST_PointN(ST_Boundary(split_terr_wall),5) ELSE NULL END'. \
        format(cfg.domain.case_schema, cfg.tables.slanted_terr_wall)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # decide which are wall points
    debug('Assigning which point is wall point')
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'is_wall1 = w_point1 = point1 OR w_point2 = point1, ' \
              'is_wall2 = w_point1 = point2 OR w_point2 = point2, ' \
              'is_wall3 = w_point1 = point3 OR w_point2 = point3, ' \
              'is_wall4 = w_point1 = point4 OR w_point2 = point4, ' \
              'is_wall5 = w_point1 = point5 OR w_point2 = point5'. \
        format(cfg.domain.case_schema, cfg.tables.slanted_terr_wall)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # assign height to each point, either from terrain height or from building (combination of terrain and building)
    # FIXME: optimeze using hash indexing and geom indexing in correct height table
    debug('Assigning height to terr_wall points')
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'z1 = (SELECT height FROM "{0}"."{2}" WHERE ST_DWithin(point1, geom, {3})' \
              '      ORDER BY ST_Distance(point1, geom) LIMIT 1),  ' \
              'z2 = (SELECT height FROM "{0}"."{2}" WHERE ST_DWithin(point2, geom, {3})' \
              '       ORDER BY ST_Distance(point2, geom) LIMIT 1),  ' \
              'z3 = (SELECT height FROM "{0}"."{2}" WHERE ST_DWithin(point3, geom, {3})' \
              '       ORDER BY ST_Distance(point3, geom) LIMIT 1),  ' \
              'z4 = (SELECT height FROM "{0}"."{2}" WHERE ST_DWithin(point4, geom, {3})' \
              '       ORDER BY ST_Distance(point4, geom) LIMIT 1),  ' \
              'z5 = (SELECT height FROM "{0}"."{2}" WHERE ST_DWithin(point5, geom, {3})' \
              '       ORDER BY ST_Distance(point5, geom) LIMIT 1)  ' \
              '' \
        .format(cfg.domain.case_schema, cfg.tables.slanted_terr_wall, cfg.tables.height_terr_corrected, cfg.slanted_pars.dist2edge)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # create 3d polygons to check
    debug('Create 3d terr_height polygons')
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'geom3d = ST_SetSRID(ST_ConvexHull(ST_Collect(ARRAY[' \
              'ST_MakePoint(ST_X(point1), ST_Y(point1), z1), ST_MakePoint(ST_X(point2), ST_Y(point2), z2), ' \
              'ST_MakePoint(ST_X(point3), ST_Y(point3), z3), ST_MakePoint(ST_X(point4), ST_Y(point4), z4), ' \
              'ST_MakePoint(ST_X(point5), ST_Y(point5), z5)' \
              '])), %s) '.format(cfg.domain.case_schema, cfg.tables.slanted_terr_wall)
    cur.execute(sqltext, (cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

def create_slanted_terrain(cfg, connection, cur):
    """ create slanted terrain """
    progress('Creating slanted terrain')
    if cfg.has_buildings:
        create_slanted_walls_terrain(cfg, connection, cur)

    debug('Creating table {}', cfg.tables.slanted_terrain)
    # FIXME: get lid inside this table from slanted_terr_wall
    sqltext = 'DROP TABLE IF EXISTS "{0}"."{1}";' \
              'CREATE TABLE "{0}"."{1}" ' \
              '(id INTEGER NOT NULL, lid INTEGER, i INTEGER NOT NULL, j INTEGER NOT NULL, n_edges INTEGER, ' \
              ' p1 geometry("POINTZ", %s), p2 geometry("POINTZ", %s), p3 geometry("POINTZ", %s), ' \
              ' p4 geometry("POINTZ", %s), p5 geometry("POINTZ", %s)) '.format(cfg.domain.case_schema, cfg.tables.slanted_terrain)
    cur.execute(sqltext, (cfg.srid_palm, cfg.srid_palm, cfg.srid_palm, cfg.srid_palm, cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

    if cfg.has_buildings:
        debug('Insert buildings surrounding into slanted terrain table')
        sqltext = 'INSERT INTO "{0}"."{1}" ' \
                  '(id, lid, i, j, n_edges, p1, p2, p3, p4, p5)' \
                  'SELECT id, CAST(NULL AS integer) AS lid, i, j, CAST(NULL AS integer) AS n_edges, ' \
                  '       ST_SetSRID(ST_MakePoint(ST_X(point1), ST_Y(point1), z1), %s) AS p1, ' \
                  '       ST_SetSRID(ST_MakePoint(ST_X(point2), ST_Y(point2), z2), %s) AS p2, ' \
                  '       ST_SetSRID(ST_MakePoint(ST_X(point3), ST_Y(point3), z3), %s) AS p3, ' \
                  '       ST_SetSRID(ST_MakePoint(ST_X(point4), ST_Y(point4), z4), %s) AS p4, ' \
                  '       ST_SetSRID(ST_MakePoint(ST_X(point5), ST_Y(point5), z5), %s) AS p5 ' \
                  'FROM "{0}"."{2}"' \
                  .format(cfg.domain.case_schema, cfg.tables.slanted_terrain, cfg.tables.slanted_terr_wall)
        cur.execute(sqltext, (cfg.srid_palm,cfg.srid_palm,cfg.srid_palm,cfg.srid_palm,cfg.srid_palm,))
        sql_debug(connection)
        connection.commit()

    debug('Calculated number of edges in slanted terrain table')
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'n_edges = CASE WHEN p1 IS NOT NULL THEN 1 ELSE 0 END + ' \
              '          CASE WHEN p2 IS NOT NULL THEN 1 ELSE 0 END + ' \
              '          CASE WHEN p3 IS NOT NULL THEN 1 ELSE 0 END + ' \
              '          CASE WHEN p4 IS NOT NULL THEN 1 ELSE 0 END + ' \
              '          CASE WHEN p5 IS NOT NULL THEN 1 ELSE 0 END '\
              .format(cfg.domain.case_schema, cfg.tables.slanted_terrain)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    if cfg.has_buildings:
        debug('Insert the rest of polygons into table {}', cfg.tables.slanted_terrain)
        sqltext = 'WITH nb AS (SELECT geom FROM "{0}"."{3}") ' \
                  'INSERT INTO "{0}"."{1}" (id, lid, i, j, n_edges) ' \
                  'SELECT id, g.lid, g.i, g.j, 4 ' \
                  'FROM "{0}"."{2}" AS g ' \
                  'WHERE (g.id NOT IN (SELECT id FROM "{0}"."{1}") AND ' \
                  '       g.id NOT IN (SELECT id FROM "{0}"."{2}" AS gg, nb ' \
                  '                   WHERE ST_Intersects(ST_SetSRID(ST_Point(gg.xcen,gg.ycen), %s), nb.geom))) ' \
                  .format(cfg.domain.case_schema, cfg.tables.slanted_terrain, cfg.tables.grid,
                          cfg.tables.build_new)
        cur.execute(sqltext, (cfg.srid_palm,))
        sql_debug(connection)
        connection.commit()
    else:
        debug('Insert the rest of polygons into table {}', cfg.tables.slanted_terrain)
        sqltext = 'INSERT INTO "{0}"."{1}" (id, lid, i, j, n_edges) ' \
                  'SELECT id, g.lid, g.i, g.j, 4 ' \
                  'FROM "{0}"."{2}" AS g ' \
                  'WHERE g.id NOT IN (SELECT id FROM "{0}"."{1}") ' \
                  .format(cfg.domain.case_schema, cfg.tables.slanted_terrain, cfg.tables.grid)
        cur.execute(sqltext, (cfg.srid_palm,))
        sql_debug(connection)
        connection.commit()

    sqltext = 'UPDATE "{0}"."{1}" AS st SET ' \
              ' p1 = ST_SetSRID(ST_MakePoint(ST_X(th1.geom), ST_Y(th1.geom), th1.height),%s)' \
              'FROM "{0}"."{3}" AS th1 ' \
              'WHERE st.i = th1.i AND st.j = th1.j ' \
              '      AND p1 IS NULL AND n_edges >= 1 AND NOT th1.iswall' \
              .format(cfg.domain.case_schema, cfg.tables.slanted_terrain, cfg.tables.grid,
                      cfg.tables.height_terr_corrected)
    cur.execute(sqltext, (cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

    sqltext = 'UPDATE "{0}"."{1}" AS st SET ' \
              ' p2 = ST_SetSRID(ST_MakePoint(ST_X(th2.geom), ST_Y(th2.geom), th2.height),%s)' \
              'FROM "{0}"."{3}" AS th2 ' \
              'WHERE st.i+1 = th2.i AND st.j = th2.j ' \
              '      AND p2 IS NULL AND n_edges >= 2  AND NOT th2.iswall' \
              .format(cfg.domain.case_schema, cfg.tables.slanted_terrain, cfg.tables.grid,
                      cfg.tables.height_terr_corrected)
    cur.execute(sqltext, (cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

    sqltext = 'UPDATE "{0}"."{1}" AS st SET ' \
              ' p3 = ST_SetSRID(ST_MakePoint(ST_X(th3.geom), ST_Y(th3.geom), th3.height),%s)' \
              'FROM "{0}"."{3}" AS th3 ' \
              'WHERE st.i+1 = th3.i AND st.j+1 = th3.j ' \
              '      AND p3 IS NULL AND n_edges >= 3  AND NOT th3.iswall' \
              .format(cfg.domain.case_schema, cfg.tables.slanted_terrain, cfg.tables.grid,
                      cfg.tables.height_terr_corrected)
    cur.execute(sqltext, (cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

    sqltext = 'UPDATE "{0}"."{1}" AS st SET ' \
              ' p4 = ST_SetSRID(ST_MakePoint(ST_X(th4.geom), ST_Y(th4.geom), th4.height),%s)' \
              'FROM "{0}"."{3}" AS th4 ' \
              'WHERE st.i = th4.i AND st.j+1 = th4.j ' \
              '      AND p4 IS NULL  AND n_edges >= 4  AND NOT th4.iswall' \
              .format(cfg.domain.case_schema, cfg.tables.slanted_terrain, cfg.tables.grid,
                      cfg.tables.height_terr_corrected)
    cur.execute(sqltext, (cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

    # sqltext = 'UPDATE "{0}"."{1}" AS st SET ' \
    #           ' p1 = ST_SetSRID(ST_MakePoint(ST_X(th1.geom), ST_Y(th1.geom), th1.height),%s)' \
    #           ' p2 = ST_SetSRID(ST_MakePoint(ST_X(th2.geom), ST_Y(th2.geom), th2.height),%s) AS p2,   ' \
    #           ' p3 = ST_SetSRID(ST_MakePoint(ST_X(th3.geom), ST_Y(th3.geom), th3.height),%s) AS p3,   ' \
    #           ' p4 = ST_SetSRID(ST_MakePoint(ST_X(th4.geom), ST_Y(th4.geom), th4.height),%s) AS p4, ' \
    #           'LEFT JOIN "{0}"."{3}" AS th1 ON st.i   = th1.i AND st.j = th1.j ' \
    #           'LEFT JOIN "{0}"."{3}" AS th2 ON st.i+1 = th2.i AND st.j = th2.j ' \
    #           'LEFT JOIN "{0}"."{3}" AS th3 ON st.i+1 = th3.i AND st.j+1 = th3.j ' \
    #           'LEFT JOIN "{0}"."{3}" AS th4 ON st.i   = th4.i AND st.j+1 = th4.j ' \
    #           'WHERE p1 IS NULL' \
    #           .format(cfg.domain.case_schema, cfg.tables.slanted_terrain, cfg.tables.grid,
    #                   cfg.tables.height_terr_corrected)
    # cur.execute(sqltext, (cfg.srid_palm,cfg.srid_palm,cfg.srid_palm,cfg.srid_palm,))
    # sql_debug(connection)
    # connection.commit()

    debug('Calculated number of edges in slanted terrain table')
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'n_edges = CASE WHEN p1 IS NOT NULL THEN 1 ELSE 0 END + ' \
              '          CASE WHEN p2 IS NOT NULL THEN 1 ELSE 0 END + ' \
              '          CASE WHEN p3 IS NOT NULL THEN 1 ELSE 0 END + ' \
              '          CASE WHEN p4 IS NOT NULL THEN 1 ELSE 0 END + ' \
              '          CASE WHEN p5 IS NOT NULL THEN 1 ELSE 0 END '\
              .format(cfg.domain.case_schema, cfg.tables.slanted_terrain)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # debug('Removing terrain where number of points is lower than 4')
    # sqltext = 'DELETE  FROM "{0}"."{1}" ' \
    #           'WHERE n_edges < 4'\
    #           .format(cfg.domain.case_schema, cfg.tables.slanted_terrain)
    # cur.execute(sqltext)
    # sql_debug(connection)
    # connection.commit()

    debug('Building 3d polygon for slanted terrain')
    sqltext = 'ALTER TABLE "{0}"."{1}" ' \
              'ADD COLUMN IF NOT EXISTS geom3d geometry("POLYGONZ", %s)'\
              .format(cfg.domain.case_schema, cfg.tables.slanted_terrain)
    cur.execute(sqltext, (cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'geom3d = ST_SetSRID(ST_ConvexHull(ST_Collect(ARRAY[' \
              '         p1, p2, p3, p4, p5])), %s)   '.format(cfg.domain.case_schema, cfg.tables.slanted_terrain)
    cur.execute(sqltext, (cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

    debug('Updating lid')
    sqltext = 'UPDATE "{0}"."{1}" AS st SET ' \
              'lid = (SELECT l.lid FROM "{0}"."{2}" AS l ' \
              '       WHERE NOT l.type BETWEEN 900 AND 999' \
              '       ORDER BY ST_Distance(geom3d, geom)  LIMIT 1' \
              '      )' \
              'WHERE lid IS NULL '\
              .format(cfg.domain.case_schema, cfg.tables.slanted_terrain, cfg.tables.landcover)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

def create_slanted_walls(cfg, connection, cur):
    """ Process slanted walls """
    progress('Creating slanted walls')
    # Create intersection between line and grid
    progress('Creating intersect between wall line and grid')
    cur.execute('DROP TABLE IF EXISTS "{0}"."{1}"'.format(cfg.domain.case_schema, cfg.tables.slanted_wall))
    sql_debug(connection)
    connection.commit()

    sqltext = 'CREATE TABLE "{0}"."{3}" AS ' \
              'WITH outer_walls AS (SELECT (ST_Dump(geom)).geom AS geom FROM "{0}"."{2}") ' \
              'SELECT id, i, j, CAST(NULL AS integer) AS wid, g.geom, w.geom as wall_geom, ST_Intersection(w.geom, ST_Boundary(g.geom)) as intersec ' \
              'FROM "{0}"."{1}" as g ' \
              'JOIN outer_walls as W ON ST_Intersects(g.geom, w.geom) ' \
        .format(cfg.domain.case_schema, cfg.tables.grid, cfg.tables.walls_outer, cfg.tables.slanted_wall)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()


    # sqltext = 'CREATE UNIQUE INDEX {1}_i_j on "{0}"."{1}" (i asc, j asc)' \
    #     .format(cfg.domain.case_schema, cfg.tables.slanted_wall)
    # cur.execute(sqltext)
    # sql_debug(connection)


    # Update places where multiple intersection occurred
    # TODO: FIXME: UGLY HACK, NO there should be no point that needs this, hopefully
    debug('Fixing multiple intersections 1')
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'intersec = CASE WHEN  ST_Distance(geom, ST_StartPoint(wall_geom)) > 0.20 THEN ' \
              '                      ST_Collect(ST_StartPoint(ST_GeometryN(ST_Intersection(geom, wall_geom),1)), ' \
              '                                 ST_EndPoint(  ST_GeometryN(ST_Intersection(geom, wall_geom),2)) ' \
              '                                 ) ' \
              '     ELSE ST_Collect(ST_StartPoint(ST_GeometryN(ST_Intersection(geom, wall_geom),2)), ' \
              '                     ST_EndPoint(  ST_GeometryN(ST_Intersection(geom, wall_geom),1))' \
              '                     ) ' \
              '     END ' \
              'WHERE ST_NumGeometries(ST_Intersection(geom, wall_geom)) = 2 AND ST_NPoints(intersec) > 2'.\
              format(cfg.domain.case_schema, cfg.tables.slanted_wall, cfg.tables.walls_outer, cfg.tables.grid)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # TODO: FIXME: UGLY HACK
    debug('Fixing multiple intersections 2')
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'intersec = ST_Collect(ST_StartPoint(ST_GeometryN(ST_Intersection(geom, wall_geom),2)), ' \
              '                      ST_EndPoint(  ST_GeometryN(ST_Intersection(geom, wall_geom),1)) ' \
              '                                    ) ' \
              'WHERE ST_NumGeometries(ST_Intersection(geom, wall_geom)) = 3 AND ST_NPoints(intersec) > 2'.format(cfg.domain.case_schema, cfg.tables.slanted_wall, cfg.tables.walls_outer, cfg.tables.grid)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # TODO: check and drop points that are inside structure, e.g. in corners etc

    # Add columns with points and splits
    debug('Adding extra columns to table {}', cfg.tables.slanted_wall)
    sqltext = 'ALTER TABLE "{0}"."{1}" ' \
              'ADD COLUMN IF NOT EXISTS {4} geometry({2}, %s), ' \
              'ADD COLUMN IF NOT EXISTS {5} geometry({2}, %s), ' \
              'ADD COLUMN IF NOT EXISTS {6} geometry({2}, %s), ' \
              'ADD COLUMN IF NOT EXISTS {7} geometry({2}, %s), ' \
              'ADD COLUMN IF NOT EXISTS {8} geometry({3}, %s), ' \
              'ADD COLUMN IF NOT EXISTS {9} geometry({3}, %s)' \
              .format(cfg.domain.case_schema, cfg.tables.slanted_wall,
                      'POLYGON', 'POINT', 'split', 'split_terr_wall', 'split1', 'split2', 'point1', 'point2')
    cur.execute(sqltext, (cfg.srid_palm,cfg.srid_palm,cfg.srid_palm,cfg.srid_palm,cfg.srid_palm,cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

    # Set this columns
    debug('Updating slanted walls (splits, points)')
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'split1 = ST_SetSRID(ST_GeometryN(ST_Split(geom,ST_LineFromMultiPoint(intersec)),1),%s), ' \
              'split2 = ST_SetSRID(ST_GeometryN(ST_Split(geom,ST_LineFromMultiPoint(intersec)),2),%s), ' \
              'point1 = ST_SetSRID(ST_GeometryN(intersec, 1), %s), ' \
              'point2 = ST_SetSRID(ST_GeometryN(intersec, 2), %s) ' \
              'WHERE ST_NPoints(intersec) > 1' \
        .format(cfg.domain.case_schema, cfg.tables.slanted_wall, str('ST_MultiPoint'))
    cur.execute(sqltext, (cfg.srid_palm,cfg.srid_palm,cfg.srid_palm,cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

    sqltext = 'CREATE INDEX point1_geom_index ON "{0}"."{1}" USING gist(point1) ;' \
              'CREATE INDEX point2_geom_index ON "{0}"."{1}" USING gist(point2)'.format(cfg.domain.case_schema, cfg.tables.slanted_wall)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # TODO: FIXME: UGLY HACK
    # remove the row that where split1 or split2 IS NULL
    sqltext = 'DELETE FROM "{0}"."{1}" ' \
              'WHERE split1 IS NULL OR split2 IS NULL'.format(cfg.domain.case_schema, cfg.tables.slanted_wall)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # Calculate Normal vectors
    debug('Calculating normal vectors')
    sqltext = 'ALTER TABLE "{0}"."{1}" ' \
              'ADD COLUMN IF NOT EXISTS cent  geometry("POINT", %s), ' \
              'ADD COLUMN IF NOT EXISTS norm1 geometry("POINT", %s), ' \
              'ADD COLUMN IF NOT EXISTS norm2 geometry("POINT", %s), ' \
              'ADD COLUMN IF NOT EXISTS norm  geometry("POINT", %s),' \
              'ADD COLUMN IF NOT EXISTS z1 double precision,' \
              'ADD COLUMN IF NOT EXISTS z1b double precision,' \
              'ADD COLUMN IF NOT EXISTS z2 double precision, ' \
              'ADD COLUMN IF NOT EXISTS z2b double precision, ' \
              'ADD COLUMN IF NOT EXISTS geom3d geometry("POLYGONZ", %s)'. \
              format(cfg.domain.case_schema, cfg.tables.slanted_wall)
    cur.execute(sqltext, (cfg.srid_palm,cfg.srid_palm,cfg.srid_palm,cfg.srid_palm,cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

    # Update Normal vectors
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'cent = ST_SetSRID(ST_Centroid(ST_Collect(ARRAY[point1, point2])), %s)' \
            .format(cfg.domain.case_schema, cfg.tables.slanted_wall)
    cur.execute(sqltext,(cfg.srid_palm,  ))
    sql_debug(connection)
    connection.commit()

    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'norm1 = ST_SetSRID(ST_Point(ST_X(cent) - (ST_Y(point1)-ST_Y(point2))/SQRT((ST_Y(point1)-ST_Y(point2))^2 + (ST_X(point1)-ST_X(point2))^2), ' \
              '                            ST_Y(cent) + (ST_X(point1)-ST_X(point2))/SQRT((ST_Y(point1)-ST_Y(point2))^2 + (ST_X(point1)-ST_X(point2))^2)), %s), ' \
              'norm2 = ST_SetSRID(ST_Point(ST_X(cent) + (ST_Y(point1)-ST_Y(point2))/SQRT((ST_Y(point1)-ST_Y(point2))^2 + (ST_X(point1)-ST_X(point2))^2), ' \
              '                            ST_Y(cent) - (ST_X(point1)-ST_X(point2))/SQRT((ST_Y(point1)-ST_Y(point2))^2 + (ST_X(point1)-ST_X(point2))^2)), %s)' \
            .format(cfg.domain.case_schema, cfg.tables.slanted_wall)
    cur.execute(sqltext,(cfg.srid_palm, cfg.srid_palm, ))
    sql_debug(connection)
    connection.commit()

    # Now decide which norm vector is the right one
    # sqltext = 'UPDATE "{0}"."{1}" AS sw SET norm = ' \
    #           'CASE WHEN ST_Distance((SELECT geom FROM "{0}"."{2}" WHERE type BETWEEN {3} AND {4} ORDER BY ST_Distance(geom, sw.geom) LIMIT 1), norm1) > ' \
    #           '          ST_Distance((SELECT geom FROM "{0}"."{2}" WHERE type BETWEEN {3} AND {4} ORDER BY ST_Distance(geom, sw.geom) LIMIT 1), norm2) ' \
    #           'THEN norm1 ' \
    #           'ELSE norm2 ' \
    #           'END'\
    #     .format(cfg.domain.case_schema, cfg.tables.slanted_wall, cfg.tables.landcover, cfg.type_range.building_min, cfg.type_range.building_max)
    # cur.execute(sqltext)
    # sql_debug(connection)
    # connection.commit()

    # sqltext = 'UPDATE "{0}"."{1}" AS sw SET norm = ' \
    #           'CASE WHEN ST_Distance((SELECT geom FROM "{0}"."{2}" WHERE type BETWEEN {3} AND {4} ORDER BY ST_Distance(geom, sw.cent) LIMIT 1), norm1) > ' \
    #           '          ST_Distance((SELECT geom FROM "{0}"."{2}" WHERE type BETWEEN {3} AND {4} ORDER BY ST_Distance(geom, sw.cent) LIMIT 1), norm2) ' \
    #           'THEN norm1 ' \
    #           'ELSE norm2 ' \
    #           'END'\
    #     .format(cfg.domain.case_schema, cfg.tables.slanted_wall, cfg.tables.landcover, cfg.type_range.building_min, cfg.type_range.building_max)
    # cur.execute(sqltext)
    # sql_debug(connection)
    # connection.commit()

    # debug('Deciding which normal vector is correct one (inside or outside of building)')
    # sqltext = 'UPDATE "{0}"."{1}" AS sw SET norm = ' \
    #           'CASE WHEN ST_Distance((SELECT geom FROM "{0}"."{2}" ORDER BY ST_Distance(geom, sw.cent) LIMIT 1), norm1) > ' \
    #           '          ST_Distance((SELECT geom FROM "{0}"."{2}" ORDER BY ST_Distance(geom, sw.cent) LIMIT 1), norm2) ' \
    #           'THEN norm1 ' \
    #           'ELSE norm2 ' \
    #           'END'\
    #     .format(cfg.domain.case_schema, cfg.tables.slanted_wall, cfg.tables.build_new)
    # cur.execute(sqltext)
    # sql_debug(connection)
    # connection.commit()

    # debug('Deciding which normal vector is correct one (inside or outside of building)')
    # sqltext = 'UPDATE "{0}"."{1}" AS sw SET norm = ' \
    #           'CASE WHEN ST_Distance((SELECT geom FROM "{0}"."{2}" ORDER BY ST_Distance(geom, sw.cent) LIMIT 1), ' \
    #           '                      ST_SetSRID(ST_MakePoint(ST_X(sw.cent)+{3}*(-ST_X(sw.cent)+ST_X(sw.norm1)), ' \
    #           '                                              ST_Y(sw.cent)+{3}*(-ST_Y(sw.cent)+ST_Y(sw.norm1))), %s))' \
    #           '          > ' \
    #           '          ST_Distance((SELECT geom FROM "{0}"."{2}" ORDER BY ST_Distance(geom, sw.cent) LIMIT 1), ' \
    #           '                      ST_SetSRID(ST_MakePoint(ST_X(sw.cent)+{3}*(-ST_X(sw.cent)+ST_X(sw.norm2)), ' \
    #           '                                              ST_Y(sw.cent)+{3}*(-ST_Y(sw.cent)+ST_Y(sw.norm2))), %s)) ' \
    #           'THEN norm1 ' \
    #           'ELSE norm2 ' \
    #           'END'\
    #     .format(cfg.domain.case_schema, cfg.tables.slanted_wall, cfg.tables.build_new, 0.50)
    # cur.execute(sqltext, (cfg.srid_palm, cfg.srid_palm,))
    # sql_debug(connection)
    # connection.commit()

    # # Now decide which split is the one
    # debug('Deciding which split take')
    # sqltext = 'UPDATE "{0}"."{1}" AS sw SET split = ' \
    #           'CASE WHEN ST_Distance(norm, split1) > ' \
    #           '          ST_Distance(norm, split2) ' \
    #           'THEN split1 ' \
    #           'ELSE split2 ' \
    #           'END'\
    #     .format(cfg.domain.case_schema, cfg.tables.slanted_wall)
    # cur.execute(sqltext)
    # sql_debug(connection)
    # connection.commit()

    debug('Deciding which split take')
    sqltext = 'WITH lb AS (SELECT geom FROM "{0}"."{2}")' \
              'UPDATE "{0}"."{1}" AS sw SET split = ' \
              'CASE WHEN ST_Area(ST_Intersection(lb.geom, split1)) / ST_Area(split1) > ' \
              '          ST_Area(ST_Intersection(lb.geom, split2)) / ST_Area(split2) ' \
              'THEN split1 ' \
              'ELSE split2 ' \
              'END ' \
              'FROM lb ' \
              'WHERE ST_Intersects(sw.geom,lb.geom)'\
        .format(cfg.domain.case_schema, cfg.tables.slanted_wall, cfg.tables.build_new)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # Now decide which split is the one
    debug('Deciding which normal vector is correct one (inside or outside of building)')
    sqltext = 'UPDATE "{0}"."{1}" AS sw SET norm = ' \
              'CASE WHEN ST_Distance(split, norm1) > ' \
              '          ST_Distance(split, norm2) ' \
              'THEN norm1 ' \
              'ELSE norm2 ' \
              'END'\
        .format(cfg.domain.case_schema, cfg.tables.slanted_wall)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # Now assign split that is for connection terrain and building
    debug('Deciding which split between terrain and building to take')
    sqltext = 'UPDATE "{0}"."{1}" AS sw SET split_terr_wall = ' \
              'CASE WHEN ST_Distance(norm, split1) < ' \
              '          ST_Distance(norm, split2) ' \
              'THEN split1 ' \
              'ELSE split2 ' \
              'END'\
        .format(cfg.domain.case_schema, cfg.tables.slanted_wall)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # Delete unnecessary columns point1, point2, norm1, norm2, split1, split2
    if cfg.slanted_pars.clean_up:
        sqltext = 'ALTER TABLE "{0}"."{1}" ' \
                  'DROP COLUMN norm1, ' \
                  'DROP COLUMN norm2, ' \
                  'DROP COLUMN split1, ' \
                  'DROP COLUMN split2, ' \
                  'DROP COLUMN geom' \
                  .format(cfg.domain.case_schema, cfg.tables.slanted_wall)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

    # Process walls heights
    debug('Calculation of z1 (top), z1b (bottom), z2, z2b in the slanted walls from buildings heights')
    max_dist = cfg.slanted_pars.wall_build_height_max_dist
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'z1 = FLOOR((SELECT height FROM "{0}"."{3}" hc' \
              '             WHERE ST_DWithin(point1, hc.geom, {4}) AND height IS NOT NULL ' \
              '             ORDER BY ST_Distance(point1, hc.geom) ' \
              '             LIMIT 1) / {2}) * {2}, ' \
              'z2 = FLOOR((SELECT height FROM "{0}"."{3}" hc' \
              '             WHERE ST_DWithin(point2, hc.geom, {4}) AND height IS NOT NULL ' \
              '             ORDER BY ST_Distance(point2, hc.geom) ' \
              '             LIMIT 1) / {2}) * {2}, ' \
              'z1b = FLOOR((SELECT height_bottom FROM "{0}"."{3}" hc' \
              '             WHERE ST_DWithin(point1, hc.geom, {4}) AND height IS NOT NULL ' \
              '             ORDER BY ST_Distance(point1, hc.geom) ' \
              '             LIMIT 1) / {2}) * {2}, ' \
              'z2b = FLOOR((SELECT height_bottom FROM "{0}"."{3}" hc' \
              '             WHERE ST_DWithin(point2, hc.geom, {4}) AND height IS NOT NULL ' \
              '             ORDER BY ST_Distance(point2, hc.geom) ' \
              '             LIMIT 1) / {2}) * {2} ' \
              ''.format(cfg.domain.case_schema, cfg.tables.slanted_wall, cfg.domain.dz,
                        cfg.tables.height_corrected, max_dist)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # process rest of the wall where height and height_bottom is missing
    debug('Correcting missing heights')
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'z1 = FLOOR((SELECT height FROM "{0}"."{3}" hc' \
              '             WHERE height IS NOT NULL ' \
              '             ORDER BY ST_Distance(point1, hc.geom) ' \
              '             LIMIT 1) / {2}) * {2}' \
              'WHERE z1 IS NULL' \
              .format(cfg.domain.case_schema, cfg.tables.slanted_wall, cfg.domain.dz, cfg.tables.height_corrected)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'z1b = FLOOR((SELECT height_bottom FROM "{0}"."{3}" hc' \
              '             WHERE height_bottom IS NOT NULL ' \
              '             ORDER BY ST_Distance(point1, hc.geom) ' \
              '             LIMIT 1) / {2}) * {2}' \
              'WHERE z1b IS NULL' \
              .format(cfg.domain.case_schema, cfg.tables.slanted_wall, cfg.domain.dz, cfg.tables.height_corrected)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'z2 = FLOOR((SELECT height FROM "{0}"."{3}" hc' \
              '             WHERE height IS NOT NULL ' \
              '             ORDER BY ST_Distance(point2, hc.geom) ' \
              '             LIMIT 1) / {2}) * {2}' \
              'WHERE z2 IS NULL' \
              .format(cfg.domain.case_schema, cfg.tables.slanted_wall, cfg.domain.dz, cfg.tables.height_corrected)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'z2b = FLOOR((SELECT height_bottom FROM "{0}"."{3}" hc' \
              '             WHERE height_bottom IS NOT NULL ' \
              '             ORDER BY ST_Distance(point2, hc.geom) ' \
              '             LIMIT 1) / {2}) * {2}' \
              'WHERE z2b IS NULL' \
              .format(cfg.domain.case_schema, cfg.tables.slanted_wall, cfg.domain.dz, cfg.tables.height_corrected)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # Create 3d polygon for debugging and checking
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'geom3d = ST_SetSRID(ST_MakePolygon(ST_MakeLine(ARRAY' \
              '[' \
              ' ST_MakePoint(ST_X(point1), ST_Y(point1), z1b), ' \
              ' ST_MakePoint(ST_X(point1), ST_Y(point1), z1 ), ' \
              ' ST_MakePoint(ST_X(point2), ST_Y(point2), z2 ), ' \
              ' ST_MakePoint(ST_X(point2), ST_Y(point2), z1b), ' \
              ' ST_MakePoint(ST_X(point1), ST_Y(point1), z1b)' \
              ']' \
              ')), %s)'\
              .format(cfg.domain.case_schema, cfg.tables.slanted_wall)
    cur.execute(sqltext,(cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

    # update wid
    debug('Updating wid')
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'wid = (SELECT wid FROM "{0}"."{2}" AS w' \
              '       ORDER BY ST_Distance(w.geom, geom3d) ' \
              '       LIMIT 1)'\
              .format(cfg.domain.case_schema, cfg.tables.slanted_wall, cfg.tables.walls)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    debug('Creating supplementary table {}', cfg.tables.slanted_wall_points)
    sqltext = 'DROP TABLE IF EXISTS "{0}"."{1}";' \
              'CREATE TABLE "{0}"."{1}" AS ' \
              '(SELECT ST_SetSRID(ST_MakePoint(ST_X(point1), ST_Y(point1)), %s) AS w_point, z1 AS z FROM "{0}"."{2}" ' \
              ' UNION ALL  ' \
              ' SELECT ST_SetSRID(ST_MakePoint(ST_X(point2), ST_Y(point2)), %s) AS w_point, z2 AS z FROM "{0}"."{2}"' \
              ' ); ' \
              'CREATE INDEX sw_point_geom_index ON "{0}"."{1}" USING gist(w_point)'\
              .format(cfg.domain.case_schema, cfg.tables.slanted_wall_points, cfg.tables.slanted_wall)
    cur.execute(sqltext, (cfg.srid_palm, cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

def create_slated_roof(cfg, connection, cur):
    """ Create slanted roof """
    progress('Creating slanted roof')
    sqltext = 'DROP TABLE IF EXISTS "{0}"."{1}"'.format(cfg.domain.case_schema, cfg.tables.slanted_roof)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    sqltext = 'CREATE TABLE "{0}"."{1}" AS SELECT ' \
              'id, i, j, CAST(NULL AS integer) AS rid, split FROM "{0}"."{2}"'\
              .format(cfg.domain.case_schema, cfg.tables.slanted_roof, cfg.tables.slanted_wall)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # add remaining buildings grid
    # sqltext = 'INSERT INTO "{0}"."{1}" SELECT ' \
    #           'bg.id, bg.i, bg.j, bg.geom ' \
    #           'FROM "{0}"."{2}" AS bg ' \
    #           'LEFT JOIN "{0}"."{1}" AS sr ON sr.id = bg.id ' \
    #           'WHERE sr.id IS NULL'.format(cfg.domain.case_schema, cfg.tables.slanted_roof,
    #                                     cfg.tables.buildings_grid)
    # cur.execute(sqltext)
    # sql_debug(connection)
    # connection.commit()

    # add rest of the roofs
    sqltext = 'WITH nb AS (SELECT geom FROM "{0}"."{3}")' \
              'INSERT INTO "{0}"."{1}" AS sr SELECT ' \
              'g.id, g.i, g.j, NULL, g.geom ' \
              'FROM "{0}"."{2}" AS g, nb ' \
              'WHERE (ST_Intersects(ST_SetSRID(ST_Point(g.xcen,g.ycen), %s), nb.geom) AND ' \
              '       g.id NOT IN (SELECT id FROM "{0}"."{1}")) ' \
              .format(cfg.domain.case_schema, cfg.tables.slanted_roof, cfg.tables.grid,
                      cfg.tables.build_new)
    cur.execute(sqltext, (cfg.srid_palm, ))
    sql_debug(connection)
    connection.commit()

    # put geom index on split polygon
    sqltext = 'CREATE INDEX slanted_roof_geom_idx ON "{0}"."{1}" USING gist(split)'.format(
        cfg.domain.case_schema, cfg.tables.slanted_roof)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # Add edge points in roof geoms, max 5 points
    sqltext = 'ALTER TABLE "{0}"."{1}" ' \
              'ADD COLUMN IF NOT EXISTS p1 geometry("POINT", %s), ' \
              'ADD COLUMN IF NOT EXISTS p2 geometry("POINT", %s), ' \
              'ADD COLUMN IF NOT EXISTS p3 geometry("POINT", %s), ' \
              'ADD COLUMN IF NOT EXISTS p4 geometry("POINT", %s), ' \
              'ADD COLUMN IF NOT EXISTS p5 geometry("POINT", %s), ' \
              'ADD COLUMN IF NOT EXISTS p6 geometry("POINT", %s), ' \
              'ADD COLUMN IF NOT EXISTS norm geometry("POINTZ", %s), ' \
              'ADD COLUMN IF NOT EXISTS cent geometry("POINTZ", %s), ' \
              'ADD COLUMN IF NOT EXISTS geom3d geometry("POLYGONZ", %s), ' \
              'ADD COLUMN IF NOT EXISTS z1 double precision, ' \
              'ADD COLUMN IF NOT EXISTS z2 double precision, ' \
              'ADD COLUMN IF NOT EXISTS z3 double precision, ' \
              'ADD COLUMN IF NOT EXISTS z4 double precision, ' \
              'ADD COLUMN IF NOT EXISTS z5 double precision, ' \
              'ADD COLUMN IF NOT EXISTS z6 double precision, ' \
              'ADD COLUMN IF NOT EXISTS n_edges integer'. \
              format(cfg.domain.case_schema, cfg.tables.slanted_roof)
    cur.execute(sqltext, (cfg.srid_palm,cfg.srid_palm,cfg.srid_palm,
                          cfg.srid_palm,cfg.srid_palm,cfg.srid_palm,
                          cfg.srid_palm,cfg.srid_palm,cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

    # fill the points
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'n_edges = ST_NPoints(split)-1, ' \
              'p1 = CASE WHEN ST_NPoints(split) > 1 THEN ST_PointN(ST_Boundary(split),1) ELSE NULL END, ' \
              'p2 = CASE WHEN ST_NPoints(split) > 1 THEN ST_PointN(ST_Boundary(split),2) ELSE NULL END, ' \
              'p3 = CASE WHEN ST_NPoints(split) > 1 THEN ST_PointN(ST_Boundary(split),3) ELSE NULL END, ' \
              'p4 = CASE WHEN ST_NPoints(split) > 1 THEN ST_PointN(ST_Boundary(split),4) ELSE NULL END, ' \
              'p5 = CASE WHEN ST_NPoints(split) > 1 THEN ST_PointN(ST_Boundary(split),5) ELSE NULL END'.\
              format(cfg.domain.case_schema, cfg.tables.slanted_roof)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # TODO: Add p6 point in case of z1 != z2 from walls
    # sqltext = 'UPDATE "{0}"."{1}" AS sr SET ' \
    #           'n_edges = ST_NPoints(sr.split)-1+1, ' \
    #           'p6 = CASE WHEN sw.z1 > sw.z2 THEN ST_SetSRID(ST_MakePoint(ST_X(sw.point1), ST_Y(sw.point1)), %s) ' \
    #           '         ELSE ST_SetSRID(ST_MakePoint(ST_X(sw.point2), ST_Y(sw.point2)), %s) END, ' \
    #           'z6 = CASE WHEN sw.z1 > sw.z2 THEN sw.z2 ELSE sw.z1 END ' \
    #           'FROM "{0}"."{2}" AS sw ' \
    #           'WHERE sw.id = sr.id AND sw.z1 != sw.z2'.\
    #           format(cfg.domain.case_schema, cfg.tables.slanted_roof, cfg.tables.slanted_wall)
    # cur.execute(sqltext, (cfg.srid_palm, cfg.srid_palm,))
    # sql_debug(connection)
    # connection.commit()

    # add indexes on p1 .. p6
    sqltext = 'CREATE INDEX p1_geom_index ON "{0}"."{1}" USING gist(p1); ' \
              'CREATE INDEX p2_geom_index ON "{0}"."{1}" USING gist(p2); ' \
              'CREATE INDEX p3_geom_index ON "{0}"."{1}" USING gist(p3); ' \
              'CREATE INDEX p4_geom_index ON "{0}"."{1}" USING gist(p4); ' \
              'CREATE INDEX p5_geom_index ON "{0}"."{1}" USING gist(p5); ' \
              'CREATE INDEX p6_geom_index ON "{0}"."{1}" USING gist(p6); '. \
              format(cfg.domain.case_schema, cfg.tables.slanted_roof)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()


    #### SOLUTION WITH USING DIRECTLY RASTER
    # sqltext = 'DROP TABLE IF EXISTS nb_temp;' \
    #           'CREATE TEMP TABLE nb_temp AS ' \
    #           'SELECT (ST_PixelAsPoints(b.rast)).geom as geom, (ST_PixelAsPoints(b.rast)).val as val ' \
    #               '            FROM "{0}"."{1}" AS b'.format(cfg.domain.case_schema, cfg.tables.buildings_height)
    # cur.execute(sqltext)
    # sql_debug(connection)
    # connection.commit()
    #
    # sqltext = 'CREATE INDEX nb_geom_index ON nb_temp USING gist(geom)'
    # cur.execute(sqltext)
    # sql_debug(connection)
    # connection.commit()
    #
    # # Add z1, z2 to slanted walls
    # debug('Calculation of z1, z2 in the slanted walls from buildings heights')
    # sqltext = 'UPDATE "{0}"."{1}" SET ' \
    #           'z1 = FLOOR((SELECT ST_Value(rast, point1) FROM "{0}"."{3}" WHERE ST_Intersects(rast, point1)) / {2}), ' \
    #           'z2 = FLOOR((SELECT ST_Value(rast, point2) FROM "{0}"."{3}" WHERE ST_Intersects(rast, point2)) / {2})' \
    #           ''.format(cfg.domain.case_schema, cfg.tables.slanted_wall, cfg.domain.dz, cfg.tables.buildings_height)
    # cur.execute(sqltext)
    # sql_debug(connection)
    # connection.commit()
    #
    # # fill the remaining heights z1, z2
    # sqltext = 'UPDATE "{0}"."{1}" SET ' \
    #           'z1 = FLOOR((SELECT nb.val FROM nb_temp AS nb ORDER BY ST_Distance(nb.geom, point1) LIMIT 1 ) / {2})' \
    #           'WHERE z1 IS NULL'.format(cfg.domain.case_schema, cfg.tables.slanted_wall, cfg.domain.dz, cfg.tables.buildings_height)
    # cur.execute(sqltext)
    # sql_debug(connection)
    # connection.commit()
    #
    # sqltext = 'UPDATE "{0}"."{1}" SET ' \
    #           'z2 = FLOOR((SELECT nb.val FROM nb_temp AS nb ORDER BY ST_Distance(nb.geom, point2) LIMIT 1 ) / {2})' \
    #           'WHERE z2 IS NULL'.format(cfg.domain.case_schema, cfg.tables.slanted_wall, cfg.domain.dz, cfg.tables.buildings_height)
    # cur.execute(sqltext)
    # sql_debug(connection)
    # connection.commit()

    #### SOLUTION USING PREPROCESSED CORRECTED BUILDING HEIGHTS


    ### USING TEMP TABLE FROM RASTER
    # # fill height
    # sqltext = 'UPDATE "{0}"."{1}" SET ' \
    #           'z1 = (CASE WHEN n_edges >= 1 THEN (SELECT ST_Value(rast, p1) FROM "{0}"."{2}" WHERE ST_Intersects(rast, p1)) ELSE NULL END), '\
    #           'z2 = (CASE WHEN n_edges >= 2 THEN (SELECT ST_Value(rast, p2) FROM "{0}"."{2}" WHERE ST_Intersects(rast, p2)) ELSE NULL END), ' \
    #           'z3 = (CASE WHEN n_edges >= 3 THEN (SELECT ST_Value(rast, p3) FROM "{0}"."{2}" WHERE ST_Intersects(rast, p3)) ELSE NULL END), ' \
    #           'z4 = (CASE WHEN n_edges >= 4 THEN (SELECT ST_Value(rast, p4) FROM "{0}"."{2}" WHERE ST_Intersects(rast, p4)) ELSE NULL END), ' \
    #           'z5 = (CASE WHEN n_edges >= 5 THEN (SELECT ST_Value(rast, p5) FROM "{0}"."{2}" WHERE ST_Intersects(rast, p5)) ELSE NULL END)'.\
    #           format(cfg.domain.case_schema, cfg.tables.slanted_roof, cfg.tables.buildings_height)
    # cur.execute(sqltext)
    # sql_debug(connection)
    # connection.commit()
    #
    # # fill remaining heights
    # # TODO: definitelly optimize, temporally table nb, add spacial index
    # for ie in [1, 2, 3, 4, 5]:
    #     debug('Filling missing z{} from spit average', ie)
    #     sqltext = 'UPDATE "{0}"."{1}" AS sr ' \
    #               'SET z{3} = (SELECT nb.val FROM nb_temp AS nb  ORDER BY ST_Distance(nb.geom, sr.p{3}) LIMIT 1) ' \
    #               ' ' \
    #               'WHERE z{3} IS NULL AND {3} <= n_edges'\
    #               .format(cfg.domain.case_schema, cfg.tables.slanted_roof, cfg.tables.buildings_height, ie)
    #     cur.execute(sqltext)
    #     sql_debug(connection)
    #     connection.commit()

    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'z1 = (CASE WHEN n_edges >= 1 THEN (SELECT height FROM "{0}"."{2}" hc' \
              '             WHERE ST_DWithin(p1, hc.geom, {3}) AND height IS NOT NULL ' \
              '             ORDER BY ST_Distance(p1, hc.geom) ' \
              '             LIMIT 1) ELSE NULL END), '\
              'z2 = (CASE WHEN n_edges >= 2 THEN (SELECT height FROM "{0}"."{2}" hc' \
              '             WHERE ST_DWithin(p2, hc.geom, {3}) AND height IS NOT NULL ' \
              '             ORDER BY ST_Distance(p2, hc.geom) ' \
              '             LIMIT 1) ELSE NULL END), ' \
              'z3 = (CASE WHEN n_edges >= 3 THEN (SELECT height FROM "{0}"."{2}" hc' \
              '             WHERE ST_DWithin(p3, hc.geom, {3}) AND height IS NOT NULL ' \
              '             ORDER BY ST_Distance(p3, hc.geom) ' \
              '             LIMIT 1) ELSE NULL END), ' \
              'z4 = (CASE WHEN n_edges >= 4 THEN (SELECT height FROM "{0}"."{2}" hc' \
              '             WHERE ST_DWithin(p4, hc.geom, {3}) AND height IS NOT NULL ' \
              '             ORDER BY ST_Distance(p4, hc.geom) ' \
              '             LIMIT 1) ELSE NULL END), ' \
              'z5 = (CASE WHEN n_edges >= 5 THEN (SELECT height FROM "{0}"."{2}" hc' \
              '             WHERE ST_DWithin(p5, hc.geom, {3}) AND height IS NOT NULL ' \
              '             ORDER BY ST_Distance(p5, hc.geom) ' \
              '             LIMIT 1) ELSE NULL END)'.\
              format(cfg.domain.case_schema, cfg.tables.slanted_roof, cfg.tables.height_corrected, cfg.slanted_pars.wall_build_height_max_dist)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # fill the missing ones
    debug('Filling missing height in slanted roof')
    for pi in range(1,6):
        verbose('Filling missing {} height in slanted roof tables', pi)
        sqltext = 'UPDATE "{0}"."{1}" AS sr SET ' \
                  'z{4} = (SELECT z FROM "{0}"."{2}" AS sw ' \
                  '        WHERE ST_DWithin(w_point, p{4}, {3}) ' \
                  '        ORDER BY p{4} <-> sw.w_point LIMIT 1)' \
                  'WHERE z{4} IS NULL AND p{4} IS NOT NULL'\
                  .format(cfg.domain.case_schema, cfg.tables.slanted_roof, cfg.tables.slanted_wall_points,
                          4.0 * cfg.domain.dx, pi)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

    # modify height at the roof edge, set to floor(height)
    # FIXME FLOOR will not work if dz != 1
    # FIXME optimize using spatial index
    roofs_dist2edge = cfg.slanted_pars.roofs_dist2edge
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'z1 = CASE WHEN (SELECT ST_Distance(p1, geom) FROM "{0}"."{3}" ORDER BY ST_Distance(p1, geom) LIMIT 1) < {4} THEN FLOOR(z1 / {2}) * {2} ELSE z1 END, '\
              'z2 = CASE WHEN (SELECT ST_Distance(p2, geom) FROM "{0}"."{3}" ORDER BY ST_Distance(p2, geom) LIMIT 1) < {4} THEN FLOOR(z2 / {2}) * {2} ELSE z2 END, ' \
              'z3 = CASE WHEN (SELECT ST_Distance(p3, geom) FROM "{0}"."{3}" ORDER BY ST_Distance(p3, geom) LIMIT 1) < {4} THEN FLOOR(z3 / {2}) * {2} ELSE z3 END, ' \
              'z4 = CASE WHEN (SELECT ST_Distance(p4, geom) FROM "{0}"."{3}" ORDER BY ST_Distance(p4, geom) LIMIT 1) < {4} THEN FLOOR(z4 / {2}) * {2} ELSE z4 END, ' \
              'z5 = CASE WHEN (SELECT ST_Distance(p5, geom) FROM "{0}"."{3}" ORDER BY ST_Distance(p5, geom) LIMIT 1) < {4} THEN FLOOR(z5 / {2}) * {2} ELSE z5 END '.\
              format(cfg.domain.case_schema, cfg.tables.slanted_roof, cfg.domain.dz, cfg.tables.walls_outer, roofs_dist2edge)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()
    # dist2edge = np.sqrt(2.0 * cfg.domain.dx) * 1.1
    roofs_dist2edge = cfg.slanted_pars.roofs_dist2edge
    for pi in range(1,6):
        verbose('Correcting point {} height in slanted roof tables', pi)
        sqltext = 'UPDATE "{0}"."{1}" AS sr SET ' \
                  'z{4} = (SELECT z FROM "{0}"."{2}" AS sw ' \
                  '        WHERE ST_DWithin(w_point, p{4}, {3}) ' \
                  '        ORDER BY p{4} <-> sw.w_point LIMIT 1) ' \
                  'WHERE (SELECT p{4} <-> sw.w_point  ' \
                  '       FROM "{0}"."{2}" AS sw ' \
                  '       WHERE ST_DWithin(w_point, p{4}, {3}) ORDER BY p{4} <-> sw.w_point LIMIT 1 )' \
                  '       < {3} ' \
                  'AND (SELECT z FROM "{0}"."{2}" AS sw ' \
                  '     WHERE ST_DWithin(w_point, p{4}, {3}) ' \
                  '     ORDER BY p{4} <-> sw.w_point LIMIT 1) > z{4}'\
                  .format(cfg.domain.case_schema, cfg.tables.slanted_roof, cfg.tables.slanted_wall_points, roofs_dist2edge, pi)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

    # Further correct --> WITH all wall points aggregated
    # If dist(point, wall geom) < roofs_dist2edge THEN zi = (ST_Z(points) ORDER BY ST_DISTANCE(pi, points))
    # Or do that properly and add wall point into building height and then do some filtering there.

    #FIXME
    # # Correct height in near distance to edge
    # for ie in [1, 2, 3, 4, 5]:
    #     debug('Correcting height of p{} that are near edges', ie)
    #     sqltext = 'WITH p1c AS ' \
    #               '    (SELECT sr.i,sr.j, MAX(sw1.z1) AS zz1, MAX(sw2.z2) AS zz2 ' \
    #               '     FROM "{0}"."{1}" AS sr ' \
    #               '     LEFT JOIN "{0}"."{2}" AS sw1 ON ST_DWithin(p{3}, sw1.point1, {4} ) AND ST_Distance(p{3}, sw1.wall_geom) > 0.001  ' \
    #               '     LEFT JOIN "{0}"."{2}" AS sw2 ON ST_DWithin(p{3}, sw2.point2, {4} ) AND ST_Distance(p{3}, sw2.wall_geom) > 0.001 ' \
    #               '     WHERE ((ST_DWithin(sr.p{3}, sw1.point1, {4} ) AND NOT ST_DWithin(sr.p{3}, sw1.point1, 0.001 )) OR ' \
    #               '            (ST_DWithin(sr.p{3}, sw2.point2, {4} ) AND NOT ST_DWithin(sr.p{3}, sw1.point2, 0.001 ))) AND ' \
    #               '            (sr.z{3} < sw1.z1 AND sr.z{3} < sw2.z2) ' \
    #               '     GROUP BY sr.i, sr.j ' \
    #               '     )' \
    #               'UPDATE "{0}"."{1}" AS sr ' \
    #               'SET z{3} = CASE WHEN zz1 > zz2 THEN zz1 ELSE zz2 END  ' \
    #               'FROM p1c ' \
    #               'WHERE p1c.i = sr.i AND p1c.j = sr.j '.\
    #               format(cfg.domain.case_schema, cfg.tables.slanted_roof, cfg.tables.slanted_wall, ie, dist2edge)
    #     cur.execute(sqltext)
    #     sql_debug(connection)
    #     connection.commit()

    # find center and normal vector
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'cent = ST_SetSRID(ST_MakePoint(ST_X(ST_Centroid(split)), ' \
              '                               ST_Y(ST_Centroid(split)), ' \
              '                               (COALESCE(z1, 0.0) + COALESCE(z2, 0.0) + COALESCE(z3, 0.0) + COALESCE(z4, 0.0) + COALESCE(z5, 0.0))/n_edges), %s) ' \
              ''.\
              format(cfg.domain.case_schema, cfg.tables.slanted_roof)
    cur.execute(sqltext, (cfg.srid_palm, ))
    sql_debug(connection)
    connection.commit()

    # sqltext = 'UPDATE "{0}"."{1}" SET ' \
    #           'norm = ST_SetSRID(ST_MakePoint(' \
    #           'ST_X(cent) + SIN((SELECT ST_NearestValue(st_aspect, cent) FROM "{0}"."{2}" WHERE ST_Intersects(st_aspect, cent))*PI()/180)*' \
    #           'SIN((SELECT ST_NearestValue(st_slope, cent) FROM "{0}"."{3}" WHERE ST_Intersects(st_slope, cent))*PI()/180) , ' \
    #           'ST_Y(cent) + COS((SELECT ST_NearestValue(st_aspect, cent) FROM "{0}"."{2}" WHERE ST_Intersects(st_aspect, cent))*PI()/180)*' \
    #           'SIN((SELECT ST_NearestValue(st_slope, cent) FROM "{0}"."{3}" WHERE ST_Intersects(st_slope, cent))*PI()/180) ,  ' \
    #           'ST_Z(cent) + COS((SELECT ST_NearestValue(st_slope, cent)  FROM "{0}"."{3}" WHERE ST_Intersects(st_slope, cent))*PI()/180)' \
    #           '), %s) ' \
    #           ''.\
    #           format(cfg.domain.case_schema, cfg.tables.slanted_roof, cfg.tables.aspect, cfg.tables.slope)
    # cur.execute(sqltext, (cfg.srid_palm, ))
    # sql_debug(connection)
    # connection.commit()

    # Create 3D polygon of the roof
    debug('Creating 3d polygon')
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'geom3d = ST_SetSRID(ST_ConvexHull(ST_Collect(ARRAY[' \
              'ST_MakePoint(ST_X(p1), ST_Y(p1), z1), ST_MakePoint(ST_X(p2), ST_Y(p2), z2), ST_MakePoint(ST_X(p3), ST_Y(p3), z3), ' \
              'ST_MakePoint(ST_X(p4), ST_Y(p4), z4), ST_MakePoint(ST_X(p5), ST_Y(p5), z5)' \
              '])), %s) ' \
              'WHERE ST_NPoints(ST_Collect(ARRAY[' \
              '         ST_MakePoint(ST_X(p1), ST_Y(p1), z1), ST_MakePoint(ST_X(p2), ST_Y(p2), z2), ST_MakePoint(ST_X(p3), ST_Y(p3), z3), ' \
              '         ST_MakePoint(ST_X(p4), ST_Y(p4), z4), ST_MakePoint(ST_X(p5), ST_Y(p5), z5)' \
              '])) > 2'.format(cfg.domain.case_schema, cfg.tables.slanted_roof)
    cur.execute(sqltext, (cfg.srid_palm, ))
    sql_debug(connection)
    connection.commit()

    # # Delete Unnecessary columns
    # sqltext = 'ALTER TABLE "{0}"."{1}" ' \
    #           'DROP COLUMN p1, DROP COLUMN p2, DROP COLUMN p3, DROP COLUMN p4, DROP COLUMN p5, ' \
    #           'DROP COLUMN z1, DROP COLUMN z2, DROP COLUMN z3, DROP COLUMN z4, DROP COLUMN z5, ' \
    #           'DROP COLUMN split' \
    #           .format(cfg.domain.case_schema, cfg.tables.slanted_roof)
    # cur.execute(sqltext)
    # sql_debug(connection)
    # connection.commit()

    debug('Updating rid')
    # TODO: Add ST_Dwithin
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'rid = (SELECT rid FROM "{0}"."{2}" AS r' \
              '       ORDER BY ST_Distance(r.geom, geom3d) ' \
              '       LIMIT 1)'\
              .format(cfg.domain.case_schema, cfg.tables.slanted_roof, cfg.tables.roofs)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

def create_grid_slanted_walls(cfg, connection, cur):
    """ Grid slanted walls """
    progress('Processing slanted walls into gridded structure')
    # CREATE TABLE of id, gridded geoms ... in each cycle fill it with geoms
    sqltext = 'DROP TABLE IF EXISTS "{0}"."{1}"'.format(cfg.domain.case_schema, cfg.tables.slanted_wall_gridded)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    sqltext = 'CREATE TABLE "{0}"."{1}" (' \
              'id integer, wid integer, ' \
              'geom geometry("POLYGONZ", %s), ' \
              'norm geometry("POINTZ", %s)' \
              ')'.format(cfg.domain.case_schema, cfg.tables.slanted_wall_gridded)
    cur.execute(sqltext, (cfg.srid_palm,cfg.srid_palm, ))
    sql_debug(connection)
    connection.commit()

    # CREATE TEMP TABLE, destroy after cycle
    # FIND MAX HEIGHT -> k_max
    sqltext = 'SELECT MAX(z1), MAX(z2) FROM "{0}"."{1}"'\
              .format(cfg.domain.case_schema, cfg.tables.slanted_wall)
    cur.execute(sqltext)
    z_max = cur.fetchall()
    z_max = [x for x in z_max[0]]
    k_max = int(max(z_max) / cfg.domain.dz) + 1
    sql_debug(connection)
    connection.commit()

    for k in range(1,k_max):
        sqltext = 'DROP TABLE IF EXISTS "{0}"."{1}"'.format(cfg.domain.case_schema, cfg.tables.slanted_wall_gridded_temp)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()


        sqltext = 'CREATE TABLE "{0}"."{1}" AS SELECT id, wid, ' \
                  'ST_X(point1) AS x1, ST_Y(point1) AS y1, ' \
                  'ST_X(point2) AS x2, ST_Y(point2) AS y2,' \
                  'z1, z2, z1b, z2b, norm  ' \
                  'FROM "{0}"."{2}"'\
                  .format(cfg.domain.case_schema, cfg.tables.slanted_wall_gridded_temp, cfg.tables.slanted_wall)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        sqltext = 'ALTER TABLE "{0}"."{1}" ' \
                  '{2} t_T  {3}, {2} t_L {3},  {2} x1L {3}, {2} y1L  {3}, {2} z1L  {3}, ' \
                  '{2} x1T  {3}, {2} y1T {3},  {2} z1T {3}, {2} x2L  {3}, {2} y2L  {3}, {2} z2L  {3}, ' \
                  '{2} x2T  {3}, {2} y2T {3},  {2} z2T {3}, {2} x12L {3}, {2} y12L {3}, {2} z12L {3}, ' \
                  '{2} x12T {3}, {2} y12T {3}, {2} z12T {3}, ' \
                  '{2} z_cent {3}, {2} n_vert integer, ' \
                  '{2} p1L {4}, {2} p1T {4}, {2} p2L {4}, {2} p2T {4}, {2} p12L {4}, {2} p12T {4} ' \
                  ''.format(cfg.domain.case_schema, cfg.tables.slanted_wall_gridded_temp,
                    'ADD COLUMN IF NOT EXISTS', 'double precision', 'geometry(POINTZ, %s)')
        cur.execute(sqltext, (cfg.srid_palm,cfg.srid_palm,cfg.srid_palm,
                              cfg.srid_palm,cfg.srid_palm,cfg.srid_palm,))
        sql_debug(connection)
        connection.commit()

        sqltext = 'UPDATE "{0}"."{1}" SET ' \
                  't_T = CASE WHEN z1 >= {2} AND z2 >= {2} THEN 0.0 ' \
                      'WHEN z1 <= {2} AND z2 >= {2} THEN ({2} - z1) / (z2 - z1)' \
                      'WHEN z1 >= {2} AND z2 <= {2} THEN ({2} - z1) / (z2 - z1)' \
                      'ELSE 0.0 END,' \
                  't_L = CASE WHEN z1 >= {2} AND z2 >= {2} THEN 0.0 ' \
                      'WHEN z1 >= ({2} - {3}) AND z2 >= ({2} - {3}) THEN 0.0' \
                      'WHEN z1 <= ({2} - {3}) AND z2 >= ({2} - {3}) THEN ({2} - {3} - z1) / (z2 - z1)' \
                      'WHEN z1 >= ({2} - {3}) AND z2 <= ({2} - {3}) THEN ({2} - {3} - z1) / (z2 - z1)' \
                      'ELSE 0.0 END' \
                  .format(cfg.domain.case_schema, cfg.tables.slanted_wall_gridded_temp, k*cfg.domain.dz, cfg.domain.dz)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        sqltext = 'UPDATE "{0}"."{1}" SET ' \
                  'x1L = CASE WHEN z1 < ({2} - {3}) OR z1b > ({2}) THEN NULL ELSE x1 END, ' \
                  'y1L = CASE WHEN z1 < ({2} - {3}) OR z1b > ({2}) THEN NULL ELSE y1 END, ' \
                  'z1L = CASE WHEN z1 < ({2} - {3}) OR z1b > ({2}) THEN NULL ELSE ({2} - {3}) END, ' \
                  'x1T = CASE WHEN z1 < ({2} - {3}) OR z1b > ({2}) THEN NULL ELSE x1 END, ' \
                  'y1T = CASE WHEN z1 < ({2} - {3}) OR z1b > ({2}) THEN NULL ELSE y1 END, ' \
                  'z1T = CASE WHEN z1 < ({2} - {3}) OR z1b > ({2}) THEN NULL WHEN z1 < {2} THEN z1 ELSE {2} END, ' \
                  'x2L = CASE WHEN z2 < ({2} - {3}) OR z2b > ({2}) THEN NULL ELSE x2  END, ' \
                  'y2L = CASE WHEN z2 < ({2} - {3}) OR z2b > ({2}) THEN NULL ELSE y2  END, ' \
                  'z2L = CASE WHEN z2 < ({2} - {3}) OR z2b > ({2}) THEN NULL ELSE ({2} - {3}) END, ' \
                  'x2T = CASE WHEN z2 < ({2} - {3}) OR z2b > ({2}) THEN NULL ELSE x2 END, ' \
                  'y2T = CASE WHEN z2 < ({2} - {3}) OR z2b > ({2}) THEN NULL ELSE y2 END, ' \
                  'z2T = CASE WHEN z2 < ({2} - {3}) OR z2b > ({2}) THEN NULL WHEN z2 < {2} THEN z2 ELSE {2} END, ' \
                  'x12L =CASE WHEN (z1 < ({2} - {3}) AND z2 > ({2} - {3})) OR ' \
                  '                (z1 > ({2} - {3}) AND z2 < ({2} - {3})) ' \
                  '           THEN x1 + t_L*(x2-x1) ELSE NULL END, ' \
                  'y12L =CASE WHEN (z1 < ({2} - {3}) AND z2 > ({2} - {3})) OR ' \
                  '                (z1 > ({2} - {3}) AND z2 < ({2} - {3})) ' \
                  '           THEN y1 + t_L*(y2-y1) ELSE NULL END, ' \
                  'z12L =CASE WHEN (z1 < ({2} - {3}) AND z2 > ({2} - {3})) OR' \
                  '                (z1 > ({2} - {3}) AND z2 < ({2} - {3})) ' \
                  '           THEN ({2} - {3}) ELSE NULL END, ' \
                  'x12T =CASE WHEN (z1 < {2} AND z2 > {2}) OR (z1 > {2} AND z2 < {2}) ' \
                  '           THEN x1 + t_T*(x2-x1) ELSE NULL END, ' \
                  'y12T =CASE WHEN (z1 < {2} AND z2 > {2}) OR (z1 > {2} AND z2 < {2}) ' \
                  '           THEN y1 + t_T*(y2-y1) ELSE NULL END, ' \
                  'z12T =CASE WHEN (z1 < {2} AND z2 > {2}) OR (z1 > {2} AND z2 < {2}) THEN {2} ELSE NULL END '\
                  .format(cfg.domain.case_schema, cfg.tables.slanted_wall_gridded_temp,
                          k*cfg.domain.dz, cfg.domain.dz)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        # calculate number of vertices
        sqltext = 'UPDATE "{0}"."{1}" SET n_vert = ' \
                  'CASE WHEN z1L  IS NOT NULL THEN 1 ELSE 0 END + ' \
                  'CASE WHEN z1T  IS NOT NULL THEN 1 ELSE 0 END + ' \
                  'CASE WHEN z2L  IS NOT NULL THEN 1 ELSE 0 END + ' \
                  'CASE WHEN z2T  IS NOT NULL THEN 1 ELSE 0 END + ' \
                  'CASE WHEN z12L IS NOT NULL THEN 1 ELSE 0 END + ' \
                  'CASE WHEN z12T IS NOT NULL THEN 1 ELSE 0 END'\
                  .format(cfg.domain.case_schema, cfg.tables.slanted_wall_gridded_temp)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        # calculate z center of the wall
        sqltext = 'UPDATE "{0}"."{1}" SET z_cent = ' \
                  '(COALESCE(z1L, 0) + COALESCE(z1T,  0) + COALESCE(z2L,  0) + ' \
                  ' COALESCE(z2T, 0) + COALESCE(z12L, 0) + COALESCE(z12T, 0) ' \
                  ') / n_vert ' \
                  'WHERE n_vert > 0' \
                  .format(cfg.domain.case_schema, cfg.tables.slanted_wall_gridded_temp)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        sqltext = 'UPDATE "{0}"."{1}" SET ' \
                  'p1L = CASE WHEN x1L  IS NOT NULL THEN ST_SetSRID(ST_MakePoint(x1L, y1L, z1L), %s) ELSE NULL END, ' \
                  'p1T = CASE WHEN x1T  IS NOT NULL THEN ST_SetSRID(ST_MakePoint(x1T, y1T, z1T), %s) ELSE NULL END, ' \
                  'p2L = CASE WHEN x2L  IS NOT NULL THEN ST_SetSRID(ST_MakePoint(x2L, y2L, z2L), %s) ELSE NULL END, ' \
                  'p2T = CASE WHEN x2T  IS NOT NULL THEN ST_SetSRID(ST_MakePoint(x2T, y2T, z2T), %s) ELSE NULL END, ' \
                  'p12L= CASE WHEN x12L IS NOT NULL THEN ST_SetSRID(ST_MakePoint(x12L,y12L,z12L),%s) ELSE NULL END, ' \
                  'p12T= CASE WHEN x12T IS NOT NULL THEN ST_SetSRID(ST_MakePoint(x12T,y12T,z12T),%s) ELSE NULL END' \
                    .format(cfg.domain.case_schema, cfg.tables.slanted_wall_gridded_temp)
        cur.execute(sqltext, (cfg.srid_palm,cfg.srid_palm,cfg.srid_palm,
                              cfg.srid_palm,cfg.srid_palm,cfg.srid_palm,))
        sql_debug(connection)
        connection.commit()

        sqltext = 'INSERT INTO "{0}"."{1}" ' \
                  'SELECT id, wid, ST_SetSRID(ST_MakePolygon(ST_MakeLine(ARRAY[p1L, p1T, p12T, p2T, p2L, p12L, ' \
                  'CASE ' \
                  'WHEN p1L  IS NOT NULL THEN p1L ' \
                  'WHEN p1T  IS NOT NULL THEN p1T ' \
                  'WHEN p12T IS NOT NULL THEN p12T ' \
                  'WHEN p2T  IS NOT NULL THEN p2T ' \
                  'WHEN p2L  IS NOT NULL THEN p2L ' \
                  'WHEN p12L IS NOT NULL THEN p12L ' \
                  'END])), %s), ' \
                  'ST_SetSRID(ST_MakePoint(ST_X(norm), ST_Y(norm), z_cent), %s)' \
                  'FROM "{0}"."{2}"' \
                  'WHERE (CASE WHEN p1L  IS NULL THEN 0 ELSE 1 END + CASE WHEN p1T  IS NULL THEN 0 ELSE 1 END  + ' \
                  '       CASE WHEN p2L  IS NULL THEN 0 ELSE 1 END + CASE WHEN p2T  IS NULL THEN 0 ELSE 1 END  + ' \
                  '       CASE WHEN p12L IS NULL THEN 0 ELSE 1 END + CASE WHEN p12T IS NULL THEN 0 ELSE 1 END ' \
                  '      )>2 ' \
                    .format(cfg.domain.case_schema, cfg.tables.slanted_wall_gridded, cfg.tables.slanted_wall_gridded_temp)
        cur.execute(sqltext, (cfg.srid_palm,cfg.srid_palm,))
        sql_debug(connection)
        connection.commit()

def create_grid_slanted_terrain(cfg, connection, cur):
    """ Process slanted terrain into gridded slanted terrain """
    progress('Gridding the slanted terrain')
    sqltext = 'DROP TABLE IF EXISTS "{0}"."{1}";' \
              'CREATE TABLE "{0}"."{1}" (' \
              'id integer, lid integer,  ' \
              'geom geometry("POLYGONZ", %s), ' \
              'points geometry("MULTIPOINTZ", %s),' \
              'lines geometry("LINESTRINGZ", %s)' \
              ')'.format(cfg.domain.case_schema, cfg.tables.slanted_terrain_gridded)
    cur.execute(sqltext, (cfg.srid_palm,cfg.srid_palm,cfg.srid_palm, ))
    sql_debug(connection)
    connection.commit()

    # get terrain max height
    debug('Selecting max height')
    sqltext = 'SELECT MAX(height) FROM "{0}"."{1}"' \
              .format(cfg.domain.case_schema, cfg.tables.height_terr_corrected)
    cur.execute(sqltext)
    z_max = cur.fetchone()[0]
    k_max = int(np.ceil(z_max) / cfg.domain.dz) + 2
    sql_debug(connection)
    connection.commit()

    debug('Loop over all k-levels')
    for k in range(0, k_max):
        debug('k={}', k)
        # Create temp table
        verbose('Drop temp table')
        sqltext = 'DROP TABLE IF EXISTS "{0}"."{1}"'.format(cfg.domain.case_schema, cfg.tables.slanted_terrain_gridded_temp)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        verbose('Create temp table')
        sqltext = 'CREATE TABLE "{0}"."{1}" AS SELECT id, lid, geom3d ' \
                  'FROM "{0}"."{2}"'\
                  .format(cfg.domain.case_schema, cfg.tables.slanted_terrain_gridded_temp, cfg.tables.slanted_terrain)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        verbose('Add new columns in temp table')
        sqltext = 'ALTER TABLE "{0}"."{1}" ' \
                  '{2} npoints integer, ' \
                  '{2} z1  {3}, {2} z2  {3}, ' \
                  '{2} mezi_dolni boolean DEFAULT FALSE, {2} mezi_horni boolean DEFAULT FALSE, ' \
                  '{2} p1 {4}, {2} p2 {4}' \
                  ''.format(cfg.domain.case_schema, cfg.tables.slanted_terrain_gridded_temp,
                            'ADD COLUMN IF NOT EXISTS', 'double precision', 'geometry(POINTZ, %s)')
        cur.execute(sqltext, (cfg.srid_palm, cfg.srid_palm, ))
        sql_debug(connection)
        connection.commit()

        verbose('Loop over all points')
        for po in range(1, 6): # 6 points is maximum and the last is repeated
            verbose('Point {}', po)
            sqltext = 'ALTER TABLE "{0}"."{1}" ' \
                      '{2} p{5}_meziD {4}, {2} p{5}_meziH {4}, {2} p{5}_f {4}' \
                      ''.format(cfg.domain.case_schema, cfg.tables.slanted_terrain_gridded_temp,
                                'ADD COLUMN IF NOT EXISTS', 'double precision', 'geometry(POINTZ, %s)', po)
            cur.execute(sqltext, (cfg.srid_palm,cfg.srid_palm,cfg.srid_palm,))
            sql_debug(connection)
            connection.commit()

            sqltext = 'UPDATE "{0}"."{1}" SET ' \
                      'p1 = ST_PointN(ST_ExteriorRing(geom3d),{2}), ' \
                      'p2 = ST_PointN(ST_ExteriorRing(geom3d),{2}+1)' \
                .format(cfg.domain.case_schema, cfg.tables.slanted_terrain_gridded_temp, po)
            cur.execute(sqltext)
            sql_debug(connection)
            connection.commit()

            sqltext = 'UPDATE "{0}"."{1}" SET ' \
                      ' mezi_dolni = CASE WHEN ST_Z(p1) < ({2}-{3}) AND (ST_Z(p2) > ({2}-{3})) THEN TRUE' \
                      '                   WHEN ST_Z(p1) > ({2}-{3}) AND (ST_Z(p2) < ({2}-{3})) THEN TRUE' \
                      '                   ELSE FALSE END, ' \
                      ' mezi_horni = CASE WHEN ST_Z(p1) < {2} AND (ST_Z(p2) > {2}) THEN TRUE' \
                      '                   WHEN ST_Z(p1) > {2} AND (ST_Z(p2) < {2}) THEN TRUE ' \
                      '                   ELSE FALSE END, ' \
                      ' z2 = CASE WHEN ST_Z(p2) < ({2}-{3}) THEN NULL' \
                      '           WHEN ST_Z(p2) > {2} THEN {2}' \
                      '           ELSE ST_Z(p2) END '\
                .format(cfg.domain.case_schema, cfg.tables.slanted_terrain_gridded_temp, k*cfg.domain.dz, cfg.domain.dz)
            cur.execute(sqltext)
            sql_debug(connection)
            connection.commit()

            sqltext = 'UPDATE "{0}"."{1}" SET ' \
                      'p{2}_meziD = CASE WHEN mezi_dolni THEN ST_SetSRID(ST_MakePoint(' \
                      '                                       ST_X(p1) + ({3}-{4}-ST_Z(p1))/(ST_Z(p2)-ST_Z(p1))*(ST_X(p2)-ST_X(p1)), ' \
                      '                                       ST_Y(p1) + ({3}-{4}-ST_Z(p1))/(ST_Z(p2)-ST_Z(p1))*(ST_Y(p2)-ST_Y(p1)), ' \
                      '                                       {3}-{4}), %s) ELSE NULL END, ' \
                      'p{2}_meziH = CASE WHEN mezi_horni THEN ST_SetSRID(ST_MakePoint(' \
                      '                                       ST_X(p1) + ({3}-ST_Z(p1))/(ST_Z(p2)-ST_Z(p1))*(ST_X(p2)-ST_X(p1)), ' \
                      '                                       ST_Y(p1) + ({3}-ST_Z(p1))/(ST_Z(p2)-ST_Z(p1))*(ST_Y(p2)-ST_Y(p1)), ' \
                      '                                       {3}), %s) ELSE NULL END, ' \
                      'p{2}_f = CASE WHEN ST_Z(p2) >= ({3}-{4}) AND ST_Z(p2) <= {3} THEN ST_SetSRID(ST_MakePoint(ST_X(p2), ST_Y(p2), z2), %s)' \
                      '              ELSE NULL END'\
                    .format(cfg.domain.case_schema, cfg.tables.slanted_terrain_gridded_temp, po, k*cfg.domain.dz, cfg.domain.dz)

            cur.execute(sqltext, (cfg.srid_palm,cfg.srid_palm,cfg.srid_palm,))
            sql_debug(connection)
            connection.commit()

            # set mezi_dolni, mezi_horni to false
            sqltext = 'UPDATE "{0}"."{1}" SET mezi_dolni = FALSE, mezi_horni = FALSE'.format(cfg.domain.case_schema, cfg.tables.slanted_terrain_gridded_temp)
            cur.execute(sqltext)
            sql_debug(connection)
            connection.commit()

        # Gather point to form polygon
        # find number of points
        debug('Gathering points')
        sqltext = 'UPDATE "{0}"."{1}" SET ' \
                  'npoints = (' \
                  '{2} p1_mezid {3} + {2} p1_mezih {3} + {2} p1_f {3} + ' \
                  '{2} p2_mezid {3} + {2} p2_mezih {3} + {2} p2_f {3} + ' \
                  '{2} p3_mezid {3} + {2} p3_mezih {3} + {2} p3_f {3} + ' \
                  '{2} p4_mezid {3} + {2} p4_mezih {3} + {2} p4_f {3} + ' \
                  '{2} p5_mezid {3} + {2} p5_mezih {3} + {2} p5_f {3} ' \
                  ')'\
            .format(cfg.domain.case_schema, cfg.tables.slanted_terrain_gridded_temp, 'CASE WHEN', 'IS NOT NULL THEN 1 ELSE 0 END ')
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        debug('Insert point into points collection')
        sqltext = 'INSERT INTO "{0}"."{1}" (id, lid, points) ' \
                  'SELECT id, lid, ST_SetSRID(ST_Collect(ARRAY[' \
                  'p1_mezid, p1_mezih, p1_f, p2_mezid, p2_mezih, p2_f, p3_mezid, p3_mezih, p3_f, ' \
                  'p4_mezid, p4_mezih, p4_f, p5_mezid, p5_mezih, p5_f ' \
                  ']), %s) AS points ' \
                  'FROM "{0}"."{2}" ' \
                  'WHERE npoints > 2'.format(cfg.domain.case_schema, cfg.tables.slanted_terrain_gridded, cfg.tables.slanted_terrain_gridded_temp)
        cur.execute(sqltext, (cfg.srid_palm,))
        sql_debug(connection)
        connection.commit()

    debug('Creating polygon geom')
    sqltext = 'UPDATE "{0}"."{1}" SET geom = ST_SetSRID(ST_ConvexHull(points), %s)'\
        .format(cfg.domain.case_schema, cfg.tables.slanted_terrain_gridded)
    cur.execute(sqltext, (cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

    verbose('Updating its SRID')
    sqltext = 'UPDATE "{0}"."{1}" SET lines = ST_SetSRID(ST_Boundary(geom), %s)'\
        .format(cfg.domain.case_schema, cfg.tables.slanted_terrain_gridded)
    cur.execute(sqltext, (cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

    # DELETE DUPLICATES (dont know where they come from)
    verbose('Deleting duplicates')
    sqltext = 'ALTER TABLE "{0}"."{1}" ADD COLUMN IF NOT EXISTS rgid SERIAL'\
              .format(cfg.domain.case_schema, cfg.tables.slanted_terrain_gridded)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    sqltext = 'DELETE FROM "{0}"."{1}" ' \
              'WHERE rgid IN (' \
              '               SELECT rgid FROM (' \
              '                     SELECT rgid, ROW_NUMBER() OVER (partition BY geom ORDER BY ID) AS RowNumber ' \
              '                     FROM "{0}"."{1}") AS T ' \
              '               WHERE T.RowNumber > 1)'\
              .format(cfg.domain.case_schema, cfg.tables.slanted_terrain_gridded)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # MB HERE
    debug('Separate polygon when 3 or more vertices has z = (k+1)*dz')
    sqltext = 'ALTER TABLE "{0}"."{1}" ' \
              'ADD COLUMN IF NOT EXISTS vert1 geometry(POINTZ, %s), ' \
              'ADD COLUMN IF NOT EXISTS vert2 geometry(POINTZ, %s), ' \
              'ADD COLUMN IF NOT EXISTS vert3 geometry(POINTZ, %s), ' \
              'ADD COLUMN IF NOT EXISTS vert4 geometry(POINTZ, %s), ' \
              'ADD COLUMN IF NOT EXISTS vert5 geometry(POINTZ, %s), ' \
              'ADD COLUMN IF NOT EXISTS vert6 geometry(POINTZ, %s), ' \
              'ADD COLUMN IF NOT EXISTS vert7 geometry(POINTZ, %s), ' \
              'ADD COLUMN IF NOT EXISTS z_max double precision, ' \
              'ADD COLUMN IF NOT EXISTS z_min double precision, ' \
              'ADD COLUMN IF NOT EXISTS k_max integer, ' \
              'ADD COLUMN IF NOT EXISTS n_vert integer' \
              .format(cfg.domain.case_schema, cfg.tables.slanted_terrain_gridded)
    cur.execute(sqltext, (cfg.srid_palm,cfg.srid_palm, cfg.srid_palm,cfg.srid_palm,cfg.srid_palm,cfg.srid_palm,cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

    verbose('Create vert points')
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'vert1 = ST_PointN(ST_Boundary(geom), 1), ' \
              'vert2 = ST_PointN(ST_Boundary(geom), 2), ' \
              'vert3 = ST_PointN(ST_Boundary(geom), 3), ' \
              'vert4 = ST_PointN(ST_Boundary(geom), 4), ' \
              'vert5 = ST_PointN(ST_Boundary(geom), 5), ' \
              'vert6 = ST_PointN(ST_Boundary(geom), 6), ' \
              'vert7 = ST_PointN(ST_Boundary(geom), 7)'\
        .format(cfg.domain.case_schema, cfg.tables.slanted_terrain_gridded)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # fetch max height
    verbose('Calculate number of points')
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'n_vert = ST_NPoints(points)+1' \
        .format(cfg.domain.case_schema, cfg.tables.slanted_terrain_gridded)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # fetch max height
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'z_max = ' \
              'GREATEST(' \
              '  ST_Z(ST_GeometryN(points, 1)), ' \
              '  ST_Z(ST_GeometryN(points, 2)), ' \
              '  ST_Z(ST_GeometryN(points, 3)), ' \
              '  ST_Z(ST_GeometryN(points, 4)), ' \
              '  ST_Z(ST_GeometryN(points, 5)), ' \
              '  ST_Z(ST_GeometryN(points, 6)), ' \
              '  ST_Z(ST_GeometryN(points, 7)) ' \
              '),' \
              'z_min = ' \
              'LEAST(' \
              '  ST_Z(ST_GeometryN(points, 1)), ' \
              '  ST_Z(ST_GeometryN(points, 2)), ' \
              '  ST_Z(ST_GeometryN(points, 3)), ' \
              '  ST_Z(ST_GeometryN(points, 4)), ' \
              '  ST_Z(ST_GeometryN(points, 5)), ' \
              '  ST_Z(ST_GeometryN(points, 6)), ' \
              '  ST_Z(ST_GeometryN(points, 7)) ' \
              ') ' \
        .format(cfg.domain.case_schema, cfg.tables.slanted_terrain_gridded)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    verbose('Selecting polygons')
    verbose('Updating all planar-horizontal faces, update their k')
    verbose('Fetching max k from slanted faces')
    sqltext = 'SELECT MAX(z_max) FROM "{0}"."{1}"'.format(cfg.domain.case_schema, cfg.tables.slanted_terrain_gridded)
    cur.execute(sqltext)
    z_max = cur.fetchone()[0]
    sql_debug(connection)
    connection.commit()
    k_max = int(np.floor(z_max) + 2)

    for k in range(0, k_max):
        sqltext = 'SELECT id, rgid, lid, n_vert, ' \
                  ' ARRAY[ST_X(vert1), ST_Y(vert1), ST_Z(vert1)],  ' \
                  ' ARRAY[ST_X(vert2), ST_Y(vert2), ST_Z(vert2)], ' \
                  ' ARRAY[ST_X(vert3), ST_Y(vert3), ST_Z(vert3)], ' \
                  ' ARRAY[ST_X(vert4), ST_Y(vert4), ST_Z(vert4)], ' \
                  ' ARRAY[ST_X(vert5), ST_Y(vert5), ST_Z(vert5)], ' \
                  ' ARRAY[ST_X(vert6), ST_Y(vert6), ST_Z(vert6)], ' \
                  ' ARRAY[ST_X(vert7), ST_Y(vert7), ST_Z(vert7)]  ' \
                  'FROM "{0}"."{1}" AS s ' \
                  'WHERE (CASE WHEN ST_Z(vert1) = {3} THEN 1 ELSE 0 END + ' \
                  '       CASE WHEN ST_Z(vert2) = {3} THEN 1 ELSE 0 END + ' \
                  '       CASE WHEN ST_Z(vert3) = {3} THEN 1 ELSE 0 END + ' \
                  '       CASE WHEN ST_Z(vert4) = {3} AND n_vert > 4 THEN 1 ELSE 0 END + ' \
                  '       CASE WHEN ST_Z(vert5) = {3} AND n_vert > 5 THEN 1 ELSE 0 END + ' \
                  '       CASE WHEN ST_Z(vert6) = {3} AND n_vert > 6 THEN 1 ELSE 0 END + ' \
                  '       CASE WHEN ST_Z(vert7) = {3} AND n_vert > 7 THEN 1 ELSE 0 END) > 2 ' \
                  '      AND n_vert > 3 AND z_min != z_max ' \
                  ''\
                  .format(cfg.domain.case_schema, cfg.tables.slanted_terrain_gridded, k, (k+1) * cfg.domain.dz)
        cur.execute(sqltext)
        verts = cur.fetchall()
        sql_debug(connection)
        connection.commit()

        # TODO: Now here separate this polygon into two. One with 3 or more top points and the rest
        to_insert = []
        to_delete = []
        for vert in verts:
            p1 = []
            p2 = []
            np_verts_p1 = np.empty((7,3), dtype=object)
            np_verts_p2 = np.empty((7,3), dtype=object)
            id, rgid, lid, n_vert = vert[0], vert[1], vert[2], vert[3]
            lastidx = 4
            x_vert, y_vert, z_vert = [], [], []
            for idx in range(n_vert):
                if vert[lastidx + idx][0] is not None:
                    x_vert.append(vert[lastidx + idx][0])
                    y_vert.append(vert[lastidx + idx][1])
                    z_vert.append(vert[lastidx + idx][2])
            x_vert = np.asarray(x_vert)
            y_vert = np.asarray(y_vert)
            z_vert = np.asarray(z_vert)
            to_keep = []
            tt_max = max(z_vert)
            for it, ttt in enumerate(z_vert):
                # print(ttt)
                if it == 0:
                    ileft = len(z_vert) - 1
                    iright = it + 1
                elif it == len(z_vert) - 1:
                    iright = 0
                    ileft = it - 1
                else:
                    iright = it + 1
                    ileft = it - 1
                if z_vert[ileft] == tt_max and z_vert[it] == tt_max and z_vert[iright] == tt_max:
                    # print('drop this one')
                    p2.append(it)
                else:
                    if z_vert[it] == tt_max:
                        p1.append(it)
                        p2.append(it)
                    else:
                        p1.append(it)

            if len(p1) > 2:
                n_vert1 = len(p1) + 1
                np_verts_p1[:n_vert1 - 1, 0] = x_vert[p1]
                np_verts_p1[:n_vert1 - 1, 1] = y_vert[p1]
                np_verts_p1[:n_vert1 - 1, 2] = z_vert[p1]
                np_verts_p1[n_vert1 - 1, :] = np_verts_p1[0, :]

                to_insert.append((rgid, lid, n_vert1,
                                  np_verts_p1[0, 0], np_verts_p1[0, 1], np_verts_p1[0, 2], cfg.srid_palm,
                                  np_verts_p1[1, 0], np_verts_p1[1, 1], np_verts_p1[1, 2], cfg.srid_palm,
                                  np_verts_p1[2, 0], np_verts_p1[2, 1], np_verts_p1[2, 2], cfg.srid_palm,
                                  np_verts_p1[3, 0], np_verts_p1[3, 1], np_verts_p1[3, 2], cfg.srid_palm,
                                  np_verts_p1[4, 0], np_verts_p1[4, 1], np_verts_p1[4, 2], cfg.srid_palm,
                                  np_verts_p1[5, 0], np_verts_p1[5, 1], np_verts_p1[5, 2], cfg.srid_palm,
                                  np_verts_p1[6, 0], np_verts_p1[6, 1], np_verts_p1[6, 2], cfg.srid_palm,
                                  )
                                  )
            if len(p2) > 2:
                n_vert2 = len(p2) + 1
                np_verts_p2[:n_vert2 - 1, 0] = x_vert[p2]
                np_verts_p2[:n_vert2 - 1, 1] = y_vert[p2]
                np_verts_p2[:n_vert2 - 1, 2] = z_vert[p2]
                np_verts_p2[n_vert2 - 1, :] = np_verts_p2[0, :]
                to_insert.append((rgid, lid, n_vert2,
                                  np_verts_p2[0, 0], np_verts_p2[0, 1], np_verts_p2[0, 2], cfg.srid_palm,
                                  np_verts_p2[1, 0], np_verts_p2[1, 1], np_verts_p2[1, 2], cfg.srid_palm,
                                  np_verts_p2[2, 0], np_verts_p2[2, 1], np_verts_p2[2, 2], cfg.srid_palm,
                                  np_verts_p2[3, 0], np_verts_p2[3, 1], np_verts_p2[3, 2], cfg.srid_palm,
                                  np_verts_p2[4, 0], np_verts_p2[4, 1], np_verts_p2[4, 2], cfg.srid_palm,
                                  np_verts_p2[5, 0], np_verts_p2[5, 1], np_verts_p2[5, 2], cfg.srid_palm,
                                  np_verts_p2[6, 0], np_verts_p2[6, 1], np_verts_p2[6, 2], cfg.srid_palm,
                                  )
                                  )
            if len(p2) > 2 or len(p1) > 2:
                to_delete.append((rgid,))

        debug('Deleting all unwanted rows')
        sqltext = 'DELETE FROM "{0}"."{1}" ' \
                  'WHERE rgid = ANY(%s)'.format(cfg.domain.case_schema, cfg.tables.slanted_terrain_gridded)
        cur.execute(sqltext, (to_delete, ))
        sql_debug(connection)
        connection.commit()

        debug('Inserting all new entries into slanted faces table')
        sqltext = 'INSERT INTO "{0}"."{1}" (id, rgid, lid, n_vert, ' \
                  '                         vert1, vert2, vert3, vert4, vert5, vert6, vert7) ' \
                  'VALUES                  (NULL, %s, %s, %s, ' \
                  '                        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
                  '                        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
                  '                        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
                  '                        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
                  '                        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
                  '                        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
                  '                        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s) ' \
                  '                        )  '.format(cfg.domain.case_schema, cfg.tables.slanted_terrain_gridded)

        cur.executemany(sqltext, to_insert)
        sql_debug(connection)
        connection.commit()

        # Create polygon
        debug('Updating 3d polygon')
        sqltext = 'UPDATE "{0}"."{1}" SET geom =  ' \
                  'ST_ForceRHR(' \
                  'ST_SetSRID(ST_MakePolygon(ST_MakeLine(ARRAY[vert1, vert2, vert3, vert4, vert5, vert6, vert7,' \
                  '     CASE WHEN vert1 IS NOT NULL THEN vert1 ' \
                  '          WHEN vert2 IS NOT NULL THEN vert2 ' \
                  '          WHEN vert3 IS NOT NULL THEN vert3 ' \
                  '          WHEN vert4 IS NOT NULL THEN vert4 ' \
                  '          WHEN vert5 IS NOT NULL THEN vert5 ' \
                  '          WHEN vert6 IS NOT NULL THEN vert6 ' \
                  '          WHEN vert7 IS NOT NULL THEN vert7 END' \
                  '])), %s)) ' \
                  'WHERE geom IS NULL'.format(cfg.domain.case_schema, cfg.tables.slanted_terrain_gridded)
        cur.execute(sqltext, (cfg.srid_palm,))
        sql_debug(connection)
        connection.commit()



    # TODO: Create new serial index
    sqltext = 'ALTER TABLE "{0}"."{1}" DROP COLUMN IF EXISTS id; ' \
              'ALTER TABLE "{0}"."{1}" ADD COLUMN id SERIAL'.format(cfg.domain.case_schema, cfg.tables.slanted_terrain_gridded)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    sqltext = 'ALTER TABLE "{0}"."{1}" ADD PRIMARY KEY (id)' \
        .format(cfg.domain.case_schema, cfg.tables.slanted_terrain_gridded)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # sqltext = 'ALTER TABLE "{0}"."{1}" DROP COLUMN IF EXISTS i; ' \
    #           'ALTER TABLE "{0}"."{1}" DROP COLUMN IF EXISTS j;' \
    #           'ALTER TABLE "{0}"."{1}" DROP COLUMN IF EXISTS k;'.format(cfg.domain.case_schema, cfg.tables.slanted_terrain_gridded)
    # cur.execute(sqltext)
    # sql_debug(connection)
    # connection.commit()

    # sqltext = 'ALTER TABLE "{0}"."{1}" ADD COLUMN i integer,' \
    #           'ALTER TABLE "{0}"."{1}" ADD COLUMN j integer,' \
    #           'ALTER TABLE "{0}"."{1}" ADD COLUMN k integer,' \
    #           'ALTER TABLE "{0}"."{1}" ADD COLUMN center geometry(POINTZ, %s)'.format(cfg.domain.case_schema, cfg.tables.slanted_terrain_gridded)
    # cur.execute(sqltext, (cfg.srid_palm, ))
    # sql_debug(connection)
    # connection.commit()
    #
    # debug('Updating faces centers')
    # sqltext = 'UPDATE "{0}"."{1}" SET ' \
    #           'center = ST_SetSRID(ST_MakePoint(ST_X(ST_Centroid(geom)), ST_Y(ST_Centroid(geom)), ' \
    #           '     (COALESCE(ST_Z(ST_PointN(ST_ExteriorRing(geom),1)), 0.0) + COALESCE(ST_Z(ST_PointN(ST_ExteriorRing(geom),2)), 0.0) + ' \
    #           '      COALESCE(ST_Z(ST_PointN(ST_ExteriorRing(geom),3)), 0.0) + COALESCE(ST_Z(ST_PointN(ST_ExteriorRing(geom),4)), 0.0) + ' \
    #           '      COALESCE(ST_Z(ST_PointN(ST_ExteriorRing(geom),5)), 0.0) + COALESCE(ST_Z(ST_PointN(ST_ExteriorRing(geom),6)), 0.0) + ' \
    #           '      COALESCE(ST_Z(ST_PointN(ST_ExteriorRing(geom),7)), 0.0)) / ST_NPoints(geom))' \
    #           '         , %s)'\
    #     .format(cfg.domain.case_schema, cfg.tables.slanted_terrain_gridded)
    # cur.execute(sqltext, (cfg.srid_palm, ))
    # sql_debug(connection)
    # connection.commit()
    #
    # debug('Calculating face i,j,k')
    # sqltext = 'UPDATE "{0}"."{1}" SET ' \
    #           'i = FLOOR((ST_X(center) - {3}) / {2}) , ' \
    #           'j = FLOOR((ST_Y(center) - {4}) / {2}), ' \
    #           'k = FLOOR(ST_Z(center) / {2})' \
    #           ''.format(cfg.domain.case_schema, cfg.tables.slanted_terrain_gridded, cfg.domain.dx,
    #                     cfg.domain.origin_x, cfg.domain.origin_y)
    # cur.execute(sqltext)
    # sql_debug(connection)
    # connection.commit()



    if cfg.slanted_pars.clean_up:
        sqltext = 'ALTER TABLE "{0}"."{1}" DROP COLUMN rgid, ' \
                  '                        DROP COLUMN rijk'\
                  .format(cfg.domain.case_schema, cfg.tables.slanted_terrain_gridded)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

def create_grid_slanted_roof(cfg, connection, cur):
    """ Processing gridding slanted roof """
    # TODO: Add debug, verbose reports
    progress('Processing, gridding slanted roof')

    sqltext = 'DROP TABLE IF EXISTS "{0}"."{1}"; ' \
              'CREATE TABLE "{0}"."{1}" (' \
              'id integer, rid integer,  ' \
              'geom geometry("POLYGONZ", %s), ' \
              'points geometry("MULTIPOINTZ", %s),' \
              'lines geometry("LINESTRINGZ", %s)' \
              ')'.format(cfg.domain.case_schema, cfg.tables.slanted_roof_gridded)
    cur.execute(sqltext, (cfg.srid_palm,cfg.srid_palm,cfg.srid_palm, ))
    sql_debug(connection)
    connection.commit()

    # CREATE TEMP TABLE, destroy after cycle
    # FIND MAX HEIGHT -> k_max
    sqltext = 'SELECT MAX(height) FROM "{0}"."{1}"'\
              .format(cfg.domain.case_schema, cfg.tables.height_corrected)
    cur.execute(sqltext)
    z_max = cur.fetchall()
    z_max = [x for x in z_max[0]]
    k_max = int(np.ceil(max(z_max)) / cfg.domain.dz) + 2
    sql_debug(connection)
    connection.commit()

    for k in range(1, k_max):
        verbose('Gridding slanted roof: {}', k)
        # Create temp table
        sqltext = 'DROP TABLE IF EXISTS "{0}"."{1}"'.format(cfg.domain.case_schema, cfg.tables.slanted_roof_gridded_temp)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        sqltext = 'CREATE TABLE "{0}"."{1}" AS SELECT id, rid, geom3d ' \
                  'FROM "{0}"."{2}"'\
                  .format(cfg.domain.case_schema, cfg.tables.slanted_roof_gridded_temp, cfg.tables.slanted_roof)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        sqltext = 'ALTER TABLE "{0}"."{1}" ' \
                  '{2} npoints integer, ' \
                  '{2} z1  {3}, {2} z2  {3}, ' \
                  '{2} mezi_dolni boolean DEFAULT FALSE, {2} mezi_horni boolean DEFAULT FALSE, ' \
                  '{2} p1 {4}, {2} p2 {4}' \
                  ''.format(cfg.domain.case_schema, cfg.tables.slanted_roof_gridded_temp,
                            'ADD COLUMN IF NOT EXISTS', 'double precision', 'geometry(POINTZ, %s)')
        cur.execute(sqltext, (cfg.srid_palm, cfg.srid_palm, ))
        sql_debug(connection)
        connection.commit()

        for po in range(1, 6): # 6 points is maximum and the last is repeated
            sqltext = 'ALTER TABLE "{0}"."{1}" ' \
                      '{2} p{5}_meziD {4}, {2} p{5}_meziH {4}, {2} p{5}_f {4}' \
                      ''.format(cfg.domain.case_schema, cfg.tables.slanted_roof_gridded_temp,
                                'ADD COLUMN IF NOT EXISTS', 'double precision', 'geometry(POINTZ, %s)', po)
            cur.execute(sqltext, (cfg.srid_palm,cfg.srid_palm,cfg.srid_palm,))
            sql_debug(connection)
            connection.commit()

            sqltext = 'UPDATE "{0}"."{1}" SET ' \
                      'p1 = ST_PointN(ST_ExteriorRing(geom3d),{2}), ' \
                      'p2 = ST_PointN(ST_ExteriorRing(geom3d),{2}+1)' \
                .format(cfg.domain.case_schema, cfg.tables.slanted_roof_gridded_temp, po)
            cur.execute(sqltext)
            sql_debug(connection)
            connection.commit()

            sqltext = 'UPDATE "{0}"."{1}" SET ' \
                      ' mezi_dolni = CASE WHEN ST_Z(p1) < ({2}-{3}) AND (ST_Z(p2) > ({2}-{3})) THEN TRUE' \
                      '                   WHEN ST_Z(p1) > ({2}-{3}) AND (ST_Z(p2) < ({2}-{3})) THEN TRUE' \
                      '                   ELSE FALSE END, ' \
                      ' mezi_horni = CASE WHEN ST_Z(p1) < {2} AND (ST_Z(p2) > {2}) THEN TRUE' \
                      '                   WHEN ST_Z(p1) > {2} AND (ST_Z(p2) < {2}) THEN TRUE ' \
                      '                   ELSE FALSE END, ' \
                      ' z2 = CASE WHEN ST_Z(p2) < ({2}-{3}) THEN NULL' \
                      '           WHEN ST_Z(p2) > {2} THEN {2}' \
                      '           ELSE ST_Z(p2) END '\
                .format(cfg.domain.case_schema, cfg.tables.slanted_roof_gridded_temp, k*cfg.domain.dz, cfg.domain.dz)
            cur.execute(sqltext)
            sql_debug(connection)
            connection.commit()

            sqltext = 'UPDATE "{0}"."{1}" SET ' \
                      'p{2}_meziD = CASE WHEN mezi_dolni THEN ST_SetSRID(ST_MakePoint(' \
                      '                                       ST_X(p1) + ({3}-{4}-ST_Z(p1))/(ST_Z(p2)-ST_Z(p1))*(ST_X(p2)-ST_X(p1)), ' \
                      '                                       ST_Y(p1) + ({3}-{4}-ST_Z(p1))/(ST_Z(p2)-ST_Z(p1))*(ST_Y(p2)-ST_Y(p1)), ' \
                      '                                       {3}-{4}), %s) ELSE NULL END, ' \
                      'p{2}_meziH = CASE WHEN mezi_horni THEN ST_SetSRID(ST_MakePoint(' \
                      '                                       ST_X(p1) + ({3}-ST_Z(p1))/(ST_Z(p2)-ST_Z(p1))*(ST_X(p2)-ST_X(p1)), ' \
                      '                                       ST_Y(p1) + ({3}-ST_Z(p1))/(ST_Z(p2)-ST_Z(p1))*(ST_Y(p2)-ST_Y(p1)), ' \
                      '                                       {3}), %s) ELSE NULL END, ' \
                      'p{2}_f = CASE WHEN ST_Z(p2) >= ({3}-{4}) AND ST_Z(p2) <= {3} THEN ST_SetSRID(ST_MakePoint(ST_X(p2), ST_Y(p2), z2), %s)' \
                      '              ELSE NULL END'\
                      .format(cfg.domain.case_schema, cfg.tables.slanted_roof_gridded_temp, po, k*cfg.domain.dz, cfg.domain.dz)

            cur.execute(sqltext, (cfg.srid_palm,cfg.srid_palm,cfg.srid_palm,))
            sql_debug(connection)
            connection.commit()

            # set mezi_dolni, mezi_horni to false
            sqltext = 'UPDATE "{0}"."{1}" SET mezi_dolni = FALSE, mezi_horni = FALSE'.format(cfg.domain.case_schema, cfg.tables.slanted_roof_gridded_temp)
            cur.execute(sqltext)
            sql_debug(connection)
            connection.commit()

        # Gather point to form polygon
        # find number of points
        sqltext = 'UPDATE "{0}"."{1}" SET ' \
                  'npoints = (' \
                  '{2} p1_mezid {3} + {2} p1_mezih {3} + {2} p1_f {3} + ' \
                  '{2} p2_mezid {3} + {2} p2_mezih {3} + {2} p2_f {3} + ' \
                  '{2} p3_mezid {3} + {2} p3_mezih {3} + {2} p3_f {3} + ' \
                  '{2} p4_mezid {3} + {2} p4_mezih {3} + {2} p4_f {3} + ' \
                  '{2} p5_mezid {3} + {2} p5_mezih {3} + {2} p5_f {3} ' \
                  ')'\
            .format(cfg.domain.case_schema, cfg.tables.slanted_roof_gridded_temp, 'CASE WHEN', 'IS NOT NULL THEN 1 ELSE 0 END ')
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        sqltext = 'INSERT INTO "{0}"."{1}" (id, rid, points)' \
                  'SELECT id, rid, ST_SetSRID(ST_Collect(ARRAY[' \
                  'p1_mezid, p1_mezih, p1_f, p2_mezid, p2_mezih, p2_f, p3_mezid, p3_mezih, p3_f, ' \
                  'p4_mezid, p4_mezih, p4_f, p5_mezid, p5_mezih, p5_f' \
                  ']), %s) AS points ' \
                  'FROM "{0}"."{2}" ' \
                  'WHERE npoints > 2'.format(cfg.domain.case_schema, cfg.tables.slanted_roof_gridded, cfg.tables.slanted_roof_gridded_temp)
        cur.execute(sqltext, (cfg.srid_palm,))
        sql_debug(connection)
        connection.commit()

    sqltext = 'UPDATE "{0}"."{1}" SET geom = ST_SetSRID(ST_ConvexHull(points), %s)'\
        .format(cfg.domain.case_schema, cfg.tables.slanted_roof_gridded)
    cur.execute(sqltext, (cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

    sqltext = 'UPDATE "{0}"."{1}" SET lines = ST_SetSRID(ST_Boundary(geom), %s)'\
        .format(cfg.domain.case_schema, cfg.tables.slanted_roof_gridded)
    cur.execute(sqltext, (cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

    # DELETE DUPLICATES (dont know where they come from)
    verbose('Deleting duplicates')
    sqltext = 'ALTER TABLE "{0}"."{1}" ADD COLUMN IF NOT EXISTS rgid SERIAL'\
              .format(cfg.domain.case_schema, cfg.tables.slanted_roof_gridded)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    sqltext = 'DELETE FROM "{0}"."{1}" ' \
              'WHERE rgid IN (' \
              '               SELECT rgid FROM (' \
              '                     SELECT rgid, ROW_NUMBER() OVER (partition BY geom ORDER BY ID) AS RowNumber ' \
              '                     FROM "{0}"."{1}") AS T ' \
              '               WHERE T.RowNumber > 1)'\
              .format(cfg.domain.case_schema, cfg.tables.slanted_roof_gridded)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    debug('Separate polygon when 3 or more vertices has z = (k+1)*dz or z = k * dz')
    sqltext = 'ALTER TABLE "{0}"."{1}" ' \
              'ADD COLUMN IF NOT EXISTS vert1 geometry(POINTZ, %s), ' \
              'ADD COLUMN IF NOT EXISTS vert2 geometry(POINTZ, %s), ' \
              'ADD COLUMN IF NOT EXISTS vert3 geometry(POINTZ, %s), ' \
              'ADD COLUMN IF NOT EXISTS vert4 geometry(POINTZ, %s), ' \
              'ADD COLUMN IF NOT EXISTS vert5 geometry(POINTZ, %s), ' \
              'ADD COLUMN IF NOT EXISTS vert6 geometry(POINTZ, %s), ' \
              'ADD COLUMN IF NOT EXISTS vert7 geometry(POINTZ, %s), ' \
              'ADD COLUMN IF NOT EXISTS z_max double precision, ' \
              'ADD COLUMN IF NOT EXISTS z_min double precision, ' \
              'ADD COLUMN IF NOT EXISTS k_max integer, ' \
              'ADD COLUMN IF NOT EXISTS n_vert integer' \
              .format(cfg.domain.case_schema, cfg.tables.slanted_roof_gridded)
    cur.execute(sqltext, (cfg.srid_palm, cfg.srid_palm, cfg.srid_palm,cfg.srid_palm,cfg.srid_palm,cfg.srid_palm,cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

    verbose('Create vert points')
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'vert1 = ST_PointN(ST_Boundary(geom), 1), ' \
              'vert2 = ST_PointN(ST_Boundary(geom), 2), ' \
              'vert3 = ST_PointN(ST_Boundary(geom), 3), ' \
              'vert4 = ST_PointN(ST_Boundary(geom), 4), ' \
              'vert5 = ST_PointN(ST_Boundary(geom), 5), ' \
              'vert6 = ST_PointN(ST_Boundary(geom), 6), ' \
              'vert7 = ST_PointN(ST_Boundary(geom), 7)'\
        .format(cfg.domain.case_schema, cfg.tables.slanted_roof_gridded)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    verbose('Calculate number of points')
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'n_vert = ST_NPoints(points)+1' \
        .format(cfg.domain.case_schema, cfg.tables.slanted_roof_gridded)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # fetch max height
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'z_max = ' \
              'GREATEST(' \
              '  ST_Z(ST_GeometryN(points, 1)), ' \
              '  ST_Z(ST_GeometryN(points, 2)), ' \
              '  ST_Z(ST_GeometryN(points, 3)), ' \
              '  ST_Z(ST_GeometryN(points, 4)), ' \
              '  ST_Z(ST_GeometryN(points, 5)), ' \
              '  ST_Z(ST_GeometryN(points, 6)), ' \
              '  ST_Z(ST_GeometryN(points, 7)) ' \
              '),' \
              'z_min = ' \
              'LEAST(' \
              '  ST_Z(ST_GeometryN(points, 1)), ' \
              '  ST_Z(ST_GeometryN(points, 2)), ' \
              '  ST_Z(ST_GeometryN(points, 3)), ' \
              '  ST_Z(ST_GeometryN(points, 4)), ' \
              '  ST_Z(ST_GeometryN(points, 5)), ' \
              '  ST_Z(ST_GeometryN(points, 6)), ' \
              '  ST_Z(ST_GeometryN(points, 7)) ' \
              ') ' \
        .format(cfg.domain.case_schema, cfg.tables.slanted_roof_gridded)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    verbose('Selecting polygons')
    verbose('Updating all planar-horizontal faces, update their k')
    verbose('Fetching max k from slanted faces')
    sqltext = 'SELECT MAX(z_max) FROM "{0}"."{1}"'.format(cfg.domain.case_schema, cfg.tables.slanted_roof_gridded)
    cur.execute(sqltext)
    z_max = cur.fetchone()[0]
    sql_debug(connection)
    connection.commit()
    k_max = int(np.floor(z_max) + 2)
    # kk = 0 check z = k * dz, kk = 1 check z = (k+1) * dz
    for kk in [1]:
        for k in range(0, k_max):
            extra_verbose('loop over k: {}', k)
            sqltext = 'SELECT id, rgid, rid, n_vert, ' \
                      ' ARRAY[ST_X(vert1), ST_Y(vert1), ST_Z(vert1)],  ' \
                      ' ARRAY[ST_X(vert2), ST_Y(vert2), ST_Z(vert2)], ' \
                      ' ARRAY[ST_X(vert3), ST_Y(vert3), ST_Z(vert3)], ' \
                      ' ARRAY[ST_X(vert4), ST_Y(vert4), ST_Z(vert4)], ' \
                      ' ARRAY[ST_X(vert5), ST_Y(vert5), ST_Z(vert5)], ' \
                      ' ARRAY[ST_X(vert6), ST_Y(vert6), ST_Z(vert6)], ' \
                      ' ARRAY[ST_X(vert7), ST_Y(vert7), ST_Z(vert7)]  ' \
                      'FROM "{0}"."{1}" AS s ' \
                      'WHERE (CASE WHEN ST_Z(vert1) = {3} THEN 1 ELSE 0 END + ' \
                      '       CASE WHEN ST_Z(vert2) = {3} THEN 1 ELSE 0 END + ' \
                      '       CASE WHEN ST_Z(vert3) = {3} THEN 1 ELSE 0 END + ' \
                      '       CASE WHEN ST_Z(vert4) = {3} AND n_vert > 4 THEN 1 ELSE 0 END + ' \
                      '       CASE WHEN ST_Z(vert5) = {3} AND n_vert > 5 THEN 1 ELSE 0 END + ' \
                      '       CASE WHEN ST_Z(vert6) = {3} AND n_vert > 6 THEN 1 ELSE 0 END + ' \
                      '       CASE WHEN ST_Z(vert7) = {3} AND n_vert > 7 THEN 1 ELSE 0 END) > 2 ' \
                      '      AND n_vert > 3 AND z_min != z_max ' \
                      ''\
                      .format(cfg.domain.case_schema, cfg.tables.slanted_roof_gridded, k, (k+kk) * cfg.domain.dz)
            cur.execute(sqltext)
            verts = cur.fetchall()
            sql_debug(connection)
            connection.commit()

            # TODO: Now here separate this polygon into two. One with 3 or more top points and the rest
            to_insert = []
            to_delete = []
            for vert in verts:
                p1 = []
                p2 = []
                np_verts_p1 = np.empty((7,3), dtype=object)
                np_verts_p2 = np.empty((7,3), dtype=object)
                id, rgid, rid, n_vert = vert[0], vert[1], vert[2], vert[3]
                lastidx = 4
                x_vert, y_vert, z_vert = [], [], []
                for idx in range(n_vert):
                    if vert[lastidx + idx][0] is not None:
                        x_vert.append(vert[lastidx + idx][0])
                        y_vert.append(vert[lastidx + idx][1])
                        z_vert.append(vert[lastidx + idx][2])
                x_vert = np.asarray(x_vert)
                y_vert = np.asarray(y_vert)
                z_vert = np.asarray(z_vert)
                to_keep = []
                tt_max = max(z_vert) if kk == 1 else min(z_vert)
                for it, ttt in enumerate(z_vert):
                    # print(ttt)
                    if it == 0:
                        ileft = len(z_vert) - 1
                        iright = it + 1
                    elif it == len(z_vert) - 1:
                        iright = 0
                        ileft = it - 1
                    else:
                        iright = it + 1
                        ileft = it - 1
                    if z_vert[ileft] == tt_max and z_vert[it] == tt_max and z_vert[iright] == tt_max:
                        # print('drop this one')
                        p2.append(it)
                    else:
                        if z_vert[it] == tt_max:
                            p1.append(it)
                            p2.append(it)
                        else:
                            p1.append(it)

                if len(p1) > 2:
                    n_vert1 = len(p1)
                    np_verts_p1[:n_vert1, 0] = x_vert[p1]
                    np_verts_p1[:n_vert1, 1] = y_vert[p1]
                    np_verts_p1[:n_vert1, 2] = z_vert[p1]
                    if not (x_vert[0] == x_vert[n_vert1-1] and y_vert[0] == y_vert[n_vert1-1] and z_vert[0] == z_vert[n_vert1-1]):
                        n_vert1 = len(p1) + 1
                        np_verts_p1[n_vert1 - 1, :] = np_verts_p1[0, :]

                    to_insert.append((rgid, rid, n_vert1,
                                      np_verts_p1[0, 0], np_verts_p1[0, 1], np_verts_p1[0, 2], cfg.srid_palm,
                                      np_verts_p1[1, 0], np_verts_p1[1, 1], np_verts_p1[1, 2], cfg.srid_palm,
                                      np_verts_p1[2, 0], np_verts_p1[2, 1], np_verts_p1[2, 2], cfg.srid_palm,
                                      np_verts_p1[3, 0], np_verts_p1[3, 1], np_verts_p1[3, 2], cfg.srid_palm,
                                      np_verts_p1[4, 0], np_verts_p1[4, 1], np_verts_p1[4, 2], cfg.srid_palm,
                                      np_verts_p1[5, 0], np_verts_p1[5, 1], np_verts_p1[5, 2], cfg.srid_palm,
                                      np_verts_p1[6, 0], np_verts_p1[6, 1], np_verts_p1[6, 2], cfg.srid_palm,
                                      )
                                      )
                if len(p2) > 2:
                    n_vert2 = len(p2) + 1
                    np_verts_p2[:n_vert2 - 1, 0] = x_vert[p2]
                    np_verts_p2[:n_vert2 - 1, 1] = y_vert[p2]
                    np_verts_p2[:n_vert2 - 1, 2] = z_vert[p2]
                    np_verts_p2[n_vert2 - 1, :] = np_verts_p2[0, :]
                    to_insert.append((rgid, rid, n_vert2,
                                      np_verts_p2[0, 0], np_verts_p2[0, 1], np_verts_p2[0, 2], cfg.srid_palm,
                                      np_verts_p2[1, 0], np_verts_p2[1, 1], np_verts_p2[1, 2], cfg.srid_palm,
                                      np_verts_p2[2, 0], np_verts_p2[2, 1], np_verts_p2[2, 2], cfg.srid_palm,
                                      np_verts_p2[3, 0], np_verts_p2[3, 1], np_verts_p2[3, 2], cfg.srid_palm,
                                      np_verts_p2[4, 0], np_verts_p2[4, 1], np_verts_p2[4, 2], cfg.srid_palm,
                                      np_verts_p2[5, 0], np_verts_p2[5, 1], np_verts_p2[5, 2], cfg.srid_palm,
                                      np_verts_p2[6, 0], np_verts_p2[6, 1], np_verts_p2[6, 2], cfg.srid_palm,
                                      )
                                      )
                if len(p2) > 2 or len(p1) > 2:
                    to_delete.append((rgid,))

            debug('Deleting all unwanted rows')
            sqltext = 'DELETE FROM "{0}"."{1}" ' \
                      'WHERE rgid = ANY(%s)'.format(cfg.domain.case_schema, cfg.tables.slanted_roof_gridded)
            cur.execute(sqltext, (to_delete, ))
            sql_debug(connection)
            connection.commit()

            debug('Inserting all new entries into slanted faces table')
            sqltext = 'INSERT INTO "{0}"."{1}" (id, rgid, rid, n_vert, ' \
                      '                         vert1, vert2, vert3, vert4, vert5, vert6, vert7) ' \
                      'VALUES                  (NULL, %s, %s, %s, ' \
                      '                        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
                      '                        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
                      '                        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
                      '                        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
                      '                        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
                      '                        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
                      '                        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s) ' \
                      '                        )  '.format(cfg.domain.case_schema, cfg.tables.slanted_roof_gridded)

            cur.executemany(sqltext, to_insert)
            sql_debug(connection)
            connection.commit()

            # Create polygon
            debug('Updating 3d polygon')
            sqltext = 'UPDATE "{0}"."{1}" SET geom =  ' \
                      'ST_ForceRHR(' \
                      'ST_SetSRID(ST_MakePolygon(ST_MakeLine(ARRAY[vert1, vert2, vert3, vert4, vert5, vert6, vert7,' \
                      '     CASE WHEN vert1 IS NOT NULL THEN vert1 ' \
                      '          WHEN vert2 IS NOT NULL THEN vert2 ' \
                      '          WHEN vert3 IS NOT NULL THEN vert3 ' \
                      '          WHEN vert4 IS NOT NULL THEN vert4 ' \
                      '          WHEN vert5 IS NOT NULL THEN vert5 ' \
                      '          WHEN vert6 IS NOT NULL THEN vert6 ' \
                      '          WHEN vert7 IS NOT NULL THEN vert7 END' \
                      '])), %s)) ' \
                      'WHERE geom IS NULL'.format(cfg.domain.case_schema, cfg.tables.slanted_roof_gridded)
            cur.execute(sqltext, (cfg.srid_palm,))
            sql_debug(connection)
            connection.commit()



    # TODO: Create new serial index
    sqltext = 'ALTER TABLE "{0}"."{1}" DROP COLUMN IF EXISTS id; ' \
              'ALTER TABLE "{0}"."{1}" ADD COLUMN id SERIAL'.format(cfg.domain.case_schema, cfg.tables.slanted_roof_gridded)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    sqltext = 'ALTER TABLE "{0}"."{1}" ADD PRIMARY KEY (id)' \
        .format(cfg.domain.case_schema, cfg.tables.slanted_roof_gridded)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    if cfg.slanted_pars.clean_up:
        sqltext = 'ALTER TABLE "{0}"."{1}" DROP COLUMN rgid'\
                  .format(cfg.domain.case_schema, cfg.tables.slanted_roof_gridded)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

def merge_walls_terrain(cfg, connection, cur):
    """ Merging roof edges with walls"""
    progress('Merging grid faces of terrain and walls')

    debug('Finding duplicates')
    # find all aggregated faces, create new table with aggregated ones
    sqltext = 'SELECT i,j,k FROM "{0}"."{1}" ' \
              'WHERE NOT isroof ' \
              'GROUP BY i,j,k ' \
              'HAVING COUNT(*)>1' \
              'ORDER BY i,j,k'.format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext)
    duplicits = cur.fetchall()
    sql_debug(connection)
    connection.commit()
    i_all = [x[0] for x in duplicits]
    j_all = [x[1] for x in duplicits]
    k_all = [x[2] for x in duplicits]

    # SELECT ALL wid
    debug('Finding all wids')
    sqltext = 'SELECT wid FROM "{0}"."{1}" WHERE iswall AND ('\
              .format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    sqltext_ijk = '('
    for i,j,k in zip(i_all, j_all, k_all):
    # for i, j, k in zip(i_all[:2], j_all[:2], k_all[:2]):
        sqltext_ijk += '(i = {0} AND j = {1} AND k = {2}) OR '.format(i,j,k)
    sqltext_ijk = sqltext_ijk[:-3] + ')' # dont end with OR
    sqltext = sqltext + sqltext_ijk + ' ) ORDER BY i,j,k'
    cur.execute(sqltext)
    wid = cur.fetchall()
    sql_debug(connection)
    connection.commit()

    if len(i_all) != len(wid):
        error('Number of items in i_all is not the same as in wid')
        exit(1)

    # Get all wall coordinates
    debug('Selecting all wall coordinates')
    sqltext = 'SELECT ' \
              ' ARRAY[ST_X(vert1), ST_Y(vert1), ST_Z(vert1)],  ' \
              ' ARRAY[ST_X(vert2), ST_Y(vert2), ST_Z(vert2)], ' \
              ' ARRAY[ST_X(vert3), ST_Y(vert3), ST_Z(vert3)], ' \
              ' ARRAY[ST_X(vert4), ST_Y(vert4), ST_Z(vert4)], ' \
              ' ARRAY[ST_X(vert5), ST_Y(vert5), ST_Z(vert5)], ' \
              ' ARRAY[ST_X(vert6), ST_Y(vert6), ST_Z(vert6)], ' \
              ' ARRAY[ST_X(vert7), ST_Y(vert7), ST_Z(vert7)] ' \
              'FROM "{0}"."{1}" ' \
              'WHERE {2} AND iswall ' \
              'ORDER BY i,j,k' \
        .format(cfg.domain.case_schema, cfg.tables.slanted_faces, sqltext_ijk)
    cur.execute(sqltext)
    wall_coord = cur.fetchall()
    sql_debug(connection)
    connection.commit()

    # transform into numpy array
    verbose('Transforming coordinates into numpy array')
    wall_coord_np = np.zeros((len(wall_coord), 7, 3), dtype=object)
    wall_vert = np.zeros(len(wall_coord), dtype=np.int)
    for idx, coord in enumerate(wall_coord):
        wall_vert[idx] = len([wc[0] for wc in wall_coord[idx] if wc[0] is not None]) - 1
        wall_coord_np[idx, 0, :] = wall_coord[idx][0]
        wall_coord_np[idx, 1, :] = wall_coord[idx][1]
        wall_coord_np[idx, 2, :] = wall_coord[idx][2]
        if wall_vert[idx] >= 3:
            wall_coord_np[idx, 3, :] = wall_coord[idx][3]
        if wall_vert[idx] >= 4:
            wall_coord_np[idx, 4, :] = wall_coord[idx][4]
        if wall_vert[idx] >= 5:
            wall_coord_np[idx, 5, :] = wall_coord[idx][5]
        if wall_vert[idx] >= 6:
            wall_coord_np[idx, 6, :] = wall_coord[idx][6]


    # Get all terrain coordinates
    debug('Selecting all terrain coordinates')
    sqltext = 'SELECT ' \
              ' ARRAY[ST_X(vert1), ST_Y(vert1), ST_Z(vert1)],  ' \
              ' ARRAY[ST_X(vert2), ST_Y(vert2), ST_Z(vert2)], ' \
              ' ARRAY[ST_X(vert3), ST_Y(vert3), ST_Z(vert3)], ' \
              ' ARRAY[ST_X(vert4), ST_Y(vert4), ST_Z(vert4)], ' \
              ' ARRAY[ST_X(vert5), ST_Y(vert5), ST_Z(vert5)], ' \
              ' ARRAY[ST_X(vert6), ST_Y(vert6), ST_Z(vert6)], ' \
              ' ARRAY[ST_X(vert7), ST_Y(vert7), ST_Z(vert7)] ' \
              'FROM "{0}"."{1}" ' \
              'WHERE {2} AND isterr ' \
              'ORDER BY i,j,k' \
              .format(cfg.domain.case_schema, cfg.tables.slanted_faces, sqltext_ijk)
    cur.execute(sqltext)
    terr_coord = cur.fetchall()
    sql_debug(connection)
    connection.commit()

    # transform into numpy array
    verbose('Transforming into numpy array')
    terr_coord_np = np.zeros((len(terr_coord), 7, 3), dtype=object)
    terr_vert = np.zeros(len(terr_coord), dtype=np.int)
    for idx, coord in enumerate(terr_coord):
        terr_vert[idx] = len([wc[0] for wc in terr_coord[idx] if wc[0] is not None]) - 1
        terr_coord_np[idx, 0, :] = terr_coord[idx][0]
        terr_coord_np[idx, 1, :] = terr_coord[idx][1]
        terr_coord_np[idx, 2, :] = terr_coord[idx][2]
        if terr_vert[idx] >= 3:
            terr_coord_np[idx, 3, :] = terr_coord[idx][3]
        if terr_vert[idx] >= 4:
            terr_coord_np[idx, 4, :] = terr_coord[idx][4]
        if terr_vert[idx] >= 5:
            terr_coord_np[idx, 5, :] = terr_coord[idx][5]
        if terr_vert[idx] >= 6:
            terr_coord_np[idx, 6, :] = terr_coord[idx][6]

    del terr_coord
    del wall_coord

    # Join them and find correct coordinates
    debug('Merging duplicates')
    vert_final = []
    for idx, i in enumerate(i_all):
        # extra_verbose('Merging i,j,k, wall_coords, terr_coords, wall_vert, terr_vert, {}, {}, {}, {}, {}, {}, {}',
        #               i_all[idx], j_all[idx], k_all[idx],
        #               wall_coord_np[idx], terr_coord_np[idx],
        #               wall_vert[idx], terr_vert[idx] )
        v_f = merge_local_wall_terrain(i_all[idx], j_all[idx], k_all[idx],
                                       wall_coord_np[idx], terr_coord_np[idx],
                                       wall_vert[idx], terr_vert[idx])
        vert_final.append(v_f)

    # COLLECT ALL INTO TUPLE
    debug('Final collecting of new vertices into one tuple')
    to_insert = []
    for idx, verts in enumerate(vert_final):
        to_insert.append((i_all[idx], j_all[idx], k_all[idx], wid[idx],
                          verts[0, 0], verts[0, 1], verts[0, 2], cfg.srid_palm,
                          verts[1, 0], verts[1, 1], verts[1, 2], cfg.srid_palm,
                          verts[2, 0], verts[2, 1], verts[2, 2], cfg.srid_palm,
                          verts[3, 0], verts[3, 1], verts[3, 2], cfg.srid_palm,
                          verts[4, 0], verts[4, 1], verts[4, 2], cfg.srid_palm,
                          verts[5, 0], verts[5, 1], verts[5, 2], cfg.srid_palm,
                          verts[6, 0], verts[6, 1], verts[6, 2], cfg.srid_palm,)
                         )
        # print(to_insert[idx])

    # DELETE old ones
    # Remove original faces
    debug('Deleting original faces')
    sqltext = 'DELETE FROM "{0}"."{1}" ' \
              'WHERE {2}'.format(cfg.domain.case_schema, cfg.tables.slanted_faces, sqltext_ijk)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # INSERT ALL NEW ENTRIES
    debug('Inserting all new entries into slanted faces table')
    sqltext = 'INSERT INTO "{0}"."{1}" (i, j, k, wid, iswall, vert1, vert2, vert3, vert4, vert5, vert6, vert7) ' \
              'VALUES                  (%s, %s, %s, %s, TRUE, ' \
              '                        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
              '                        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
              '                        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
              '                        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
              '                        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
              '                        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
              '                        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s) ' \
              '                        )  '.format(cfg.domain.case_schema, cfg.tables.slanted_faces)

    cur.executemany(sqltext, to_insert)
    sql_debug(connection)
    connection.commit()


    # Create polygon
    debug('Updating 3d polygon')
    sqltext = 'UPDATE "{0}"."{1}" SET geom =  ' \
              'ST_ForceRHR(' \
              'ST_SetSRID(ST_MakePolygon(ST_MakeLine(ARRAY[vert1, vert2, vert3, vert4, vert5, vert6, vert7,' \
              '     CASE WHEN vert1 IS NOT NULL THEN vert1 ' \
              '          WHEN vert2 IS NOT NULL THEN vert2 ' \
              '          WHEN vert3 IS NOT NULL THEN vert3 ' \
              '          WHEN vert4 IS NOT NULL THEN vert4 ' \
              '          WHEN vert5 IS NOT NULL THEN vert5 ' \
              '          WHEN vert6 IS NOT NULL THEN vert6 ' \
              '          WHEN vert7 IS NOT NULL THEN vert7 END' \
              '])), %s)) ' \
              'WHERE geom IS NULL'.format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext, (cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

    debug('Reupdate vertices')
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'vert1 = ST_PointN(ST_Boundary(geom), 1), ' \
              'vert2 = ST_PointN(ST_Boundary(geom), 2), ' \
              'vert3 = ST_PointN(ST_Boundary(geom), 3), ' \
              'vert4 = ST_PointN(ST_Boundary(geom), 4), ' \
              'vert5 = ST_PointN(ST_Boundary(geom), 5), ' \
              'vert6 = ST_PointN(ST_Boundary(geom), 6), ' \
              'vert7 = ST_PointN(ST_Boundary(geom), 7)' \
        .format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    debug('Update center of faces, just in case')
    debug('Updating center')
    sqltext = 'UPDATE "{0}"."{1}" SET center = ' \
              'ST_SetSRID(ST_MakePoint(ST_X(ST_Centroid(geom)), ST_Y(ST_Centroid(geom)), ' \
              '                       (ST_Z(vert1) + ST_Z(vert2) + ' \
              '                          CASE WHEN ST_NPoints(geom) > 3 THEN ST_Z(ST_PointN(ST_ExteriorRing(geom),3)) ELSE 0.0 END + ' \
              '                          CASE WHEN ST_NPoints(geom) > 4 THEN ST_Z(ST_PointN(ST_ExteriorRing(geom),4)) ELSE 0.0 END + ' \
              '                          CASE WHEN ST_NPoints(geom) > 5 THEN ST_Z(ST_PointN(ST_ExteriorRing(geom),5)) ELSE 0.0 END + ' \
              '                          CASE WHEN ST_NPoints(geom) > 6 THEN ST_Z(ST_PointN(ST_ExteriorRing(geom),6)) ELSE 0.0 END) / (ST_NPoints(geom)-1))' \
              '         , %s) ' \
              'WHERE center IS NULL'.format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext, (cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

    ### OLD ONE
    # sqltext = 'SELECT i,j,k FROM "{0}"."{1}" ' \
    #           'WHERE isterr OR iswall ' \
    #           'GROUP BY i,j,k ' \
    #           'HAVING COUNT(*)>1'.format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    # cur.execute(sqltext)
    # duplicits = cur.fetchall()
    # sql_debug(connection)
    # connection.commit()
    #
    # debug('Processing duplicates')
    # for i, j, k in duplicits:
    #     verbose('\tNow processing duplicate [k,j,i] = [{},{},{}]', k,j,i)
    #     sqltext = 'SELECT wid ' \
    #               'FROM "{0}"."{1}" ' \
    #               'WHERE i = {2} AND j = {3} AND k = {4} AND iswall' \
    #         .format(cfg.domain.case_schema, cfg.tables.slanted_faces, i, j, k)
    #     cur.execute(sqltext)
    #     wid = cur.fetchone()[0]
    #     sql_debug(connection)
    #     connection.commit()
    #
    #     sqltext = 'SELECT ' \
    #               ' ARRAY[ST_X(vert1), ST_Y(vert1), ST_Z(vert1)],  ' \
    #               ' ARRAY[ST_X(vert2), ST_Y(vert2), ST_Z(vert2)], ' \
    #               ' ARRAY[ST_X(vert3), ST_Y(vert3), ST_Z(vert3)], ' \
    #               ' ARRAY[ST_X(vert4), ST_Y(vert4), ST_Z(vert4)], ' \
    #               ' ARRAY[ST_X(vert5), ST_Y(vert5), ST_Z(vert5)], ' \
    #               ' ARRAY[ST_X(vert6), ST_Y(vert6), ST_Z(vert6)], ' \
    #               ' ARRAY[ST_X(vert7), ST_Y(vert7), ST_Z(vert7)] ' \
    #               'FROM "{0}"."{1}" ' \
    #               'WHERE i = {2} AND j = {3} AND k = {4} AND iswall'\
    #               .format(cfg.domain.case_schema, cfg.tables.slanted_faces, i, j, k)
    #     cur.execute(sqltext)
    #     wall_coord = list(cur.fetchall()[0])
    #     wall_coord1 = []
    #     for wall_c in wall_coord:
    #         if not wall_c[0] is None:
    #             wall_coord1.append(wall_c)
    #     wall_coord = np.asarray(wall_coord1[:-1])
    #     sql_debug(connection)
    #     connection.commit()
    #
    #     sqltext = 'SELECT ' \
    #               ' ARRAY[ST_X(vert1), ST_Y(vert1), ST_Z(vert1)],  ' \
    #               ' ARRAY[ST_X(vert2), ST_Y(vert2), ST_Z(vert2)], ' \
    #               ' ARRAY[ST_X(vert3), ST_Y(vert3), ST_Z(vert3)], ' \
    #               ' ARRAY[ST_X(vert4), ST_Y(vert4), ST_Z(vert4)], ' \
    #               ' ARRAY[ST_X(vert5), ST_Y(vert5), ST_Z(vert5)], ' \
    #               ' ARRAY[ST_X(vert6), ST_Y(vert6), ST_Z(vert6)], ' \
    #               ' ARRAY[ST_X(vert7), ST_Y(vert7), ST_Z(vert7)] ' \
    #               'FROM "{0}"."{1}" ' \
    #               'WHERE i = {2} AND j = {3} AND k = {4} AND isterr'\
    #               .format(cfg.domain.case_schema, cfg.tables.slanted_faces, i, j, k)
    #     cur.execute(sqltext)
    #     terr_coord = list(cur.fetchall()[0])
    #     terr_coord1 = []
    #     for terr_c in terr_coord:
    #         if not terr_c[0] is None:
    #             terr_coord1.append(terr_c)
    #     terr_coord = np.asarray(terr_coord1[:-1])
    #     sql_debug(connection)
    #     connection.commit()
    #
    #     # Remove original faces
    #     verbose('\tDeleting original faces')
    #     sqltext = 'DELETE FROM "{0}"."{1}" ' \
    #               'WHERE i = {2} AND j = {3} AND k = {4}'.format(cfg.domain.case_schema, cfg.tables.slanted_faces, i, j, k)
    #     cur.execute(sqltext)
    #     sql_debug(connection)
    #     connection.commit()
    #
    #     # FIND which terrain points are wall point and find appropriate ones
    #     # only top ones
    #     wall_coord_t = wall_coord[wall_coord[:, 2] >= np.max(wall_coord[:, 2])]
    #     trtw = []
    #     tr_new = []
    #     for tri, tr in enumerate(terr_coord):
    #         tr_new.append(tr)
    #         trtw0 = [tri, -1, False]
    #         found = False
    #         for twi, tw in enumerate(wall_coord_t):
    #             if tw[0] == tr[0] and tw[1] == tr[1]:
    #                 trtw0 = [tri, twi, True]
    #                 trtw.append(trtw0)
    #                 if tw[2] > tr[2]:
    #                     trtw.append(trtw0)
    #                     tr_new.append(tw)
    #                 found = True
    #                 break
    #         if not found:
    #             trtw.append(trtw0)
    #
    #     tr_new.append(tr_new[0])
    #     vertices = np.asarray(tr_new)
    #
    #     n_vert = vertices.shape[0]
    #
    #     cent = np.mean(vertices, 0)
    #
    #     # calculate normal vector of some triangles
    #     vi = 0
    #     vj = 1
    #     A = cent - vertices[vi, :]
    #     B = cent - vertices[vj, :]
    #     # n = np.array([A[1] * B[2] - A[2] * B[1], A[2] * B[0] - A[0] * B[2], A[0] * B[1] - A[1] * B[0]])
    #     #
    #
    #     # b = vertices[:,2].T
    #     # Ab = np.concatenate((vertices[:,:2], np.ones((1,n_vert)).T),axis=1)
    #     # Ab[:, 0] = Ab[:, 0] # - cfg.domain.origin_x
    #     # Ab[:, 1] = Ab[:, 1] #- cfg.domain.origin_y
    #     # fit, residual, rnk, s = lstsq(Ab, b)
    #     # n = fit / norm(fit)
    #
    #     svd = np.linalg.svd(vertices.T - np.mean(vertices.T, axis=1, keepdims=True))
    #     left = svd[0]
    #     left[:, -1]
    #
    #     n = left[:, -1] / norm(left[:, -1])
    #
    #     q = np.cross(A, n)
    #
    #     angle = np.zeros(n_vert)
    #     for ir in range(n_vert):
    #         r = cent - vertices[ir, :]
    #         t = np.dot(n, np.cross(r, A))
    #         u = np.dot(n, np.cross(r, q))
    #         angle[ir] = np.arctan2(u, t) / np.pi * 180
    #         angle[ir] = angle[ir] if angle[ir] > 0 else 360.0 + angle[ir]
    #
    #     order = angle.argsort()
    #     order_all = np.append(order, np.arange(n_vert, 7))
    #
    #     # fig = plt.figure()
    #     # from mpl_toolkits import mplot3d
    #     # ax = plt.axes(projection='3d')
    #     # plt.title('i,j,k: {},{},{}'.format(i, j, k))
    #     # ax.plot3D(vertices[:, 0], vertices[:, 1], vertices[:, 2], 'o')
    #     # ax.plot3D(vertices[:, 0], vertices[:, 1], vertices[:, 2], '-b')
    #     # ax.plot3D(vertices[np.append(order, order[0]), 0],
    #     #           vertices[np.append(order, order[0]), 1],
    #     #           vertices[np.append(order, order[0]), 2], 'k-')
    #     # plt.show()
    #     vert_order = np.zeros((8, 3), dtype=object)
    #     vert_order[:] = None
    #     vert_order[:n_vert] = vertices[order]
    #
    #     verbose('\tInserting new merged face')
    #     sqltext = 'INSERT INTO "{0}"."{1}" (i, j, k, wid, iswall, isroof, isterr,' \
    #               '                         vert1, vert2, vert3, vert4, vert5, vert6, vert7) ' \
    #               'SELECT {2}, {3}, {4}, {5}, TRUE, FALSE, FALSE, ' \
    #               '       ST_SetSRID(ST_MakePoint(%s,%s,%s), %s), ' \
    #               '       ST_SetSRID(ST_MakePoint(%s,%s,%s), %s), ' \
    #               '       ST_SetSRID(ST_MakePoint(%s,%s,%s), %s), ' \
    #               '       ST_SetSRID(ST_MakePoint(%s,%s,%s), %s), ' \
    #               '       ST_SetSRID(ST_MakePoint(%s,%s,%s), %s), ' \
    #               '       ST_SetSRID(ST_MakePoint(%s,%s,%s), %s), ' \
    #               '       ST_SetSRID(ST_MakePoint(%s,%s,%s), %s)' \
    #         .format(cfg.domain.case_schema, cfg.tables.slanted_faces, i, j, k, wid, *order_all + 1)
    #     cur.execute(sqltext, (vert_order[0,0], vert_order[0,1], vert_order[0,2], cfg.srid_palm,
    #                           vert_order[1,0], vert_order[1,1], vert_order[1,2], cfg.srid_palm,
    #                           vert_order[2,0], vert_order[2,1], vert_order[2,2], cfg.srid_palm,
    #                           vert_order[3,0], vert_order[3,1], vert_order[3,2], cfg.srid_palm,
    #                           vert_order[4,0], vert_order[4,1], vert_order[4,2], cfg.srid_palm,
    #                           vert_order[5,0], vert_order[5,1], vert_order[5,2], cfg.srid_palm,
    #                           vert_order[6,0], vert_order[6,1], vert_order[6,2], cfg.srid_palm,))
    #     sql_debug(connection)
    #     connection.commit()
    #
    #     verbose('Update 3d polygon')
    #     # Create polygon
    #     sqltext = 'UPDATE "{0}"."{1}" SET geom =  ' \
    #               'ST_ForceRHR(' \
    #               'ST_SetSRID(ST_MakePolygon(ST_MakeLine(ARRAY[vert1, vert2, vert3, vert4, vert5, vert6, vert7,' \
    #               '     CASE WHEN vert1 IS NOT NULL THEN vert1 ' \
    #               '          WHEN vert2 IS NOT NULL THEN vert2 ' \
    #               '          WHEN vert3 IS NOT NULL THEN vert3 ' \
    #               '          WHEN vert4 IS NOT NULL THEN vert4 ' \
    #               '          WHEN vert5 IS NOT NULL THEN vert5 ' \
    #               '          WHEN vert6 IS NOT NULL THEN vert6 ' \
    #               '          WHEN vert7 IS NOT NULL THEN vert7 END' \
    #               '])), %s)) ' \
    #               'WHERE i = {2} AND j = {3} AND k = {4}'.format(cfg.domain.case_schema, cfg.tables.slanted_faces, i, j, k)
    #     cur.execute(sqltext, (cfg.srid_palm,))
    #     sql_debug(connection)
    #     connection.commit()
    #
    #     sqltext = 'UPDATE "{0}"."{1}" SET ' \
    #               'vert1 = ST_PointN(ST_Boundary(geom), 1), ' \
    #               'vert2 = ST_PointN(ST_Boundary(geom), 2), ' \
    #               'vert3 = ST_PointN(ST_Boundary(geom), 3), ' \
    #               'vert4 = ST_PointN(ST_Boundary(geom), 4), ' \
    #               'vert5 = ST_PointN(ST_Boundary(geom), 5), ' \
    #               'vert6 = ST_PointN(ST_Boundary(geom), 6), ' \
    #               'vert7 = ST_PointN(ST_Boundary(geom), 7)' \
    #         .format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    #     cur.execute(sqltext)
    #     sql_debug(connection)
    #     connection.commit()
    #
    #     sqltext = ' UPDATE "{0}"."{1}" SET center = ST_SetSRID(ST_MakePoint(ST_X(ST_Centroid(geom)), ST_Y(ST_Centroid(geom)), ' \
    #               '     (COALESCE(ST_Z(ST_PointN(ST_ExteriorRing(geom),1)), 0.0) + COALESCE(ST_Z(ST_PointN(ST_ExteriorRing(geom),2)), 0.0) + ' \
    #               '      COALESCE(ST_Z(ST_PointN(ST_ExteriorRing(geom),3)), 0.0) + COALESCE(ST_Z(ST_PointN(ST_ExteriorRing(geom),4)), 0.0) + ' \
    #               '      COALESCE(ST_Z(ST_PointN(ST_ExteriorRing(geom),5)), 0.0) + COALESCE(ST_Z(ST_PointN(ST_ExteriorRing(geom),6)), 0.0) + ' \
    #               '      COALESCE(ST_Z(ST_PointN(ST_ExteriorRing(geom),7)), 0.0)) / ST_NPoints(geom))' \
    #               '         , %s)' \
    #               'WHERE i = {2} AND j = {3} AND k = {4}'.format(cfg.domain.case_schema, cfg.tables.slanted_faces, i, j, k)
    #     cur.execute(sqltext, (cfg.srid_palm,))
    #     sql_debug(connection)
    #     connection.commit()

def normal_vector_triangle(p1, p2, p3):
    N = np.cross(p2-p1, p3-p1)
    return N[::-1] / np.sqrt(np.sum(N ** 2))


def merge_walls_roofs(cfg, connection, cur):
    """ Function that merge wall and roof occupying same grid box"""
    progress('Merging grid faces of walls and roofs')
    # find all aggregated faces, create new table with aggregated ones
    debug('Finding duplicates')
    sqltext = 'SELECT i,j,k FROM "{0}"."{1}" ' \
              'WHERE NOT isterr ' \
              'GROUP BY i,j,k ' \
              'HAVING COUNT(*)>1' \
              'ORDER BY i,j,k'.format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext)
    duplicits = cur.fetchall()
    sql_debug(connection)
    connection.commit()
    i_all = [x[0] for x in duplicits]
    j_all = [x[1] for x in duplicits]
    k_all = [x[2] for x in duplicits]

    if len(i_all) == 0:
        return 1

    # SELECT ALL rid
    debug('Finding all rids')
    sqltext = 'SELECT rid FROM "{0}"."{1}" WHERE isroof AND ('.format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    sqltext_ijk = '('
    for i,j,k in zip(i_all, j_all, k_all):
    # for i, j, k in zip([291], [136], [17]):
        sqltext_ijk += '(i = {0} AND j = {1} AND k = {2}) OR '.format(i,j,k)
    sqltext_ijk = sqltext_ijk[:-3] + ')'
    sqltext = sqltext + sqltext_ijk + ' ) ORDER BY i,j,k'  # dont end with OR
    cur.execute(sqltext)
    rid = cur.fetchall()
    sql_debug(connection)
    connection.commit()

    # SELECT ALL Vertices
    debug('Selecting all coordinates')
    sqltext = 'SELECT COUNT(*), ARRAY_AGG(ST_X(points)), ARRAY_AGG(ST_Y(points)), ARRAY_AGG(ST_Z(points)), i, j, k ' \
              'FROM (SELECT points, i, j, k' \
              '      FROM (SELECT (ST_DumpPoints(geom)).geom AS points, i, j, k FROM "{0}"."{1}" ' \
              '            WHERE NOT isterr AND {2}) AS s ' \
              '      GROUP BY points, i, j, k ) AS ss ' \
              'GROUP BY i,j,k ' \
              'ORDER BY i,j,k'.format(cfg.domain.case_schema, cfg.tables.slanted_faces, sqltext_ijk)
    cur.execute(sqltext)
    verts = cur.fetchall()
    sql_debug(connection)
    connection.commit()

    verbose('Creating counts')
    count = [x[0] for x in verts]

    verbose('Creating x verts')
    x_vert = [x[1] for x in verts]
    x_vert_np = np.empty((len(count), 7), dtype=object)
    for i, j in enumerate(x_vert):
        x_vert_np[i][0:len(j)] = j
    del x_vert

    verbose('Creating y verts')
    y_vert = [x[2] for x in verts]
    y_vert_np = np.empty((len(count), 7), dtype=object)
    for i, j in enumerate(y_vert):
        y_vert_np[i][0:len(j)] = j
    del y_vert

    verbose('Creating z verts')
    z_vert = [x[3] for x in verts]
    z_vert_np = np.empty((len(count), 7), dtype=object)
    for i, j in enumerate(z_vert):
        z_vert_np[i][0:len(j)] = j
    del z_vert

    verbose('Deleting verts')
    del verts

    debug('Delete the ones that does not lie at grid lines')
    x_vert_np_n, y_vert_np_n, z_vert_np_n = np.empty((len(count), 7), dtype=object), np.empty((len(count), 7), dtype=object), np.empty((len(count), 7), dtype=object)
    for i in range(len(count)):
        q = 0
        for j in range(count[i]):
            xdist = np.abs(
                (x_vert_np[i, j] - cfg.domain.origin_x) / cfg.domain.dx
                - np.round((x_vert_np[i, j] - cfg.domain.origin_x) / cfg.domain.dx))

            ydist = np.abs(
                (y_vert_np[i, j] - cfg.domain.origin_y) / cfg.domain.dy
                - np.round((y_vert_np[i, j] - cfg.domain.origin_y) / cfg.domain.dy))
            if xdist > 1e-10 and ydist > 1e-10:
                pass
            else:
                x_vert_np_n[i, q] = x_vert_np[i, j]
                y_vert_np_n[i, q] = y_vert_np[i, j]
                z_vert_np_n[i, q] = z_vert_np[i, j]
                q += 1
        count[i] = q

    x_vert_np, y_vert_np, z_vert_np = x_vert_np_n, y_vert_np_n, z_vert_np_n

    debug('Merging all coordinates into final form')
    order_all = []
    for idx, c in enumerate(count):
        # extra_verbose('Merging i,j,k, x coords, y coords, z coords, count, {}, {}, {}, {}, {}, {}, {}',
        #               i_all[idx], j_all[idx], k_all[idx],
        #               x_vert_np[idx], y_vert_np[idx], y_vert_np[idx], count[idx])
        or_all = merge_local_wall_roof(i_all[idx], j_all[idx], k_all[idx], x_vert_np[idx], y_vert_np[idx], z_vert_np[idx], count[idx], cfg)
        order_all.append(or_all)
    order_all = np.asarray(order_all)

    # COLLECT ALL INTO TUPLE
    debug('Final collecting of new vertices into one tuple')
    to_insert = []
    for idx, order in enumerate(order_all):
        to_insert.append((i_all[idx], j_all[idx], k_all[idx], rid[idx],
                          x_vert_np[idx, order[0]], y_vert_np[idx, order[0]], z_vert_np[idx, order[0]], cfg.srid_palm,
                          x_vert_np[idx, order[1]], y_vert_np[idx, order[1]], z_vert_np[idx, order[1]], cfg.srid_palm,
                          x_vert_np[idx, order[2]], y_vert_np[idx, order[2]], z_vert_np[idx, order[2]], cfg.srid_palm,
                          x_vert_np[idx, order[3]], y_vert_np[idx, order[3]], z_vert_np[idx, order[3]], cfg.srid_palm,
                          x_vert_np[idx, order[4]], y_vert_np[idx, order[4]], z_vert_np[idx, order[4]], cfg.srid_palm,
                          x_vert_np[idx, order[5]], y_vert_np[idx, order[5]], z_vert_np[idx, order[5]], cfg.srid_palm,
                          x_vert_np[idx, order[6]], y_vert_np[idx, order[6]], z_vert_np[idx, order[6]], cfg.srid_palm,)
                         )
        # print(to_insert[idx])


    # DELETE previous entry
    debug('Deleting original faces')
    sqltext = 'DELETE FROM "{0}"."{1}" ' \
              'WHERE {2}'.format(cfg.domain.case_schema, cfg.tables.slanted_faces, sqltext_ijk)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # INSERT ALL NEW ENTRIES
    debug('Inserting all new entries into slanted faces table')
    sqltext = 'INSERT INTO "{0}"."{1}" (i, j, k, rid, isroof, vert1, vert2, vert3, vert4, vert5, vert6, vert7) ' \
              'VALUES                  (%s, %s, %s, %s, TRUE, ' \
              '                        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
              '                        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
              '                        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
              '                        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
              '                        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
              '                        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
              '                        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s) ' \
              '                        )  '.format(cfg.domain.case_schema, cfg.tables.slanted_faces)

    cur.executemany(sqltext, to_insert)
    sql_debug(connection)
    connection.commit()

    # Create polygon
    debug('Updating 3d polygon')
    sqltext = 'UPDATE "{0}"."{1}" SET geom =  ' \
              'ST_ForceRHR(' \
              'ST_SetSRID(ST_MakePolygon(ST_MakeLine(ARRAY[vert1, vert2, vert3, vert4, vert5, vert6, vert7,' \
              '     CASE WHEN vert1 IS NOT NULL THEN vert1 ' \
              '          WHEN vert2 IS NOT NULL THEN vert2 ' \
              '          WHEN vert3 IS NOT NULL THEN vert3 ' \
              '          WHEN vert4 IS NOT NULL THEN vert4 ' \
              '          WHEN vert5 IS NOT NULL THEN vert5 ' \
              '          WHEN vert6 IS NOT NULL THEN vert6 ' \
              '          WHEN vert7 IS NOT NULL THEN vert7 END' \
              '])), %s)) ' \
              'WHERE geom IS NULL'.format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext, (cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

    debug('Reupdate vertices')
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'vert1 = ST_PointN(ST_Boundary(geom), 1), ' \
              'vert2 = ST_PointN(ST_Boundary(geom), 2), ' \
              'vert3 = ST_PointN(ST_Boundary(geom), 3), ' \
              'vert4 = ST_PointN(ST_Boundary(geom), 4), ' \
              'vert5 = ST_PointN(ST_Boundary(geom), 5), ' \
              'vert6 = ST_PointN(ST_Boundary(geom), 6), ' \
              'vert7 = ST_PointN(ST_Boundary(geom), 7)' \
        .format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    debug('Update center of faces, just in case')
    debug('Updating center')
    sqltext = 'UPDATE "{0}"."{1}" SET center = ' \
              'ST_SetSRID(ST_MakePoint(ST_X(ST_Centroid(geom)), ST_Y(ST_Centroid(geom)), ' \
              '                       (ST_Z(vert1) + ST_Z(vert2) + ' \
              '                          CASE WHEN ST_NPoints(geom) > 3 THEN ST_Z(ST_PointN(ST_ExteriorRing(geom),3)) ELSE 0.0 END + ' \
              '                          CASE WHEN ST_NPoints(geom) > 4 THEN ST_Z(ST_PointN(ST_ExteriorRing(geom),4)) ELSE 0.0 END + ' \
              '                          CASE WHEN ST_NPoints(geom) > 5 THEN ST_Z(ST_PointN(ST_ExteriorRing(geom),5)) ELSE 0.0 END + ' \
              '                          CASE WHEN ST_NPoints(geom) > 6 THEN ST_Z(ST_PointN(ST_ExteriorRing(geom),6)) ELSE 0.0 END) / (ST_NPoints(geom) - 1))' \
              '         , %s) ' \
              'WHERE center IS NULL'.format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext, (cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

    sqltext = 'WITH ijk AS (SELECT i,j,k ' \
          '                 FROM "{0}"."{1}"  ' \
              '             WHERE isterr ' \
              '             GROUP BY i,j,k ' \
              '             HAVING COUNT(*) = 2 ' \
              '             ORDER BY i,j,k ) ' \
              'SELECT ST_Z(center), id, s.i, s.j, s.k  ' \
              'FROM "{0}"."{1}" AS s ' \
              'RIGHT JOIN ijk ON ijk.i = s.i AND ijk.j = s.j AND ijk.k = s.k ' \
              'WHERE isterr ' \
              'ORDER BY s.i, s.j, s.k'.format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext)
    zijk = cur.fetchall()
    sql_debug(connection)
    connection.commit()

    id2del = []
    for ijk in range(0, len(zijk), 2):
        z1, z2 = zijk[ijk][0], zijk[ijk+1][0]
        id1, id2 = zijk[ijk][1], zijk[ijk+1][1]
        id2del.append((id1 if z1 < z2 else id2,))

    debug('Deleting all unwanted rows')
    sqltext = 'DELETE FROM "{0}"."{1}" ' \
              'WHERE id = ANY(%s)'.format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext, (id2del,))
    sql_debug(connection)
    connection.commit()

    #
    # for i, j, k in duplicits:
    #     verbose('Processing [i,j,k] = [{},{},{}]', i,j,k)
    #     sqltext = 'SELECT rid FROM "{0}"."{1}" ' \
    #               'WHERE i = {2} AND j = {3} AND k = {4} AND isroof'\
    #               .format(cfg.domain.case_schema, cfg.tables.slanted_faces, i, j, k)
    #     cur.execute(sqltext)
    #     rid = cur.fetchone()[0]
    #     sql_debug(connection)
    #     connection.commit()
    #
    #     # delete old one and create new one
    #     sqltext = 'SELECT ST_Collect(ARRAY_AGG(points)) ' \
    #               '     FROM (' \
    #               '           SELECT points FROM ( ' \
    #               '                               SELECT i, j, k, iswall, (ST_DumpPoints(geom)).geom AS points ' \
    #               '                               FROM "{0}"."{1}"' \
    #               '                               WHERE i = {2} AND j = {3} AND k = {4}) AS s ' \
    #               '     GROUP BY points) AS ss'.format(cfg.domain.case_schema, cfg.tables.slanted_faces, i, j, k)
    #     cur.execute(sqltext, (cfg.srid_palm, ))
    #     polygon = cur.fetchone()[0]
    #     sql_debug(connection)
    #     connection.commit()
    #
    #     sqltext = 'DELETE FROM "{0}"."{1}" ' \
    #               'WHERE i = {2} AND j = {3} AND k = {4}'.format(cfg.domain.case_schema, cfg.tables.slanted_faces, i, j, k)
    #     cur.execute(sqltext)
    #     sql_debug(connection)
    #     connection.commit()
    #
    #     # reorder the data to form polygon
    #     sqltext = 'SELECT ST_X(ST_GeometryN(%s,1)), ST_X(ST_GeometryN(%s,2)), ST_X(ST_GeometryN(%s,3)), ' \
    #               '       ST_X(ST_GeometryN(%s,4)), ST_X(ST_GeometryN(%s,5)), ST_X(ST_GeometryN(%s,6)),' \
    #               '       ST_X(ST_GeometryN(%s,7)) '. \
    #               format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    #     cur.execute(sqltext, (polygon, polygon, polygon, polygon, polygon, polygon, polygon, ))
    #     x_vert = cur.fetchall()
    #     sql_debug(connection)
    #     connection.commit()
    #
    #     sqltext = 'SELECT ST_Y(ST_GeometryN(%s,1)), ST_Y(ST_GeometryN(%s,2)), ST_Y(ST_GeometryN(%s,3)), ' \
    #               '       ST_Y(ST_GeometryN(%s,4)), ST_Y(ST_GeometryN(%s,5)), ST_Y(ST_GeometryN(%s,6)),' \
    #               '       ST_Y(ST_GeometryN(%s,7)) '. \
    #               format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    #     cur.execute(sqltext, (polygon, polygon, polygon, polygon, polygon, polygon, polygon, ))
    #     y_vert = cur.fetchall()
    #     sql_debug(connection)
    #     connection.commit()
    #
    #     sqltext = 'SELECT ST_Z(ST_GeometryN(%s,1)), ST_Z(ST_GeometryN(%s,2)), ST_Z(ST_GeometryN(%s,3)), ' \
    #               '       ST_Z(ST_GeometryN(%s,4)), ST_Z(ST_GeometryN(%s,5)), ST_Z(ST_GeometryN(%s,6)),' \
    #               '       ST_Z(ST_GeometryN(%s,7)) '. \
    #               format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    #     cur.execute(sqltext, (polygon, polygon, polygon, polygon, polygon, polygon, polygon, ))
    #     z_vert = cur.fetchall()
    #     sql_debug(connection)
    #     connection.commit()
    #
    #     x_vert = np.asarray([np.nan if x is None else x for x in x_vert[0]])
    #     y_vert = np.asarray([np.nan if x is None else x for x in y_vert[0]])
    #     z_vert = np.asarray([np.nan if x is None else x for x in z_vert[0]])
    #
    #     n_vert = np.sum(~np.isnan(x_vert))
    #
    #     vertices = np.zeros((n_vert,3))
    #     vertices[:, 0] = x_vert[:n_vert]
    #     vertices[:, 1] = y_vert[:n_vert]
    #     vertices[:, 2] = z_vert[:n_vert]
    #
    #
    #     cent = np.mean(vertices,0)
    #
    #     # calculate normal vector of some triangles
    #     vi = 0
    #     vj = 1
    #     A = cent - vertices[vi, :]
    #     B = cent - vertices[vj, :]
    #     # n = np.array([A[1] * B[2] - A[2] * B[1], A[2] * B[0] - A[0] * B[2], A[0] * B[1] - A[1] * B[0]])
    #     #
    #
    #     # b = vertices[:,2].T
    #     # Ab = np.concatenate((vertices[:,:2], np.ones((1,n_vert)).T),axis=1)
    #     # Ab[:, 0] = Ab[:, 0] # - cfg.domain.origin_x
    #     # Ab[:, 1] = Ab[:, 1] #- cfg.domain.origin_y
    #     # fit, residual, rnk, s = lstsq(Ab, b)
    #     # n = fit / norm(fit)
    #
    #     svd = np.linalg.svd(vertices.T - np.mean(vertices.T , axis=1, keepdims=True))
    #     left = svd[0]
    #     left[:, -1]
    #
    #     n = left[:, -1] / norm(left[:, -1])
    #
    #     q = np.cross(A, n)
    #
    #     angle = np.zeros(n_vert)
    #     for ir in range(n_vert):
    #         r = cent - vertices[ir, :]
    #         t = np.dot(n, np.cross(r, A))
    #         u = np.dot(n, np.cross(r, q))
    #         angle[ir] = np.arctan2(u, t) / np.pi * 180
    #         angle[ir] = angle[ir] if angle[ir] > 0 else 360.0 + angle[ir]
    #
    #     order = angle.argsort()
    #     order_all = np.append(order, np.arange(n_vert, 7))
    #     #
    #     # fig = plt.figure()
    #     # from mpl_toolkits import mplot3d
    #     # ax = plt.axes(projection='3d')
    #     # plt.title('i,j,k: {},{},{}'.format(i, j, k))
    #     # ax.plot3D(vertices[:, 0], vertices[:, 1], vertices[:, 2], 'o')
    #     # ax.plot3D([cent[0], cent[0] + n[0]], [cent[1], cent[1] + n[1]], [cent[2], cent[2] + n[2]], 'k-')
    #     # ax.plot3D([cent[0], cent[0] + q[0]], [cent[1], cent[1] + q[1]], [cent[2], cent[2] + q[2]], 'k:')
    #     # ax.plot3D([cent[0], cent[0] + A[0]], [cent[1], cent[1] + A[1]], [cent[2], cent[2] + A[2]], 'k--')
    #     # ax.plot3D(vertices[np.append(order, order[0]), 0],
    #     #           vertices[np.append(order, order[0]), 1],
    #     #           vertices[np.append(order, order[0]), 2], 'k-')
    #     # # xlim = ax.get_xlim()
    #     # # ylim = ax.get_ylim()
    #     # # X, Y = np.meshgrid(np.arange(xlim[0], xlim[1], 0.1),
    #     # #                    np.arange(ylim[0], ylim[1], 0.1))
    #     # # Z = np.zeros(X.shape)
    #     # # for r in range(X.shape[0]):
    #     # #     for c in range(X.shape[1]):
    #     # #         Z[r, c] = n[0] * X[r, c] + n[1] * Y[r, c] + n[2] * Z[r ,c]
    #     # # ax.plot_wireframe(X, Y, Z, color='k')
    #     # plt.show()
    #
    #     sqltext = 'INSERT INTO "{0}"."{1}" (i, j, k, rid, vert1, vert2, vert3, vert4, vert5, vert6, vert7) ' \
    #               'SELECT {2}, {3}, {4}, {5}, ' \
    #               'ST_GeometryN(%s,{6}), ST_GeometryN(%s,{7}), ST_GeometryN(%s,{8}), ' \
    #               'ST_GeometryN(%s,{9}), ST_GeometryN(%s,{10}), ST_GeometryN(%s,{11}), ST_GeometryN(%s,{12}) '\
    #         .format(cfg.domain.case_schema, cfg.tables.slanted_faces, i, j, k, rid, *order_all+1)
    #     cur.execute(sqltext, (polygon, polygon, polygon, polygon, polygon, polygon, polygon, ))
    #     sql_debug(connection)
    #     connection.commit()

        # # Create polygon
        # sqltext = 'UPDATE "{0}"."{1}" SET geom =  ' \
        #           'ST_ForceRHR(' \
        #           'ST_SetSRID(ST_MakePolygon(ST_MakeLine(ARRAY[vert1, vert2, vert3, vert4, vert5, vert6, vert7,' \
        #           '     CASE WHEN vert1 IS NOT NULL THEN vert1 ' \
        #           '          WHEN vert2 IS NOT NULL THEN vert2 ' \
        #           '          WHEN vert3 IS NOT NULL THEN vert3 ' \
        #           '          WHEN vert4 IS NOT NULL THEN vert4 ' \
        #           '          WHEN vert5 IS NOT NULL THEN vert5 ' \
        #           '          WHEN vert6 IS NOT NULL THEN vert6 ' \
        #           '          WHEN vert7 IS NOT NULL THEN vert7 END' \
        #           '])), %s)) ' \
        #           'WHERE i = {2} AND j = {3} AND k = {4}'.format(cfg.domain.case_schema, cfg.tables.slanted_faces, i, j, k)
        # cur.execute(sqltext, (cfg.srid_palm, ))
        # sql_debug(connection)
        # connection.commit()
        #
        # sqltext = 'UPDATE "{0}"."{1}" SET ' \
        #           'vert1 = ST_PointN(ST_Boundary(geom), 1), ' \
        #           'vert2 = ST_PointN(ST_Boundary(geom), 2), ' \
        #           'vert3 = ST_PointN(ST_Boundary(geom), 3), ' \
        #           'vert4 = ST_PointN(ST_Boundary(geom), 4), ' \
        #           'vert5 = ST_PointN(ST_Boundary(geom), 5), ' \
        #           'vert6 = ST_PointN(ST_Boundary(geom), 6), ' \
        #           'vert7 = ST_PointN(ST_Boundary(geom), 7)' \
        #     .format(cfg.domain.case_schema, cfg.tables.slanted_faces)
        # cur.execute(sqltext)
        # sql_debug(connection)
        # connection.commit()
        #
        # sqltext = ' UPDATE "{0}"."{1}" SET center = ST_SetSRID(ST_MakePoint(ST_X(ST_Centroid(geom)), ST_Y(ST_Centroid(geom)), ' \
        #           '     (COALESCE(ST_Z(ST_PointN(ST_ExteriorRing(geom),1)), 0.0) + COALESCE(ST_Z(ST_PointN(ST_ExteriorRing(geom),2)), 0.0) + ' \
        #           '      COALESCE(ST_Z(ST_PointN(ST_ExteriorRing(geom),3)), 0.0) + COALESCE(ST_Z(ST_PointN(ST_ExteriorRing(geom),4)), 0.0) + ' \
        #           '      COALESCE(ST_Z(ST_PointN(ST_ExteriorRing(geom),5)), 0.0) + COALESCE(ST_Z(ST_PointN(ST_ExteriorRing(geom),6)), 0.0) + ' \
        #           '      COALESCE(ST_Z(ST_PointN(ST_ExteriorRing(geom),7)), 0.0)) / ST_NPoints(geom))' \
        #           '         , %s)' \
        #           'WHERE i = {2} AND j = {3} AND k = {4}'.format(cfg.domain.case_schema, cfg.tables.slanted_faces, i, j, k)
        # cur.execute(sqltext, (cfg.srid_palm,))
        # sql_debug(connection)
        # connection.commit()

def merge_local_wall_roof(i, j, k, x_vert, y_vert, z_vert, n_vert, cfg):
    """ merge vertices on local machine """
    vertices = np.zeros((n_vert, 3))
    vertices[:, 0] = np.asarray(x_vert[:n_vert]) - cfg.domain.origin_x
    vertices[:, 1] = np.asarray(y_vert[:n_vert]) - cfg.domain.origin_y
    vertices[:, 2] = np.asarray(z_vert[:n_vert])

    # cent = np.mean(vertices, 0)
    # cent = np.mean(np.unique(vertices, axis=0), 0)
    cent = np.mean(np.unique(vertices[:, :2], axis=0), 0)

    # calculate normal vector of some triangles
    # A = cent - vertices[1, :]
    #
    # svd = np.linalg.svd(vertices.T - np.mean(vertices.T, axis=1, keepdims=True))
    # left = svd[0]
    # # left[:, -1]
    #
    # n = left[:, -1] / norm(left[:, -1])
    #
    # q = np.cross(A, n)
    #
    # angle = np.zeros(n_vert)
    # for ir in range(n_vert):
    #     r = cent - vertices[ir, :]
    #     t = np.dot(n, np.cross(r, A))
    #     u = np.dot(n, np.cross(r, q))
    #     angle[ir] = np.arctan2(u, t) / np.pi * 180
    #     angle[ir] = angle[ir] if angle[ir] > 0 else 360.0 + angle[ir]
    #
    # order = angle.argsort()

    # cx, cy = list_of_xy_coords.mean(0)

    angles = np.arctan2(vertices[:, 0] - cent[0], vertices[:, 1] - cent[1])
    order = np.argsort(angles)


    # vertices[:, 0] = x_vert[:n_vert]
    # vertices[:, 1] = y_vert[:n_vert]
    # vertices[:, 2] = z_vert[:n_vert]
    # if i == 17 and j == 19 and k == 14:
    #     print(vertices, order)
    #     print(vertices[order])
    #     # exit(1)
    for ir in range(n_vert):
        ir_curr = ir
        ir_prev = ir - 1 if ir > 0 else n_vert - 1
        ir_next = ir + 1 if ir < n_vert - 1 else 0
        curr = order[ir]
        prev = order[ir_prev]
        next = order[ir_next]
        # check if following point lies above each other
        if vertices[curr, 0] == vertices[next, 0] and vertices[curr, 1] == vertices[next, 1]:
        #     # if current and previos has the same height, keep, else change curr and next:
            if not vertices[curr, 2] == vertices[prev, 2] and vertices[next, 2] == vertices[prev, 2]:
                order[ir_curr], order[ir_next] = order[ir_next], order[ir_curr]
            # elif vertices[next, 2] == vertices[prev, 2]:
            #     order[ir_curr], order[ir_next] = order[ir_next], order[ir_curr]
        # if vertices[curr, 0] == vertices[prev, 0] and vertices[curr, 1] == vertices[prev, 1]:
        #     # if current and previos has the same height, keep, else change curr and next:
        #     if not vertices[curr, 2] == vertices[next, 2] and vertices[next, 2] == vertices[prev, 2]:
        #         order[ir_curr], order[ir_prev] = order[ir_prev], order[ir_curr]
            # elif vertices[next, 2] == vertices[prev, 2]:
            #     order[ir_curr], order[ir_next] = order[ir_next], order[ir_curr]
    # if i == 17 and j == 19 and k == 14:
    #     print(vertices, order)
    #     print(vertices[order])
    #     exit(1)

    # fig = plt.figure()
    # from mpl_toolkits import mplot3d
    # ax = plt.axes(projection='3d')
    # plt.title('i,j,k: {},{},{}'.format(i, j, k))
    # ax.plot3D(vertices[:, 0], vertices[:, 1], vertices[:, 2], 'o')
    # ax.plot3D([cent[0], cent[0] + n[0]], [cent[1], cent[1] + n[1]], [cent[2], cent[2] + n[2]], 'k-')
    # ax.plot3D([cent[0], cent[0] + q[0]], [cent[1], cent[1] + q[1]], [cent[2], cent[2] + q[2]], 'k:')
    # ax.plot3D([cent[0], cent[0] + A[0]], [cent[1], cent[1] + A[1]], [cent[2], cent[2] + A[2]], 'k--')
    # ax.plot3D(vertices[np.append(order, order[0]), 0],
    #           vertices[np.append(order, order[0]), 1],
    #           vertices[np.append(order, order[0]), 2], 'k-')
    # ax.set_xlabel('x coordinate')
    # ax.set_ylabel('y coordinate')
    # ax.set_zlabel('z coordinate')
    # # xlim = ax.get_xlim()
    # # ylim = ax.get_ylim()
    # # X, Y = np.meshgrid(np.arange(xlim[0], xlim[1], 0.1),
    # #                    np.arange(ylim[0], ylim[1], 0.1))
    # # Z = np.zeros(X.shape)
    # # for r in range(X.shape[0]):
    # #     for c in range(X.shape[1]):
    # #         Z[r, c] = n[0] * X[r, c] + n[1] * Y[r, c] + n[2] * Z[r ,c]
    # # ax.plot_wireframe(X, Y, Z, color='k')
    # plt.show()

    order_all = np.append(order, np.arange(n_vert, 7))
    return order_all

def merge_local_wall_terrain(i, j, k, wall_coord, terr_coord, wall_vert, terr_vert ):
    """ merging faces between terrain and wall """
    wall_coord = wall_coord[:wall_vert, :]
    terr_coord = terr_coord[:terr_vert, :]

    wall_coord_t = wall_coord[wall_coord[:, 2] >= np.max(wall_coord[:, 2])]
    # only_top = True if wall_coord_t[0,2] == wall_coord_t[1,2] else False
    # loop for detection when points are above each other
    trtw = []
    tr_new = []
    for tri, tr in enumerate(terr_coord):
        tr_new.append(tr)
        trtw0 = [tri, -1, False]
        found = False
        for twi, tw in enumerate(wall_coord_t):
            if tw[0] == tr[0] and tw[1] == tr[1]:
                trtw0 = [tri, twi, True]
                trtw.append(trtw0)
                if tw[2] > tr[2]:
                    trtw.append(trtw0)
                    terr_coord[tri, 2] = tw[2]
                    # tr_new.append(tw)
                found = True
                break
        if not found:
            trtw.append(trtw0)

    tr_new.append(tr_new[0])
    vertices = np.asarray(tr_new).astype(np.float)

    n_vert = vertices.shape[0]

    cent = np.mean(vertices, 0)

    # calculate normal vector of some triangles
    vi = 0
    vj = 1
    A = cent - vertices[vi, :]

    svd = np.linalg.svd(vertices.T - np.mean(vertices.T, axis=1, keepdims=True))
    left = svd[0]

    n = left[:, -1] / norm(left[:, -1])

    q = np.cross(A, n)

    angle = np.zeros(n_vert)
    for ir in range(n_vert):
        r = cent - vertices[ir, :]
        t = np.dot(n, np.cross(r, A))
        u = np.dot(n, np.cross(r, q))
        angle[ir] = np.arctan2(u, t) / np.pi * 180
        angle[ir] = angle[ir] if angle[ir] > 0 else 360.0 + angle[ir]

    order = angle.argsort()
    order_all = np.append(order, np.arange(n_vert, 7))

    # fig = plt.figure()
    # from mpl_toolkits import mplot3d
    # ax = plt.axes(projection='3d')
    # plt.title('i,j,k: {},{},{}'.format(i, j, k))
    # ax.plot3D(vertices[:, 0], vertices[:, 1], vertices[:, 2], 'o')
    # ax.plot3D(vertices[:, 0], vertices[:, 1], vertices[:, 2], '-b')
    # ax.plot3D(vertices[np.append(order, order[0]), 0],
    #           vertices[np.append(order, order[0]), 1],
    #           vertices[np.append(order, order[0]), 2], 'k-')
    # ax.set_xlabel('x coordinate')
    # ax.set_ylabel('y coordinate')
    # ax.set_zlabel('z coordinate')
    # plt.show()

    vert_order = np.zeros((8, 3), dtype=object)
    vert_order[:] = None
    vert_order[:n_vert] = vertices[order]

    return vert_order

def initialize_slanted_faces(cfg, connection, cur):
    """ Initialization of slanted faces from terrain, wall, roof faces """
    # create final structure of slanted faces
    progress('Processing merging of final structure')
    sqltext = 'DROP TABLE IF EXISTS "{0}"."{1}"; ' \
              'CREATE TABLE "{0}"."{1}" ( ' \
              ' i integer, j integer, k integer, wid integer,' \
              ' rid integer, lid integer, ' \
              ' iswall boolean, ' \
              ' isroof boolean, ' \
              ' isterr boolean, ' \
              ' n_vert integer, ' \
              ' center geometry(POINTZ, %s), ' \
              ' geom geometry("POLYGONZ", %s),' \
              ' vert1 geometry(POINTZ, %s), ' \
              ' vert2 geometry(POINTZ, %s), ' \
              ' vert3 geometry(POINTZ, %s), ' \
              ' vert4 geometry(POINTZ, %s), ' \
              ' vert5 geometry(POINTZ, %s), ' \
              ' vert6 geometry(POINTZ, %s), ' \
              ' vert7 geometry(POINTZ, %s), ' \
              ' norm geometry("POINTZ", %s), ' \
              ' vert1i integer, vert2i integer, vert3i integer, ' \
              ' vert4i integer, vert5i integer, vert6i integer, vert7i integer  ' \
              '' \
              ')'.format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext, (cfg.srid_palm,cfg.srid_palm,cfg.srid_palm,cfg.srid_palm,cfg.srid_palm,
                          cfg.srid_palm,cfg.srid_palm,cfg.srid_palm,cfg.srid_palm,cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

    if cfg.has_buildings:
        debug('Inserting slanted gridded roof into slanted faces')
        sqltext = 'INSERT INTO "{0}"."{1}" (geom, rid, wid, lid, iswall, isroof, isterr) ' \
                  'SELECT ST_ForceRHR(geom), rid, NULL, NULL, FALSE, TRUE, FALSE ' \
                  'FROM "{0}"."{2}"'\
            .format(cfg.domain.case_schema, cfg.tables.slanted_faces, cfg.tables.slanted_roof_gridded)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        debug('Inserting slanted gridded wall into slanted faces')
        sqltext = 'INSERT INTO "{0}"."{1}" (geom, rid, wid, lid, iswall, isroof, isterr, norm) ' \
                  'SELECT ST_ForceRHR(geom), NULL, wid, NULL, TRUE, FALSE, FALSE, norm  ' \
                  'FROM "{0}"."{2}"'\
            .format(cfg.domain.case_schema, cfg.tables.slanted_faces, cfg.tables.slanted_wall_gridded)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

    debug('Inserting slanted gridded terrain into slanted faces')
    sqltext = 'INSERT INTO "{0}"."{1}" (geom, rid, wid, lid, iswall, isroof, isterr) ' \
              'SELECT ST_ForceRHR(geom), NULL, NULL, lid, FALSE, FALSE, TRUE  ' \
              'FROM "{0}"."{2}"'\
        .format(cfg.domain.case_schema, cfg.tables.slanted_faces, cfg.tables.slanted_terrain_gridded)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    debug('Updating faces centers')
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'center = ST_SetSRID(ST_MakePoint(ST_X(ST_Centroid(geom)), ST_Y(ST_Centroid(geom)), ' \
              '     (COALESCE(ST_Z(ST_PointN(ST_ExteriorRing(geom),1)), 0.0) + COALESCE(ST_Z(ST_PointN(ST_ExteriorRing(geom),2)), 0.0) + ' \
              '      COALESCE(ST_Z(ST_PointN(ST_ExteriorRing(geom),3)), 0.0) + COALESCE(ST_Z(ST_PointN(ST_ExteriorRing(geom),4)), 0.0) + ' \
              '      COALESCE(ST_Z(ST_PointN(ST_ExteriorRing(geom),5)), 0.0) + COALESCE(ST_Z(ST_PointN(ST_ExteriorRing(geom),6)), 0.0) + ' \
              '      COALESCE(ST_Z(ST_PointN(ST_ExteriorRing(geom),7)), 0.0)) / ST_NPoints(geom))' \
              '         , %s)'\
        .format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext, (cfg.srid_palm, ))
    sql_debug(connection)
    connection.commit()

    debug('Calculating face i,j,k')
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'i = FLOOR((ST_X(center) - {3}) / {2}) , ' \
              'j = FLOOR((ST_Y(center) - {4}) / {2}), ' \
              'k = FLOOR(ST_Z(center) / {2})' \
              ''.format(cfg.domain.case_schema, cfg.tables.slanted_faces, cfg.domain.dx,
                        cfg.domain.origin_x, cfg.domain.origin_y)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    debug('Add id into table slanted faces table')
    cur.execute('ALTER TABLE "{0}"."{1}" ADD COLUMN id SERIAL'
                .format(cfg.domain.case_schema, cfg.tables.slanted_faces))
    sql_debug(connection)
    connection.commit()

    verbose('Delete all wall that are under terr_wall faces')
    sqltext = 'WITH ij_terr AS (SELECT i,j,k FROM "{0}"."{1}" WHERE isterr) ' \
              'SELECT s.id ' \
              'FROM "{0}"."{1}" AS s ' \
              'RIGHT JOIN ij_terr AS ts ON ts.i = s.i AND ts.j = s.j ' \
              'WHERE iswall AND s.k < ts.k' \
              .format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext)
    ids = cur.fetchall()
    sql_debug(connection)
    connection.commit()
    ids = [x[0] for x in ids]

    sqltext = 'DELETE FROM "{0}"."{1}" ' \
              'WHERE id IN {2}'\
              .format(cfg.domain.case_schema, cfg.tables.slanted_faces, tuple(ids))
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()


    # verbose('Deleting duplicates')
    # sqltext = 'ALTER TABLE "{0}"."{1}" ADD COLUMN IF NOT EXISTS rijk SERIAL'\
    #           .format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    # cur.execute(sqltext)
    # sql_debug(connection)
    # connection.commit()
    #
    # sqltext = 'DELETE FROM "{0}"."{1}" ' \
    #           'WHERE rijk IN (' \
    #           '               SELECT rijk FROM (' \
    #           '                     SELECT rijk, ROW_NUMBER() OVER (partition BY i,j,k ) AS RowNumber ' \
    #           '                     FROM "{0}"."{1}"' \
    #           '                     WHERE isroof) AS T ' \
    #           '               WHERE T.RowNumber > 1)'\
    #           .format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    # cur.execute(sqltext)
    # sql_debug(connection)
    # connection.commit()

    # sqltext = 'ALTER TABLE "{0}"."{1}" DROP COLUMN rijk'\
    #           .format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    # cur.execute(sqltext)
    # sql_debug(connection)
    # connection.commit()
    #
    #
    # verbose('Deleting duplicates')
    # sqltext = 'ALTER TABLE "{0}"."{1}" ADD COLUMN IF NOT EXISTS rijk SERIAL'\
    #           .format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    # cur.execute(sqltext)
    # sql_debug(connection)
    # connection.commit()
    #
    # sqltext = 'DELETE FROM "{0}"."{1}" ' \
    #           'WHERE rijk IN (' \
    #           '               SELECT rijk FROM (' \
    #           '                     SELECT rijk, ROW_NUMBER() OVER (partition BY i,j,k ) AS RowNumber ' \
    #           '                     FROM "{0}"."{1}"' \
    #           '                     WHERE isterr) AS T ' \
    #           '               WHERE T.RowNumber > 1)'\
    #           .format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    # cur.execute(sqltext)
    # sql_debug(connection)
    # connection.commit()
    #
    # sqltext = 'ALTER TABLE "{0}"."{1}" DROP COLUMN rijk'\
    #           .format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    # cur.execute(sqltext)
    # sql_debug(connection)
    # connection.commit()

    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'vert1 = ST_PointN(ST_Boundary(geom), 1), ' \
              'vert2 = ST_PointN(ST_Boundary(geom), 2), ' \
              'vert3 = ST_PointN(ST_Boundary(geom), 3), ' \
              'vert4 = ST_PointN(ST_Boundary(geom), 4), ' \
              'vert5 = ST_PointN(ST_Boundary(geom), 5), ' \
              'vert6 = ST_PointN(ST_Boundary(geom), 6), ' \
              'vert7 = ST_PointN(ST_Boundary(geom), 7)'\
        .format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    debug('Updating n_vert')
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'n_vert = CASE WHEN vert7 IS NOT NULL THEN 7 ' \
              '              WHEN vert6 IS NOT NULL THEN 6 ' \
              '              WHEN vert5 IS NOT NULL THEN 5 ' \
              '              WHEN vert4 IS NOT NULL THEN 4 ' \
              '              WHEN vert3 IS NOT NULL THEN 3 ' \
              '              WHEN vert2 IS NOT NULL THEN 2 ' \
              '              WHEN vert1 IS NOT NULL THEN 1 END '.format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    debug('Updating all planar-horizontal faces, update their k')
    verbose('Fetching max k from slanted faces')
    sqltext = 'SELECT MAX(k) FROM "{0}"."{1}"'.format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext)
    k_max = cur.fetchone()[0]
    sql_debug(connection)
    connection.commit()

    for k in range(k_max+1):
        verbose('\tDelete all wall that has the same height in all vertices k={}', k)
        sqltext = 'DELETE FROM "{0}"."{1}" ' \
                  'WHERE id IN ' \
                  '  (SELECT id FROM "{0}"."{1}" ' \
                  '   WHERE       (CASE WHEN ST_Z(vert1) = {3} THEN 1 ELSE 0 END + ' \
                  '                CASE WHEN ST_Z(vert2) = {3} THEN 1 ELSE 0 END + ' \
                  '                CASE WHEN ST_Z(vert3) = {3} THEN 1 ELSE 0 END + ' \
                  '                CASE WHEN ST_Z(vert4) = {3} AND n_vert > 4 THEN 1 ELSE 0 END + ' \
                  '                CASE WHEN ST_Z(vert5) = {3} AND n_vert > 5 THEN 1 ELSE 0 END + ' \
                  '                CASE WHEN ST_Z(vert6) = {3} AND n_vert > 6 THEN 1 ELSE 0 END + ' \
                  '                CASE WHEN ST_Z(vert7) = {3} AND n_vert > 7 THEN 1 ELSE 0 END) = n_vert-1 ' \
                  '               AND iswall)'\
                  .format(cfg.domain.case_schema, cfg.tables.slanted_faces, k, (k+1) * cfg.domain.dz)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

    for k in range(k_max+1):
        verbose('\tUpdating k={}', k)
        sqltext = 'SELECT id, i, j, k, n_vert, rid, wid, lid, isterr, iswall, isroof, ' \
                  ' ARRAY[ST_X(vert1), ST_Y(vert1), ST_Z(vert1)],  ' \
                  ' ARRAY[ST_X(vert2), ST_Y(vert2), ST_Z(vert2)], ' \
                  ' ARRAY[ST_X(vert3), ST_Y(vert3), ST_Z(vert3)], ' \
                  ' ARRAY[ST_X(vert4), ST_Y(vert4), ST_Z(vert4)], ' \
                  ' ARRAY[ST_X(vert5), ST_Y(vert5), ST_Z(vert5)], ' \
                  ' ARRAY[ST_X(vert6), ST_Y(vert6), ST_Z(vert6)], ' \
                  ' ARRAY[ST_X(vert7), ST_Y(vert7), ST_Z(vert7)]  ' \
                  'FROM "{0}"."{1}" AS s ' \
                  'WHERE (CASE WHEN ST_Z(vert1) = {3} THEN 1 ELSE 0 END + ' \
                  '       CASE WHEN ST_Z(vert2) = {3} THEN 1 ELSE 0 END + ' \
                  '       CASE WHEN ST_Z(vert3) = {3} THEN 1 ELSE 0 END + ' \
                  '       CASE WHEN ST_Z(vert4) = {3} AND n_vert > 4 THEN 1 ELSE 0 END + ' \
                  '       CASE WHEN ST_Z(vert5) = {3} AND n_vert > 5 THEN 1 ELSE 0 END + ' \
                  '       CASE WHEN ST_Z(vert6) = {3} AND n_vert > 6 THEN 1 ELSE 0 END + ' \
                  '       CASE WHEN ST_Z(vert7) = {3} AND n_vert > 7 THEN 1 ELSE 0 END) > 2 ' \
                  '      AND n_vert > 4 AND k = {2} ' \
                  'ORDER BY i,j,k'\
                  .format(cfg.domain.case_schema, cfg.tables.slanted_faces, k, (k+1) * cfg.domain.dz)
        cur.execute(sqltext)
        verts = cur.fetchall()
        sql_debug(connection)
        connection.commit()

        to_insert = []
        to_delete = []
        for vert in verts:
            np_verts = np.empty((7,3), dtype=object)
            id, i, j, k, n_vert = vert[0], vert[1], vert[2], vert[3], vert[4]
            rid, wid, lid, isterr, iswall, isroof = vert[5], vert[6], vert[7], vert[8], vert[9], vert[10]
            lastidx = 11
            x_vert, y_vert, z_vert = [], [], []
            for idx in range(n_vert-1):
                if vert[lastidx+idx][0] is not None:
                    x_vert.append(vert[lastidx+idx][0])
                    y_vert.append(vert[lastidx+idx][1])
                    z_vert.append(vert[lastidx+idx][2])
            x_vert = np.asarray(x_vert)
            y_vert = np.asarray(y_vert)
            z_vert = np.asarray(z_vert)
            to_keep = []
            tt_max = max(z_vert)
            for it, ttt in enumerate(z_vert):
                # print(ttt)
                if it == 0:
                    ileft = len(z_vert) - 1
                    iright = it + 1
                elif it == len(z_vert) - 1:
                    iright = 0
                    ileft = it - 1
                else:
                    iright = it + 1
                    ileft = it - 1
                # print(ileft, it, iright, )
                # print(z_vert[ileft], z_vert[it], z_vert[iright])
                if z_vert[ileft] == tt_max and z_vert[it] == tt_max and z_vert[iright] == tt_max:
                    # print('drop this one')
                    pass
                else:
                    to_keep.append(it)
            n_vert = len(to_keep) + 1
            np_verts[:n_vert - 1, 0] = x_vert[to_keep]
            np_verts[:n_vert - 1, 1] = y_vert[to_keep]
            np_verts[:n_vert - 1, 2] = z_vert[to_keep]
            np_verts[n_vert - 1, :] = np_verts[0, :]
            to_insert.append((id, i, j, k, rid, wid, lid, iswall, isroof, isterr,
                              np_verts[0, 0], np_verts[0, 1], np_verts[0, 2], cfg.srid_palm,
                              np_verts[1, 0], np_verts[1, 1], np_verts[1, 2], cfg.srid_palm,
                              np_verts[2, 0], np_verts[2, 1], np_verts[2, 2], cfg.srid_palm,
                              np_verts[3, 0], np_verts[3, 1], np_verts[3, 2], cfg.srid_palm,
                              np_verts[4, 0], np_verts[4, 1], np_verts[4, 2], cfg.srid_palm,
                              np_verts[5, 0], np_verts[5, 1], np_verts[5, 2], cfg.srid_palm,
                              np_verts[6, 0], np_verts[6, 1], np_verts[6, 2], cfg.srid_palm,
                              )
                             )
            to_delete.append((id,))

        debug('Deleting all unwanted rows')
        sqltext = 'DELETE FROM "{0}"."{1}" ' \
                  'WHERE id = ANY(%s)'.format(cfg.domain.case_schema, cfg.tables.slanted_faces)
        cur.execute(sqltext, (to_delete, ))
        sql_debug(connection)
        connection.commit()

        debug('Inserting all new entries into slanted faces table')
        sqltext = 'INSERT INTO "{0}"."{1}" (id, i, j, k, rid, wid, lid, iswall, isroof, isterr, ' \
                  '                         vert1, vert2, vert3, vert4, vert5, vert6, vert7) ' \
                  'VALUES                  (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, ' \
                  '                        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
                  '                        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
                  '                        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
                  '                        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
                  '                        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
                  '                        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
                  '                        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s) ' \
                  '                        )  '.format(cfg.domain.case_schema, cfg.tables.slanted_faces)

        cur.executemany(sqltext, to_insert)
        sql_debug(connection)
        connection.commit()

        # Create polygon
        debug('Updating 3d polygon')
        sqltext = 'UPDATE "{0}"."{1}" SET geom =  ' \
                  'ST_ForceRHR(' \
                  'ST_SetSRID(ST_MakePolygon(ST_MakeLine(ARRAY[vert1, vert2, vert3, vert4, vert5, vert6, vert7,' \
                  '     CASE WHEN vert1 IS NOT NULL THEN vert1 ' \
                  '          WHEN vert2 IS NOT NULL THEN vert2 ' \
                  '          WHEN vert3 IS NOT NULL THEN vert3 ' \
                  '          WHEN vert4 IS NOT NULL THEN vert4 ' \
                  '          WHEN vert5 IS NOT NULL THEN vert5 ' \
                  '          WHEN vert6 IS NOT NULL THEN vert6 ' \
                  '          WHEN vert7 IS NOT NULL THEN vert7 END' \
                  '])), %s)) ' \
                  'WHERE geom IS NULL'.format(cfg.domain.case_schema, cfg.tables.slanted_faces)
        cur.execute(sqltext, (cfg.srid_palm,))
        sql_debug(connection)
        connection.commit()


    # debug('Update center of faces, just in case')
    # debug('Updating center')
    # sqltext = 'UPDATE "{0}"."{1}" SET center = ' \
    #           'ST_SetSRID(ST_MakePoint(ST_X(ST_Centroid(geom)), ST_Y(ST_Centroid(geom)), ' \
    #           '                       (ST_Z(vert1) + ST_Z(vert2) + ' \
    #           '                          CASE WHEN ST_NPoints(geom) > 3 THEN ST_Z(ST_PointN(ST_ExteriorRing(geom),3)) ELSE 0.0 END + ' \
    #           '                          CASE WHEN ST_NPoints(geom) > 4 THEN ST_Z(ST_PointN(ST_ExteriorRing(geom),4)) ELSE 0.0 END + ' \
    #           '                          CASE WHEN ST_NPoints(geom) > 5 THEN ST_Z(ST_PointN(ST_ExteriorRing(geom),5)) ELSE 0.0 END + ' \
    #           '                          CASE WHEN ST_NPoints(geom) > 6 THEN ST_Z(ST_PointN(ST_ExteriorRing(geom),6)) ELSE 0.0 END) / (ST_NPoints(geom) - 1))' \
    #           '         , %s)'.format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    # cur.execute(sqltext, (cfg.srid_palm,))
    # sql_debug(connection)
    # connection.commit()

    debug('Updating faces centers')
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'center = ST_SetSRID(ST_MakePoint(ST_X(ST_Centroid(geom)), ST_Y(ST_Centroid(geom)), ' \
              '     (COALESCE(ST_Z(ST_PointN(ST_ExteriorRing(geom),1)), 0.0) + COALESCE(ST_Z(ST_PointN(ST_ExteriorRing(geom),2)), 0.0) + ' \
              '      COALESCE(ST_Z(ST_PointN(ST_ExteriorRing(geom),3)), 0.0) + COALESCE(ST_Z(ST_PointN(ST_ExteriorRing(geom),4)), 0.0) + ' \
              '      COALESCE(ST_Z(ST_PointN(ST_ExteriorRing(geom),5)), 0.0) + COALESCE(ST_Z(ST_PointN(ST_ExteriorRing(geom),6)), 0.0) + ' \
              '      COALESCE(ST_Z(ST_PointN(ST_ExteriorRing(geom),7)), 0.0) + COALESCE(ST_Z(ST_PointN(ST_ExteriorRing(geom),8)), 0.0)) / ST_NPoints(geom))' \
              '         , %s)'\
        .format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext, (cfg.srid_palm, ))
    sql_debug(connection)
    connection.commit()

    debug('Update i,j,k according to new centers')
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'i = FLOOR((ST_X(center) - {3}) / {2}), ' \
              'j = FLOOR((ST_Y(center) - {4}) / {2}), ' \
              'k = FLOOR((ST_Z(center) / {2}) + 1)' \
              ''.format(cfg.domain.case_schema, cfg.tables.slanted_faces, cfg.domain.dx,
                        cfg.domain.origin_x, cfg.domain.origin_y)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    for k in range(k_max + 1):
        verbose('\tUpdating k={} (to k+1) in faces where all vertices lies in the same horizontal plane = k*dz [{}]',
                k, k*cfg.domain.dz)
        sqltext = 'UPDATE "{0}"."{1}" SET k = ' \
                  'k + 1 ' \
                  'WHERE k = {2} AND (CASE WHEN ST_Z(vert1) = {3} THEN 1 ELSE 0 END + ' \
                  '                   CASE WHEN ST_Z(vert2) = {3} THEN 1 ELSE 0 END + ' \
                  '                   CASE WHEN ST_Z(vert3) = {3} THEN 1 ELSE 0 END + ' \
                  '                   CASE WHEN ST_Z(vert4) = {3} AND n_vert > 4 THEN 1 ELSE 0 END + ' \
                  '                   CASE WHEN ST_Z(vert5) = {3} AND n_vert > 5 THEN 1 ELSE 0 END + ' \
                  '                   CASE WHEN ST_Z(vert6) = {3} AND n_vert > 6 THEN 1 ELSE 0 END + ' \
                  '                   CASE WHEN ST_Z(vert7) = {3} AND n_vert > 7 THEN 1 ELSE 0 END) >= n_vert-1'\
                  .format(cfg.domain.case_schema, cfg.tables.slanted_faces, k, k*cfg.domain.dz)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

    # DELETE THE ONES THAT ARE IN GROUP BY, HAS ALL Z vert the same equal (k-1) * dz
    sqltext = 'WITH ijk AS (SELECT i,j,k ' \
          '                 FROM "{0}"."{1}"  ' \
              '             WHERE NOT iswall  ' \
              '             GROUP BY i,j,k ' \
              '             HAVING COUNT(*) = 2 ' \
              '             ORDER BY i,j,k ) ' \
              'SELECT ST_Z(center), id, s.i, s.j, s.k  ' \
              'FROM "{0}"."{1}" AS s ' \
              'RIGHT JOIN ijk ON ijk.i = s.i AND ijk.j = s.j AND ijk.k = s.k ' \
              'WHERE NOT iswall ' \
              'ORDER BY s.i, s.j, s.k'.format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext)
    zijk = cur.fetchall()
    sql_debug(connection)
    connection.commit()

    id2del = []
    for ijk in range(0, len(zijk), 2):
        z1, z2 = zijk[ijk][0], zijk[ijk+1][0]
        id1, id2 = zijk[ijk][1], zijk[ijk+1][1]
        id2del.append((id1 if z1 < z2 else id2,))

    debug('Deleting all unwanted rows')
    sqltext = 'DELETE FROM "{0}"."{1}" ' \
              'WHERE id = ANY(%s)'.format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext, (id2del,))
    sql_debug(connection)
    connection.commit()

def normal_vector_trinagulation(cfg, connection, cur):
    """ calculate face normal vector using triangulation"""
    sqltext  = 'ALTER TABLE "{0}"."{1}" ' \
               'ADD COLUMN IF NOT EXISTS norm_line geometry("LINESTRINGZ", %s),' \
               ' ' \
               'ADD COLUMN IF NOT EXISTS normx double precision, ' \
               'ADD COLUMN IF NOT EXISTS normy double precision,' \
               'ADD COLUMN IF NOT EXISTS normz double precision, ' \
               'ADD COLUMN IF NOT EXISTS area double precision, ' \
               '' \
               'ADD COLUMN IF NOT EXISTS n1x double precision, ' \
               'ADD COLUMN IF NOT EXISTS n1y double precision,' \
               'ADD COLUMN IF NOT EXISTS n1z double precision, ' \
               'ADD COLUMN IF NOT EXISTS area1 double precision, ' \
               '' \
               'ADD COLUMN IF NOT EXISTS n2x double precision, ' \
               'ADD COLUMN IF NOT EXISTS n2y double precision,' \
               'ADD COLUMN IF NOT EXISTS n2z double precision, ' \
               'ADD COLUMN IF NOT EXISTS area2 double precision, ' \
               '' \
               'ADD COLUMN IF NOT EXISTS n3x double precision, ' \
               'ADD COLUMN IF NOT EXISTS n3y double precision,' \
               'ADD COLUMN IF NOT EXISTS n3z double precision, ' \
               'ADD COLUMN IF NOT EXISTS area3 double precision, ' \
               '' \
               'ADD COLUMN IF NOT EXISTS n4x double precision, ' \
               'ADD COLUMN IF NOT EXISTS n4y double precision,' \
               'ADD COLUMN IF NOT EXISTS n4z double precision, ' \
               'ADD COLUMN IF NOT EXISTS area4 double precision, ' \
               '' \
               'ADD COLUMN IF NOT EXISTS n5x double precision, ' \
               'ADD COLUMN IF NOT EXISTS n5y double precision,' \
               'ADD COLUMN IF NOT EXISTS n5z double precision, ' \
               'ADD COLUMN IF NOT EXISTS area5 double precision, ' \
               '' \
               'ADD COLUMN IF NOT EXISTS n6x double precision, ' \
               'ADD COLUMN IF NOT EXISTS n6y double precision,' \
               'ADD COLUMN IF NOT EXISTS n6z double precision, ' \
               'ADD COLUMN IF NOT EXISTS area6 double precision, ' \
               '' \
               'ADD COLUMN IF NOT EXISTS n7x double precision, ' \
               'ADD COLUMN IF NOT EXISTS n7y double precision,' \
               'ADD COLUMN IF NOT EXISTS n7z double precision, ' \
               'ADD COLUMN IF NOT EXISTS area7 double precision ' \
               ''.format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext, (cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

    for ni in range(1,7):
        debug('Trinagulating normal vector, {}', ni)
        sqltext = 'UPDATE "{0}"."{1}" SET (n{2}x, n{2}y, n{2}z) = (' \
                  '' \
                  'CASE WHEN {2} < n_vert THEN ((ST_Y(vert{3})-ST_Y(center))*(ST_Z(vert{2})-ST_Z(center)) - (ST_Z(vert{3})-ST_Z(center))*(ST_Y(vert{2})-ST_Y(center)))' \
                  ' ELSE NULL END, ' \
                  'CASE WHEN {2} < n_vert THEN ((ST_Z(vert{3})-ST_Z(center))*(ST_X(vert{2})-ST_X(center)) - (ST_X(vert{3})-ST_X(center))*(ST_Z(vert{2})-ST_Z(center)))' \
                  ' ELSE NULL END, ' \
                  'CASE WHEN {2} < n_vert THEN ((ST_X(vert{3})-ST_X(center))*(ST_Y(vert{2})-ST_Y(center)) - (ST_Y(vert{3})-ST_Y(center))*(ST_X(vert{2})-ST_X(center)))' \
                  ' ELSE NULL END ' \
                  ')' \
                  ''.format(cfg.domain.case_schema, cfg.tables.slanted_faces, ni, ni+1)
        cur.execute(sqltext, (cfg.srid_palm,))
        sql_debug(connection)
        connection.commit()

        sqltext = 'UPDATE "{0}"."{1}" SET area{2} = ' \
                  'CASE WHEN {2} < n_vert THEN SQRT((n{2}x)^2 + (n{2}y)^2 + (n{2}z)^2) / 2.0' \
                  ' ELSE NULL END'.format(cfg.domain.case_schema, cfg.tables.slanted_faces, ni)
        cur.execute(sqltext, (cfg.srid_palm,))
        sql_debug(connection)
        connection.commit()

        # normalize vector
        # sqltext = 'UPDATE "{0}"."{1}" SET (n{2}x, n{2}y, n{2}z) = ' \
        #           ' (CASE WHEN {2} < n_vert THEN n{2}x/SQRT(n{2}x^2 + n{2}y^2 + n{2}z^2) ELSE NULL END, ' \
        #           '  CASE WHEN {2} < n_vert THEN n{2}y/SQRT(n{2}x^2 + n{2}y^2 + n{2}z^2) ELSE NULL END, ' \
        #           '  CASE WHEN {2} < n_vert THEN n{2}z/SQRT(n{2}x^2 + n{2}y^2 + n{2}z^2) ELSE NULL END)'.format(cfg.domain.case_schema, cfg.tables.slanted_faces, ni)
        # cur.execute(sqltext)
        # sql_debug(connection)
        # connection.commit()


    # get overall area
    debug('Calculation of overall face area')
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'area = area1 + COALESCE(area2, 0.0) + COALESCE(area3, 0.0) + COALESCE(area4, 0.0) + ' \
              '               COALESCE(area5, 0.0) + COALESCE(area6, 0.0) + COALESCE(area7, 0.0)' \
              .format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # get average normal vector
    debug('Normal vector for faces with zero area')
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
             'normx = 0, ' \
             'normy = 0, ' \
             'normz = 1 ' \
             'WHERE area = 0.0 AND norm IS NULL'\
             .format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    debug('Calculating average unit normal vector')
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'normx = (n1x*area1 + COALESCE(n2x, 0.0)*COALESCE(area2, 0.0) + COALESCE(n3x, 0.0)*COALESCE(area3, 0.0) + ' \
              '                     COALESCE(n4x, 0.0)*COALESCE(area4, 0.0) + COALESCE(n5x, 0.0)*COALESCE(area5, 0.0) + ' \
              '                     COALESCE(n6x, 0.0)*COALESCE(area6, 0.0) + COALESCE(n7x, 0.0)*COALESCE(area7, 0.0) ' \
              '         )/area, ' \
              'normy = (n1y*area1 + COALESCE(n2y, 0.0)*COALESCE(area2, 0.0) + COALESCE(n3y, 0.0)*COALESCE(area3, 0.0) + ' \
              '                     COALESCE(n4y, 0.0)*COALESCE(area4, 0.0) + COALESCE(n5y, 0.0)*COALESCE(area5, 0.0) + ' \
              '                     COALESCE(n6y, 0.0)*COALESCE(area6, 0.0) + COALESCE(n7y, 0.0)*COALESCE(area7, 0.0) ' \
              '         )/area,' \
              'normz = (n1z*area1 + COALESCE(n2z, 0.0)*COALESCE(area2, 0.0) + COALESCE(n3z, 0.0)*COALESCE(area3, 0.0) + ' \
              '                     COALESCE(n4z, 0.0)*COALESCE(area4, 0.0) + COALESCE(n5z, 0.0)*COALESCE(area5, 0.0) + ' \
              '                     COALESCE(n6z, 0.0)*COALESCE(area6, 0.0) + COALESCE(n7z, 0.0)*COALESCE(area7, 0.0) ' \
              '         )/area ' \
              'WHERE norm IS NULL AND area > 0.0' \
              ''.format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    verbose('Calculate normals from vertical walls norm')
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'normx = ST_X(norm) - ST_X(center), ' \
              'normy = ST_Y(norm) - ST_Y(center), ' \
              'normz = ST_Z(norm) - ST_Z(center) ' \
              'WHERE norm IS NOT NULL'\
              .format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    debug('Normal vector for faces with zero area')
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
             'normx = 0, ' \
             'normy = 0, ' \
             'normz = 1 ' \
             'WHERE norm IS NULL AND normx^2 + normy^2 + normz^2 = 0'\
             .format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    verbose('Correct normal vector, z component for vertical faces')
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'normz = 0.0' \
              'WHERE ABS(normz/area) < 1e-5 AND area > 0.0'\
              .format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()


    # normalize vector
    debug('Normalizing normal vector')
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'normx = normx / SQRT(normx^2 + normy^2 + normz^2), ' \
              'normy = normy / SQRT(normx^2 + normy^2 + normz^2), ' \
              'normz = normz / SQRT(normx^2 + normy^2 + normz^2)' \
              ''.format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # DELETE unnecessary column
    # if cfg.slanted_pars.clean_up:
    # debug('Deleting extra columns')
    # sqltext = 'ALTER TABLE "{0}"."{1}" ' \
    #           '{2} n1x, {2} n2x, {2} n3x, {2} n4x, {2} n5x, {2} n6x, {2} n7x, ' \
    #           '{2} n1y, {2} n2y, {2} n3y, {2} n4y, {2} n5y, {2} n6y, {2} n7y, ' \
    #           '{2} n1z, {2} n2z, {2} n3z, {2} n4z, {2} n5z, {2} n6z, {2} n7z, ' \
    #           '{2} area1, {2} area2, {2} area3, {2} area4, {2} area5, ' \
    #           '{2} area6, {2} area7'\
    #           .format(cfg.domain.case_schema, cfg.tables.slanted_faces, 'DROP COLUMN IF EXISTS')
    # cur.execute(sqltext)
    # sql_debug(connection)
    # connection.commit()

    # Create point and line for normal vector
    debug('Calculation of normal point')
    sqltext = 'UPDATE "{0}"."{1}" SET norm = ' \
              'ST_SetSRID(ST_MakePoint(ST_X(center)+normx, ST_Y(center)+normy, ST_Z(center)+normz),%s) ' \
              'WHERE norm IS NULL'\
              .format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext, (cfg.srid_palm, ))
    sql_debug(connection)
    connection.commit()

    debug('Calculation of normal line')
    sqltext = 'UPDATE "{0}"."{1}" SET norm_line = ' \
              'ST_SetSRID(ST_MakeLine(center, norm), %s)'\
              .format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext, (cfg.srid_palm, ))
    sql_debug(connection)
    connection.commit()

def create_integer_vertices(cfg, connection, cur):
    """ create coordinates of vertices using 4 integer and 1 float
        i,j,k position of the line where vertex sits
        dir direction in which inside orthogonal face vertex is placed
            0 +k, 1 -k, 2 +j, 3 -j, 4 +i, 5 -i
        length normalized lenght between inside point and slanted face vertex
    """
    progress('Creating new coordinates of the slanted faces')
    for ni in range(1,8):
        debug('Processing {}. coordinates', ni)
        verbose('\tAdding new columns')
        sqltext = 'ALTER TABLE "{0}"."{1}" ' \
                  'ADD COLUMN IF NOT EXISTS ii{2} integer, ' \
                  'ADD COLUMN IF NOT EXISTS jj{2} integer, ' \
                  'ADD COLUMN IF NOT EXISTS kk{2} integer, ' \
                  'ADD COLUMN IF NOT EXISTS dir{2} integer, ' \
                  'ADD COLUMN IF NOT EXISTS len{2} double precision, ' \
                  'ADD COLUMN IF NOT EXISTS iline boolean, ' \
                  'ADD COLUMN IF NOT EXISTS jline boolean, ' \
                  'ADD COLUMN IF NOT EXISTS kline boolean,' \
                  'ADD COLUMN IF NOT EXISTS nx{2}_l double precision, ' \
                  'ADD COLUMN IF NOT EXISTS ny{2}_l double precision, ' \
                  'ADD COLUMN IF NOT EXISTS nz{2}_l double precision ' \
                  .format(cfg.domain.case_schema, cfg.tables.slanted_faces, ni)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        debug('Calculate for each vertex its own normal vector from adjacent vertices')
        # TODO: optimize
        sqltext = 'WITH c AS ( ' \
                  '   SELECT '\
                  '       vert{2} AS vert_m, '\
                  '       CASE WHEN ({2} + 1) < n_vert THEN vert{3} ELSE vert1 END AS vert_r, '\
                  '       CASE WHEN ({2} - 1) > 0 THEN vert{4} ELSE ' \
                  '            CASE WHEN n_vert = 4 THEN vert3 WHEN n_vert = 5 THEN vert4 WHEN n_vert = 6 THEN vert5 WHEN n_vert = 7 THEN vert6  END '\
                  '       END AS vert_l, '\
                  '       id '\
                  '    FROM "{0}"."{1}") ' \
                  'UPDATE "{0}"."{1}" AS s SET ( nx{2}_l, ny{2}_l, nz{2}_l ) = ' \
                  '(CASE WHEN {2} < n_vert THEN ((ST_Y(vert_l)-ST_Y(vert_m))*(ST_Z(vert_r)-ST_Z(vert_m)) - (ST_Z(vert_l)-ST_Z(vert_m))*(ST_Y(vert_r)-ST_Y(vert_m)))' \
                  '  ELSE NULL END, '\
                   'CASE WHEN {2} < n_vert THEN ((ST_Z(vert_l)-ST_Z(vert_m))*(ST_X(vert_r)-ST_X(vert_m)) - (ST_X(vert_l)-ST_X(vert_m))*(ST_Z(vert_r)-ST_Z(vert_m)))'\
                  '  ELSE NULL END, '\
                  ' CASE WHEN {2} < n_vert THEN ((ST_X(vert_l)-ST_X(vert_m))*(ST_Y(vert_r)-ST_Y(vert_m)) - (ST_Y(vert_l)-ST_Y(vert_m))*(ST_X(vert_r)-ST_X(vert_m)))'\
                  '  ELSE NULL END) '\
                  'FROM c WHERE c.id = s.id AND NOT iswall'\
                   .format(cfg.domain.case_schema, cfg.tables.slanted_faces, ni, min([ni + 1, 7]), max([ni - 1, 1]) )
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        sqltext = ('UPDATE "{0}"."{1}" AS s SET ( nx{2}_l, ny{2}_l, nz{2}_l ) = '
                   '(nx1_l, ny1_l, nz1_l) '
                   'WHERE {2} = n_vert '
                   .format(cfg.domain.case_schema, cfg.tables.slanted_faces, ni))
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        sqltext = ('UPDATE "{0}"."{1}" AS s SET ( nx{2}_l, ny{2}_l, nz{2}_l ) = '
                   '(normx, normy, normz) '
                   'WHERE iswall '
                   .format(cfg.domain.case_schema, cfg.tables.slanted_faces, ni))
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        verbose('\tUpdating ii, jj, kk')
        sqltext = 'UPDATE "{0}"."{1}" SET ' \
                  'ii{2} = FLOOR((ST_X(vert{2}) - {6}) / {3}), ' \
                  'jj{2} = FLOOR((ST_Y(vert{2}) - {7}) / {4}), ' \
                  'kk{2} = FLOOR(ST_Z(vert{2}) / {5}) '\
                  .format(cfg.domain.case_schema, cfg.tables.slanted_faces, ni,
                          cfg.domain.dx, cfg.domain.dy, cfg.domain.dz,
                          cfg.domain.origin_x, cfg.domain.origin_y)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        verbose('\tUpdating iline, jline, kline')
        sqltext = 'UPDATE "{0}"."{1}" SET ' \
                  'iline = CASE WHEN ABS(ST_X(vert{2}) - ii{2} * {3} - {6}) > 1e-6 THEN TRUE ELSE FALSE END, ' \
                  'jline = CASE WHEN ABS(ST_Y(vert{2}) - jj{2} * {4} - {7}) > 1e-6 THEN TRUE ELSE FALSE END, ' \
                  'kline = CASE WHEN ABS(ST_Z(vert{2}) - kk{2} * {5}) > 1e-6       THEN TRUE ELSE FALSE END '\
                  .format(cfg.domain.case_schema, cfg.tables.slanted_faces, ni,
                          cfg.domain.dx, cfg.domain.dy, cfg.domain.dz,
                          cfg.domain.origin_x, cfg.domain.origin_y)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        verbose('\tUpdating dir')
        sqltext = 'UPDATE "{0}"."{1}" SET dir{2} = ' \
                  'CASE WHEN iline THEN CASE WHEN nx{2}_l > 0.0 THEN 4 ' \
                  '                          WHEN nx{2}_l < 0.0 THEN 5 ' \
                  '                          ELSE 6 END ' \
                  '     WHEN jline THEN CASE WHEN ny{2}_l > 0.0 THEN 2 ' \
                  '                          WHEN ny{2}_l < 0.0 THEN 3' \
                  '                          ELSE 6 END ' \
                  '     WHEN kline THEN CASE WHEN nz{2}_l > 0.0 THEN 0 ' \
                  '                          WHEN nz{2}_l < 0.0 THEN 1 ' \
                  '                          ELSE 6 END ' \
                  '     ELSE 6 END'\
                  .format(cfg.domain.case_schema, cfg.tables.slanted_faces, ni)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        verbose('\tUpdating ii,jj,kk')
        sqltext = 'UPDATE "{0}"."{1}" SET ' \
                  'ii{2} = CASE WHEN dir{2} = 5 THEN ii{2}+1 ELSE ii{2} END, ' \
                  'jj{2} = CASE WHEN dir{2} = 3 THEN jj{2}+1 ELSE jj{2} END, ' \
                  'kk{2} = CASE WHEN dir{2} = 1 THEN kk{2}+1 ELSE kk{2} END ' \
            .format(cfg.domain.case_schema, cfg.tables.slanted_faces, ni)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

        verbose('\tUpdating len')
        sqltext = 'UPDATE "{0}"."{1}" SET len{2} = ' \
                  'SQRT(' \
                  '     ((ST_X(vert{2}) - {6} - ii{2} * {3}) / {3} ) ^ 2 + ' \
                  '     ((ST_Y(vert{2}) - {7} - jj{2} * {4}) / {4} ) ^ 2 + ' \
                  '     ((ST_Z(vert{2})       - kk{2} * {5}) / {5} ) ^ 2' \
                  '     )' \
                  .format(cfg.domain.case_schema, cfg.tables.slanted_faces, ni,
                          cfg.domain.dx, cfg.domain.dy, cfg.domain.dz,
                          cfg.domain.origin_x, cfg.domain.origin_y)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()
        # 'CASE WHEN dir{2} = 1 THEN (ST_X(vert{2}) - {6} - ii{2} * {3}) / {3} ' \
        # '     WHEN dir{2} = 0 THEN (ST_X(vert{2}) - {6} - ii{2} * {3}) / {3}' \
        # '     WHEN dir{2} = 3 THEN (ST_Y(vert{2}) - {7} - jj{2} * {4}) / {4} ' \
        # '     WHEN dir{2} = 2 THEN (ST_Y(vert{2}) - {7} - jj{2} * {4}) / {4}' \
        # '     WHEN dir{2} = 5 THEN (ST_Z(vert{2})       - kk{2} * {5}) / {5} ' \
        # '     WHEN dir{2} = 4 THEN (ST_Z(vert{2})       - kk{2} * {5}) / {6}' \
        # 'END '\

def create_vertices_indexes(cfg, connection, cur):
    """ join all vertices and create back index to faces """
    progress('Creating vertices table and back indexing')
    sqltext = 'DROP TABLE IF EXISTS "{0}"."{1}" CASCADE '.format(cfg.domain.case_schema, cfg.tables.vertices)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    debug('Create table of individual vertices 2, based on kjidir index')
    sqltext = 'CREATE TABLE "{0}"."{1}" AS ' \
              'SELECT ROW_NUMBER() OVER () AS id, k, j, i, dir, len, point ' \
              'FROM (' \
              '      SELECT kk1 AS k, jj1 AS j, ii1 AS i, dir1 AS dir, len1 AS len, vert1 AS point FROM "{0}"."{2}" AS v2_1 WHERE ii1 IS NOT NULL' \
              '      UNION ALL ' \
              '      SELECT kk2 AS k, jj2 AS j, ii2 AS i, dir2 AS dir, len2 AS len, vert2 AS point FROM "{0}"."{2}" AS v2_2 WHERE ii2 IS NOT NULL' \
              '      UNION ALL ' \
              '      SELECT kk3 AS k, jj3 AS j, ii3 AS i, dir3 AS dir, len3 AS len, vert3 AS point FROM "{0}"."{2}" AS v2_3 WHERE ii3 IS NOT NULL' \
              '      UNION ALL ' \
              '      SELECT kk4 AS k, jj4 AS j, ii4 AS i, dir4 AS dir, len4 AS len, vert4 AS point FROM "{0}"."{2}" AS v2_4 WHERE ii4 IS NOT NULL' \
              '      UNION ALL ' \
              '      SELECT kk5 AS k, jj5 AS j, ii5 AS i, dir5 AS dir, len5 AS len, vert5 AS point FROM "{0}"."{2}" AS v2_5 WHERE ii5 IS NOT NULL' \
              '      UNION ALL ' \
              '      SELECT kk6 AS k, jj6 AS j, ii6 AS i, dir6 AS dir, len6 AS len, vert6 AS point FROM "{0}"."{2}" AS v2_6 WHERE ii6 IS NOT NULL' \
              '      UNION ALL ' \
              '      SELECT kk7 AS k, jj7 AS j, ii7 AS i, dir7 AS dir, len7 AS len, vert7 AS point FROM "{0}"."{2}" AS v2_7 WHERE ii7 IS NOT NULL' \
              '      ) AS s ' \
              'GROUP BY k,j,i,dir,len,point'.format(cfg.domain.case_schema, cfg.tables.vertices, cfg.tables.slanted_faces)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    verbose('\tAdding kjidir index')
    sqltext = 'CREATE index IF NOT EXISTS vert2_kjidir_idx ON "{0}"."{1}" (k,j,i,dir)' \
        .format(cfg.domain.case_schema, cfg.tables.vertices)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    for i in range(1, 8):
        verbose('Creating kjidir {} index on slanted faces table', i)
        sqltext = 'CREATE INDEX IF NOT EXISTS vert2_kjidur_{2}_geom_idx ON "{0}"."{1}" (kk{2},jj{2},ii{2},dir{2})' \
            .format(cfg.domain.case_schema, cfg.tables.slanted_faces, i)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

    verbose('\tAdding integer index')
    sqltext = 'CREATE index IF NOT EXISTS index_id ON "{0}"."{1}" (id)' \
        .format(cfg.domain.case_schema, cfg.tables.vertices)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # join vert with vertices table
    for i in range(1, 8):
        debug('Joining vertices n. {}', i)
        sqltext = 'UPDATE "{0}"."{2}" SET ' \
                  'vert{3}i = v.id FROM "{0}"."{1}" AS v ' \
                  'WHERE ii{3} = v.i AND jj{3} = v.j AND kk{3} = v.k AND dir{3} = v.dir'\
                  .format(cfg.domain.case_schema, cfg.tables.vertices, cfg.tables.slanted_faces, i)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()
        debug('End joining vertices n. {}', i)

    # ADD ID column
    debug('Renewing id column on slanted faces table')
    sqltext = 'ALTER TABLE "{0}"."{1}" DROP COLUMN IF EXISTS id; ' \
              'ALTER TABLE "{0}"."{1}" ADD COLUMN id SERIAL'.format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()


    # debug('Creating vertices tables from individual vertices')
    # sqltext = 'CREATE TABLE "{0}"."{1}" AS ' \
    #             'SELECT ROW_NUMBER() OVER () AS id, point  ' \
    #             'FROM ( ' \
    #             '    SELECT v1.vert1 AS point FROM "{0}"."{2}" AS v1  WHERE v1.vert1 IS NOT NULL ' \
    #             '    UNION ALL ' \
    #             '    SELECT v2.vert2 AS point FROM "{0}"."{2}" AS v2  WHERE v2.vert2 IS NOT NULL ' \
    #             '    UNION ALL ' \
    #             '    SELECT v3.vert3 AS point FROM "{0}"."{2}" AS v3  WHERE v3.vert3 IS NOT NULL ' \
    #             '    UNION ALL ' \
    #             '    SELECT v4.vert4 AS point FROM "{0}"."{2}" AS v4  WHERE v4.vert4 IS NOT NULL ' \
    #             '    UNION ALL ' \
    #             '    SELECT v5.vert5 AS point FROM "{0}"."{2}" AS v5  WHERE v5.vert5 IS NOT NULL ' \
    #             '    UNION ALL ' \
    #             '    SELECT v6.vert6 AS point FROM "{0}"."{2}" AS v6  WHERE v6.vert6 IS NOT NULL ' \
    #             '    UNION ALL ' \
    #             '    SELECT v7.vert7 AS point FROM "{0}"."{2}" AS v7  WHERE v7.vert7 IS NOT NULL ' \
    #             ') AS s ' \
    #             'GROUP BY point'.format(cfg.domain.case_schema, cfg.tables.vertices, cfg.tables.slanted_faces)
    # cur.execute(sqltext)
    # sql_debug(connection)
    # connection.commit()

    # verbose('Add extra columns for new coordinate system')
    # sqltext = 'ALTER TABLE "{0}"."{1}" ' \
    #           'ADD COLUMN IF NOT EXISTS ii integer, ' \
    #           'ADD COLUMN IF NOT EXISTS jj integer, ' \
    #           'ADD COLUMN IF NOT EXISTS kk integer, ' \
    #           'ADD COLUMN IF NOT EXISTS dir integer, ' \
    #           'ADD COLUMN IF NOT EXISTS len double precision '\
    #           .format(cfg.domain.case_schema, cfg.tables.vertices)
    # cur.execute(sqltext)
    # sql_debug(connection)
    # connection.commit()

    # add geom indexes
    # TODO: Merge point based on ii,jj,kk,dir INDEX
    # verbose('\tAdding HASH indexes')
    # sqltext = 'CREATE index IF NOT EXISTS hash_geom_idx ON "{0}"."{1}" USING HASH(point)' \
    #     .format(cfg.domain.case_schema, cfg.tables.vertices)
    # cur.execute(sqltext)
    # sql_debug(connection)
    # connection.commit()


    # for i in range(1, 8):
    #     sqltext = 'CREATE INDEX vert_hash_{2}_geom_idx ON "{0}"."{1}" USING HASH(vert{2})' \
    #         .format(cfg.domain.case_schema, cfg.tables.slanted_faces, i)
    #     cur.execute(sqltext)
    #     sql_debug(connection)
    #     connection.commit()
    #
    # verbose('\tAdding integer index')
    # sqltext = 'CREATE index IF NOT EXISTS index_id ON "{0}"."{1}" (id)' \
    #     .format(cfg.domain.case_schema, cfg.tables.vertices)
    # cur.execute(sqltext)
    # sql_debug(connection)
    # connection.commit()

    # for i in range(1, 8):
    #     sqltext = 'CREATE index IF NOT EXISTS index_vert_{2} ON "{0}"."{1}" (vert{2}i)' \
    #         .format(cfg.domain.case_schema, cfg.tables.slanted_faces, i)
    #     cur.execute(sqltext)
    #     sql_debug(connection)
    #     connection.commit()
    #
    # # join vert with vertices table
    # for i in range(1, 8):
    #     debug('Joining vertices n. {}', i)
    #     sqltext = 'UPDATE "{0}"."{2}" SET ' \
    #               'vert{3}i = v.id FROM "{0}"."{1}" AS v ' \
    #               'WHERE vert{3} = v.point'.format(cfg.domain.case_schema, cfg.tables.vertices, cfg.tables.slanted_faces, i)
    #     cur.execute(sqltext)
    #     sql_debug(connection)
    #     connection.commit()
    #     debug('End joining vertices n. {}', i)
    #
    # debug('Updating ii, jj, kk, dir, len in table: {}', cfg.tables.vertices)
    # for i in range(1, 8):
    #     verbose('\tProcessing n. {}', i)
    #     sqltext = 'UPDATE "{0}"."{1}" AS v SET ' \
    #               '(ii, jj, kk, dir, len) = ' \
    #               '(SELECT ii{3}, jj{3}, kk{3}, dir{3}, len{3} ' \
    #               ' FROM "{0}"."{2}" AS s ' \
    #               ' WHERE v.id = s.vert{3}i ' \
    #               ' LIMIT 1' \
    #               ') ' \
    #               'WHERE ii IS NULL'\
    #               .format(cfg.domain.case_schema, cfg.tables.vertices, cfg.tables.slanted_faces, i)
    #     cur.execute(sqltext)
    #     sql_debug(connection)
    #     connection.commit()
    #
    # # ADD ID column
    # sqltext = 'ALTER TABLE "{0}"."{1}" DROP COLUMN IF EXISTS id; ' \
    #           'ALTER TABLE "{0}"."{1}" ADD COLUMN id SERIAL'.format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    # cur.execute(sqltext)
    # sql_debug(connection)
    # connection.commit()

def check_for_vertex_singularities(cfg, connection, cur):
    """ Check if there some vertices that lies in gridboxes corners and can have multiple vertices

        Each interface Air / Solid must have 1 Vertex, In case of singularity, there will be two point with the same z,y,x (kk,jj,ii) but different dir
    """
    # dirs = np.array([[[ 0, 1, 0], [ 1, 0, 0], [ 0, 0, 1]],
    #                  [[ 0, 0,-1], [ 1, 0, 0], [ 0, 1, 0]],
    #                  [[ 0,-1, 0], [ 1, 0, 0], [ 0, 0,-1]],
    #                  [[ 0, 0, 1], [ 1, 0, 0], [ 0,-1, 0]],
    #                  [[ 0, 1, 0], [-1, 0, 0], [ 0, 0, 1]],
    #                  [[ 0, 0,-1], [-1, 0, 0], [ 0, 1, 0]],
    #                  [[ 0,-1, 0], [-1, 0, 0], [ 0, 0,-1]],
    #                  [[ 0, 0, 1], [-1, 0, 0], [ 0,-1, 0]],
    #                  ])
    # exit(1)
    dirs_vect = np.array([[ 1,  0,  1],
                          [-1,  0, -1],
                          [ 0,  1,  0],
                          [ 0, -1,  0],
                          [ 0,  0,  1],
                          [ 0,  0, -1]])

    vert2dirs = np.array([[4, 0, 2],
                          [2, 0, 5],
                          [5, 0, 3],
                          [3, 0, 4],
                          [4, 1, 2],
                          [2, 1, 5],
                          [5, 1, 3],
                          [3, 1, 4],])

    dir_corns = np.array([[1, 4, 3],
                          [2, 5, 0],
                          [3, 6, 1],
                          [0, 7, 2],
                          [5, 0, 7],
                          [6, 1, 4],
                          [7, 2, 5],
                          [4, 3, 6],
                          ])

    corners = np.array([[ 0, 0, 0],
                        [ 0, 0, 1],
                        [ 0, 1, 1],
                        [ 0, 1, 0],
                        [ 1, 0, 0],
                        [ 1, 0, 1],
                        [ 1, 1, 1],
                        [ 1, 1, 0]])
    progress('Checking for singularities')
    debug('Identification of polygons with singularities')
    sqltext = 'SELECT id, k, j, i, n_vert, ' \
              '       wid, rid, lid, isterr, iswall, isroof, ' \
              '       len1, len2, len3, len4, len5, len6, len7, ' \
              '       dir1, dir2, dir3, dir4, dir5, dir6, dir7, ' \
              '       normz, normy, normx, area, center, ' \
              '       ii1, ii2, ii3, ii4, ii5, ii6, ii7, ' \
              '       jj1, jj2, jj3, jj4, jj5, jj6, jj7, ' \
              '       kk1, kk2, kk3, kk4, kk5, kk6, kk7, ' \
              'ARRAY[ST_Z(vert1), ST_Y(vert1), ST_X(vert1)], ' \
              'ARRAY[ST_Z(vert2), ST_Y(vert2), ST_X(vert2)], ' \
              'ARRAY[ST_Z(vert3), ST_Y(vert3), ST_X(vert3)], ' \
              'ARRAY[ST_Z(vert4), ST_Y(vert4), ST_X(vert4)], ' \
              'ARRAY[ST_Z(vert5), ST_Y(vert5), ST_X(vert5)], ' \
              'ARRAY[ST_Z(vert6), ST_Y(vert6), ST_X(vert6)], ' \
              'ARRAY[ST_Z(vert7), ST_Y(vert7), ST_X(vert7)] ' \
              ' ' \
              'FROM "{0}"."{1}" ' \
              'WHERE (len1=0 OR len2=0 OR len3=0 OR len4=0 OR len5=0 OR len6=0 OR' \
              '       dir1=6 OR dir2=6 OR dir3=6 OR (dir4=6 AND ii4 IS NOT NULL) OR (dir5=6 AND ii5 IS NOT NULL) OR (dir6=6 AND ii6 IS NOT NULL)) ' \
              ''\
              .format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext)
    singulars = cur.fetchall()
    sql_debug(connection)
    connection.commit()

    to_delete = []
    to_insert = []

    for singular in singulars:
        ids, k, j, i, n_vert = singular[0], singular[1], singular[2], singular[3], singular[4]
        # if k == 15 and j == 19 and i == 25: # [25,19,15]
        #     break
        verbose('Processing singular id: {}', ids)
        wid, rid, lid, isterr, iswall, isroof = singular[5], singular[6], singular[7],  singular[8],  singular[9],  singular[10]
        lens = [irun for irun in singular[11:18]]
        dirs = [irun for irun in singular[18:25]]
        norm = [singular[25], singular[26], singular[27]]
        area = singular[28]
        center = singular[29]
        ii = [irun for irun in singular[30:37]]
        jj = [irun for irun in singular[37:44]]
        kk = [irun for irun in singular[44:51]]
        lastidx = 51
        x_vert, y_vert, z_vert = [], [], []
        for idx in range(n_vert):
            if singular[lastidx + idx][0] is not None:
                z_vert.append(singular[lastidx + idx][0])
                y_vert.append(singular[lastidx + idx][1])
                x_vert.append(singular[lastidx + idx][2])

        x_vert_n, y_vert_n, z_vert_n, ii_n, jj_n, kk_n, dirs_n, lens_n = [], [], [], [], [], [], [], []
        for irun in range(len(x_vert) - 1):
            if x_vert[irun] == x_vert[irun + 1] and y_vert[irun] == y_vert[irun + 1] and z_vert[irun] == z_vert[irun + 1]:
                verbose('duplicate {}', irun)
            else:
                x_vert_n.append(x_vert[irun])
                y_vert_n.append(y_vert[irun])
                z_vert_n.append(z_vert[irun])
                ii_n.append(ii[irun])
                jj_n.append(jj[irun])
                kk_n.append(kk[irun])
                dirs_n.append(dirs[irun])
                lens_n.append(lens[irun])

        irun = len(x_vert)-1
        if not (x_vert[irun] == x_vert[0] and y_vert[irun] == y_vert[0] and z_vert[irun] == z_vert[0]):
            x_vert_n.append(x_vert[irun])
            y_vert_n.append(y_vert[irun])
            z_vert_n.append(z_vert[irun])
            ii_n.append(ii[irun])
            jj_n.append(jj[irun])
            kk_n.append(kk[irun])
            dirs_n.append(dirs[irun])
            lens_n.append(lens[irun])

        x_vert, y_vert, z_vert, ii, jj, kk, dirs, lens = x_vert_n, y_vert_n, z_vert_n, ii_n, jj_n, kk_n, dirs_n, lens_n

        x_vert = np.asarray(x_vert)
        y_vert = np.asarray(y_vert)
        z_vert = np.asarray(z_vert)
        # ii = ((x_vert - cfg.domain.origin_x) / cfg.domain.dx).astype(np.int)
        # jj = ((y_vert - cfg.domain.origin_y) / cfg.domain.dy).astype(np.int)
        # kk = ((z_vert) / cfg.domain.dz).astype(np.int) + 1   # because k = 0 is flat ground
        k_temp = k - 1
        corner_id = -1 * np.ones(7, dtype=np.int)
        for idx in range(len(x_vert)):
            if   (ii[idx] == i)     & (jj[idx] == j)     & (kk[idx] == k_temp)     & (lens[idx] == 0):
                corner_id[idx] = 0
            elif (ii[idx] == i + 1) & (jj[idx] == j)     & (kk[idx] == k_temp)     & (lens[idx] == 0):
                corner_id[idx] = 1
            elif (ii[idx] == i + 1) & (jj[idx] == j + 1) & (kk[idx] == k_temp)     & (lens[idx] == 0):
                corner_id[idx] = 2
            elif (ii[idx] == i)     & (jj[idx] == j + 1) & (kk[idx] == k_temp)     & (lens[idx] == 0):
                corner_id[idx] = 3
            elif (ii[idx] == i)     & (jj[idx] == j)     & (kk[idx] == k_temp + 1) & (lens[idx] == 0):
                corner_id[idx] = 4
            elif (ii[idx] == i + 1) & (jj[idx] == j)     & (kk[idx] == k_temp + 1) & (lens[idx] == 0):
                corner_id[idx] = 5
            elif (ii[idx] == i + 1) & (jj[idx] == j + 1) & (kk[idx] == k_temp + 1) & (lens[idx] == 0):
                corner_id[idx] = 6
            elif (ii[idx] == i)     & (jj[idx] == j + 1) & (kk[idx] == k_temp + 1) & (lens[idx] == 0):
                corner_id[idx] = 7

        # now the suspicious points has corner_id != -1 (in range 0 .. 7)

        # sa_corner = np.zeros(8, dtype=np.bool)
        # adj_corners = corners + np.array([k_temp, j, i])

        # assign SOLID/AIR corner, solid = True, Air = False
        sa_corner = np.ones(8, dtype=np.bool)
        adj_corners = corners + np.array([k_temp, j, i])
        for cidx in range(8):
            if norm[0] > 0:
                less_z = adj_corners[cidx, 0] * cfg.domain.dz >= z_vert
            elif norm[0] < 0:
                less_z = adj_corners[cidx, 0] * cfg.domain.dz < z_vert
            else:
                less_z = adj_corners[cidx, 0] * cfg.domain.dz != z_vert
            if norm[1] > 0:
                less_y = adj_corners[cidx, 1] * cfg.domain.dy >= y_vert - cfg.domain.origin_y
            elif norm[1] < 0:
                less_y = adj_corners[cidx, 1] * cfg.domain.dy <= y_vert - cfg.domain.origin_y
            else:
                less_y = adj_corners[cidx, 1] * cfg.domain.dy != y_vert - cfg.domain.origin_y
            if norm[2] > 0:
                less_x = adj_corners[cidx, 2] * cfg.domain.dx >= x_vert - cfg.domain.origin_x
            elif norm[2] < 0:
                less_x = adj_corners[cidx, 2] * cfg.domain.dx <= x_vert - cfg.domain.origin_x
            else:
                less_x = adj_corners[cidx, 2] * cfg.domain.dx != x_vert - cfg.domain.origin_x

            less_all = np.any(less_z & less_y & less_x)
            if norm[0] == 1:
                less_all = np.any(less_all | less_z)

            if less_all:
                # print('corner: ', cidx, ' is Air')
                sa_corner[cidx] = False

        # adjust singularities
        for idx in range(len(x_vert)):
            if corner_id[idx] != -1:
                # print('corner: ', corner_id[idx], ' is Solid')
                sa_corner[corner_id[idx]] = True

        # more detailed search is some obscure cases
        for cidx in range(8):
            for idx in range(len(x_vert)):
                midx = idx # middle point
                lidx = idx - 1 if idx > 0 else len(x_vert) - 1
                ridx = idx + 1 if idx < len(x_vert) - 1 else 0
                p1, p2, p3 = np.array([x_vert[lidx] - cfg.domain.origin_x, y_vert[lidx] - cfg.domain.origin_y, z_vert[lidx]]), \
                             np.array([x_vert[midx] - cfg.domain.origin_x, y_vert[midx] - cfg.domain.origin_y, z_vert[midx]]), \
                             np.array([x_vert[ridx] - cfg.domain.origin_x, y_vert[ridx] - cfg.domain.origin_y, z_vert[ridx]])

                # due to ordering RHS forcing
                norm_temp = normal_vector_triangle(p3, p2, p1)
                # print(idx)
                # print(p1, p2, p3)
                # print(norm_temp)
                if corner_id[idx] != -1:
                    continue
                air_x, air_y, air_z = False, False, False
                sol_x, sol_y, sol_z = False, False, False
                if np.abs(x_vert[idx] - cfg.domain.origin_x - ii[idx] * cfg.domain.dx) > 1.0e-6 and norm_temp[2] != 0:
                    # lies on i-line, select which corner and look if air or solid -> create direction
                    if norm_temp[2] < 0:
                        air_x = adj_corners[cidx, 2] * cfg.domain.dx < x_vert[idx] - cfg.domain.origin_x
                        sol_x = adj_corners[cidx, 2] * cfg.domain.dx > x_vert[idx] - cfg.domain.origin_x
                    elif norm_temp[2] > 0:
                        air_x = adj_corners[cidx, 2] * cfg.domain.dx > x_vert[idx] - cfg.domain.origin_x
                        sol_x = adj_corners[cidx, 2] * cfg.domain.dx < x_vert[idx] - cfg.domain.origin_x
                    # also point must line on correct i/line
                    air_y = sol_y = adj_corners[cidx, 1] * cfg.domain.dy == y_vert[idx] - cfg.domain.origin_y
                    air_z = sol_z = adj_corners[cidx, 0] * cfg.domain.dz == z_vert[idx]

                elif np.abs(y_vert[idx] - cfg.domain.origin_y - jj[idx] * cfg.domain.dy) > 1.0e-6 and norm_temp[1] != 0:
                    # lies on j-line, select which corner and look if air or solid -> create direction
                    if norm_temp[1] < 0:
                        air_y = adj_corners[cidx, 1] * cfg.domain.dy < y_vert[idx] - cfg.domain.origin_y
                        sol_y = adj_corners[cidx, 1] * cfg.domain.dy > y_vert[idx] - cfg.domain.origin_y
                    elif norm_temp[1] > 0:
                        air_y = adj_corners[cidx, 1] * cfg.domain.dy > y_vert[idx] - cfg.domain.origin_y
                        sol_y = adj_corners[cidx, 1] * cfg.domain.dy < y_vert[idx] - cfg.domain.origin_y
                    air_x = sol_x = adj_corners[cidx, 2] * cfg.domain.dx == x_vert[idx] - cfg.domain.origin_x
                    air_z = sol_z = adj_corners[cidx, 0] * cfg.domain.dz == z_vert[idx]

                elif np.abs(z_vert[idx] - kk[idx] * cfg.domain.dz) > 1.0e-6 and norm_temp[0] != 0:
                    # lies on k-line, select which corner and look if air or solid -> create direction
                    if norm_temp[0] < 0:
                        air_z = adj_corners[cidx, 0] * cfg.domain.dz < z_vert[idx]
                        sol_z = adj_corners[cidx, 0] * cfg.domain.dz > z_vert[idx]
                    elif norm_temp[0] > 0:
                        air_z = adj_corners[cidx, 0] * cfg.domain.dz > z_vert[idx]
                        sol_z = adj_corners[cidx, 0] * cfg.domain.dz < z_vert[idx]
                    air_x = sol_x = adj_corners[cidx, 2] * cfg.domain.dx == x_vert[idx] - cfg.domain.origin_x
                    air_y = sol_y = adj_corners[cidx, 1] * cfg.domain.dy == y_vert[idx] - cfg.domain.origin_y

                if norm_temp[0] == 1:
                    air_x = sol_x = adj_corners[cidx, 2] * cfg.domain.dx == x_vert[idx] - cfg.domain.origin_x
                    air_y = sol_y = adj_corners[cidx, 1] * cfg.domain.dy == y_vert[idx] - cfg.domain.origin_y
                    air_z = adj_corners[cidx, 0] * cfg.domain.dz >= z_vert[idx]
                    sol_z = adj_corners[cidx, 0] * cfg.domain.dz <= z_vert[idx]

                air_all = air_x & air_y & air_z
                sol_all = sol_x & sol_y & sol_z
                # print(idx, air_all)
                if air_all:
                    # print(cidx)
                    sa_corner[cidx] = False
                if sol_all:
                    sa_corner[cidx] = True

        # adjust singularities
        for idx in range(len(x_vert)):
            if corner_id[idx] != -1:
                # print('corner: ', corner_id[idx], ' is Solid')
                sa_corner[corner_id[idx]] = True

        # Due to approach that downward facing faces are forbiden, modify all point bellow.
        # If upper corner is solid and corner bellow is air, mark as solid
        if sa_corner[4]:
            sa_corner[0] = True
        if sa_corner[5]:
            sa_corner[1] = True
        if sa_corner[6]:
            sa_corner[2] = True
        if sa_corner[7]:
            sa_corner[3] = True


        # Loop over points, suspicious points are skipped
        #  Loop over direction and check if

        ii_new = [None for irun in range(10)]
        jj_new = [None for irun in range(10)]
        kk_new = [None for irun in range(10)]
        dir_new = [None for irun in range(10)]
        len_new = [None for irun in range(10)]
        x_vert_new = [None for irun in range(10)]
        y_vert_new = [None for irun in range(10)]
        z_vert_new = [None for irun in range(10)]
        ivert = 0
        for idx in range(len(x_vert)):
            if corner_id[idx] == -1:
                x_vert_new[ivert] = x_vert[idx]
                y_vert_new[ivert] = y_vert[idx]
                z_vert_new[ivert] = z_vert[idx]

                if 1==1: #dirs[idx] == 6:
                    if np.abs(x_vert[idx] - cfg.domain.origin_x - ii[idx] * cfg.domain.dx) > 1.0e-6:
                        # lies on i-line, select which corner and look if air or solid -> create direction
                        f1 = False # index if corner was found
                        for cidx in range(8):
                            if (ii[idx]   == adj_corners[cidx, 2]) & (jj[idx] == adj_corners[cidx, 1]) & (kk[idx] == adj_corners[cidx, 0]):
                                f1 = True
                                break
                        cidx1 = cidx
                        f2 = False
                        for cidx in range(8):
                            if (ii[idx] + 1 == adj_corners[cidx, 2]) & (jj[idx] == adj_corners[cidx, 1]) & (kk[idx] == adj_corners[cidx, 0]):
                                f2 = True
                                break
                        cidx2 = cidx
                        f3 = False
                        for cidx in range(8):
                            if (ii[idx] - 1 == adj_corners[cidx, 2]) & (jj[idx] == adj_corners[cidx, 1]) & (kk[idx] == adj_corners[cidx, 0]):
                                f3 = True
                                break
                        cidx3 = cidx
                        if sa_corner[cidx1] and f1:
                            # corner is SOLID, take
                            ccidx = cidx1
                        elif sa_corner[cidx2] and f2:
                            # corner 2 is SOLID take corner 2
                            ccidx = cidx2
                        elif sa_corner[cidx3] and f3:
                            # corner 3 is SOLID, take
                            ccidx = cidx3
                        else:
                            verbose('Some issues with corner, [{},{},{}], [{},{},{}]', ii[idx], jj[idx], kk[idx], i, j, k)
                            exit(1)
                        ii_new[ivert] = int(adj_corners[ccidx, 2])
                        jj_new[ivert] = int(adj_corners[ccidx, 1])
                        kk_new[ivert] = int(adj_corners[ccidx, 0])
                        len_new[ivert] = np.abs(x_vert[idx] - cfg.domain.origin_x - ii_new[ivert] * cfg.domain.dx) / cfg.domain.dx
                        if x_vert[idx] - cfg.domain.origin_x > ii_new[ivert] * cfg.domain.dx:
                            dir_new[ivert] = 4
                        else:
                            dir_new[ivert] = 5
                    elif np.abs(y_vert[idx] - cfg.domain.origin_y - jj[idx] * cfg.domain.dy) > 1.0e-6:
                        # lies on j-line, select which corner and look if air or solid -> create direction
                        f1 = False
                        for cidx in range(8):
                            if (ii[idx] == adj_corners[cidx, 2]) & (jj[idx] == adj_corners[cidx, 1]) & (kk[idx] == adj_corners[cidx, 0]):
                                f1 = True
                                break
                        cidx1 = cidx
                        f2 = False
                        for cidx in range(8):
                            if (ii[idx] == adj_corners[cidx, 2]) & (jj[idx] + 1 == adj_corners[cidx, 1]) & (kk[idx] == adj_corners[cidx, 0]):
                                f2 = True
                                break
                        cidx2 = cidx
                        f3 = False
                        for cidx in range(8):
                            if (ii[idx] == adj_corners[cidx, 2]) & (jj[idx] - 1 == adj_corners[cidx, 1]) & (kk[idx] == adj_corners[cidx, 0]):
                                f3 = True
                                break
                        cidx3 = cidx
                        if sa_corner[cidx1] and f1:
                            # corner is SOLID, take
                            ccidx = cidx1
                        elif sa_corner[cidx2] and f2:
                            # corner 2 is SOLID take corner 2
                            ccidx = cidx2
                        elif sa_corner[cidx3] and f3:
                            # corner 3 is SOLID, take
                            ccidx = cidx3
                        else:
                            verbose('Some issues with corner, [{},{},{}], [{},{},{}]', ii[idx], jj[idx], kk[idx], i, j, k)
                            exit(1)
                        ii_new[ivert] = int(adj_corners[ccidx, 2])
                        jj_new[ivert] = int(adj_corners[ccidx, 1])
                        kk_new[ivert] = int(adj_corners[ccidx, 0])
                        len_new[ivert] = np.abs(y_vert[idx] - cfg.domain.origin_y - jj_new[ivert] * cfg.domain.dy) / cfg.domain.dy
                        if y_vert[idx] - cfg.domain.origin_y > jj_new[ivert] * cfg.domain.dy:
                            dir_new[ivert] = 2
                        else:
                            dir_new[ivert] = 3
                    elif z_vert[idx] - kk[idx] * cfg.domain.dz > 1.0e-6:
                        # lies on k-line, select which corner and look if air or solid -> create direction
                        for cidx in range(8):
                            if (ii[idx] == adj_corners[cidx, 2]) & (jj[idx] == adj_corners[cidx, 1]) & (kk[idx] == adj_corners[cidx, 0]):
                                break
                        cidx1 = cidx
                        for cidx in range(8):
                            if (ii[idx] == adj_corners[cidx, 2]) & (jj[idx] == adj_corners[cidx, 1]) & (kk[idx] + 1 == adj_corners[cidx, 0]):
                                break
                        cidx2 = cidx
                        if sa_corner[cidx1]:
                            # corner is SOLID, take
                            ccidx = cidx1
                        else:
                            # corner 1 is AIR take corner 2
                            ccidx = cidx2
                        ii_new[ivert] = int(adj_corners[ccidx, 2])
                        jj_new[ivert] = int(adj_corners[ccidx, 1])
                        kk_new[ivert] = int(adj_corners[ccidx, 0])
                        len_new[ivert] = np.abs(z_vert[idx] - kk_new[ivert] * cfg.domain.dz) / cfg.domain.dz
                        if z_vert[idx] > kk_new[ivert] * cfg.domain.dz:
                            dir_new[ivert] = 0
                        else:
                            dir_new[ivert] = 1
                else:
                    ii_new[ivert] = int(ii[idx])
                    jj_new[ivert] = int(jj[idx])
                    kk_new[ivert] = int(kk[idx])
                    dir_new[ivert] = int(dirs[idx])
                    len_new[ivert] = lens[idx]
                ivert += 1
                continue
            interfaces = 0
            append_dirs = []
            for idir in range(3):
                dirr = dir_corns[corner_id[idx], idir]

                # is air cell in dir?
                if not sa_corner[dirr]:
                    # print('there is Air cell', idx, corner_id[idx], dirr)
                    interfaces += 1
                    append_dirs.append(vert2dirs[corner_id[idx], idir])

            # add new vertex if interface > 1. case when interface == 1 is the already there
            for fidx in range(interfaces):
                ii_new[ivert] = int(ii[idx])
                jj_new[ivert] = int(jj[idx])
                kk_new[ivert] = int(kk[idx])
                x_vert_new[ivert] = x_vert[idx]
                y_vert_new[ivert] = y_vert[idx]
                z_vert_new[ivert] = z_vert[idx]
                dir_new[ivert] = int(append_dirs[fidx])
                len_new[ivert] = 0.0
                ivert += 1

        n_vert = ivert

        # remove duplicates that has the same direction, but one has len0
        # or check if point is pointing to solid corner
        drop_idx = []
        for ivert in range(n_vert):
            if ii_new[ivert] is None:
                drop_idx.append(ivert)
                continue
            if ivert == 0:
                ileft = n_vert - 1
            elif ivert == n_vert - 1:
                ileft = ivert - 1
            else:
                ileft = ivert - 1

            # check left
            if ii_new[ivert] == ii_new[ileft] and jj_new[ivert] == jj_new[ileft] and \
               kk_new[ivert] == kk_new[ileft] and dir_new[ivert] == dir_new[ileft]:
                # drop the one with zeros index
                if len_new[ivert] == 0.0:
                    drop_idx.append(ivert)
                    # continue
                else:
                    drop_idx.append(ileft)
                    # continue

            # check where point is pointing
            # idir = dir_new[ivert]
            # dirr = dirs_vect[idir]
            # ii_l, jj_l, kk_l = ii_new[ivert], jj_new[ivert], kk_new[ivert]
            # ii_t, jj_t, kk_t = ii_l + dirr[2], jj_l + dirr[1], kk_l + dirr[0]
            # cidx = np.argwhere(adj_corner == np.array([kk_t, jj_t, ii_t]))
            # if cidx.size > 0:
            #     verbose('Found corner that points to solid corner')
                # drop_idx.append(ivert)


        x_vert_n, y_vert_n, z_vert_n, ii_n, jj_n, kk_n, dir_n, len_n = [None for irun in range(10)], [None for irun in range(10)], [None for irun in range(10)], [None for irun in range(10)], [None for irun in range(10)], [None for irun in range(10)], [None for irun in range(10)], [None for irun in range(10)]
        irun = 0
        for ivert in range(n_vert):
            if not ivert in drop_idx:
                x_vert_n[irun] = x_vert_new[ivert]
                y_vert_n[irun] = y_vert_new[ivert]
                z_vert_n[irun] = z_vert_new[ivert]
                ii_n[irun] = ii_new[ivert]
                jj_n[irun] = jj_new[ivert]
                kk_n[irun] = kk_new[ivert]
                dir_n[irun] = dir_new[ivert]
                len_n[irun] = len_new[ivert]
                irun += 1



        x_vert_new, y_vert_new, z_vert_new, ii_new, jj_new, kk_new, dir_new, len_new = x_vert_n, y_vert_n, z_vert_n, ii_n, jj_n, kk_n, dir_n, len_n
        x_vert_new[irun] = x_vert_new[0]
        y_vert_new[irun] = y_vert_new[0]
        z_vert_new[irun] = z_vert_new[0]
        ii_new[irun] = ii_new[0]
        jj_new[irun] = jj_new[0]
        kk_new[irun] = kk_new[0]
        dir_new[irun] = dir_new[0]
        len_new[irun] = len_new[0]

        n_vert = min([irun+1, 7])
        # now delete old entry and create new entry
        to_delete.append(ids)
        to_insert.append((ids, k, j, i, n_vert,
                          wid, rid, lid,
                          isterr, iswall, isroof,
                          ii_new[0], ii_new[1], ii_new[2], ii_new[3], ii_new[4], ii_new[5], ii_new[6],
                          jj_new[0], jj_new[1], jj_new[2], jj_new[3], jj_new[4], jj_new[5], jj_new[6],
                          kk_new[0], kk_new[1], kk_new[2], kk_new[3], kk_new[4], kk_new[5], kk_new[6],
                          len_new[0], len_new[1], len_new[2], len_new[3], len_new[4], len_new[5], len_new[6],
                          dir_new[0], dir_new[1], dir_new[2], dir_new[3], dir_new[4], dir_new[5], dir_new[6],
                          norm[0], norm[1], norm[2], area,
                          center,
                          x_vert_new[0], y_vert_new[0], z_vert_new[0], cfg.srid_palm,
                          x_vert_new[1], y_vert_new[1], z_vert_new[1], cfg.srid_palm,
                          x_vert_new[2], y_vert_new[2], z_vert_new[2], cfg.srid_palm,
                          x_vert_new[3], y_vert_new[3], z_vert_new[3], cfg.srid_palm,
                          x_vert_new[4], y_vert_new[4], z_vert_new[4], cfg.srid_palm,
                          x_vert_new[5], y_vert_new[5], z_vert_new[5], cfg.srid_palm,
                          x_vert_new[6], y_vert_new[6], z_vert_new[6], cfg.srid_palm, ))

    debug('Deleting all unwanted rows')
    sqltext = 'DELETE FROM "{0}"."{1}" ' \
              'WHERE id = ANY(%s)'.format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext, (to_delete,))
    sql_debug(connection)
    connection.commit()

    debug('Inserting all new entries into slanted faces table')
    sqltext = 'INSERT INTO "{0}"."{1}"  ' \
              '   (id, k, j, i, n_vert, ' \
              '    wid, rid, lid, ' \
              '    isterr, iswall, isroof, ' \
              '    ii1, ii2, ii3, ii4, ii5, ii6, ii7, ' \
              '    jj1, jj2, jj3, jj4, jj5, jj6, jj7, ' \
              '    kk1, kk2, kk3, kk4, kk5, kk6, kk7, ' \
              '    len1, len2, len3, len4, len5, len6, len7, ' \
              '    dir1, dir2, dir3, dir4, dir5, dir6, dir7, ' \
              '    normz, normy, normx, area, ' \
              '    center,  ' \
              '    vert1, vert2, vert3, vert4, vert5, vert6, vert7) ' \
              'VALUES (%s, %s, %s, %s, %s, ' \
              '        %s, %s, %s, ' \
              '        %s, %s, %s, ' \
              '        %s, %s, %s, %s, %s, %s, %s, ' \
              '        %s, %s, %s, %s, %s, %s, %s, ' \
              '        %s, %s, %s, %s, %s, %s, %s, ' \
              '        %s, %s, %s, %s, %s, %s, %s, ' \
              '        %s, %s, %s, %s, %s, %s, %s, ' \
              '        %s, %s, %s, %s, ' \
              '        %s, ' \
              '        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
              '        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
              '        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
              '        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
              '        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
              '        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
              '        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s) ' \
              '        )  '.format(cfg.domain.case_schema, cfg.tables.slanted_faces)

    cur.executemany(sqltext, to_insert)
    sql_debug(connection)
    connection.commit()

    debug('Updating 3d polygon')
    sqltext = 'UPDATE "{0}"."{1}" SET geom =  ' \
              'ST_ForceRHR(' \
              'ST_SetSRID(ST_MakePolygon(ST_MakeLine(ARRAY[vert1, vert2, vert3, vert4, vert5, vert6, vert7,' \
              '     CASE WHEN vert1 IS NOT NULL THEN vert1 ' \
              '          WHEN vert2 IS NOT NULL THEN vert2 ' \
              '          WHEN vert3 IS NOT NULL THEN vert3 ' \
              '          WHEN vert4 IS NOT NULL THEN vert4 ' \
              '          WHEN vert5 IS NOT NULL THEN vert5 ' \
              '          WHEN vert6 IS NOT NULL THEN vert6 ' \
              '          WHEN vert7 IS NOT NULL THEN vert7 END' \
              '])), %s)) ' \
              'WHERE geom IS NULL'.format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext, (cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

def check_for_vertex_singularities2(cfg, connection, cur):
    """ """
    dirs_singular = np.array([[[+1, 0, 0], [+1,-1, 0], [+1, 0,-1], [+1,-1,-1]],
                              [[ 0, 0, 0], [ 0,-1, 0], [ 0, 0,-1], [ 0,-1,-1]],
                              [[ 0, 0, 0], [ 0, 0,-1], [+1, 0, 0], [+1, 0,-1]],
                              [[ 0,-1, 0], [ 0,-1,-1], [+1,-1, 0], [+1,-1,-1]],
                              [[ 0, 0, 0], [ 0,-1, 0], [+1, 0, 0], [+1,-1, 0]],
                              [[ 0, 0,-1], [ 0,-1,-1], [+1, 0,-1], [1, -1,-1]]
                             ])
    vert2dirs = np.array([[4, 0, 2],
                          [2, 0, 5],
                          [5, 0, 3],
                          [3, 0, 4],
                          [4, 1, 2],
                          [2, 1, 5],
                          [5, 1, 3],
                          [3, 1, 4],])
    dir_corns = np.array([[1, 4, 3],
                          [2, 5, 0],
                          [3, 6, 1],
                          [0, 7, 2],
                          [5, 0, 7],
                          [6, 1, 4],
                          [7, 2, 5],
                          [4, 3, 6],
                          ])

    corners = np.array([[ 0, 0, 0],
                        [ 0, 0, 1],
                        [ 0, 1, 1],
                        [ 0, 1, 0],
                        [ 1, 0, 0],
                        [ 1, 0, 1],
                        [ 1, 1, 1],
                        [ 1, 1, 0]])

    sqltext = 'CREATE INDEX IF NOT EXISTS slanted_faces_kji_idx ON "{0}"."{1}" (k,j,i); ' \
              .format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()
    # exit(1)
    progress('Checking if vertices has 4 adjacent faces')
    sqltext = 'SELECT id, i,j,k, normx, normy, normz, ' \
              '       wid, rid, lid, isterr, iswall, isroof,' \
              '       ii1, ii2, ii3, ii4, ii5, ii6, ' \
              '       jj1, jj2, jj3, jj4, jj5, jj6, ' \
              '       kk1, kk2, kk3, kk4, kk5, kk6, ' \
              '       dir1, dir2, dir3, dir4, dir5, dir6, ' \
              '       CASE WHEN len1 = 0              THEN 1 ELSE NULL END, ' \
              '       CASE WHEN len2 = 0              THEN 2 ELSE NULL END, ' \
              '       CASE WHEN len3 = 0              THEN 3 ELSE NULL END, ' \
              '       CASE WHEN len4 = 0 AND n_vert>4 THEN 4 ELSE NULL END, ' \
              '       CASE WHEN len5 = 0 AND n_vert>5 THEN 5 ELSE NULL END, ' \
              '       CASE WHEN len6 = 0 AND n_vert>6 THEN 6 ELSE NULL END,  ' \
              '       ARRAY[ST_Z(vert1), ST_Y(vert1), ST_X(vert1)], ' \
              '       ARRAY[ST_Z(vert2), ST_Y(vert2), ST_X(vert2)], ' \
              '       ARRAY[ST_Z(vert3), ST_Y(vert3), ST_X(vert3)], ' \
              '       ARRAY[ST_Z(vert4), ST_Y(vert4), ST_X(vert4)], ' \
              '       ARRAY[ST_Z(vert5), ST_Y(vert5), ST_X(vert5)], ' \
              '       ARRAY[ST_Z(vert6), ST_Y(vert6), ST_X(vert6)]' \
              '       '  \
              '  FROM "{0}"."{1}" ' \
              '  WHERE (' \
              '  len1 = 0 OR' \
              '  len2 = 0 OR' \
              '  len3 = 0 OR' \
              '  len4 = 0 OR' \
              '  len5 = 0 OR' \
              '  len6 = 0) ' \
              ''.format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext)
    singulars = cur.fetchall()
    sql_debug(connection)
    connection.commit()

    # max id of slanted faces
    sqltext = 'SELECT MAX(id) FROM "{0}"."{1}"'.format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext)
    max_id = cur.fetchone()[0]
    sql_debug(connection)
    connection.commit()

    to_insert = []
    polygons = []
    polygons_props = []
    polygons2 = {}
    for singular in singulars:
        id, i, j, k, normx, normy, normz = singular[0], singular[1], singular[2], singular[3], singular[4], singular[5], singular[6]
        wid, rid, lid, isterr, iswall, isroof = singular[7], singular[8], singular[9], singular[10], singular[11], singular[12]
        ii0 = [irun for irun in singular[13:19]]
        jj0 = [irun for irun in singular[19:25]]
        kk0 = [irun for irun in singular[25:31]]
        dirs0 = [irun for irun in singular[31:37]]
        sidx = [irun for irun in singular[37:43]]
        lastidx=43
        x_vert, y_vert, z_vert = [], [], []
        ii, jj, kk, dirs = [], [], [], []
        for idx in range(6):
            if sidx[idx] is not None:
                z_vert.append(singular[lastidx + idx][0])
                y_vert.append(singular[lastidx + idx][1])
                x_vert.append(singular[lastidx + idx][2])
                ii.append(ii0[idx])
                jj.append(jj0[idx])
                kk.append(kk0[idx])
                dirs.append(dirs0[idx])

        # check surrounding face of the vertex, there has to be 4 of them
        sqltext = 'SELECT COUNT(*) FROM "{0}"."{1}" ' \
                  'WHERE i={2} AND j={3} AND k={4} '

        for vidx in range(len(ii)):
            counter = 0
            for ndirs in range(4):
                dk, dj, di = dirs_singular[dirs[vidx], ndirs]
                cur.execute(sqltext.format(cfg.domain.case_schema, cfg.tables.slanted_faces, ii[vidx]+di, jj[vidx]+dj, kk[vidx]+dk))
                count = cur.fetchone()
                sql_debug(connection)
                connection.commit()
                if count[0] == 0:
                    # new polygon, if not exists must be created
                    if not [ii[vidx]+di, jj[vidx]+dj, kk[vidx]+dk] in polygons:
                        # print('new polygon', ii[vidx]+di, jj[vidx]+dj, kk[vidx]+dk)

                        polygons.append([ii[vidx]+di, jj[vidx]+dj, kk[vidx]+dk])
                        polygons_props.append([wid, rid, lid, isterr, iswall, isroof, normx, normy, normz])
                        pidx = polygons.index([ii[vidx]+di, jj[vidx]+dj, kk[vidx]+dk])
                        polygons2[pidx] = []
                        polygons2[pidx].append([z_vert[vidx], y_vert[vidx], x_vert[vidx], kk[vidx], jj[vidx], ii[vidx]])
                    else:
                        # print('existing polygon', ii[vidx]+di, jj[vidx]+dj, kk[vidx]+dk)
                        pidx = polygons.index([ii[vidx]+di, jj[vidx]+dj, kk[vidx]+dk])
                        # check if vertex in not there already
                        if not [z_vert[vidx], y_vert[vidx], x_vert[vidx], kk[vidx], jj[vidx], ii[vidx]] in polygons2[pidx]:
                            polygons2[pidx].append([z_vert[vidx], y_vert[vidx], x_vert[vidx], kk[vidx], jj[vidx], ii[vidx]])
                # counter += count[0]
                # print(count, ii[vidx]+di, jj[vidx]+dj, kk[vidx]+dk)

    n_new_polygons = len(polygons)
    for npol in range(n_new_polygons):
        x_vert_new, y_vert_new, z_vert_new, ii_new, jj_new, kk_new, dir_new, len_new = [None for irun in range(7)], [None for irun in range(7)], [None for irun in range(7)], [None for irun in range(7)], [None for irun in range(7)], [None for irun in range(7)], [None for irun in range(7)], [None for irun in range(7)]
        n_vert = len(polygons2[npol])
        i,j,k = polygons[npol]
        wid, rid, lid, isterr, iswall, isroof, normx, normy, normz = polygons_props[npol]
        # print(max_id+npol+1, k,j,i, n_vert,)
        for nv in range(n_vert):
            x_vert_new[nv] = polygons2[npol][nv][2]
            y_vert_new[nv] = polygons2[npol][nv][1]
            z_vert_new[nv] = polygons2[npol][nv][0]
            ii_new[nv]     = polygons2[npol][nv][5]
            jj_new[nv]     = polygons2[npol][nv][4]
            kk_new[nv]     = polygons2[npol][nv][3]

        # Sort the point, no duplicity
        sa_corner = np.zeros(8, dtype=np.bool)
        k_temp = k - 1
        corner_id = -1 * np.ones(7, dtype=np.int)
        for idx in range(n_vert):
            if   (ii_new[idx] == i)     & (jj_new[idx] == j)     & (kk_new[idx] == k_temp):
                corner_id[idx] = 0
            elif (ii_new[idx] == i + 1) & (jj_new[idx] == j)     & (kk_new[idx] == k_temp):
                corner_id[idx] = 1
            elif (ii_new[idx] == i + 1) & (jj_new[idx] == j + 1) & (kk_new[idx] == k_temp):
                corner_id[idx] = 2
            elif (ii_new[idx] == i)     & (jj_new[idx] == j + 1) & (kk_new[idx] == k_temp):
                corner_id[idx] = 3
            elif (ii_new[idx] == i)     & (jj_new[idx] == j)     & (kk_new[idx] == k_temp + 1):
                corner_id[idx] = 4
            elif (ii_new[idx] == i + 1) & (jj_new[idx] == j)     & (kk_new[idx] == k_temp + 1):
                corner_id[idx] = 5
            elif (ii_new[idx] == i + 1) & (jj_new[idx] == j + 1) & (kk_new[idx] == k_temp + 1):
                corner_id[idx] = 6
            elif (ii_new[idx] == i)     & (jj_new[idx] == j + 1) & (kk_new[idx] == k_temp + 1):
                corner_id[idx] = 7
        adj_corners = corners + np.array([k_temp, j, i])
        for idx in range(n_vert):
            if corner_id[idx] != -1:
                # print('corner: ', corner_id[idx], ' is Solid')
                sa_corner[corner_id[idx]] = True

        ii, jj, kk, dirs, lens, x_vert, y_vert, z_vert = [None for irun in range(7)], [None for irun in range(7)], [None for irun in range(7)], [None for irun in range(7)], [None for irun in range(7)], [None for irun in range(7)], [None for irun in range(7)], [None for irun in range(7)]
        ivert = 0
        for idx in range(n_vert):
            interfaces = 0
            for idir in range(3):
                dirr = dir_corns[corner_id[idx], idir]
                # is air cell in dir?
                if not sa_corner[dirr]:
                    interfaces += 1
                    ii[ivert] = int(ii_new[idx])
                    jj[ivert] = int(jj_new[idx])
                    kk[ivert] = int(kk_new[idx])
                    x_vert[ivert] = x_vert_new[idx]
                    y_vert[ivert] = y_vert_new[idx]
                    z_vert[ivert] = z_vert_new[idx]
                    dirs[ivert] = int(vert2dirs[corner_id[idx], idir])
                    lens[ivert] = 0.0
                    ivert += 1

        n_vert = ivert + 1
        x_vert[ivert] = x_vert[0]
        y_vert[ivert] = y_vert[0]
        z_vert[ivert] = z_vert[0]
        ii[ivert] = ii[0]
        jj[ivert] = jj[0]
        kk[ivert] = kk[0]
        dirs[ivert] = dirs[0]
        lens[ivert] = lens[0]

        to_insert.append((max_id+npol+1, int(k), int(j), int(i), int(n_vert),
                          wid, rid, lid,
                          isterr, iswall, isroof,
                          ii[0], ii[1], ii[2], ii[3], ii[4], ii[5], ii[6],
                          jj[0], jj[1], jj[2], jj[3], jj[4], jj[5], jj[6],
                          kk[0], kk[1], kk[2], kk[3], kk[4], kk[5], kk[6],
                          lens[0], lens[1], lens[2], lens[3], lens[4], lens[5], lens[6],
                          dirs[0], dirs[1], dirs[2], dirs[3], dirs[4], dirs[5], dirs[6],
                          normz, normy, normx, 0.0,
                          x_vert[0], y_vert[0], z_vert[0], cfg.srid_palm,
                          x_vert[1], y_vert[1], z_vert[1], cfg.srid_palm,
                          x_vert[2], y_vert[2], z_vert[2], cfg.srid_palm,
                          x_vert[3], y_vert[3], z_vert[3], cfg.srid_palm,
                          x_vert[4], y_vert[4], z_vert[4], cfg.srid_palm,
                          x_vert[5], y_vert[5], z_vert[5], cfg.srid_palm,
                          x_vert[6], y_vert[6], z_vert[6], cfg.srid_palm, ))
    debug('Inserting all new entries into slanted faces table')
    sqltext = 'INSERT INTO "{0}"."{1}"  ' \
              '   (id, k, j, i, n_vert, ' \
              '    wid, rid, lid, ' \
              '    isterr, iswall, isroof, ' \
              '    ii1, ii2, ii3, ii4, ii5, ii6, ii7, ' \
              '    jj1, jj2, jj3, jj4, jj5, jj6, jj7, ' \
              '    kk1, kk2, kk3, kk4, kk5, kk6, kk7, ' \
              '    len1, len2, len3, len4, len5, len6, len7, ' \
              '    dir1, dir2, dir3, dir4, dir5, dir6, dir7, ' \
              '    normz, normy, normx, area, ' \
              '    vert1, vert2, vert3, vert4, vert5, vert6, vert7) ' \
              'VALUES (%s, %s, %s, %s, %s, ' \
              '        %s, %s, %s, ' \
              '        %s, %s, %s, ' \
              '        %s, %s, %s, %s, %s, %s, %s, ' \
              '        %s, %s, %s, %s, %s, %s, %s, ' \
              '        %s, %s, %s, %s, %s, %s, %s, ' \
              '        %s, %s, %s, %s, %s, %s, %s, ' \
              '        %s, %s, %s, %s, %s, %s, %s, ' \
              '        %s, %s, %s, %s, ' \
              '        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
              '        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
              '        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
              '        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
              '        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
              '        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s), ' \
              '        ST_SetSRID(ST_MakePoint(%s, %s, %s), %s) ' \
              '        )  '.format(cfg.domain.case_schema, cfg.tables.slanted_faces)

    cur.executemany(sqltext, to_insert)
    sql_debug(connection)
    connection.commit()

    debug('Updating empty centers')
    sqltext = 'UPDATE "{0}"."{1}" SET center = ' \
              'ST_SetSRID(ST_MakePoint(' \
              ' (ST_X(vert1) + ST_X(vert2) + ' \
              '       CASE WHEN n_vert > 3 THEN ST_X(vert3) ELSE 0.0 END + ' \
              '       CASE WHEN n_vert > 4 THEN ST_X(vert4) ELSE 0.0 END + ' \
              '       CASE WHEN n_vert > 5 THEN ST_X(vert5) ELSE 0.0 END + ' \
              '       CASE WHEN n_vert > 6 THEN ST_X(vert6) ELSE 0.0 END) / (n_vert - 1), ' \
              ' (ST_Y(vert1) + ST_Y(vert2) + ' \
              '       CASE WHEN n_vert > 3 THEN ST_Y(vert3) ELSE 0.0 END + ' \
              '       CASE WHEN n_vert > 4 THEN ST_Y(vert4) ELSE 0.0 END + ' \
              '       CASE WHEN n_vert > 5 THEN ST_Y(vert5) ELSE 0.0 END + ' \
              '       CASE WHEN n_vert > 6 THEN ST_Y(vert6) ELSE 0.0 END) / (n_vert - 1),' \
              ' (ST_Z(vert1) + ST_Z(vert2) + ' \
              '       CASE WHEN n_vert > 3 THEN ST_Z(vert3) ELSE 0.0 END + ' \
              '       CASE WHEN n_vert > 4 THEN ST_Z(vert4) ELSE 0.0 END + ' \
              '       CASE WHEN n_vert > 5 THEN ST_Z(vert5) ELSE 0.0 END + ' \
              '       CASE WHEN n_vert > 6 THEN ST_Z(vert6) ELSE 0.0 END) / (n_vert - 1)' \
              '), %s) ' \
              'WHERE center IS NULL '.format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext, (cfg.srid_palm,))
    sql_debug(connection)
    connection.commit()

def slanted_surface_init(cfg, connection, cur):
    """ initialize slanted surfaces """
    preprocess_terrain_height(cfg, connection, cur)

    if cfg.has_buildings:
        preprocess_building_height(cfg, connection, cur)

        create_slanted_walls(cfg, connection, cur)

        calculate_aspect_slope(cfg, connection, cur)

    create_slanted_terrain(cfg, connection, cur)

    if cfg.has_buildings:
        create_slated_roof(cfg, connection, cur)

    # process slanted face into gridded form
    if cfg.has_buildings:
        create_grid_slanted_walls(cfg, connection, cur)

    create_grid_slanted_terrain(cfg, connection, cur)

    if cfg.has_buildings:
        create_grid_slanted_roof(cfg, connection, cur)

    initialize_slanted_faces(cfg, connection, cur)

    if cfg.has_buildings:
        merge_walls_terrain(cfg, connection, cur)
        merge_walls_roofs(cfg, connection, cur)


    # TODO: check if there triple duplicates wall/terrain/roof
    debug('Add index on slanted geom')
    sqltext = 'CREATE INDEX slanted_face_geom_index ON "{0}"."{1}" USING gist(geom)'\
              .format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # obtain number of vertices
    sqltext = 'UPDATE "{0}"."{1}" SET ' \
              'n_vert = CASE WHEN vert7 IS NOT NULL THEN 7 ' \
              '              WHEN vert6 IS NOT NULL THEN 6 ' \
              '              WHEN vert5 IS NOT NULL THEN 5 ' \
              '              WHEN vert4 IS NOT NULL THEN 4 ' \
              '              WHEN vert3 IS NOT NULL THEN 3 ' \
              '              WHEN vert2 IS NOT NULL THEN 2 ' \
              '              WHEN vert1 IS NOT NULL THEN 1 END '.format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # TODO: remove faces that are "under" terrain, there are walls that did not merge with terrain, PALM should take care of them
    sqltext = 'DELETE FROM "{0}"."{1}" ' \
              'WHERE iswall AND (' \
              '                  ST_Z(vert1) = 0 OR ' \
              '                  ST_Z(vert2) = 0 OR ' \
              '                  ST_Z(vert3) = 0 OR ' \
              '                  ST_Z(vert4) = 0 OR' \
              '                  ST_Z(vert5) = 0 OR' \
              '                  ST_Z(vert6) = 0 OR' \
              '                  ST_Z(vert7) = 0' \
              ')'.format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # calculate normal vector using triangulation
    normal_vector_trinagulation(cfg, connection, cur)

    create_integer_vertices(cfg, connection, cur)

    # TODO: find duplicates, there are some points with the same x,y,z

    check_for_vertex_singularities(cfg, connection, cur)
    #
    # check_for_vertex_singularities2(cfg, connection, cur)

    # calculate vertices
    create_vertices_indexes(cfg, connection, cur)

    # put and indexies at the lid, rid, gid
    sqltext = 'CREATE INDEX IF NOT EXISTS slanted_faces_rid_idx ON "{0}"."{1}" (rid); ' \
              'CREATE INDEX IF NOT EXISTS slanted_faces_lid_idx ON "{0}"."{1}" (lid); ' \
              'CREATE INDEX IF NOT EXISTS slanted_faces_wid_idx ON "{0}"."{1}" (wid)'\
              .format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # calculate normal vector using triangulation, MB just in case, do it twice
    normal_vector_trinagulation(cfg, connection, cur)

    # create vtk faces
    if cfg.slanted_pars.do_vtk:
        create_slanted_vtk(cfg, connection, cur)

    # final touch, remove extra vertex in face
    debug('Removing extra vertex')
    verbose('Decreasing n_vert by 1')
    sqltext = 'UPDATE "{0}"."{1}" SET n_vert = n_vert - 1'.format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()
    for iv in range(1,8):
        verbose('Removing extra vertices, {}', iv)
        sqltext = 'UPDATE "{0}"."{1}" SET vert{2}i = NULL ' \
                  'WHERE n_vert < {2}'.format(cfg.domain.case_schema, cfg.tables.slanted_faces, iv)
        cur.execute(sqltext)
        sql_debug(connection)
        connection.commit()

def slanted_write_nc(ncfile, cfg, connection, cur):
    """ DESC """
    # now dump point to Netcdf4
    progress('Creating nc file')
    cur.execute('SELECT COUNT(*) FROM "{0}"."{1}"'.format(cfg.domain.case_schema, cfg.tables.vertices))
    num_vert = cur.fetchone()[0]
    sql_debug(connection)
    connection.commit()

    empty_vert = cfg.slanted_pars.empty_vert # empty_vert is for missing points

    debug('Building dimension')
    nc_create_dimension(ncfile, 'cct_num_vert', num_vert)

    ncfile.setncattr('LOD', 2)
    ncfile.setncattr('empty_vert', empty_vert)

    debug('Building vertices')
    vn = 'cct_vertices'
    nc_vertices = ncfile.createVariable(vn, 'f8', ('dim_3d', 'cct_num_vert',))
    nc_write_attribute(ncfile, vn, 'long_name', '[z,y,x] coordinates of vertices')
    sqltext = 'SELECT CAST(ST_Z(point) AS FLOAT), ' \
              '       CAST(ST_Y(point) AS FLOAT), ' \
              '       CAST(ST_X(point) AS FLOAT) ' \
              'FROM "{0}"."{1}" ORDER BY id'.format(cfg.domain.case_schema, cfg.tables.vertices)
    cur.execute(sqltext)
    vert_points = cur.fetchall()
    sql_debug(connection)
    connection.commit()
    nc_vertices[0, :] = [x[0] for x in vert_points]
    nc_vertices[1, :] = [x[1] for x in vert_points]
    nc_vertices[2, :] = [x[2] for x in vert_points]
    nc_vertices[2, :] += - cfg.domain.origin_x
    nc_vertices[1, :] += - cfg.domain.origin_y

    del vert_points

    debug('Building vertices 2')
    vn = 'cct_vertex_coords'
    nc_verticex_coords = ncfile.createVariable(vn, 'i4', ('cct_dim_vertex_coords', 'cct_num_vert', ))
    nc_write_attribute(ncfile, vn, 'long_name', 'Type of vertices coordinates for radiation module. Order: [k;j;i;dir], kji index of corner, dir direction: [0.+z, 1.-z, 2.+y, 3.-y, 4.+x, 5.-x]')

    sqltext = 'SELECT k, j, i, dir FROM "{0}"."{1}" ORDER BY id' \
              .format(cfg.domain.case_schema, cfg.tables.vertices)
    cur.execute(sqltext)
    vert2_points = cur.fetchall()
    sql_debug(connection)
    connection.commit()
    nc_verticex_coords[0, :] = [x[0] for x in vert2_points]
    nc_verticex_coords[1, :] = [x[1] for x in vert2_points]
    nc_verticex_coords[2, :] = [x[2] for x in vert2_points]
    nc_verticex_coords[3, :] = [x[3] for x in vert2_points]
    del vert2_points

    vn = 'cct_vertex_shifts'
    nc_verticex_shifts = ncfile.createVariable(vn, 'f8', ('cct_dim_vertex_shifts', 'cct_num_vert', ))
    nc_write_attribute(ncfile, vn, 'long_name', 'Type of vertices coordinates for radiation module. Order: [len], len distance [corner; vertex]')

    sqltext = 'SELECT len FROM "{0}"."{1}" ORDER BY id' \
              .format(cfg.domain.case_schema, cfg.tables.vertices)
    cur.execute(sqltext)
    vert2_points = cur.fetchall()
    sql_debug(connection)
    connection.commit()
    nc_verticex_shifts[:] = [x[0] for x in vert2_points]
    del vert2_points

    cur.execute('SELECT COUNT(*) FROM "{0}"."{1}"'.format(cfg.domain.case_schema, cfg.tables.slanted_faces))
    num_faces = cur.fetchone()[0]
    sql_debug(connection)
    connection.commit()
    nc_create_dimension(ncfile, 'cct_num_faces', num_faces)

    vn = 'cct_vertices_per_face'
    faces = ncfile.createVariable(vn, 'i4', ('cct_max_num_vertices_per_face', 'cct_num_faces', ))
    nc_write_attribute(ncfile, vn, 'long_name', 'list of vertices for each face, ordered by Right-Hand-Rule')

    vn = 'cct_face_center'
    centers = ncfile.createVariable(vn, 'f8', ('dim_3d', 'cct_num_faces',))
    nc_write_attribute(ncfile, vn, 'long_name', 'positions of slanted faces centers, [z;y;x]')

    vn = 'cct_3d_grid_indices'
    kji_locs = ncfile.createVariable(vn, 'i4', ('dim_3d', 'cct_num_faces',))
    nc_write_attribute(ncfile, vn, 'long_name', '[kji] locations of slanted faces')

    vn = 'cct_offsets'
    offs     = ncfile.createVariable(vn, 'i4', ('dim_3d', 'cct_num_faces',))
    nc_write_attribute(ncfile, vn, 'long_name', '[kji] offsets for slanted faces, relates each surfaces to its building (need for properties)')

    vn = 'cct_face_normal_vector'
    normals = ncfile.createVariable(vn, 'f8', ('dim_3d', 'cct_num_faces',))
    nc_write_attribute(ncfile, vn, 'long_name', '[z;y;x] components of normalized normal vector')

    vn = 'cct_face_area'
    area = ncfile.createVariable(vn, 'f8', ('cct_num_faces',))
    nc_write_attribute(ncfile, vn, 'long_name', 'area of slanted face surface')

    vn = 'cct_num_vertices_per_face'
    num_edges = ncfile.createVariable(vn, 'i4', ('cct_num_faces',))
    nc_write_attribute(ncfile, vn, 'long_name', 'number of vertices in each face')

    vn = 'cct_surface_type_classification'
    types = ncfile.createVariable(vn, 'i4', ('cct_num_faces',))
    nc_write_attribute(ncfile, vn, 'long_name', 'index for separating lsm=0, usm wall=1, usm roof=2 surfaces')

    sqltext = 'SELECT vert1i, vert2i, vert3i, vert4i, vert5i, vert6i, vert7i, n_vert, ' \
              '       CASE WHEN isterr THEN 0 WHEN iswall THEN 1 ELSE 2 END  ' \
              'FROM "{0}"."{1}" ORDER BY k,j,i'\
              .format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext)
    vert_points = cur.fetchall()
    sql_debug(connection)
    connection.commit()
    for i in range(7):
        faces[i,:] = [empty_vert if x[i] is None else x[i] for x in vert_points]
    num_verti = np.asarray([x[7] for x in vert_points])
    num_edges[:] = num_verti
    num_types = np.asarray([x[8] for x in vert_points])
    types[:] = num_types

    vn = 'cct_vegetation_type_classification'
    vn_type = 'i4'
    vt_fill = cfg.fill_values[vn_type]
    veg_types = ncfile.createVariable(vn, vn_type, ('cct_num_faces',), fill_value=vt_fill)
    nc_write_attribute(ncfile, vn, 'long_name', 'PALM PIDS vegetation type')
    sqltext = 'SELECT l.type - {3}' \
              'FROM "{0}"."{1}" AS s ' \
              'LEFT OUTER JOIN "{0}"."{2}" AS l ON l.lid = s.lid AND l.type BETWEEN {3} AND {4}' \
              'ORDER BY k,j,i'\
              .format(cfg.domain.case_schema, cfg.tables.slanted_faces, cfg.tables.landcover,
                      cfg.type_range.vegetation_min, cfg.type_range.vegetation_max,  )
    cur.execute(sqltext)
    cct_type = cur.fetchall()
    sql_debug(connection)
    connection.commit()
    veg_types[:] = [vt_fill if x[0] is None else x[0] for x in cct_type]

    vn = 'cct_pavement_type_classification'
    vn_type = 'i4'
    vt_fill = cfg.fill_values[vn_type]
    pav_types = ncfile.createVariable(vn, vn_type, ('cct_num_faces',), fill_value=vt_fill)
    nc_write_attribute(ncfile, vn, 'long_name', 'PALM PIDS pavement type')
    sqltext = 'SELECT l.type - {3}' \
              'FROM "{0}"."{1}" AS s ' \
              'LEFT OUTER JOIN "{0}"."{2}" AS l ON l.lid = s.lid AND l.type BETWEEN {3} AND {4}' \
              'ORDER BY k,j,i'\
              .format(cfg.domain.case_schema, cfg.tables.slanted_faces, cfg.tables.landcover,
                      cfg.type_range.pavement_min, cfg.type_range.pavement_max,  )
    cur.execute(sqltext)
    cct_type = cur.fetchall()
    sql_debug(connection)
    connection.commit()
    pav_types[:] = [vt_fill if x[0] is None else x[0] for x in cct_type]

    vn = 'cct_water_type_classification'
    vn_type = 'i4'
    vt_fill = cfg.fill_values[vn_type]
    wat_types = ncfile.createVariable(vn, vn_type, ('cct_num_faces',), fill_value=vt_fill)
    nc_write_attribute(ncfile, vn, 'long_name', 'PALM PIDS water type')
    sqltext = 'SELECT l.type - {3}' \
              'FROM "{0}"."{1}" AS s ' \
              'LEFT OUTER JOIN "{0}"."{2}" AS l ON l.lid = s.lid AND l.type BETWEEN {3} AND {4}' \
              'ORDER BY k,j,i'\
              .format(cfg.domain.case_schema, cfg.tables.slanted_faces, cfg.tables.landcover,
                      cfg.type_range.water_min, cfg.type_range.water_max,  )
    cur.execute(sqltext)
    cct_type = cur.fetchall()
    sql_debug(connection)
    connection.commit()
    wat_types[:] = [vt_fill if x[0] is None else x[0] for x in cct_type]

    if cfg.has_buildings:
        vn = 'cct_building_type_classification'
        vn_type = 'i4'
        vt_fill = cfg.fill_values[vn_type]
        build_types = ncfile.createVariable(vn, vn_type, ('cct_num_faces',), fill_value=vt_fill)
        nc_write_attribute(ncfile, vn, 'long_name', 'PALM PIDS buildings type')
        sqltext = 'SELECT CASE WHEN iswall THEN lw.type - {6} ' \
                  '            WHEN isroof THEN lr.type - {6} ' \
                  '            ELSE NULL END ' \
                  'FROM "{0}"."{1}" AS s ' \
                  'LEFT OUTER JOIN "{0}"."{2}" AS w ON w.wid = s.wid ' \
                  'LEFT OUTER JOIN "{0}"."{3}" AS lw ON lw.lid = w.lid AND lw.type BETWEEN {6} AND {7} ' \
                  'LEFT OUTER JOIN "{0}"."{4}" AS r ON r.rid = s.rid ' \
                  'LEFT OUTER JOIN "{0}"."{5}" AS lr ON lr.lid = r.lid AND lr.type BETWEEN {6} AND {7} ' \
                  'ORDER BY k,j,i'\
                  .format(cfg.domain.case_schema, cfg.tables.slanted_faces, cfg.tables.walls,
                          cfg.tables.landcover, cfg.tables.roofs, cfg.tables.landcover,
                          cfg.type_range.building_min, cfg.type_range.building_max,  )
        cur.execute(sqltext)
        cct_type = cur.fetchall()
        sql_debug(connection)
        connection.commit()
        build_types[:] = [vt_fill if x[0] is None else x[0] for x in cct_type]

        vn = 'cct_building_id_classification'
        vn_type = 'i4'
        vt_fill = cfg.fill_values[vn_type]
        build_ids = ncfile.createVariable(vn, vn_type, ('cct_num_faces',), fill_value=vt_fill)
        nc_write_attribute(ncfile, vn, 'long_name', 'PALM PIDS buildings id')
        sqltext = 'SELECT CASE WHEN iswall THEN lw.lid ' \
                  '            WHEN isroof THEN lr.lid ' \
                  '            ELSE NULL END ' \
                  'FROM "{0}"."{1}" AS s ' \
                  'LEFT OUTER JOIN "{0}"."{2}" AS w ON w.wid = s.wid ' \
                  'LEFT OUTER JOIN "{0}"."{3}" AS lw ON lw.lid = w.lid AND lw.type BETWEEN {6} AND {7} ' \
                  'LEFT OUTER JOIN "{0}"."{4}" AS r ON r.rid = s.rid ' \
                  'LEFT OUTER JOIN "{0}"."{5}" AS lr ON lr.lid = r.lid AND lr.type BETWEEN {6} AND {7} ' \
                  'ORDER BY k,j,i'\
                  .format(cfg.domain.case_schema, cfg.tables.slanted_faces, cfg.tables.walls,
                          cfg.tables.landcover, cfg.tables.roofs, cfg.tables.landcover,
                          cfg.type_range.building_min, cfg.type_range.building_max,  )
        cur.execute(sqltext)
        cct_type = cur.fetchall()
        sql_debug(connection)
        connection.commit()
        build_ids[:] = [vt_fill if x[0] is None else x[0] for x in cct_type]

    sqltext = 'SELECT ST_X(center), ST_Y(center), ST_Z(center) FROM "{0}"."{1}" ORDER BY k,j,i'\
              .format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext)
    center_points = cur.fetchall()
    sql_debug(connection)
    connection.commit()
    centers[2, :] = [x[0] - cfg.domain.origin_x for x in center_points]
    centers[1, :] = [x[1] - cfg.domain.origin_y for x in center_points]
    centers[0, :] = [x[2] for x in center_points]
    del center_points

    debug('Selecting k,j,i and offk, offj, offi from slanted faces into static driver')
    if cfg.has_buildings:
        sqltext = 'SELECT * FROM ( ' \
                  '    SELECT gg.nz AS koff, gg.j AS joff, gg.i AS ioff, ' \
                  '           s.k AS k,     s.j AS j,     s.i AS i ' \
                  '    FROM "{0}"."{1}" AS s ' \
                  '    JOIN LATERAL (SELECT i,j,nz FROM "{0}"."{2}" AS g ' \
                  '                  WHERE s.lid = g.lid AND ST_DWithin(g.geom, s.center, {5}) ' \
                  '                  ORDER BY ST_Distance(g.geom, s.center) ' \
                  '                  LIMIT 1) AS gg ON TRUE ' \
                  '	   WHERE isterr ' \
                  '    UNION ALL' \
                  '    SELECT bb.k AS koff, bb.j AS joff, bb.i AS ioff, ' \
                  '           s.k AS k,     s.j AS j,     s.i AS i ' \
                  '           FROM "{0}"."{1}" AS s ' \
                  '           JOIN LATERAL (SELECT i,j,k FROM "{0}"."{4}" AS b ' \
                  '                         WHERE ST_DWithin(b.geom, s.center, {5}) ' \
                  '                         ORDER BY ST_Distance(b.geom, s.center) ' \
                  '                         LIMIT 1 ' \
                  '                         ) AS bb ON TRUE ' \
                  '    WHERE iswall' \
                  '    UNION ALL' \
                  '    SELECT bb.k AS koff, bb.j AS joff, bb.i AS ioff, ' \
                  '           s.k AS k,     s.j AS j,     s.i AS i ' \
                  '           FROM "{0}"."{1}" AS s ' \
                  '           JOIN LATERAL (SELECT i,j,k FROM "{0}"."{4}" AS b ' \
                  '                         WHERE ST_DWithin(b.geom, s.center, {5}) ' \
                  '                         ORDER BY ST_Distance(b.geom, s.geom) ' \
                  '                         LIMIT 1 ' \
                  '                         ) AS bb ON TRUE ' \
                  '    WHERE isroof ' \
                  ') AS js ' \
                  'ORDER BY k, j, i ' \
            .format(cfg.domain.case_schema, cfg.tables.slanted_faces, cfg.tables.grid,
                    cfg.tables.landcover, cfg.tables.buildings_grid, cfg.slanted_pars.off_dist * cfg.domain.dx)
    else:
        sqltext = 'SELECT * FROM ( ' \
                  '    SELECT gg.nz AS koff, gg.j AS joff, gg.i AS ioff, ' \
                  '           s.k AS k,     s.j AS j,     s.i AS i ' \
                  '    FROM "{0}"."{1}" AS s ' \
                  '    JOIN LATERAL (SELECT i, j, nz FROM "{0}"."{2}" AS g ' \
                  '                  WHERE s.lid = g.lid AND ST_DWithin(g.geom, s.center, {4}) ' \
                  '                  ORDER BY ST_Distance(g.geom, s.center) ' \
                  '                  LIMIT 1) AS gg ON TRUE ' \
                  '	   WHERE isterr ' \
                  ') AS js ' \
                  'ORDER BY k, j, i '\
                  .format(cfg.domain.case_schema, cfg.tables.slanted_faces, cfg.tables.grid,
                          cfg.tables.landcover, cfg.slanted_pars.off_dist * cfg.domain.dx)
    cur.execute(sqltext)
    kji_offs_points = cur.fetchall()
    sql_debug(connection)
    connection.commit()
    kji_locs[0, :] = [x[3] for x in kji_offs_points]
    kji_locs[1, :] = [x[4] for x in kji_offs_points]
    kji_locs[2, :] = [x[5] for x in kji_offs_points]
    offs[0, :] = [x[0] - x[3] for x in kji_offs_points]
    offs[1, :] = [x[1] - x[4] for x in kji_offs_points]
    offs[2, :] = [x[2] - x[5] for x in kji_offs_points]

    del kji_offs_points

    sqltext = 'SELECT normz, normy, normx FROM "{0}"."{1}" ORDER BY k,j,i'\
              .format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext)
    normals_vecs = cur.fetchall()
    sql_debug(connection)
    connection.commit()
    for i in range(3):
        normals[i, :] = [x[i] for x in normals_vecs]

    sqltext = 'SELECT area FROM "{0}"."{1}" ORDER BY k,j,i'\
              .format(cfg.domain.case_schema, cfg.tables.slanted_faces)
    cur.execute(sqltext)
    area_vecs = cur.fetchall()
    sql_debug(connection)
    connection.commit()
    area[:] = [x[0] for x in area_vecs]




    # TODO: complete parameters
    # nc_sl_pars = ncfile.createVariable('slanted_params','f4', ('num_pars', 'num_faces',))
    # nc_sl_pars[:] = -9999.0
    """SELECT s.albedo, s.capacity_surf
    FROM dejvice_01_slanted_code.slanted_faces AS f
    LEFT OUTER JOIN dejvice_01_slanted_code.landcover AS l ON l.lid = f.lid
    LEFT OUTER JOIN dejvice_01_slanted_code.walls AS w ON w.wid = f.wid
    LEFT OUTER JOIN dejvice_01_slanted_code.roofs AS r ON r.rid = f.rid
    LEFT OUTER JOIN dejvice_01_slanted_code.surface_params AS s ON s.code = l.katland OR s.code = r.katroof OR s.code = w.stenakatd
    """

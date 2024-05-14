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

from math import ceil
from config.logger import *

def process_trees(cfg, connection, cur):
    """ Process point trees into grid structure """
    debug('Processing 3D grided structure of the trees')
    change_log_level(cfg.logs.level_trees)

    debug('create new table trees_grid')
    sqltext = 'CREATE TABLE "{0}"."{1}" AS SELECT id, i, j from "{0}"."{2}" with data '\
              .format(cfg.domain.case_schema, cfg.tables.trees_grid, cfg.tables.grid)
    cur.execute(sqltext)
    sql_debug(connection)
    verbose('Add primary key and unique index')
    sqltext = 'ALTER TABLE "{0}"."{1}" add primary key (id)'\
              .format(cfg.domain.case_schema, cfg.tables.trees_grid)
    cur.execute(sqltext)
    sqltext = 'create unique index {1}_i_j on "{0}"."{1}" (i asc, j asc)'\
              .format(cfg.domain.case_schema, cfg.tables.trees_grid)
    cur.execute(sqltext)
    sql_debug(connection)
    debug('Get max height of trees')
    sqltext = 'select max(vysstr) from "{0}"."{1}"'.format(cfg.domain.case_schema, cfg.tables.trees)
    cur.execute(sqltext)
    sql_debug(connection)
    thm = cur.fetchone()[0]
    nzlad = ceil(thm / cfg.domain.dz) + 1
    debug('nzlad: {}', nzlad)
    sqltext = 'ALTER TABLE "{0}"."{1}" '.format(cfg.domain.case_schema, cfg.tables.trees_grid)
    for i in range(nzlad):
        sqltext += ', ' if i > 0 else ''
        sqltext += 'add "lad_{0}" double precision default 0, '.format(i)
        sqltext += 'add "bad_{0}" double precision default 0 '.format(i)
    cur.execute(sqltext)
    sql_debug(connection)
    debug('Call palm_tree_grid routine {}', [cfg.domain.case_schema, cfg.tables.trees, cfg.tables.trees_grid,
                                    cfg.tables.grid, cfg.tables.buildings_grid,
                                    cfg.domain.dz, cfg.trees.nhv, cfg.trees.nump, cfg.trees.lad_reduction,
                                    cfg.trees.bad_coef, cfg.trees.ext_coef, cfg.logs.level_trees])
    cur.callproc('palm_tree_grid', [cfg.domain.case_schema, cfg.tables.trees, cfg.tables.trees_grid,
                                    cfg.tables.grid, cfg.tables.buildings_grid,
                                    cfg.domain.dz, cfg.trees.nhv, cfg.trees.nump, cfg.trees.lad_reduction,
                                    cfg.trees.bad_coef, cfg.trees.ext_coef, cfg.logs.level_trees])
    sql_debug(connection)
    connection.commit()
    sqltext = 'ALTER TABLE "{0}"."{1}" OWNER TO {2}'\
              .format(cfg.domain.case_schema, cfg.tables.trees_grid, cfg.pg_owner)
    cur.execute(sqltext)
    sql_debug(connection)
    restore_log_level(cfg)
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

"""
This script is designed to preprocess tree leaf area density.
Inputs are satellite images of LAI and tree height.
"""

import psycopg2
from argparse import ArgumentParser
import getpass
import os
import sys
import pandas as pd
from argparse import ArgumentParser
import getopt
from config.logger import *
from config.config import load_config, cfg
from sqlalchemy import create_engine
from centerline.geometry import Centerline
from shapely import wkb
import geopandas as gpd

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

# create connection to the postgresql server
if cfg.pg_password is None or cfg.pg_password == '':
    pg_password = getpass.getpass()
connection = psycopg2.connect(database=cfg.database, host=cfg.pg_host,
                              port=cfg.pg_port, user=cfg.pg_user, password=cfg.pg_password)
connection.set_client_encoding('UTF8')
cur = connection.cursor()
sql_debug(connection)

progress('Script start')

"""
--create table aprague_slurm.lai as table inputs_prague_local.lai;
--drop table aprague_slurm.canopy_height;
--create table aprague_slurm.canopy_height as table inputs_prague_local.canopy_height;
--create table aprague_slurm.canopy_height as select rid, ST_Transform(rast, 32633) as rast from inputs_prague_local.canopy_height;
--select * from aprague_slurm.canopy_height;


drop table if exists aprague_slurm.lad_grid;
CREATE TABLE aprague_slurm.lad_grid AS 
SELECT g.id, g.i, g.j, 0.0 as lai, 0.0 as c_height, g.geom as geom, g.xcen as xcen, g.ycen as ycen
FROM aprague_slurm.grid g
--where g.id = 35817
--limit 100
;


WITH lai as (
select 
	lg.id as id, AVG(r.lai) as lai
from aprague_slurm.lad_grid lg
JOIN LATERAL ( SELECT (ST_Intersection(rast, lg.geom)).val AS lai 
			   FROM aprague_slurm.lai WHERE ST_Intersects(rast, lg.geom)) r on true 
group by lg.id			   
)
update aprague_slurm.lad_grid lg
set lai = lai.lai
from lai
where lai.id = lg.id;

WITH ch as (
select 
	lg.id as id, AVG(r.height) as height
from aprague_slurm.lad_grid lg
JOIN LATERAL ( SELECT (ST_Intersection(rast, lg.geom)).val AS height 
			   FROM aprague_slurm.canopy_height WHERE ST_Intersects(rast, lg.geom)) r on true 
group by lg.id			   
)
update aprague_slurm.lad_grid lg
set c_height = ch.height
from ch
where ch.id = lg.id;
"""
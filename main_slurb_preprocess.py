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
This script is designed to preprocess open source data into slurb readable format.
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

progress('Create centerline')
cur.execute('drop table if exists "{}"."{}"'.format(cfg.domain.case_schema, cfg.tables.centerline))
connection.commit()

engine = create_engine(f"postgresql+psycopg2://{cfg.pg_user}:{cfg.pg_password}@{cfg.pg_host}:{cfg.pg_port}/{cfg.database}")

debug('Read landcover shapefiles into pandas dataframe')
sql = """SELECT geom as geom FROM "{0}".landcover where type = 202;""".format(cfg.domain.case_schema, cfg.tables.ladcover)
# sql = """SELECT ST_Union(geom) as geom FROM "{0}".landcover where type = 202;"""
parcels = pd.read_sql(sql, engine)

debug('Transform into geopandas')
parcels.geom = wkb.loads(parcels.geom, hex=True)

debug('Loop over all shapefiles')
centerlines = []
for index, row in parcels.iterrows():
    verbose('Looping over {}/{}', index, parcels.shape[0])
    geom = row['geom']
    try:
        centerline = Centerline(geom) #parcels.geom.iloc[0]) #, **attributes)
        centerlines.append(centerline.geometry)
    except:
        pass

debug('Upload centerline back into servers')
gdf = gpd.GeoDataFrame(geometry=centerlines) #[centerline.geometry])
gdf.to_postgis(name=cfg.tables.centerline,
               schema=cfg.domain.case_schema,
               con=engine,
               index=False)

debug('Creation of centerline is done')

progress('Process centerline')
debug('Simplify centerline and make some corrections')
sqltext = """
DROP TABLE IF EXISTS "{0}".centerline_simplified;
CREATE TABLE "{0}".centerline_simplified as 
SELECT 
	geometry as geom_original,
	ST_Simplify(ST_LineMerge(geometry), 0.5) as geom
FROM "{0}".centerline_sample_segments;
"""

debug('Change ownership')
sqltext = """
ALTER TABLE "{0}".centerline_simplified OWNER TO palm;	
"""

debug('Cut centerline into segments')
sqltext = """ 
DROP TABLE IF EXISTS "{0}".centerline_simplified_segments;
CREATE TABLE "{0}".centerline_simplified_segments as 
SELECT 
	(ST_Dump(geom)).*
FROM "{0}".centerline_simplified;

alter table "{0}"."centerline_simplified_segments" owner to palm
"""

debug('Create perpendicular segments')
sqltext = """
DROP TABLE IF EXISTS "{0}".centerline_simplified_segments_hairy;
CREATE TABLE "{0}".centerline_simplified_segments_hairy as 
WITH
	geodata AS (SELECT row_number() over() AS id, geom FROM "{0}".centerline_simplified_segments WHERE ST_Length(geom)>40),
	linecut AS (SELECT id, ST_LineSubstring(d.geom, substart, CASE WHEN subend > 1 THEN 1 ELSE subend END) geom
				FROM (SELECT id, geom, ST_Length(((geom)::geometry)) len, 15 sublen FROM geodata) AS d
					CROSS JOIN LATERAL (SELECT i,  (sublen * i)/len AS substart, (sublen * (i+1))/len AS subend
										FROM generate_series(0, floor( d.len / sublen)::integer) AS t(i) 
										WHERE (sublen * i)/len <> 1.0) 
				AS d2),
	rotate AS (SELECT id, (ST_Rotate(ST_Collect(geom), -pi()/2, ST_Centroid(geom))) geom FROM linecut GROUP BY id, geom),
	tbld AS (SELECT id, (ST_Dump(geom)).geom geom FROM rotate),
    bl AS (SELECT (ST_Dump(ST_MakeLine(ST_StartPoint(geom), ST_EndPoint(geom)))).geom AS geom FROM tbld)
SELECT ST_MakeLine(ST_TRANSLATE(a, sin(az2) * len, cos(az2) * len), ST_Centroid(ST_Collect(a, b))) as geom_1, 
	   ST_MakeLine(ST_Centroid(ST_Collect(a, b)), ST_TRANSLATE(b,sin(az1) * len, cos(az1) * len)) as geom_2,
	   ST_Centroid(ST_Collect(a, b)) as geom_center
  FROM (
    SELECT a, b, ST_Azimuth(a,b) AS az1, ST_Azimuth(b, a) AS az2, 50 AS len
      FROM (
        SELECT ST_StartPoint(geom) AS a, ST_EndPoint(geom) AS b
          FROM bl
    ) AS sub
) AS sub2
;

create index if not exists hairy_geom_idx on "{0}".centerline_simplified_segments_hairy using gist(geom_1);
create index if not exists hairy_geom_idx on "{0}".centerline_simplified_segments_hairy using gist(geom_2);
alter table "{0}".centerline_simplified_segments_hairy add column if not exists gid serial;

alter table "{0}".centerline_simplified_segments_hairy owner to palm;

"""

debug('Calculate intersection points with buildings for street canyon width calculation')
sqltext = """
DROP TABLE IF EXISTS "{0}".centerline_simplified_intersections;
CREATE TABLE "{0}".centerline_simplified_intersections as 
select 
	cs.*,
	l_1.*,
	l_2.*,
	ST_Distance(l_1.point_1, l_2.point_2),
	ST_Collect(array[cs.geom_center, cs.geom_1, cs.geom_2, l_1.point_1, l_2.point_2]) as geom_col
from "{0}".centerline_simplified_segments_hairy cs
	join lateral (select (ST_Dump(ST_Intersection(cs.geom_1, ST_boundary(l.geom)))).geom as point_1
				  from "{0}".landcover l
				  where l.type between 900 and 999 
				        and ST_Intersects(cs.geom_1, l.geom)
				  order by ST_Distance(cs.geom_center, (ST_Dump(ST_Intersection(cs.geom_1, ST_boundary(l.geom)))).geom) asc
				  limit 1) as l_1 on true
	join lateral (select (ST_Dump(ST_Intersection(cs.geom_2, ST_boundary(l.geom)))).geom as point_2
				  from "{0}".landcover l
				  where l.type between 900 and 999 
				        and ST_Intersects(cs.geom_2, l.geom)
				  order by ST_Distance(cs.geom_center, (ST_Dump(ST_Intersection(cs.geom_2, ST_boundary(l.geom)))).geom) asc
				  limit 1) as l_2 on true;
alter table "{0}".centerline_simplified_intersections owner to palm;			  
"""

progress('Calculate street canyon width and orientation')
sqltext = """
ROP TABLE IF EXISTS "{0}".centerline_results;
CREATE TABLE "{0}".centerline_results as 
select
	gid,
	ST_Azimuth(point_1, point_2) * 180 / 3.14 + 90.0 as orientation,
	ST_Distance(point_1, point_2) as width,
	b1.val as val_1,
	b2.val as val_2,
	geom_center as geom
	--point_1 as point_1,
	--point_2 as point_2
from "{0}".centerline_simplified_intersections
	--left join lateral (select ST_NearestValue(rast, ST_Transform(point_1, 3035)) as val from "{0}".imported_buildings where st_intersects(rast, ST_Transform(point_1, 3035))) b1 on true
	--left join lateral (select ST_NearestValue(rast, ST_Transform(point_2, 3035)) as val from "{0}".imported_buildings where st_intersects(rast, ST_Transform(point_2, 3035))) b2 on true
	left join lateral (select bp.val as val from build_pixel bp where ST_DWithin(bp.geom, point_1, 100.0) order by ST_Distance(bp.geom, point_1) asc limit 1) b1 on true
	left join lateral (select bp.val as val from build_pixel bp where ST_DWithin(bp.geom, point_2, 100.0) order by ST_Distance(bp.geom, point_2) asc limit 1) b2 on true;

alter table "{0}".centerline_results owner to palm;
"""

progress('Calculate building front area')
sqltext = """
drop table if exists outer_wall_full;
create temp table outer_wall_full as 
select 
	ST_ForceRHR((ST_Dump(ST_UNION(ST_Buffer(l.geom, 0.001)))).geom) as geom
FROM (select geom, type from "{0}".landcover where type between 900 and 999) AS l
	where l.type between 900 and 999; 

alter table outer_wall_full add column id serial;

drop table if exists build_transf;
create temp table build_transf as 
select rid, ST_Transform(rast, 32633) as rast
from "{0}".imported_buildings;


DROP TABLE IF EXISTS "{0}".outer_wall;
CREATE TABLE "{0}".outer_wall AS
with outer_wall as (
	SELECT 
		id, 
		ST_Boundary(geom) AS geom 
	FROM outer_wall_full
)
SELECT 
	id, 
	ST_LineSubstring(d.geom, substart, CASE WHEN subend > 1 THEN 1 ELSE subend END) geom,
	ST_Centroid(ST_LineSubstring(d.geom, substart, CASE WHEN subend > 1 THEN 1 ELSE subend END)) as point,
	null :: numeric as building_height
	--20.0 as building_height
FROM (SELECT id, geom, ST_Length(((geom)::geometry)) len, 15 sublen FROM outer_wall) AS d
	CROSS JOIN LATERAL (SELECT i,  (sublen * i)/len AS substart, (sublen * (i+1))/len AS subend
						 FROM generate_series(0, floor( d.len / sublen)::integer) AS t(i) 
						 WHERE (sublen * i)/len <> 1.0) AS d2;
			  
update "{0}".outer_wall 
set building_height = (SELECT ST_NearestValue(rast, point) AS height 
                       FROM build_transf
					   WHERE ST_Intersects(rast, point)
					   LIMIT 1);

update "{0}".outer_wall 
set building_height = (select bp.val as val 
			 		   from build_pixel bp
					   where ST_DWithin(bp.geom, point, 100.0) 
					   order by ST_Distance(bp.geom, point) asc 
					   limit 1)
where building_Height is null;

update "{0}".outer_wall 
set building_height = 20
where building_Height is null;

DROP TABLE IF EXISTS "{0}".building_area;
CREATE TABLE "{0}".building_area AS
select 
	id,
	sum(building_height * ST_Length(geom)) as wall_area,
	null :: numeric as roof_area,
	null :: geometry as geom
from "{0}".outer_wall
group by id;

update "{0}".building_area ba
set (roof_area, geom) = 
(select ST_Area(geom), geom
from outer_wall_full owf
where ba.id = owf.id);

alter table inputs_prague_local.building_area add primary key (id)
"""

progress('Cleaning up')
debug('Delete uncessary tables if configured')

"""


"""
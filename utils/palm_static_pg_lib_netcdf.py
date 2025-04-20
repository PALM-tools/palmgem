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

from netCDF4 import Dataset
from config.logger import *
import sys
from datetime import datetime
import numpy as np
from math import ceil
from utils.visualization import variable_visualization
from pathlib import Path

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

def nc_write_global_attributes_slurb(ncfile, cfg):
    """ Write attributes for SLURb driver """
    ncfile = nc_write_global_attributes(ncfile, cfg)

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

def create_slurb_dims(ncfile, cfg, connection, cur):
    """ Create SLURb dimensions """
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

    debug('Create rest of the dims')
    for name in cfg.ndims_slurb._settings.keys():
        d = cfg.ndims_slurb[name]
        nc_create_dimension(ncfile, name, d)
        temp = ncfile.createVariable(name, 'i', name)
        temp[:] = np.asarray(list(range(d)))

def create_slurb_vars(ncfile, cfg, connection, cur):
    """ Create all necessary slurb related field """
    # Calculate slurb grid
    debug('Create slurb grid')
    sqltext = """
    drop table if exists "{0}"."{1}";
    create table "{0}"."{1}" as 
    select 
        g.id, 
        g.i, 
        g.j, 
        g.xcen,
        g.ycen,
        g.geom, 
        0.0 as building_plan_area_fraction,
        0.0 as urban_fraction,
        273.15 as deep_soil_temperature,
        {3} :: double precision as building_height,
        273.15 as building_indoor_temperature,
        1 as building_type,
        1 as pavement_type,
        null :: double precision as building_frontal_area_fraction,
        null :: double precision as street_canyon_aspect_ratio,
        null :: double precision as street_canyon_orientation
    from "{0}"."{2}" g
    """.format(cfg.domain.case_schema, cfg.tables.grid_slurb, cfg.tables.grid,
               cfg.default_building_height)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    debug('Calculate building fraction')
    sqltext = (('UPDATE "{0}"."{1}" g SET '
                'building_plan_area_fraction = least(0.988, coalesce(s.sum_area)) '
                'FROM ( '
                '       SELECT g.id AS gid, SUM(ST_Area(ST_Intersection(g.geom, l.geom))) / {5} AS sum_area '
                '       FROM "{0}"."{1}" g'
                '       JOIN "{0}"."{2}" l ON ST_Intersects(l.geom, g.geom) '
                '       WHERE l.type BETWEEN {3} AND {4} '
                '       GROUP BY g.id '
                ') AS s '
                'WHERE g.id = s.gid;').format(cfg.domain.case_schema, cfg.tables.grid_slurb, cfg.tables.landcover,
                                                          900, 999, cfg.domain.dx * cfg.domain.dy))
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    debug('Calculation of urban_fraction')
    sqltext = (('UPDATE "{0}"."{1}" g SET '
                'urban_fraction = least(1.0, s.sum_area) '
                'FROM ( '
                '       SELECT g.id AS gid, SUM(ST_Area(ST_Intersection(g.geom, l.geom))) / {5} AS sum_area '
                '       FROM "{0}"."{1}" g'
                '       JOIN "{0}"."{2}" l ON ST_Intersects(l.geom, g.geom) '
                '       WHERE l.type BETWEEN {3} AND {4} '
                '          or l.type between {6} and {7}'
                '       GROUP BY g.id '
                ') AS s '
                'WHERE g.id = s.gid;').format(cfg.domain.case_schema, cfg.tables.grid_slurb, cfg.tables.landcover,
                                                          900, 999, cfg.domain.dx * cfg.domain.dy,
                                              200, 299))
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # debug('Remove grid cell with lower building fraction')
    sqltext = """
    delete from "{0}"."{1}" gs
    where gs.building_plan_area_fraction < {2}    
    """.format(cfg.domain.case_schema, cfg.tables.grid_slurb, cfg.min_plan_area)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # sqltext = """
    # update "{0}"."{1}"
    # set building_plan_area_fraction = 0.0
    # where building_plan_area_fraction < {2}
    # """.format(cfg.domain.case_schema, cfg.tables.grid_slurb, cfg.min_plan_area)
    # cur.execute(sqltext)
    # sql_debug(connection)
    # connection.commit()

    # Now include all relevant information inside this grid_slurb table

    # deep_soil_temperature
    debug('Updating deep soil temperature')
    sqltext = """
    update "{0}"."{1}" gs 
    set deep_soil_temperature = {2}
    """.format(cfg.domain.case_schema, cfg.tables.grid_slurb, cfg.deep_soil_temperature_default)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # building_height
    debug('Calculation of building height in grid cell')
    if cfg.has_buildings:
        sqltext = """
        with grid_height as (
            select 
                gs.id as id,
                a.height
            from "{0}"."{1}" gs
                join lateral (select AVG(ST_NearestValue(rast, ST_SetSRID(ST_Point(gs.xcen,gs.ycen), %s))) AS height
                              FROM "{0}"."{2}"
                              WHERE ST_Intersects(rast, gs.geom)) a on true
        )
        update "{0}"."{1}" gs 
        set building_height = gh.height
        from grid_height gh
        where gh.id = gs.id
        """.format(cfg.domain.case_schema, cfg.tables.grid_slurb, cfg.tables.buildings_height)
        cur.execute(sqltext, (cfg.srid_palm,))
        sql_debug(connection)
        connection.commit()

    # building_indoor_temperature
    debug('Updating building indoor temperature')
    sqltext = """
    update "{0}"."{1}" gs 
    set building_indoor_temperature = {2}
    """.format(cfg.domain.case_schema, cfg.tables.grid_slurb, cfg.building_indoor_temperature_default)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # building_type
    debug('Calculation of building type')
    sqltext = """
    with grid_type as (
        select 
            gs.id as id,
            a.g_type - {2} as g_type
        from "{0}"."{1}" gs
            join lateral (select 
                             PERCENTILE_CONT(0.5) WITHIN GROUP(ORDER BY l.type) as g_type
                             --l.type as g_type
                          FROM "{0}".{4} l 
                          WHERE ST_Intersects(l.geom, gs.geom)
                            and l.type between {2} and {3}) a on true
    )
    update "{0}"."{1}" gs 
    set building_type = gt.g_type
    from grid_type gt
    where gt.id = gs.id
    """.format(cfg.domain.case_schema, cfg.tables.grid_slurb, cfg.type_range.building_min, cfg.type_range.building_max,
               cfg.tables.landcover)
    cur.execute(sqltext, (cfg.srid_palm, ))
    sql_debug(connection)
    connection.commit()

    # pavement_type
    debug('Calculation of pavement_type')
    sqltext = """
    with grid_type as (
        select 
            gs.id as id,
            a.g_type - {2} as g_type
        from "{0}"."{1}" gs
            join lateral (select 
                             PERCENTILE_CONT(0.5) WITHIN GROUP(ORDER BY l.type) as g_type
                             --l.type as g_type
                          FROM "{0}".{4} l 
                          WHERE ST_Intersects(l.geom, gs.geom)
                            and l.type between {2} and {3}) a on true
    )
    update "{0}"."{1}" gs 
    set pavement_type = 1 --gt.g_type
    from grid_type gt
    where gt.id = gs.id
    """.format(cfg.domain.case_schema, cfg.tables.grid_slurb, cfg.type_range.pavement_min, cfg.type_range.pavement_max,
               cfg.tables.landcover)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # building_frontal_area_fraction
    sqltext = """
    with fraction as (
        select
            gs.id, 
            sum(st_area(st_intersection(ba.geom, gs.geom)) / ba.roof_area * (wall_area + roof_area) / {3}) as f_area
        from "{0}"."{1}" gs
            join "{0}"."{2}" ba on st_intersects(ba.geom, gs.geom) 
        group by gs.id
    )
    update "{0}"."{1}" gs 
    set building_frontal_area_fraction  = least(f.f_area, 0.99)
    from fraction f
    where f.id = gs.id
    """.format(cfg.domain.case_schema, cfg.tables.grid_slurb, cfg.tables.building_area,
               cfg.domain.dx * cfg.domain.dy)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # street_canyon_aspect_ratio
    sqltext = """
    with subquery as (
        select 
            gs.id as id,
            a.hw
        from "{0}"."{1}" gs
            join lateral (select 
                             AVG(((c.val_1 + c.val_2) / 2) / c.width) as hw
                          FROM "{0}"."{2}" c 
                          WHERE ST_Intersects(c.geom, gs.geom)) a on true
    )
    update "{0}"."{1}" gs 
    set street_canyon_aspect_ratio  = s.hw
    from subquery s
    where s.id = gs.id
    """.format(cfg.domain.case_schema, cfg.tables.grid_slurb, cfg.tables.centerline)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # Update street canyon where not filled, take from nearest
    sqltext = f"""
    update "{cfg.domain.case_schema}"."{cfg.tables.grid_slurb}" gs
    set street_canyon_aspect_ratio = (select street_canyon_aspect_ratio 
                                      from "{cfg.domain.case_schema}"."{cfg.tables.grid_slurb}" gss
                                      where gss.street_canyon_aspect_ratio is not null
                                      order by ST_Distance(gss.geom, gs.geom) 
                                      limit 1)
    where street_canyon_aspect_ratio is null
    """
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # street_canyon_orientation
    sqltext = """
    with subquery as (
        select 
            gs.id as id,
            a.orientation
        from "{0}"."{1}" gs
            join lateral (select 
                             ATAN2(AVG(SIN(orientation * PI() / 180.0)), 
                                   AVG(COS(orientation * PI() / 180.0))
                                   ) * 180.0 / PI() as orientation
                          FROM "{0}"."{2}" c 
                          WHERE ST_Intersects(c.geom, gs.geom)) a on true
    )
    update "{0}"."{1}" gs 
    set street_canyon_orientation = 
        case when s.orientation < 0.0 then s.orientation + 360.0
             when s.orientation > 360.0 then s.orientation - 360.0
             else s.orientation
             end
    from subquery s
    where s.id = gs.id
    """.format(cfg.domain.case_schema, cfg.tables.grid_slurb, cfg.tables.centerline)
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # Update street canyon where not filled, take from nearest
    sqltext = f"""
    update "{cfg.domain.case_schema}"."{cfg.tables.grid_slurb}" gs
    set street_canyon_orientation = (select street_canyon_orientation 
                                      from "{cfg.domain.case_schema}"."{cfg.tables.grid_slurb}" gss
                                      where gss.street_canyon_orientation is not null
                                      order by ST_Distance(gss.geom, gs.geom) 
                                      limit 1)
    where street_canyon_orientation is null
    """
    cur.execute(sqltext)
    sql_debug(connection)
    connection.commit()

    # Urban fraction must be filled everywhere
    debug('Filling urban fraction')
    vn = 'urban_fraction'
    vt = cfg.slurb_vars[vn].vt
    long_name = cfg.slurb_vars[vn].long_name

    sqltext = """
        select coalesce(sg.{3}, 0.0) 
        from "{0}"."{1}" g
            left join "{0}"."{2}" sg on sg.id = g.id 
        order by g.j, g.i
    """.format(cfg.domain.case_schema, cfg.tables.grid, cfg.tables.grid_slurb, vn)
    cur.execute(sqltext)
    sql_debug(connection)

    write_type_variable(ncfile, cfg, vn, long_name, vt, cfg.fill_values, 0, connection, cur)

    # Download to netcdf format and insert into driver
    progress('Insert slurb variables into netcdf file')
    for vn in cfg.slurb_vars_done:
        debug('Processing slurb variable: {}', vn)
        vt = cfg.slurb_vars[vn].vt
        long_name = cfg.slurb_vars[vn].long_name

        sqltext = """
            select sg.{3} 
            from "{0}"."{1}" g
                left join "{0}"."{2}" sg on sg.id = g.id 
            order by g.j, g.i
        """.format(cfg.domain.case_schema, cfg.tables.grid, cfg.tables.grid_slurb, vn)
        cur.execute(sqltext)
        sql_debug(connection)

        write_type_variable(ncfile, cfg, vn, long_name, vt, cfg.fill_values, 0, connection, cur)

    # Fetch mask
    sqltext = """
        select 
            case when sg.id is not null then True else False end
        from "{0}"."{1}" g
            left join "{0}"."{2}" sg on sg.id = g.id 
        order by g.j, g.i
    """.format(cfg.domain.case_schema, cfg.tables.grid, cfg.tables.grid_slurb)
    cur.execute(sqltext)
    res = cur.fetchall()
    sql_debug(connection)
    # res = [False if x[0] is None else True for x in res]
    res = [x[0] for x in res]
    slurb_mask = np.reshape(np.asarray([x for x in res], dtype='bool_'), (cfg.domain.ny, cfg.domain.nx))

    # Has xy dim
    yx_dims = [
     ['albedo_road', 0.10],
     ['albedo_roof', 0.17],
     ['albedo_wall', 0.12],
     ['albedo_window', 0.12],
     ['emiss_road', 0.90],
     ['emiss_roof', 0.90],
     ['emiss_wall', 0.8],
     ['emiss_window', 0.8],
     ['window_fraction', 0.1],
     ['z0_road', 0.002],
     ['z0_roof', 0.002],
     ['z0_urb', 0.002],
     ['z0_wall', 0.002],
     ['z0_window', 0.00001],
     ['z0h_roof', 0.002],
     ['z0h_road', 0.002],
     ['z0h_urb', 0.002],
     ['z0h_wall', 0.002],
     ['z0h_window', 0.00001],
    ]

    for vn, def_val in yx_dims:
        vt = cfg.slurb_vars[vn].vt
        long_name = cfg.slurb_vars[vn].long_name
        vn_units = cfg.slurb_vars[vn].units

        ncfile.createVariable(vn, vt, ('y', 'x'), fill_value=cfg.fill_values[vt])
        fill_var = np.ones((cfg.domain.ny, cfg.domain.nx)) * cfg.fill_values[vt]
        fill_var[slurb_mask] = def_val
        ncfile[vn][...] = fill_var
        # ncfile[vn][...] = slurb_mask * cfg.fill_values[vt]

        nc_write_attribute(ncfile, vn, 'long_name', long_name)
        nc_write_attribute(ncfile, vn, 'units', vn_units)
        nc_write_attribute(ncfile, vn, 'res_orig', cfg.domain.dz)
        nc_write_attribute(ncfile, vn, 'coordinates', 'E_UTM N_UTM lon lat')
        nc_write_attribute(ncfile, vn, 'grid_mapping', 'E_UTM N_UTM lon lat')

    # other dims
    # nroadyx
    nroadyx_dims = [
     ['c_road', 1.5e6],
     ['dz_road', 0.5],
     ['lambda_road', 0.12]
    ]

    for vn, def_val in nroadyx_dims:
        vt = cfg.slurb_vars[vn].vt
        long_name = cfg.slurb_vars[vn].long_name
        vn_units = cfg.slurb_vars[vn].units

        ncfile.createVariable(vn, vt, ('nroad_3d', 'y', 'x'), fill_value=cfg.fill_values[vt])
        fill_var = np.ones((cfg.domain.ny, cfg.domain.nx)) * cfg.fill_values[vt]
        fill_var[slurb_mask] = def_val
        for idim in range(cfg.ndims_slurb.nroad_3d):
            ncfile[vn][idim, :, :] = fill_var

        nc_write_attribute(ncfile, vn, 'long_name', long_name)
        nc_write_attribute(ncfile, vn, 'units', vn_units)
        nc_write_attribute(ncfile, vn, 'res_orig', cfg.domain.dz)
        nc_write_attribute(ncfile, vn, 'coordinates', 'E_UTM N_UTM lon lat')
        nc_write_attribute(ncfile, vn, 'grid_mapping', 'E_UTM N_UTM lon lat')

    # nroofyx
    nroofyx_dims = [
     ['c_roof', 1.52e6],
     ['dz_roof', 0.5],
     ['lambda_roof', 0.12],
    ]

    for vn, def_val in nroofyx_dims:
        vt = cfg.slurb_vars[vn].vt
        long_name = cfg.slurb_vars[vn].long_name
        vn_units = cfg.slurb_vars[vn].units

        ncfile.createVariable(vn, vt, ('nroof_3d', 'y', 'x'), fill_value=cfg.fill_values[vt])
        fill_var = np.ones((cfg.domain.ny, cfg.domain.nx)) * cfg.fill_values[vt]
        fill_var[slurb_mask] = def_val
        for idim in range(cfg.ndims_slurb.nroof_3d):
            ncfile[vn][idim, :, :] = fill_var

        nc_write_attribute(ncfile, vn, 'long_name', long_name)
        nc_write_attribute(ncfile, vn, 'units', vn_units)
        nc_write_attribute(ncfile, vn, 'res_orig', cfg.domain.dz)
        nc_write_attribute(ncfile, vn, 'coordinates', 'E_UTM N_UTM lon lat')
        nc_write_attribute(ncfile, vn, 'grid_mapping', 'E_UTM N_UTM lon lat')


    # nwallyx
    nwallyx_dims = [
     ['c_wall', 1.5e6],
     ['dz_wall', 0.5],
     ['lambda_wall', 0.9],
    ]

    for vn, def_val in nwallyx_dims:
        vt = cfg.slurb_vars[vn].vt
        long_name = cfg.slurb_vars[vn].long_name
        vn_units = cfg.slurb_vars[vn].units

        ncfile.createVariable(vn, vt, ('nwall_3d', 'y', 'x'), fill_value=cfg.fill_values[vt])
        fill_var = np.ones((cfg.domain.ny, cfg.domain.nx)) * cfg.fill_values[vt]
        fill_var[slurb_mask] = def_val
        for idim in range(cfg.ndims_slurb.nwall_3d):
            ncfile[vn][idim, :, :] = fill_var

        nc_write_attribute(ncfile, vn, 'long_name', long_name)
        nc_write_attribute(ncfile, vn, 'units', vn_units)
        nc_write_attribute(ncfile, vn, 'res_orig', cfg.domain.dz)
        nc_write_attribute(ncfile, vn, 'coordinates', 'E_UTM N_UTM lon lat')
        nc_write_attribute(ncfile, vn, 'grid_mapping', 'E_UTM N_UTM lon lat')

    # nwinyx
    nwinyx_dims = [
     ['c_window', 1.7e6],
     ['dz_window', 0.02],
     ['lambda_window', 0.4],
    ]

    for vn, def_val in nwinyx_dims:
        debug(vn)
        vt = cfg.slurb_vars[vn].vt
        long_name = cfg.slurb_vars[vn].long_name
        vn_units = cfg.slurb_vars[vn].units

        ncfile.createVariable(vn, vt, ('nwin_3d', 'y', 'x'), fill_value=cfg.fill_values[vt])
        fill_var = np.ones((cfg.domain.ny, cfg.domain.nx)) * cfg.fill_values[vt]
        fill_var[slurb_mask] = def_val
        for idim in range(cfg.ndims_slurb.nwin_3d):
            ncfile[vn][idim, :, :] = fill_var

        nc_write_attribute(ncfile, vn, 'long_name', long_name)
        nc_write_attribute(ncfile, vn, 'units', vn_units)
        nc_write_attribute(ncfile, vn, 'res_orig', cfg.domain.dz)
        nc_write_attribute(ncfile, vn, 'coordinates', 'E_UTM N_UTM lon lat')
        nc_write_attribute(ncfile, vn, 'grid_mapping', 'E_UTM N_UTM lon lat')


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


def write_surface_fractions(ncfile, cfg, connection, cur):
    """ Write surface fraction if enabled """
    vn = 'surface_fraction'
    vt = 'f8'
    debug('Processing surface fractions')
    verbose('Creating dimension')
    ncfile.createVariable(vn, vt, ('nsurface_fraction', 'y', 'x'), fill_value=cfg.fill_values[vt])
    nc_write_attribute(ncfile, vn, 'long_name', 'Surface fractions 0 vegetation, 1 pavement, 2 water')
    nc_write_attribute(ncfile, vn, 'units', '')
    nc_write_attribute(ncfile, vn, 'res_orig', cfg.domain.dz)
    nc_write_attribute(ncfile, vn, 'lod', 0)
    nc_write_attribute(ncfile, vn, 'coordinates', 'E_UTM N_UTM lon lat')
    nc_write_attribute(ncfile, vn, 'grid_mapping', 'E_UTM N_UTM lon lat')

    fractions = [[0, 'veg_fraction', 'vegetation_type', 'veg_fract_type', 'Vegetation Type'],
                 [1, 'pav_fraction', 'pavement_type', 'pav_fract_type', 'Pavement Type'],
                 [2, 'wat_fraction', 'water_type', 'wat_fract_type', 'Water Type']]
    for it, name_fr, type_name, type_fr, type_full_name in fractions:
        verbose('Fraction: {}', name_fr)

        sqltext = ('SELECT case when g.{2} > 0.0 then g.{2} else 0.0 end FROM "{0}"."{1}" g '
                   'ORDER BY g.j, g.i').format(cfg.domain.case_schema, cfg.tables.grid, name_fr)
        cur.execute(sqltext)
        sql_debug(connection)
        res = cur.fetchall()

        res = [cfg.fill_values[vt] if x[0] is None else x[0] for x in res]
        var = np.reshape(np.asarray([x for x in res], dtype=vt), (cfg.domain.ny, cfg.domain.nx))
        del res
        var = np.nan_to_num(var, copy=False, nan=cfg.fill_values[vt], posinf=cfg.fill_values[vt], neginf=cfg.fill_values[vt])
        ncfile[vn][it, ...] = var

    debug('Variable {} ({}) has been written.', vn, 'Surface fraction')




def write_pavements(ncfile, cfg, connection, cur):
    # write landcover
    vn = 'pavement_type'
    vt = 'b'
    if cfg.landcover.surface_fractions:
        sqltext = f"""
            select case when g.pav_fraction > {cfg.landcover.min_fraction} and 
                             g.pav_fract_type is not null then g.pav_fract_type - %s 
                   else null end
            from "{cfg.domain.case_schema}"."{cfg.tables.grid}" g
            order by g.j, g.i
        """
        cur.execute(sqltext, (cfg.type_range.pavement_min,))
        sql_debug(connection)
    else:
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

    if cfg.lod2:
        # Process pavement_pars
        vn = "pavement_pars"
        vt = 'f4'
        ncfile.createVariable(vn, vt, ('npavement_pars', 'y', 'x'), fill_value=cfg.fill_values[vt])
        nc_write_attribute(ncfile, vn, 'long_name', 'pavement parameters')
        nc_write_attribute(ncfile, vn, 'units', '')
        nc_write_attribute(ncfile, vn, 'res_orig', cfg.domain.dz)
        nc_write_attribute(ncfile, vn, 'coordinates', 'E_UTM N_UTM lon lat')
        nc_write_attribute(ncfile, vn, 'grid_mapping', 'E_UTM N_UTM lon lat')
        # fill individual parameters in dimension npavement_pars
        for par in cfg.pavement_pars._settings.keys():
            # par - index of npars dimmension in netcdf
            # pavement_pars[par] - calculation formula for parameter
            sqltext = 'SELECT case when l.type is not null then {0} else null end ' \
                      'from "{1}"."{2}" g    ' \
                      'left outer join "{1}"."{3}" l on l.lid = g.lid and l.type >= %s and l.type < %s ' \
                      'left outer join "{1}"."{4}" p on p.code = l."{5}" ' \
                      'order by g.j, g.i' \
                .format(cfg.pavement_pars[par], cfg.domain.case_schema, cfg.tables.grid,
                        cfg.tables.landcover,   cfg.tables.surface_params,
                        cfg.landcover_params_var)
            cur.execute(sqltext, (cfg.type_range.pavement_min, cfg.type_range.pavement_max,))
            res = cur.fetchall()
            sql_debug(connection)
            res = [cfg.fill_values[vt] if x[0] is None else x[0] for x in res]
            var = np.reshape(np.asarray(res, dtype=vt), (cfg.domain.ny, cfg.domain.nx))
            del res
            var = np.nan_to_num(var, copy=False,
                                nan=cfg.fill_values[vt], posinf=cfg.fill_values[vt], neginf=cfg.fill_values[vt])
            ncfile.variables[vn][par,...] = var
            debug('Variable {} parameter {} written.', vn, par)

            if cfg.visual_check.enabled:
                variable_visualization(var=ncfile[vn][par, ...],
                                       x=np.asarray(ncfile.variables['x']), y=np.asarray(ncfile.variables['y']),
                                       var_name=vn, par_id=par, text_id='pavement_pars', path=cfg.visual_check.path,
                                       show_plots = cfg.visual_check.show_plots)

        # Process pavement_subsurface_pars
        # create zsoil dimenzion
        vn = 'zsoil'
        nc_create_dimension(ncfile, vn, cfg.ground.nzsoil)
        temp = ncfile.createVariable(vn, 'f4', vn)
        zs = np.zeros(cfg.ground.nzsoil)
        zs[0] = cfg.ground.dz_soil[0]
        for i in range(1,cfg.ground.nzsoil):
            zs[i] = zs[i-1] + cfg.ground.dz_soil[i]

        temp[...] = zs
        # create and fill pavement_subsurface_pars
        vn = 'pavement_subsurface_pars'
        vt = 'f4'
        var = ncfile.createVariable(vn, vt, ('npavement_subsurface_pars', 'zsoil', 'y', 'x'),
                                    fill_value=cfg.fill_values[vt])
        nc_write_attribute(ncfile, vn, 'long_name', 'pavement subsurface parameters')
        nc_write_attribute(ncfile, vn, 'units', '')
        nc_write_attribute(ncfile, vn, 'res_orig', cfg.domain.dz)
        nc_write_attribute(ncfile, vn, 'coordinates', 'E_UTM N_UTM lon lat')
        nc_write_attribute(ncfile, vn, 'grid_mapping', 'E_UTM N_UTM lon lat')
        # assignment of the upper and lower layers of pavement
        lrange = [range(0, cfg.ground.nzsoil_surface), range(cfg.ground.nzsoil_surface, cfg.ground.nzsoil)]
        # fill individual parameters in dimension npavement_subsurface_pars
        for par in cfg.pavement_subsurface_pars._settings.keys():
            # par - index of npars dimmension in netcdf
            # pavement_subsurface_pars[par] - calculation formula for parameters of upper and lower layers
            for k in range(0, 2):
                # 0 - param for upper layers, 1 - param for lower layers
                sqltext = 'SELECT case when l.type is not null then {0} else null end ' \
                          'from "{1}"."{2}" g    ' \
                          'left outer join "{1}"."{3}" l on l.lid = g.lid and l.type >= %s and l.type < %s ' \
                          'left outer join "{1}"."{4}" p on p.code = l."{5}" ' \
                          'order by g.j, g.i' \
                    .format(cfg.pavement_subsurface_pars[par][k], cfg.domain.case_schema, cfg.tables.grid,
                            cfg.tables.landcover, cfg.tables.surface_params, cfg.landcover_params_var)
                cur.execute(sqltext, (cfg.type_range.pavement_min, cfg.type_range.pavement_max,))

                res = cur.fetchall()
                sql_debug(connection)
                varp = np.reshape(np.asarray([x[0] for x in res], dtype=vt), (cfg.domain.ny, cfg.domain.nx))
                del res
                varp = np.nan_to_num(varp, copy=False, nan=cfg.fill_values[vt],
                                     posinf=cfg.fill_values[vt], neginf=cfg.fill_values[vt])
                # fill parameter for corresponding soil layers
                for p in lrange[k]:
                    var[par, p, :, :] = varp
                    if cfg.visual_check.enabled:
                        variable_visualization(var=var[par, p, :, :],
                                               x=np.asarray(ncfile.variables['x']), y=np.asarray(ncfile.variables['y']),
                                               var_name=vn, par_id=par, text_id='par: {}, k: {}, p: {}'.format(par,k, p),
                                               path=cfg.visual_check.path, show_plots = cfg.visual_check.show_plots)
            debug('Variable {}, parameter {} written', vn, par)

def write_water(ncfile, cfg, connection, cur):
    # write water surfaces
    vn = "water_type"
    vt = 'b'
    if cfg.landcover.surface_fractions:
        sqltext = f"""
            select case when g.wat_fraction > {cfg.landcover.min_fraction} and 
                             g.wat_fract_type is not null then g.wat_fract_type - %s 
                   else null end
            from "{cfg.domain.case_schema}"."{cfg.tables.grid}" g
            order by g.j, g.i
        """
        cur.execute(sqltext, (cfg.type_range.water_min, ))
        sql_debug(connection)
    else:
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

    if cfg.lod2:
        # Process water_pars
        for par in cfg.water_pars._settings.keys():
            # par - index of npars dimmension in netcdf
            # pavement_pars[par] - calculation formula for parameter
            sqltext = 'SELECT case when l.type is not null then {0} else null end ' \
                       'from "{1}"."{2}" g    ' \
                      'left outer join "{1}"."{3}" l on l.lid = g.lid and l.type >= %s and l.type < %s ' \
                      'left outer join "{1}"."{4}" p on p.code = l."{5}" ' \
                      'order by g.j, g.i' \
                .format(cfg.water_pars[par], cfg.domain.case_schema, cfg.tables.grid,
                        cfg.tables.landcover, cfg.tables.surface_params, cfg.landcover_params_var)
            cur.execute(sqltext, (cfg.type_range.water_min, cfg.type_range.water_max,))
            res = cur.fetchall()
            sql_debug(connection)
            res = [cfg.fill_values[vt] if x[0] is None else x[0] for x in res]
            var = np.reshape(np.asarray(res, dtype=vt), (cfg.domain.ny, cfg.domain.nx))
            del res
            var = np.nan_to_num(var, copy=False, nan=cfg.fill_values[vt],
                                posinf=cfg.fill_values[vt], neginf=cfg.fill_values[vt])
            ncfile.variables[vn][par,...] = var
            debug('Variable {}, parameter {} written.', vn, par)
            if cfg.visual_check.enabled:
                variable_visualization(var=ncfile.variables[vn][par,...],
                                       x=np.asarray(ncfile.variables['x']), y=np.asarray(ncfile.variables['y']),
                                       var_name=vn, par_id=par, text_id='water_pars', path=cfg.visual_check.path,
                                       show_plots = cfg.visual_check.show_plots)



def write_vegetation(ncfile, cfg, connection, cur):
    # write vegetation surfaces
    vn = "vegetation_type"
    vt = 'b'
    if cfg.landcover.surface_fractions:
        sqltext = f"""
            select case when g.veg_fraction > {cfg.landcover.min_fraction} and 
                             g.veg_fract_type is not null then g.veg_fract_type - %s 
                   else null end
            from "{cfg.domain.case_schema}"."{cfg.tables.grid}" g
            order by g.j, g.i
        """
        cur.execute(sqltext, (cfg.type_range.vegetation_min, ))
        sql_debug(connection)
    else:
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

    if cfg.lod2:
        # Process vegetation_pars
        vn = "vegetation_pars"
        vt = 'f4'
        ncfile.createVariable(vn, vt, ('nvegetation_pars', 'y', 'x'), fill_value=cfg.fill_values[vt])
        nc_write_attribute(ncfile, vn, 'long_name', 'vegetation parameters')
        nc_write_attribute(ncfile, vn, 'units', '1')
        nc_write_attribute(ncfile, vn, 'source', '')
        nc_write_attribute(ncfile, vn, 'res_orig', cfg.domain.dz)
        nc_write_attribute(ncfile, vn, 'coordinates', 'E_UTM N_UTM lon lat')
        nc_write_attribute(ncfile, vn, 'grid_mapping', 'E_UTM N_UTM lon lat')
        # fill individual parameters in dimension npavement_pars

        for par in cfg.vegetation_pars._settings.keys():
            # par - index of npars dimmension in netcdf
            # pavement_pars[par] - calculation formula for parameter
            sqltext = 'SELECT case when l.type is not null then {0} else null end from "{1}"."{2}" g    ' \
                      'left outer join "{1}"."{3}" l on l.lid = g.lid and l.type >= %s and l.type < %s ' \
                      'left outer join "{1}"."{4}" p on p.code = l."{5}" ' \
                      'order by g.j, g.i' \
                .format(cfg.vegetation_pars[par], cfg.domain.case_schema, cfg.tables.grid,
                        cfg.tables.landcover, cfg.tables.surface_params,
                        cfg.landcover_params_var)
            cur.execute(sqltext, (cfg.type_range.vegetation_min, cfg.type_range.vegetation_max,))
            res = cur.fetchall()
            sql_debug(connection)
            res = [cfg.fill_values[vt] if x[0] is None else x[0] for x in res]
            var = np.reshape(np.asarray(res, dtype=vt), (cfg.domain.ny, cfg.domain.nx))
            del res
            var = np.nan_to_num(var, copy=False, nan=cfg.fill_values[vt],
                                posinf=cfg.fill_values[vt], neginf=cfg.fill_values[vt])
            ncfile.variables[vn][par, ...] = var
            debug('Variable {}, parameter {} written', vn, par)
            if cfg.visual_check.enabled:
                variable_visualization(var=ncfile.variables[vn][par, ...],
                                       x=np.asarray(ncfile.variables['x']), y=np.asarray(ncfile.variables['y']),
                                       var_name=vn, par_id=par, text_id='vegetation_pars', path=cfg.visual_check.path,
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

    # write soil moisture adjustments: It is not in PIDS, it is our own extension used by our user code
    # if cfg.lod2 and len(cfg.soil_moisture_adjust._settings) > 0:
    #     vn = 'soil_moisture_adjust'
    #     vt = 'f4'
    #     # first create calculation formula
    #     sm_text = 'case '
    #     for lc in cfg.soil_moisture_adjust._settings:
    #         sm_text = sm_text + 'when l."{0}" = {1} then {2} '.format( \
    #                   cfg.landcover_params_var, lc, cfg.soil_moisture_adjust._settings[lc])
    #     sm_text = sm_text + 'else 1 end '
    #     sqltext = 'select case when l.lid is not null then {0} else null end ' \
    #               'from "{1}"."{2}" g ' \
    #               'left outer join "{1}"."{3}" l on l.lid = g.lid and l.type >= %s and l.type < %s order by g.j, g.i' \
    #         .format(sm_text, cfg.domain.case_schema, cfg.tables.grid, cfg.tables.landcover)
    #     cur.execute(sqltext, (cfg.type_range.vegetation_min, cfg.type_range.vegetation_max,))
    #     sql_debug(connection)
    #     write_type_variable(ncfile, cfg, vn, 'Soil moisture adjust', vt, cfg.fill_values, 0, connection, cur)
    #     if cfg.visual_check.enabled:
    #         variable_visualization(var=ncfile[vn],
    #                                x=np.asarray(ncfile.variables['x']), y=np.asarray(ncfile.variables['y']),
    #                                var_name=vn, par_id='', text_id='soil_moisture_adjust', path=cfg.visual_check.path,
    #                                show_plots = cfg.visual_check.show_plots)

def write_mask_usm(ncfile, cfg, connection, cur):
    """ In case of SLURB, mask buildings with vegetation """
    progress('Masking buildings in case of SLURB')
    debug('Create mask')
    sqltext = 'select case when b.lid is not null then true else false end ' \
              'from "{0}"."{1}" g ' \
              'left outer join "{0}"."{2}" b on b.id = g.id order by g.j, g.i' \
        .format(cfg.domain.case_schema, cfg.tables.grid, cfg.tables.buildings_grid)
    cur.execute(sqltext)
    res = cur.fetchall()
    sql_debug(connection)
    # res = [False if x[0] is None else True for x in res]
    res = [x[0] for x in res]
    slurb_mask = np.reshape(np.asarray([x for x in res], dtype='bool_'), (cfg.domain.ny, cfg.domain.nx))
    veg_type = ncfile.variables['vegetation_type'][:,:]
    veg_type[slurb_mask] = 1
    ncfile.variables['vegetation_type'][:, :] = veg_type

    if cfg.landcover.surface_fractions:
        debug('Replacing surface fractions')
        veg_frac = ncfile.variables['surface_fraction'][:, :, :]
        veg_frac[0, slurb_mask] = 1.0
        veg_frac[1, slurb_mask] = 0.0 #cfg.fill_values['f4']
        veg_frac[2, slurb_mask] = 0.0 #cfg.fill_values['f4']
        ncfile.variables['surface_fraction'][:, :, ] = veg_frac

def write_buildings(ncfile, cfg, connection, cur):
    """ write building_height (buildings_2d), building_id and building_type into netcdf file
    """
    if cfg.force_lsm_only:
        return
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

    if cfg.has_3d_buildings:
        debug('Write buildings_3d into static driver')
        verbose('Prepare netcdf variable')
        vt = 'b'
        vn = 'buildings_3d'
        sqltext = 'SELECT MAX(nz) FROM "{}"."{}"'.format(cfg.domain.case_schema, cfg.tables.buildings_grid)
        cur.execute(sqltext)
        sql_debug(connection)
        max_nz = int(cur.fetchone()[0]) + 1
        # create z dimension
        nc_create_dimension(ncfile, 'z', max_nz)
        ncfile.createVariable('z', 'f4', ('z'))
        ncfile['z'][:] = np.append(0, np.arange(cfg.domain.dz / 2.0, (max_nz - 1) * cfg.domain.dz, cfg.domain.dz))
        nc_write_attribute(ncfile, 'z', 'long_name', 'z')
        nc_write_attribute(ncfile, 'z', 'standard_name', 'projection_z_coordinate')
        nc_write_attribute(ncfile, 'z', 'units', 'm')
        ncfile.createVariable(vn, vt, ('z', 'y', 'x'), fill_value=cfg.fill_values[vt])
        nc_write_attribute(ncfile, vn, 'coordinates', 'E_UTM N_UTM lon lat')
        nc_write_attribute(ncfile, vn, 'flag_meanings', 'no building, building')
        nc_write_attribute(ncfile, vn, 'long_name', 'buildings_3d')
        nc_write_attribute(ncfile, vn, 'res_orig', cfg.domain.dz)
        nc_write_attribute(ncfile, vn, 'source', '')
        nc_write_attribute(ncfile, vn, 'units', '1')
        nc_write_attribute(ncfile, vn, 'lod', 2)
        var_3d = np.zeros((max_nz, cfg.domain.ny, cfg.domain.nx), dtype='int')

        verbose('Pull info from server, max height')
        sqltext = 'SELECT b.nz FROM "{0}"."{1}" AS g ' \
                  'LEFT OUTER JOIN "{0}"."{2}" AS b ON b.id = g.id ' \
                  'ORDER BY g.j, g.i' \
                  .format(cfg.domain.case_schema, cfg.tables.grid, cfg.tables.buildings_grid)
        cur.execute(sqltext)
        sql_debug(connection)
        res = cur.fetchall()
        sql_debug(connection)
        res = [0 if x[0] is None else x[0] for x in res]
        var_top = np.reshape(np.asarray([x for x in res], dtype=vt), (cfg.domain.ny, cfg.domain.nx))
        var_top = np.nan_to_num(var_top, copy=False, nan=cfg.fill_values[vt],
                                posinf=cfg.fill_values[vt], neginf=cfg.fill_values[vt])
        del res

        verbose('Pull info from server, bottom height')
        sqltext = 'SELECT CASE WHEN b.is_bridge OR b.has_bottom THEN b.nz_min ELSE 0 END FROM "{0}"."{1}" AS g ' \
                  'LEFT OUTER JOIN "{0}"."{2}" AS b ON b.id = g.id ' \
                  'ORDER BY g.j, g.i' \
            .format(cfg.domain.case_schema, cfg.tables.grid, cfg.tables.buildings_grid)
        cur.execute(sqltext)
        sql_debug(connection)
        res = cur.fetchall()
        sql_debug(connection)
        res = [0 if x[0] is None else x[0] for x in res]
        var_bottom = np.reshape(np.asarray([x for x in res], dtype=vt), (cfg.domain.ny, cfg.domain.nx))
        var_bottom = np.nan_to_num(var_bottom, copy=False, nan=cfg.fill_values[vt],
                                   posinf=cfg.fill_values[vt], neginf=cfg.fill_values[vt])
        del res

        verbose('Prepare 3d array with boolean information about topo')
        for j in range(cfg.domain.ny):
            for i in range(cfg.domain.nx):
                var_3d[var_bottom[j, i] + 1:var_top[j, i] + 1, j, i] = 1
        var_3d[0, :, :] = var_3d[1, :, :]

        verbose('A check for buildings bottoms')
        sqltext = 'SELECT CASE WHEN (b.is_bridge OR b.has_bottom) AND b.nz = 0 THEN True ELSE False END FROM "{0}"."{1}" g ' \
                  'LEFT OUTER JOIN "{0}"."{2}" b on b.id = g.id ' \
                  'ORDER BY g.j, g.i' \
            .format(cfg.domain.case_schema, cfg.tables.grid, cfg.tables.buildings_grid)
        cur.execute(sqltext)
        sql_debug(connection)
        res = cur.fetchall()
        sql_debug(connection)
        res = [cfg.fill_values[vt] if x[0] is None else x[0] for x in res]
        var_bottom = np.reshape(np.asarray([x for x in res], dtype='bool_'), (cfg.domain.ny, cfg.domain.nx))
        del res

        for j in range(cfg.domain.ny):
            for i in range(cfg.domain.nx):
                if var_bottom[j, i]:
                    var_3d[0, j, i] = 1

        ncfile.variables[vn][:] = var_3d

        debug('Variable buildings 3d has been written')

        debug('Process types under buildings 3d')
        sqltext = 'SELECT ub.j, ub.i, ux.typed FROM "{0}"."{1}" AS ub ' \
                  'LEFT OUTER JOIN "{0}"."{2}" AS ux ON ux.gid = ub.lid_extra ' \
                  'WHERE ub.under ' \
                  'ORDER BY ub.j, ub.i ' \
            .format(cfg.domain.case_schema, cfg.tables.buildings_grid, cfg.tables.extras_shp)
        cur.execute(sqltext)
        sql_debug(connection)
        res = cur.fetchall()
        for j, i, typeu in res:
            extra_verbose('Updating [j,i] = [{},{}]', j, i)
            if cfg.type_range.pavement_min < typeu < cfg.type_range.pavement_max:
                # pavement type
                vn = 'pavement_type'
                extra_verbose('\tIt is pavement type, add [i,j] into pavement_type')
                ncfile.variables[vn][j, i] = typeu - cfg.type_range.pavement_min
                vn = "soil_type"
                ncfile.variables[vn][j, i] = cfg.ground.soil_type_default

            elif cfg.type_range.vegetation_min < typeu < cfg.type_range.vegetation_max:
                # vegetation type
                vn = 'vegetation_type'
                extra_verbose('\tIt is vegetation type, add [i,j] into vegetation_type')
                ncfile.variables[vn][j, i] = typeu - cfg.type_range.vegetation_min
                vn = "soil_type"
                ncfile.variables[vn][j, i] = cfg.ground.soil_type_default

            elif cfg.type_range.water_min < typeu < cfg.type_range.water_max:
                # water type
                vn = 'water_type'
                extra_verbose('\tIt is water type, add [i,j] into water type')
                ncfile.variables[vn][j, i] = typeu - cfg.type_range.water_min

            else:
                warning('During filling types under 3d structure. In [i,j] = [{},{}] is unknown type = {}',
                        i, j, typeu)

        # FIXME: NOT IMPLEMENTED IN PALM
        # # pavement pars
        # if cfg.landcover_params_var:
        #     # Process pavement_pars
        #     vn = "pavement_pars"
        #     vt = 'f4'
        #     for par in cfg.pavement_pars._settings.keys():
        #         sqltext = 'SELECT br.j, br.i, {0} ' \
        #                   'from "{1}"."{2}" br    ' \
        #                   'left outer join "{1}"."{3}" l on l.lid = br.lid and l.typed >= %s and l.typed < %s ' \
        #                   'left outer join "{1}"."{4}" p on p.code = l."{5}" ' \
        #                   'WHERE br.under ' \
        #                   'order by br.j, br.i ' \
        #             .format(cfg.pavement_pars[par], cfg.domain.case_schema, cfg.tables.bridge_grid,
        #                     cfg.tables.bridge_shp,   cfg.tables.surface_params,
        #                     cfg.landcover_params_var+'d')
        #         cur.execute(sqltext, (cfg.type_range.pavement_min, cfg.type_range.pavement_max,))
        #         res = cur.fetchall()
        #         sql_debug(connection)
        #         for j, i, ty in res:
        #             extra_verbose('Updating {} in [j,i] = [{},{}], instead of {} write {}',
        #                           vn, j, i, ncfile.variables[vn][par, j, i], ty)
        #             ncfile.variables[vn][par, j, i] = ty
        #     connection.commit()
        #     debug_sql(connection)

def prepared_lad_netcdf(ncfile, cfg, nzlad):
    """ Prepare dimension and variables """
    vt = 'f4'
    debug('create zlad and zbad dimensions')
    zlad = [0] + [x * cfg.domain.dz + 0.5 * cfg.domain.dz for x in range(nzlad)]
    nc_create_dimension(ncfile, 'zlad', nzlad+1)
    temp = ncfile.createVariable('zlad', vt, 'zlad')
    temp[:] = zlad[:]
    # create and write lad and bad variables
    vn = 'lad'
    var = ncfile.createVariable(vn, vt, ('zlad', 'y', 'x'), fill_value=cfg.fill_values[vt])
    nc_write_attribute(ncfile, vn, 'long_name', 'leaf area density')
    nc_write_attribute(ncfile, vn, 'units', 'm2/m3')
    nc_write_attribute(ncfile, vn, 'res_orig', cfg.domain.dz)
    nc_write_attribute(ncfile, vn, 'coordinates', 'E_UTM N_UTM lon lat')
    nc_write_attribute(ncfile, vn, 'grid_mapping', 'E_UTM N_UTM lon lat')

    vn = 'bad'
    ncfile.createVariable(vn, vt, ('zlad', 'y', 'x'), fill_value=cfg.fill_values[vt])
    nc_write_attribute(ncfile, vn, 'long_name', 'branch area density')
    nc_write_attribute(ncfile, vn, 'units', 'm3/m3')
    nc_write_attribute(ncfile, vn, 'res_orig', cfg.domain.dz)
    nc_write_attribute(ncfile, vn, 'coordinates', 'E_UTM N_UTM lon lat')
    nc_write_attribute(ncfile, vn, 'grid_mapping', 'E_UTM N_UTM lon lat')

def write_trees_grid(ncfile, cfg, connection, cur):
    """ Routine to generate trees """
    change_log_level(cfg.logs.level_trees)
    debug('get max height of the tree in the domain (relative height from ground)')
    sqltext = 'select max(vysstr) from "{0}"."{1}"'.format(cfg.domain.case_schema, cfg.tables.trees)
    cur.execute(sqltext)
    ret = cur.fetchone()
    nzlad = ceil(ret[0] / cfg.domain.dz) + 1
    sql_debug(connection)
    # construct lad array
    debug('nzlad = {}', nzlad)
    vt = 'f4'
    lad = np.zeros((nzlad, cfg.domain.ny, cfg.domain.nx), dtype=vt)
    bad = np.zeros((nzlad, cfg.domain.ny, cfg.domain.nx), dtype=vt)
    sqltext = 'select i, j'
    for l in range(nzlad):
        sqltext += ', lad_{0}, bad_{0}'.format(l)
    sqltext += ' from "{0}"."{1}"'.format(cfg.domain.case_schema, cfg.tables.trees_grid)
    cur.execute(sqltext)
    trees = cur.fetchall()
    # construct lad array
    for tree in trees:
        # call palm_tree_lad sql function for itree
        i, j = tree[0], tree[1]
        extra_verbose('tree: [j, i], [{},{}]', i, j)
        for l in range(nzlad):
            lad[l, j, i] += tree[2 * l + 2]
            bad[l, j, i] += tree[2 * l + 3]

    prepared_lad_netcdf(ncfile, cfg, nzlad)

    ncfile.variables['lad'][1:, :, :] = lad
    ncfile.variables['bad'][1:, :, :] = bad
    if cfg.visual_check.enabled:
        for k in range(nzlad):
            variable_visualization(var=ncfile['lad'][k, ...],
                                   x=np.asarray(ncfile.variables['x']), y=np.asarray(ncfile.variables['y']),
                                   var_name='lad', par_id=k, text_id='leaf_area_density', path=cfg.visual_check.path,
                                   show_plots=cfg.visual_check.show_plots)
            variable_visualization(var=ncfile['bad'][k, ...],
                                   x=np.asarray(ncfile.variables['x']), y=np.asarray(ncfile.variables['y']),
                                   var_name='bad', par_id=k, text_id='trunk_area_density', path=cfg.visual_check.path,
                                   show_plots=cfg.visual_check.show_plots)

def write_lad_grid(ncfile, cfg, connection, cur):
    """ Write into netcdf LAD and BAD field """
    sqltext = """
        select max(canopy_height) 
        from "{0}"."{1}" 
    """.format(cfg.domain.case_schema, cfg.tables.grid)
    cur.execute(sqltext)
    ret = cur.fetchone()
    nzlad = ceil(ret[0] / cfg.domain.dz) + 1
    sql_debug(connection)

    debug('nzlad = {}', nzlad)
    vt = 'f4'
    lad = np.zeros((nzlad, cfg.domain.ny, cfg.domain.nx), dtype=vt)
    bad = np.zeros((nzlad, cfg.domain.ny, cfg.domain.nx), dtype=vt)

    prepared_lad_netcdf(ncfile, cfg, nzlad)

    ncfile.variables['lad'][0, :, :] = 0.0
    ncfile.variables['lad'][0, :, :] = 0.0

    for nz in range(1, nzlad+1):
        verbose('nz: {}', nz)
        sqltext = """ 
        select 
            case when {5} between 0.0 and g.canopy_height and g.canopy_height > 0
                    then least(1.6, lai / canopy_height) * greatest(0.0, least((canopy_height - {2}) / {4}, 1.0))
                else 0.0 end
        from "{0}"."{1}" g
        order by g.j, g.i
        """.format(cfg.domain.case_schema, cfg.tables.grid,
                   (nz - 1) * cfg.domain.dz, (nz) * cfg.domain.dz, cfg.domain.dz,
                   (nz - 0.5) * cfg.domain.dz)
        cur.execute(sqltext)
        res = cur.fetchall()
        sql_debug(connection)
        lad = np.reshape(np.asarray([x[0] for x in res], dtype=vt), (cfg.domain.ny, cfg.domain.nx))
        ncfile.variables['lad'][nz, :, :] = lad
        bad = 0.1 * lad
        ncfile.variables['bad'][nz, :, :] = bad

def write_building_pars(ncfile, cfg, connection, cur):
    # write building_pars - roof parameters
    vn = 'building_pars'
    vt = 'f4'
    ncfile.createVariable(vn, vt, ('nbuilding_pars', 'y', 'x'), fill_value=cfg.fill_values[vt])
    nc_write_attribute(ncfile, vn, 'long_name', 'building parameters')
    nc_write_attribute(ncfile, vn, 'units', '')
    nc_write_attribute(ncfile, vn, 'res_orig', cfg.domain.dz)
    nc_write_attribute(ncfile, vn, 'coordinates', 'E_UTM N_UTM lon lat')
    nc_write_attribute(ncfile, vn, 'grid_mapping', 'E_UTM N_UTM lon lat')
    # fill individual parameters in dimension nbuilding_pars
    for par in cfg.building_pars._settings.keys():
        extra_verbose('Processing buildings pars, par: {}', par)
        # par - index of npars dimmension in netcdf
        # building_pars[par] - calculation formula for parameter
        # first create calculation formula
        if par in cfg.building_pars_repl._settings.keys():
            # parameter replacements for green frac
            repl = cfg.building_pars_repl[par]
            g_text = 'case '
            for gf in repl:
                g_text = g_text + 'when p.code = {0} then {1} '.format(gf[0], gf[1])
            g_text = g_text + 'else {0} end '.format(cfg.building_pars[par])
        else:
            g_text = '{0} '.format(cfg.building_pars[par])
        extra_verbose('\t{}', g_text)
        # If insulation enabled
        if cfg.insulation.enabled and par in cfg.insulation.building_pars and \
                cfg.insulation.exists[cfg.insulation.pars_fields[cfg.insulation.building_pars.index(par)]]:
            # parameter replacements for insulation layers of walls
            i = cfg.insulation.building_pars.index(par)
            i_text = 'CASE WHEN "{0}" <> 0 THEN {1} ELSE {2} END '.format(
                cfg.insulation.fields[cfg.insulation.pars_fields[i]],
                cfg.insulation.values[cfg.insulation.pars_items[i]], g_text)
        else:
            i_text = '{0} '.format(g_text)
        extra_verbose('\t{}', i_text)

        # calculation formula
        # TODO: Move into one select, wall and roof are separated in parameters
        """
        select distinct on (g.i, g.j) g.j, g.i, b.id as b_id, bw.id as bw_id, 
            g.id as g_id, w.gid as w_gid, r.gid as r_gid, pr.code as pr_code, pg.code as pg_code, pu.code as pu_code,
            case when (bw.id is not null OR b.id IS NOT NULL) then 
            --w.winfrach 
            pr.capacity_surf
            else null end  as pr_capacity_surf,
            case when (bw.id is not null OR b.id IS NOT NULL) then 
            --w.winfrach 
            pu.capacity_surf
            else null end as pu_capacity_surf,
            case when (bw.id is not null OR b.id IS NOT NULL) then 
            --w.winfrach 
            w.emisivitah
            else null end as w_emisivitah
        from "evropska_martin_03_validation_summer"."grid" g 
        left outer join "evropska_martin_03_validation_summer"."building_walls" bw on bw.id = g.id and not bw.isroof 
        left outer join "evropska_martin_03_validation_summer"."buildings" b on b.id = g.id
        left outer join "evropska_martin_03_validation_summer"."walls" w on w.gid = bw.wid 
        left outer join "evropska_martin_03_validation_summer"."roofs" r on r.gid = b.rid 
        left outer join "evropska_martin_03_validation_summer"."surface_params" pr on pr.code = r.material+200
        left outer join "evropska_martin_03_validation_summer"."surface_params" pg on pg.code = w.stenakatd 
        left outer join "evropska_martin_03_validation_summer"."surface_params" pu on pu.code = w.stenakath 
        order by g.j,g.i
        limit 100
        """
        p_text = 'case when b.id is not null then {0} else null end'.format(i_text)
        sqltext = 'SELECT DISTINCT ON (g.i, g.j) ' \
                  '{0} from "{1}"."{2}" g ' \
                  'left outer join "{1}"."{3}" b on b.id = g.id ' \
                  'left outer join "{1}"."{4}" r on r.gid = b.rid ' \
                  'left outer join "{1}"."{7}" p on p.code = cast(r.material as integer)+%s ' \
                  'left outer join "{1}"."{5}" bw on bw.id = g.id and not bw.isroof ' \
                  'left outer join "{1}"."{6}" w on w.gid = bw.wid ' \
                  'left outer join "{1}"."{7}" pg on pg.code = w.stenakatd ' \
                  'left outer join "{1}"."{7}" pu on pu.code = w.stenakath ' \
                  'order by g.j,g.i ' \
            .format(p_text, cfg.domain.case_schema, cfg.tables.grid,
                    cfg.tables.buildings_grid, cfg.tables.roofs, cfg.tables.building_walls, cfg.tables.walls, cfg.tables.surface_params, )
        extra_verbose('\t{}', sqltext)
        cur.execute(sqltext, (cfg.surf_range.roof_min,))
        res = cur.fetchall()
        sql_debug(connection)
        varr = np.reshape(np.asarray([x[0] for x in res], dtype=vt), (cfg.domain.ny, cfg.domain.nx))
        del res
        # if par in cfg.building_pars_wall._settings.keys():
        #     # process also parameter of walls
        #     # create calculation formula
        #     if cfg.insulation.enabled and par in cfg.insulation.building_pars and \
        #        cfg.insulation.exists[cfg.insulation.pars_fields[cfg.insulation.building_pars.index(par)]]:
        #         # parameter replacements for insulation layers of walls
        #         i = cfg.insulation.building_pars.index(par)
        #         i_text = 'case when "{0}" <> 0 then {1} else {2} end '.format(
        #                  cfg.insulation.fields[cfg.insulation.pars_fields[i]],
        #                  cfg.insulation.values[cfg.insulation.pars_items[i]], cfg.building_pars_wall[par])
        #     else:
        #         i_text = '{0} '.format(cfg.building_pars_wall[par])
        #     # if par in cfg.building_pars_repl._settings.keys():
        #     #     # parameter replacements for green frac
        #     #     repl = cfg.building_pars_repl[par]
        #     #     g_text = 'case '
        #     #     for gf in repl:
        #     #         g_text = g_text + 'when p.code = {0} then {1} '.format(gf[0], gf[1])
        #     #     i_text = g_text + 'else {0} end '.format(i_text)
        #
        #     p_text = 'case when bw.id is not null then {0} else null end '.format(i_text)
        #     sqltext = 'select distinct on (g.i, g.j) ' \
        #               '{0} from "{1}"."{2}" g ' \
        #               'left outer join "{1}"."{3}" bw on bw.id = g.id and not bw.isroof ' \
        #               'left outer join "{1}"."{4}" w on w.gid = bw.wid ' \
        #               'left outer join "{1}"."{5}" pg on pg.code = w.stenakatd ' \
        #               'left outer join "{1}"."{5}" pu on pu.code = w.stenakath ' \
        #               'order by g.j,g.i ' \
        #         .format(p_text, cfg.domain.case_schema, cfg.tables.grid, cfg.tables.building_walls,
        #                 cfg.tables.walls, cfg.tables.surface_params, )
        #     cur.execute(sqltext)
        #     res = cur.fetchall()
        #     sql_debug(connection)
        #     varw = np.reshape(np.asarray([x[0] for x in res], dtype=vt), (cfg.domain.ny, cfg.domain.nx))
        #     del res
        #     # combine the arrays
        #     varr = np.where(np.isnan(varw), varr, varw)
        # replace nan with fillValue and save in netcdf file
        var = np.nan_to_num(varr, copy=False, nan=cfg.fill_values[vt], posinf=cfg.fill_values[vt], neginf=cfg.fill_values[vt])
        ncfile.variables[vn][par,...] = var
        debug('Variable {}, parameter {} written.', vn, par)
        if cfg.visual_check.enabled:
            variable_visualization(var=ncfile.variables[vn][par,...],
                                   x=np.asarray(ncfile.variables['x']), y=np.asarray(ncfile.variables['y']),
                                   var_name=vn, par_id=par, text_id='building_pars', path=cfg.visual_check.path,
                                   show_plots = cfg.visual_check.show_plots)

    # check winfrac < 0.95 Honza's idea
    winfrac = np.argwhere(ncfile.variables[vn][1,...] > 0.95)
    if winfrac.size > 0:
        verbose('Modifying winfrac in buildings_pars where winfrac > 0.95, '
                'setting winfrac = 0.95 and wallfrac = 0.05')
        for j, i in winfrac:
            ncfile.variables[vn][1, j, i] = 0.95
            ncfile.variables[vn][0, j, i] = 0.05

    winfrac = np.argwhere(ncfile.variables[vn][22, ...] > 0.95)
    if winfrac.size > 0:
        verbose('Modifying winfrac in buildings_pars where winfrac > 0.95, '
                'setting winfrac = 0.95 and wallfrac = 0.05')
        for j, i in winfrac:
            ncfile.variables[vn][22, j, i] = 0.95
            ncfile.variables[vn][21, j, i] = 0.05
    del winfrac


    debug('Variable {} completely written.', vn)

def test_building_insulation(cfg, connection, cur):
    # test if columns zatepd and zateph exist in table walls and add the indication to insulation cfg parameters
    cfg.insulation._settings['exists'] = []
    debug('testing buildings insulation')
    for f in cfg.insulation.fields:
        sqltext = 'select exists (select * from information_schema.columns ' \
                  'where table_schema = %s and table_name=%s and column_name = %s)'
        cur.execute(sqltext, (cfg.domain.case_schema, cfg.tables.walls, f,))
        cfg.insulation.exists.append(cur.fetchone()[0])
    sql_debug(connection)
    connection.commit()

def write_building_surface_pars(ncfile, cfg, connection, cur, vtabs):
    # create dimension s(ns) for indexing of surfaces
    # acquire max surf_id
    sqltext = 'SELECT sid FROM "{0}"."{1}" AS s ' \
              'LEFT OUTER JOIN "{0}"."{2}" AS g ON g.id = s.gid ' \
              'ORDER BY g.j ASC, g.i ASC, s.direction ASC, s.zs DESC'.format(cfg.domain.case_schema, cfg.tables.surfaces, cfg.tables.grid)
    cur.execute(sqltext)
    var = cur.fetchall()
    sql_debug(connection)
    ns = len(var)
    nc_create_dimension(ncfile, 's', ns)
    temp = nc_create_variable(ncfile, 's', 'i', 's')
    temp[:] = np.asarray([x[0] for x in var], dtype='i')
    # create all related variables
    sqltext = 'SELECT xs, ys, zs, azimuth, zenith, lons, lats, "Es_UTM", "Ns_UTM" FROM "{0}"."{1}" AS s ' \
              'LEFT OUTER JOIN "{0}"."{2}" AS g ON g.id = s.gid ' \
              'ORDER BY g.j ASC, g.i ASC, s.direction ASC, s.zs DESC'.format(cfg.domain.case_schema, cfg.tables.surfaces, cfg.tables.grid)
    cur.execute(sqltext)
    var = cur.fetchall()
    sql_debug(connection)
    vi = ['xs', 'ys', 'zs', 'azimuth', 'zenith', 'lons', 'lats']
    vt = 'f4'
    for i in range(len(vi)):
        nc_create_dimension(ncfile, vi[i], ns)
        temp = nc_create_variable(ncfile, vi[i], vt, 's')
        res = [cfg.fill_values[vt] if x[i] is None else x[i] for x in var]
        temp[:] = np.asarray(res, dtype=vt)
    vj = ['Es_UTM', 'Ns_UTM']
    vt = 'f8'
    i = 0
    for j in range(len(vi), len(vi)+len(vj)):
        nc_create_dimension(ncfile, vj[i], ns)
        temp = nc_create_variable(ncfile, vj[i], vt, 's')
        res = [cfg.fill_values[vt] if x[j] is None else x[j] for x in var]
        temp[:] = np.asarray(res, dtype=vt)
        i += 1

    # write roof and wall surface pars
    vn = 'building_surface_pars'
    vt = 'f4'
    var = nc_create_variable(ncfile, vn, vt, ('nbuilding_surface_pars', 's'), fill_value=cfg.fill_values[vt])
    nc_write_attribute(ncfile, vn, 'long_name', 'building parameters')
    nc_write_attribute(ncfile, vn, 'units', '')
    nc_write_attribute(ncfile, vn, 'res_orig', cfg.domain.dz)
    nc_write_attribute(ncfile, vn, 'coordinates', 'Es_UTM Ns_UTM lons lats')
    # fill individual parameters in dimension nbuilding_pars
    for par in cfg.building_surface_pars._settings.keys():
        # par - index of npars dimension in netcdf
        # building_pars[par] - calculation formula for parameter

        # create calculation formula
        if cfg.tables.extras_shp in vtabs:
            p = cfg.building_surface_pars[par]
            if not isinstance(p, list):
                p = [p, p, p, p, p]
            pt = ['pr', 'pg', 'pu', 'pd', 'b'] # pr ... roof, pg ... ground, pu ... upper floors, pd ... downward facing, b .. bridges
            sqlline1 = 'LEFT OUTER JOIN "{0}"."{1}" d on s.eid = d.gid '.format(cfg.domain.case_schema, cfg.tables.extras_shp)
            sqlline2 = 'LEFT OUTER JOIN "{0}"."{1}" pd on pd.code = d.katlandd '.format(cfg.domain.case_schema, cfg.tables.surface_params)
            sqlline3 = 'LEFT OUTER JOIN "{0}"."{1}" AS be ON s.eid = be.gid '.format(cfg.domain.case_schema, cfg.tables.extras_shp)
            sqlline4 = 'LEFT OUTER JOIN "{0}"."{1}" AS b ON b.code = be.katlandu '.format(cfg.domain.case_schema, cfg.tables.surface_params)
        else:
            p = cfg.building_surface_pars[par]
            if not isinstance(p, list):
                p = [p, p, p]
            p = p[:3]
            pt = ['pr', 'pg', 'pu']
            sqlline1 = ''
            sqlline2 = ''
            sqlline3 = ''
            sqlline4 = ''
        g_text = {}

        for pk in range(len(p)):  # over roof, ground, upper levels
            if par in cfg.building_surface_pars_repl._settings.keys():
                # parameter replacements for green frac
                repl = cfg.building_surface_pars_repl[par]
                g_text[pk] = 'case '
                for gf in repl: # replace keys
                    g_text[pk] = g_text[pk] + 'when {0}.code = {1} then {2} '.format(pt[pk], gf[0], gf[1])
                g_text[pk] = g_text[pk] + 'else {0} end '.format(p[pk])
            else:
                g_text[pk] = '{0} '.format(p[pk])
        # test if additional insulation applies
        if cfg.insulation.enabled and par in cfg.insulation.building_surface_pars:
            # parameter replacements for insulation layers of walls
            i = cfg.insulation.building_surface_pars.index(par)
            for pk in range(1,len(p)): # insulation for walls, not roofs
                if pk >= 3:
                    continue # skip downward facing walls
                if cfg.insulation.exists[pk-1]:  # insulation.exists does not contain roofs
                    g_text[pk] = 'case when "{0}" <> 0 then {1} else {2} end '.format(cfg.insulation.fields[pk-1],
                                  cfg.insulation.values[cfg.insulation.surface_pars_items[i]], g_text[pk])
        # set final calculation formula
        if cfg.tables.extras_shp in vtabs:
            p_text = ' case when s.isroof then CASE WHEN s.eid IS NOT NULL THEN {6} ELSE {0} END when not s.isroof and s.ishorizontal then {4} ' \
                     '  else case when s.zs<={1}+wart.nz_min_art*{5} then {2} else {3} end end ' \
                     .format(g_text[0], cfg.ground.ground_floor_height, g_text[1], g_text[2], g_text[3], cfg.domain.dz, g_text[4])
        else:
            p_text = ' case when s.isroof then {0} ' \
                     '  else case when s.zs<={1}+wart.nz_min_art*{4} then {2} else {3} end end ' \
                .format(g_text[0], cfg.ground.ground_floor_height, g_text[1], g_text[2], cfg.domain.dz)
        # build sql select text # FIXME change r.rid (r.gid)
        sqltext = 'select {0} from "{1}"."{2}" s ' \
                  'left outer join "{1}"."{3}" r on r.gid = s.rid ' \
                  'left outer join "{1}"."{4}" pr on pr.code = cast(r.material as integer) + {7} ' \
                  'left outer join "{1}"."{5}" w on w.gid = s.wid  ' \
                  '{8} ' \
                  '{9} ' \
                  'left outer join "{1}"."{4}" pg on pg.code = w.stenakatd ' \
                  'left outer join "{1}"."{4}" pu on pu.code = w.stenakath ' \
                  '{10} ' \
                  '{11} ' \
                  'left outer join "{1}"."{6}" g on g.id = s.gid ' \
                  'left outer join (SELECT id, direction, nz_min_art FROM "{1}"."{12}" GROUP BY id, direction, nz_min_art) ' \
                  '     AS wart ON wart.id = g.id AND wart.direction = s.direction ' \
                  'order by g.j ASC, g.i ASC, s.direction ASC, s.zs DESC' \
            .format(p_text, cfg.domain.case_schema, cfg.tables.surfaces, cfg.tables.roofs, cfg.tables.surface_params,
                    cfg.tables.walls, cfg.tables.grid, cfg.surf_range.roof_min, sqlline1, sqlline3, sqlline2, sqlline4, cfg.tables.building_walls)
        extra_verbose('\t{}', sqltext)
        cur.execute(sqltext)
        res = cur.fetchall()
        sql_debug(connection)
        varr = np.asarray([x[0] for x in res], dtype=vt)
        if any(np.isnan(varr)):
            nans = np.isnan(varr)
            nans = nans[nans]
            warning('NaN value(s) in parameter {} in function write_building_surface_pars [{} out of {}, which is {} %]', par, nans.size, ns, nans.size / ns)
        del res
        var[par, :] = np.nan_to_num(varr, copy=False, nan=cfg.fill_values[vt],
                                          posinf=cfg.fill_values[vt], neginf=cfg.fill_values[vt])
        debug('Variable {}, parameter {} written.', vn, par)
    # TODO: add vertical surfaces from building_walls to the building_surface_pars

    # check winfrac > 0.95 Honza's idea
    winfrac = np.argwhere(ncfile.variables[vn][1, :] > 0.95)
    if winfrac.size > 0:
        verbose('Modifying winfrac in building_surface_pars where winfrac > 0.95, '
                'setting winfrac = 0.95 and wallfrac = 0.05')
        for s_change in winfrac:
            ncfile.variables[vn][1, s_change] = 0.95
            ncfile.variables[vn][0, s_change] = 0.05
    debug('Variable {} completely written.', vn)


# write albedo parameters
def write_albedo_pars(ncfile,cfg, connection, cur):
    # write albedo_pars - roof parameters
    vn = 'albedo_pars'
    vt = 'f4'
    var = ncfile.createVariable(vn, vt, ('nalbedo_pars', 'y', 'x'), fill_value=cfg.fill_values[vt])
    nc_write_attribute(ncfile, vn, 'long_name', 'building parameters')
    nc_write_attribute(ncfile, vn, 'units', '')
    nc_write_attribute(ncfile, vn, 'res_orig', cfg.domain.dz)
    nc_write_attribute(ncfile, vn, 'coordinates', 'E_UTM N_UTM lon lat')
    nc_write_attribute(ncfile, vn, 'grid_mapping', 'E_UTM N_UTM lon lat')
    # fill individual parameters in dimension nalbedo_pars
    for par in cfg.albedo_pars._settings.keys():
        # par - index of npars dimmension in netcdf
        # albedo_pars[par] - calculation formulas for albedo for land, roof and wall
        # first check if albedo or emissivity column exists in landcover table
        sqltext = 'select {0} from "{1}"."{2}" as l limit 1'.format(cfg.albedo_pars[par][0],
                                                                    cfg.domain.case_schema, cfg.tables.landcover)
        al4l_exists = True
        try:
            cur.execute(sqltext)
        except Exception:
            al4l_exists = False
        if al4l_exists:
            al4l = cfg.albedo_pars[par][0]
        else:
            al4l = 'null'
        sql_debug(connection)
        connection.commit()
        sqltext = 'select case when b.id is null then {0} else {1} end ' \
                  'from "{2}"."{3}" g ' \
                  'left outer join "{2}"."{4}" l on l.lid = g.lid ' \
                  'left outer join "{2}"."{5}" b on b.id = g.id ' \
                  'left outer join "{2}"."{6}" r on r.gid = b.rid ' \
                  'order by g.j,g.i '.format(al4l, cfg.albedo_pars[par][1], cfg.domain.case_schema, cfg.tables.grid,
                                             cfg.tables.landcover, cfg.tables.buildings_grid, cfg.tables.roofs,
                                             cfg.tables.surface_params, )
        cur.execute(sqltext)
        res = cur.fetchall()
        sql_debug(connection)
        varr = np.reshape(np.asarray([x[0] for x in res], dtype=vt), (cfg.domain.ny, cfg.domain.nx))
        del res
        # process parameter of walls
        sqltext = 'select distinct on (g.i, g.j) ' \
                  'case when bw.id is not null then {0} else null end ' \
                  'from "{1}"."{2}" g ' \
                  'left outer join "{1}"."{3}" bw on bw.id = g.id and not bw.isroof ' \
                  'left outer join "{1}"."{4}" w on w.gid = bw.wid ' \
                  'left outer join "{1}"."{5}" pw on pw.code = w.stenakath ' \
                  'order by g.j,g.i ' \
            .format(cfg.albedo_pars[par][2], cfg.domain.case_schema, cfg.tables.grid, cfg.tables.building_walls,
                    cfg.tables.walls, cfg.tables.surface_params, )
        cur.execute(sqltext)
        res = cur.fetchall()
        sql_debug(connection)
        varw = np.reshape(np.asarray([x[0] for x in res], dtype=vt), (cfg.domain.ny, cfg.domain.nx))
        del res
        # combine the arrays
        varr = np.where(np.isnan(varw), varr, varw)
        # replace nan with fillValue and save in netcdf file
        var[par, :, :] = np.nan_to_num(varr, copy=False, nan=cfg.fill_values[vt], posinf=cfg.fill_values[vt], neginf=cfg.fill_values[vt])
        debug('Variable {}, parameter {} written.', vn, par)
        if cfg.visual_check.enabled:
            variable_visualization(var=var[par, ...],
                                   x=np.asarray(ncfile.variables['x']), y=np.asarray(ncfile.variables['y']),
                                   var_name=vn, par_id=par, text_id='albedo_type', path=cfg.visual_check.path,
                                   show_plots = cfg.visual_check.show_plots)
    debug('Variable {} completely written.', vn)
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

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
    """ write building_height (buildings_2d), building_id and building_type into netcdf file
    """
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
        var_3d = np.zeros((max_nz, cfg.domain.ny, cfg.domain.nx), dtype=np.int)

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
        var_bottom = np.reshape(np.asarray([x for x in res], dtype=np.bool), (cfg.domain.ny, cfg.domain.nx))
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

    debug('Crop lad / bad to configured max lad / bad')
    lad[lad > cfg.trees.max_lad] = cfg.trees.max_lad
    bad[bad > cfg.trees.max_bad] = cfg.trees.max_bad

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
    var[1:, :, :] = lad
    vn = 'bad'
    var = ncfile.createVariable(vn, vt, ('zlad', 'y', 'x'), fill_value=cfg.fill_values[vt])
    nc_write_attribute(ncfile, vn, 'long_name', 'branch area density')
    nc_write_attribute(ncfile, vn, 'units', 'm3/m3')
    nc_write_attribute(ncfile, vn, 'res_orig', cfg.domain.dz)
    nc_write_attribute(ncfile, vn, 'coordinates', 'E_UTM N_UTM lon lat')
    nc_write_attribute(ncfile, vn, 'grid_mapping', 'E_UTM N_UTM lon lat')
    var[1:, :, :] = bad
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
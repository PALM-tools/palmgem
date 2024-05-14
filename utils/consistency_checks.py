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

    if cfg.has_buildings_3d:
        mask_3d = np.logical_and(
                          np.logical_and(~ncfile.variables['building_id'][:].mask,
                                         ncfile.variables['buildings_3d'][0, :, :] == 0),  # is 3d structure with empty bottom
                          ~np.logical_or(np.logical_or(~ncfile.variables['pavement_type'][:,:].mask,
                                                       ~ncfile.variables['water_type'][:,:].mask),
                                         ~ncfile.variables['vegetation_type'][:,:].mask)) # not defined pavement / water / vegetation type
        missing_values = np.sum(mask_3d)
        if missing_values > 0:
            warning('There are {} missing values that were filled with default pavement type {}',
                    missing_values, pavement_type_default)
            extra_verbose('Missing grids: {}', np.where(mask_3d))
        mask = np.logical_or(mask, mask_3d)

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

    return ncfile

def check_cct_consistency(ncfile, cfg, connection, cur):
    """ Routine to check if each cct surface has it own type (land, wall, roof),
        And if land cct surface has at least (exactly) one type (vegetation, pavement, water) defined,
        and if wall or roof has building type defined.
    """
    # TODO: make it vector optimized
    cct_surface_type_classification    = ncfile.variables['cct_surface_type_classification'][:]
    cct_vegetation_type_classification = ncfile.variables['cct_vegetation_type_classification'][:]
    cct_pavement_type_classification   = ncfile.variables['cct_pavement_type_classification'][:]
    cct_water_type_classification      = ncfile.variables['cct_water_type_classification'][:]
    cct_building_type_classification   = ncfile.variables['cct_building_type_classification'][:]
    kji_locs                           = ncfile.variables['cct_3d_grid_indices'][:]

    n_surfs = cct_surface_type_classification.size
    allok = True
    for n_surf in range(n_surfs):
        k,j,i = kji_locs[:, n_surf]
        land = True if cct_surface_type_classification[n_surf] == 0 else False
        wall = True if cct_surface_type_classification[n_surf] == 1 else False
        roof = True if cct_surface_type_classification[n_surf] == 2 else False

        vege = True if cct_vegetation_type_classification[n_surf] > 0 else False
        vege_type = cct_vegetation_type_classification[n_surf]
        pave = True if cct_pavement_type_classification[n_surf] > 0 else False
        pave_type = cct_pavement_type_classification[n_surf]
        wate = True if cct_water_type_classification[n_surf] > 0 else False
        wate_type = cct_water_type_classification[n_surf]
        build = True if cct_building_type_classification[n_surf] > 0 else False
        build_type = cct_building_type_classification[n_surf]

        if not land and not wall and not roof:
            allok = False
            warning('Surf id: {}, kji: [{},{},{}], '
                    'does not have cct_surface_type_classification defined correcly',
                    n_surf, k, j, i)

        if not vege and not pave and not wate and not build:
            allok = False
            warning('Surf id: {}, kji: [{},{},{}], '
                    'does not have any type defined.',
                    n_surf, k, j, i)

        if land and not (vege or pave or wate):
            allok = False
            warning('Surf id: {}, kji: [{},{},{}], '
                    'is defined as land but vege, water, pavement types are mismatched. '
                    'Vege: {}, Pave: {}, Water: {}',
                    n_surf, k, j, i, vege_type, pave_type, wate_type)

        if (wall or roof) and not build:
            allok = False
            warning('Surf id: {}, kji: [{},{},{}], '
                    'is defined as urban but building type is mismatched. '
                    'Build type: {}',
                    n_surf, k, j, i, build_type)
    if allok:
        verbose('CCT check finished without a problem')
    else:
        error('CCT check finished with a problem')

def check_singular(edge_faces, faces, vert_kji, vert_len, ie):
    """ Check if vertex is singular, or vertex with the same [k,j,i] coordinate is singular"""
    if vert_len[ie - 1] == 0.0:
        return True
    adj_verts = []
    for iface in edge_faces:
        verts = faces[:, iface]
        for vt in verts:
            if vt > 0 and vt not in adj_verts and vt != ie:
                # print(vt, vert_kji[:, vt-1], '||||', vert_kji[:, ie])
                if (vert_kji[:, vt - 1] == vert_kji[:, ie]).all():
                    # print(vt, vert_kji[:, vt-1])
                    if vert_len[vt - 1] == 0.0:
                        # print('Possible singular point')
                        return True
                adj_verts.append(vt)
    return False

def cct_continuity_check(ncfile, cfg):
    """
    Function to check continuity in CCT. Check if structures are water-tight.
    Also check, if all TOPO / AIR gridcell corners are consitently defined thought out all surface.
    The condition that each vertex must have exactly 4 adjacent CCT surface is check.
    (Except the ones at domain boundary and singular points - only warning is printed
    """
    na_ = np.newaxis

    vertex_shift_nvect = np.array([[1, 0, 0],
                                   [-1, 0, 0],
                                   [0, 1, 0],
                                   [0, -1, 0],
                                   [0, 0, 1],
                                   [0, 0, -1]])

    extra_verbose('NOTE: All indices printed here are counted from 0!')
    progress('Checking cct consistency')

    debug('Loading static cct data')

    nx = len(ncfile.dimensions['x'])
    ny = len(ncfile.dimensions['y'])

    faces = ncfile.variables['cct_vertices_per_face'][:, :]
    face_nvert = ncfile.variables['cct_num_vertices_per_face'][:]

    vv = ncfile.variables['cct_vertex_coords']
    vert_kji = vv[0:3, :]
    vert_dir = vv[3, :]
    vv = ncfile.variables['cct_vertex_shifts']
    vert_len = vv[:].squeeze()

    kji_locs = ncfile.variables['cct_3d_grid_indices'][:]

    nfaces = faces.shape[1]
    nvert = vert_dir.shape[0]
    verbose('Loaded {} faces and {} vertices.', nfaces, nvert)

    allok = True

    debug('Checking stable assignment of corners')
    corners = {}
    for iv in range(nvert):
        full_corner = vert_kji[:, iv]
        shift = vertex_shift_nvect[vert_dir[iv]]  # from full corner towards vertex
        free_corner = tuple(full_corner + shift)
        full_corner = tuple(full_corner)

        full, last = corners.setdefault(full_corner, (True, iv))
        if not full:
            warning('Corner (k,j,i) {} is specified as full by vertex {}, '
                    'although it was already specified as free by vertex {}.', full_corner, iv + 1, last + 1)
            allok = False

        free, last = corners.setdefault(free_corner, (False, iv))
        if free:
            warning('Corner (k,j,i) {} is specified as full by vertex {}, '
                    'although it was already specified as free by vertex {}.', free_corner, iv + 1, last + 1)
            allok = False
    if allok:
        debug('Checked {} corners as correct.', len(corners))
    del corners

    debug('Checking continuity of face edges')
    for jf in range(nfaces):
        verts = faces[:, jf] - 1
        nv = face_nvert[jf]

        if verts[0] == verts[nv - 1]:
            warning('Face {} first and last vertex are same: {}.', jf, verts[:nv])
            allok = False

        full_corners = vert_kji[:, verts].T
        shifts = vertex_shift_nvect[vert_dir[verts]]  # from full corner towards vertex
        free_corners = full_corners + shifts

        for iv in range(nv - 1):
            iv2 = iv + 1
            if all(shifts[iv, :] == shifts[iv2, :]):
                # Identical vertex direction, edges must be 1 apart in
                # a different dimension.
                edge_diff = full_corners[iv2, :] - full_corners[iv, :]
                if (edge_diff[:] * shifts[iv, :]).sum():
                    warning('Face {} has consecutive edges with equal directions whose relative shift is nonzero in that direction: '
                            '({}, {}) -> ({}, {})',
                            jf, full_corners[iv], free_corners[iv], full_corners[iv2], free_corners[iv2])
                    allok = False

                if np.count_nonzero(edge_diff) != 1:
                    warning('Face {} has consecutive edges with equal directions that are not shifted in exactly one dimension: '
                            '({}, {}) -> ({}, {})',
                            jf, full_corners[iv], free_corners[iv], full_corners[iv2], free_corners[iv2])
                    allok = False

                if np.abs(edge_diff.sum()) != 1:
                    warning('Face {} has consecutive edges with equal directions that are not shifted exactly by one: '
                            '({}, {}) -> ({}, {})',
                            jf, full_corners[iv], free_corners[iv], full_corners[iv2], free_corners[iv2])
                    allok = False

            else:
                # Different directions, edges must share either inner or outer
                # corner.
                common_full = all(full_corners[iv, :] == full_corners[iv2, :])
                common_free = all(free_corners[iv, :] == free_corners[iv2, :])
                if not common_full and not common_free:
                    warning('Face {jf} has consecutive edges with different directions that do not have a common corner: '
                            '({}, {}) -> ({}, {})',
                            jf, full_corners[iv], free_corners[iv], full_corners[iv2], free_corners[iv2])
                    allok = False
        else:
            continue
    if allok:
        debug('Face edges are continuous.')

    debug('Checking number of adjacent faces per edge')
    edges = [[] for ie in range(nvert)]
    for jf in range(nfaces):
        verts = faces[:, jf] - 1
        nv = face_nvert[jf]

        for iv in range(nv):
            if not 0 <= verts[iv] < nvert:
                error('Face {} has {}. vertex = {} which is out of range.', jf, iv, verts[iv])
            edges[verts[iv]].append(jf)
    extra_verbose('... calculating sums')
    faces_per_edge = np.array(list(map(len, edges)))
    fhist = np.bincount(faces_per_edge, minlength=5)
    verbose('Number of edges adjacent to X faces: ' + ', '.join(f'{n}: {ne}' for n, ne in enumerate(fhist)))
    corner_edges = bnd_edges = 0
    for ie, edge_faces in enumerate(edges):
        c1 = vert_kji[:, ie]
        shift = vertex_shift_nvect[vert_dir[ie]]  # from full corner towards vertex
        c2 = c1 + shift

        nbound = 0
        if (c1[2] == c2[2] == 0) or (c1[2] == c2[2] == nx):
            nbound += 1
        if (c1[1] == c2[1] == 0) or (c1[1] == c2[1] == ny):
            nbound += 1
        if (c1[0] == c2[0] == 0):
            nbound += 1

        if nbound == 0:
            if len(edge_faces) != 4:
                msg = f'Standard edge with vertex {ie} has {len(edge_faces)} adjacent faces instead of 4: {edge_faces}.'
                singular = check_singular(edge_faces, faces, vert_kji, vert_len, ie)
                if singular:
                    extra_verbose('Singular edge')
                    extra_verbose('Singular edge; ' + msg)
                else:
                    warning(msg)
                    allok = False
        elif nbound == 1:
            if len(edge_faces) != 2:
                msg = f'Edge at domain boundary with vertex {ie} has {len(edge_faces)} adjacent faces instead of 2: {edge_faces}.'
                singular = check_singular(edge_faces, faces, vert_kji, vert_len, ie)
                if singular:
                    extra_verbose('Singular edge')
                    extra_verbose('Singular edge; ' + msg)
                else:
                    warning(msg)
                    allok = False
        elif nbound == 2:
            if len(edge_faces) != 1:
                msg = f'Edge with vertex {ie} at domain edge has {len(edge_faces)} adjacent faces instead of 1: {edge_faces}.'
                singular = check_singular(edge_faces, faces, vert_kji, vert_len, ie)
                if singular:
                    extra_verbose('Singular edge')
                    extra_verbose('Singular edge; ' + msg)
                else:
                    warning(msg)
                    allok = False
        else:
            error('Unexpected number of boundaries, this should never happen.')
    if allok:
        debug('Number of adjacent faces per edge is correct everywhere.')

    verbose('Examined {} cell edges along domain edges and {} edges on domain boundaries.', corner_edges, bnd_edges)

    if allok:
        progress('The file is fully consistent.')
    else:
        error('There were consistency errors (see above).')
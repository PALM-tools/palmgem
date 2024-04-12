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
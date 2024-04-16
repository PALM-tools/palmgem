import vtk
from vtk.util import numpy_support
import matplotlib.pyplot  as plt
import numpy as np
import os
from config.logger import *

def variable_visualization(var, x, y, var_name, par_id, text_id, path, scale=15, color_name='gist_rainbow', show_plots=False): #hsv
    """ create 2d plot from var 2d matrix"""
    # TODO: add posibility to directly add X, Y = lat, lon
    if not show_plots:
        plt.ioff()
    ny, nx = var.shape
    fig = plt.figure()
    ax = plt.subplot(111)
    X, Y = np.meshgrid(x, y)

    ax_sl = ax.imshow(var, aspect='equal', cmap=plt.get_cmap(color_name), origin='lower',
                      )
    fig.colorbar(ax_sl, extend='max', orientation = 'horizontal')

    plt.title('{}: {}\n'
              '[{}], '
              'min: {:.2f}, max: {:.2f}, avg: {:.2f}, stdev: {:.2f}'.format(
                var_name, par_id, text_id,
                np.min(var), np.max(var), np.average(var), np.std(var)))

    fig.savefig(os.path.join(path,'{}_{}.png'.format(var_name, par_id)),
                dpi=400)
    debug('{}: {} was successfully plotted', var_name, par_id)
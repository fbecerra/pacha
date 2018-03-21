import matplotlib
matplotlib.use('TkAgg')
import src as pr
import numpy as np
import pylab as pl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1 import ImageGrid

import time
import h5py

def size_label(size):
  if size >= 1e-3:
    label = '%i pc' %(size*1e3)
  elif size >= 1e-4:
    label = '%.1f pc' %(size*1e3)
  elif size >= 1e-5:
    label = '%.2f pc' %(size*1e3)
  return label

def type_label(sim):
  if sim == 'nahw1r4sm1' or sim == 'nahw1r4sm2' or sim == 'nahw1r4sm3':
    label = 'sink'
  elif sim == 'nahw1r4ad1' or sim == 'nahw1r4ad2' or sim == 'nahw1r4ad3':
    label = 'core'
  return label

def nth_label(sim):
  if sim == 'nahw1r4sm3' or sim == 'nahw1r4ad3':
    label = r'$n_{\rm th} = 10^{12}\,{\rm cm}^{-3}$'
  elif sim == 'nahw1r4sm2' or sim == 'nahw1r4ad2':
    label = r'$n_{\rm th} = 10^{10}\,{\rm cm}^{-3}$'
  elif sim == 'nahw1r4sm1' or sim == 'nahw1r4ad1':
    label = r'$n_{\rm th} = 10^8\,{\rm cm}^{-3}$'
  return label

start_total = time.time()

#sims, snaps, size, nhmin, nhmax, tmin, tmax, step = ['nahw1r4sm1', 'nahw1r4sm2'], [7, 15], 3e-3, 4.7, 8, 3.7, 4.0, 0.1
#sims, snaps, size, nhmin, nhmax, tmin, tmax, step = ['nahw1r4sm2', 'nahw1r4sm3'], [10, 26], 3e-4, 7, 10, 3.7, 4.0, 0.1
#sims, snaps, size, nhmin, nhmax, tmin, tmax, step = ['nahw1r4sm1', 'nahw1r4ad1'], [7, 31], 3e-3, 4.7, 8, 3.7, 4.0, 0.1
#sims, snaps, size, nhmin, nhmax, tmin, tmax, step = ['nahw1r4sm2', 'nahw1r4ad2'], [10, 15], 3e-4, 7, 10, 3.7, 4.0, 0.1
sims, snaps, size, nhmin, nhmax, tmin, tmax, step = ['nahw1r4sm3', 'nahw1r4ad3'], [17, 32], 3e-5, 8, 12, 3.7, 4.0, 0.1
base = '/n/hernquistfs2/fbecerra/'
fields = ['nh', 'temp']
nsims = len(sims)

fig = pl.figure()

ncols, nrows = len(sims), len(fields)
left, right, bottom, top = 0.1, 0.9, 0.1, 0.9
width = (right - left)/ncols
height = (top - bottom)/nrows

for idx_sim, sim in enumerate(sims):
 
  print 'Sim: %s' %sim

  snap = snaps[idx_sim]

  path = base+sim+'/snapdir_%03i/' %snap
  file = sim+'_%03i' %snap
  snapbase = path+file
  
  MySnap = pr.snap.Snap()
  MySnap.read_header(snapbase)
  MySnap.read_fields(snapbase)
  MySnap.center_box()
  MySnap.rotate_box()

  for idx_field, field in enumerate(fields):

    MyImage = pr.image.Image()
    MyImage.calculate_image(MySnap, field, ImgWidth=size)

    row = idx_field
    col = idx_sim

    idx_ax = row * ncols + col + 1

    ax = fig.add_subplot(nrows, ncols, idx_ax)
    ax.set(aspect=1)
    # Bottom left: size
    ax.text(MyImage.xbins/20., MyImage.ybins*1./20., size_label(size), fontsize=36, color='w')
    # Bottom right: type
    ax.text(MyImage.xbins*8./10., MyImage.ybins*1./20., type_label(sim), fontsize=36, color='w')
    # Top left: n_th
    ax.text(MyImage.xbins/20., MyImage.ybins*9./10., nth_label(sim), fontsize=36, color='w')

    if field == 'nh':
      vmin, vmax = nhmin, nhmax
    else:
      vmin, vmax = tmin, tmax
    #elif col == 1:
    #  if field == 'nh':
    #    vmin, vmax = 7.1, 10
    #  else:
    #    vmin, vmax = 3.6, 4.2
    #elif col == 2:
    #  if field == 'nh':
    #    vmin, vmax = 7.9, 12
    #  else:
    #    vmin, vmax = 3.6, 4.8
    im = ax.imshow(MyImage.img, cmap = pr.utils.get_colormap(field), extent = [0, MyImage.xbins, 0, MyImage.ybins], vmin = vmin, vmax = vmax)

    try:
      nsinks = len(MySnap.new_fields['sinks']['id'])
      if nsinks > 0:
        x_sinks = MyImage.xbins / 2. + MyImage.xbins * MySnap.new_fields['sinks']['x'] / MyImage.width
        y_sinks = MyImage.ybins / 2. + MyImage.ybins * MySnap.new_fields['sinks']['y'] / MyImage.height
        for cx, cy in zip(x_sinks, y_sinks):
          if ((cx > 0) & (cx < MyImage.xbins) & (cy > 0) & (cy < MyImage.ybins)):
            ax.plot(cx, cy, 'w+', ms=5, mew=3)
    except:
      print 'Snapshot does not have sink particles'

    # Fixing limits
    ax.set_xlim(0, MyImage.xbins)
    ax.set_ylim(0, MyImage.ybins)

    ax.set_xticklabels([])
    ax.set_yticklabels([])
  
    # Colorbar
    if col == 0:
      if field == 'nh':
        ticks = np.arange(vmin, vmax + 1).astype(int)
        labelpad = -85
      elif field == 'temp':
        ticks = np.arange(vmin, vmax + 0.05, step)
        labelpad = -85
      if row == 0:
        cax = pl.axes([left+col*width+0.02, bottom+nrows*height+0.04, 2*width-0.04, 0.02])
        cb = pl.colorbar(im, cax = cax, orientation='horizontal', ticks=ticks)
        cax.set_xlabel(pr.utils.get_label(field), fontsize=30, labelpad=labelpad)
        cb.ax.tick_params(labelsize=24)
      elif row == 1:
        cax = pl.axes([left+col*width+0.02, bottom - 0.07, 2*width-0.04, 0.02])
        cb = pl.colorbar(im, cax = cax, orientation='horizontal', ticks=ticks)
        cax.set_xlabel(pr.utils.get_label(field), fontsize=30, labelpad=labelpad)
        cb.ax.tick_params(labelsize=24)

pl.subplots_adjust(left=left, right=right, bottom=bottom, top=top, wspace=0.0, hspace=0.0)
fig.set_size_inches(ncols*8, nrows*8)

fig.savefig('/n/home02/fbecerra/pacha_new/images/convergence_'+'_'.join(sims)+'.pdf', dpi=100)
pl.show()
pl.close()

end_total = time.time()
print 'Took: %f seconds' %(end_total - start_total) 

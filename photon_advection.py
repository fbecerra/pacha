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

def get_accretion_radius(mass):

  n_thres = 1e8
  temp = 8000.
  n_bondi = 4.66e13 * (temp / 6000.)**3.  / mass**2.

  if n_thres >= n_bondi:
    acc_radius = 17 * mass * (6000 / temp) * pr.constants.ASTRONOMICAL_UNIT / pr.constants.UNIT_LENGTH
  else:
    acc_radius = np.sqrt(1.49e17 / n_thres) * np.sqrt(temp / 6000) * pr.constants.ASTRONOMICAL_UNIT / pr.constants.UNIT_LENGTH

  return acc_radius

start_total = time.time()

sim = 'nahw1r4rt1'
snaps = [4, 8, 12, 15]
size = 2e-2
base = '/n/hernquistfs2/fbecerra/'
fields = ['phodens']
time0 = 1.27917e-6

fig = pl.figure()

ncols, nrows = 2, 2
left, right, bottom, top = 0.1, 0.9, 0.1, 0.9
width = (right - left)/ncols
height = (top - bottom)/nrows

for idx_snap, snap in enumerate(snaps):

  print 'Snapshot: %s' %snap

  path = base+sim+'/advection/snapdir_%03i/' %snap
  file = sim+'_%03i' %snap
  snapbase = path+file
  
  MySnap = pr.snap.Snap()
  MySnap.read_header(snapbase)
  MySnap.read_fields(snapbase)
  MySnap.center_box()
  MySnap.rotate_box()

  time_snap = (MySnap.params['time'] - time0) * pr.constants.UNIT_TIME / pr.constants.SEC_PER_YEAR

  # Center on protostar
  try:
    idx_sink = np.argmax(MySnap.new_fields['sinks']['mass'])
    x_center = MySnap.new_fields['sinks']['x'][idx_sink]
    y_center = MySnap.new_fields['sinks']['y'][idx_sink]
    z_center = MySnap.new_fields['sinks']['z'][idx_sink]

    MySnap.fields['x'] -= x_center
    MySnap.fields['y'] -= y_center
    MySnap.fields['z'] -= z_center

    MySnap.new_fields['sinks']['x'] -= x_center
    MySnap.new_fields['sinks']['y'] -= y_center
    MySnap.new_fields['sinks']['z'] -= z_center
  except:
    print 'Snapshot does not have sink particles'

  for idx_field, field in enumerate(fields):

    MyImage = pr.image.Image()
    MyImage.calculate_image(MySnap, field, ImgWidth=size)

    row = int(idx_snap/ncols)
    col = (idx_snap%ncols)

    idx_ax = row * ncols + col + 1

    ax = fig.add_subplot(nrows, ncols, idx_ax)
    ax.set(aspect=1)
    # Top left: size
    ax.text(MyImage.xbins/20., MyImage.ybins*9./10., size_label(size), fontsize=36, color='w')
    # Bottom right: time
    ax.text(MyImage.xbins*7.5/10., MyImage.ybins*1./20., '%i yr' %time_snap, fontsize=36, color='w')
    # Top left: n_th
    #ax.text(MyImage.xbins/20., MyImage.ybins*9./10., nth_label(sim), fontsize=36, color='w')

    vmin, vmax = np.min(MyImage.img), np.max(MyImage.img)
    print vmin, vmax
    if field == 'nh':
      vmin, vmax = 5.5, 7.5
    elif field == 'phodens':
      vmin, vmax = 0, 14
    im = ax.imshow(MyImage.img, cmap = pr.utils.get_colormap(field), extent = [0, MyImage.xbins, 0, MyImage.ybins], vmin = vmin, vmax = vmax)

    try:
      nsinks = len(MySnap.new_fields['sinks']['id'])
      if nsinks > 0:
        print MySnap.new_fields['sinks']['x'], MySnap.new_fields['sinks']['y'], MySnap.new_fields['sinks']['mass']
        print get_accretion_radius(MySnap.new_fields['sinks']['mass'])
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
    if field == 'nh':
      ticks = np.arange(vmin, vmax + 1).astype(int)
      labelpad = -85
    elif field == 'phodens':
      ticks = np.arange(vmin, vmax + 1, 4).astype(int)
      labelpad = -85

    if row == 0:
      cax = pl.axes([left+col*width+0.01, bottom+nrows*height+0.04, width-0.02, 0.02])
      cb = pl.colorbar(im, cax = cax, orientation='horizontal', ticks=ticks)
      cax.set_xlabel(pr.utils.get_label(field), fontsize=30, labelpad=labelpad)
      cb.ax.tick_params(labelsize=24)
    elif row == 1:
      cax = pl.axes([left+col*width+0.01, bottom - 0.07, width-0.02, 0.02])
      cb = pl.colorbar(im, cax = cax, orientation='horizontal', ticks=ticks)
      cax.set_xlabel(pr.utils.get_label(field), fontsize=30, labelpad=labelpad)
      cb.ax.tick_params(labelsize=24)

pl.subplots_adjust(left=left, right=right, bottom=bottom, top=top, wspace=0.0, hspace=0.0)
fig.set_size_inches(ncols*8, nrows*8)

fig.savefig('/n/home02/fbecerra/pacha_new/images/RT_advection_'+field+'.pdf', dpi=100)
pl.show()
pl.close()

end_total = time.time()
print 'Took: %f seconds' %(end_total - start_total) 

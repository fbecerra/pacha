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
  if size > 1e3:
    label = '%i kpc' %(size/1e3)
  elif size < 1:
    label = '%f pc' %size
  else:
    label = '%i pc' %size
  return label

start_total = time.time()

#sim = 'nahw1'
#snap = 137
sim = 'nahw1r4ad3'
snap = 26
base = '/n/hernquistfs2/fbecerra/'
#sizes = [20000, 2000, 200, 20]
#fields = ['nh', 'temp']
sizes = [0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001]
fields = ['nh']
nsizes = len(sizes)
nfields = len(fields)

path = base+sim+'/snapdir_%03i/' %snap
file = sim+'_%03i' %snap
snapbase = path+file

MySnap = pr.snap.Snap()
MySnap.read_header(snapbase)
MySnap.read_fields(snapbase)
MySnap.center_box()
MySnap.rotate_box()

fig = pl.figure()

ncols, nrows = nfields, nsizes
left, right, bottom, top = 0.1, 0.9, 0.1, 0.9
width = (right - left)/ncols
height = (top - bottom)/nrows

for idx_size, size in enumerate(sizes):
  for idx_field, field in enumerate(fields):

    print 'Size: %i, field: %s' %(size, field)

    MyImage = pr.image.Image()
    MyImage.calculate_image(MySnap, field, ImgWidth=size)

    idx_ax = nfields * idx_size + idx_field + 1

    ax = fig.add_subplot(nrows, ncols, idx_ax)
    ax.set(aspect=1)
    ax.text(MyImage.xbins/20., MyImage.ybins*9./10., size_label(size), fontsize=36, color='w')

    vmin, vmax = np.min(MyImage.img), np.max(MyImage.img)
    #if idx_size == 0:
    #  if idx_field == 0:
    #    vmin, vmax = -4, 4
    #  else:
    #    vmin, vmax = 2.4, 4.0
    #elif idx_size == 1:
    #  if idx_field == 0:
    #    vmin, vmax = -2, 5
    #  else:
    #    vmin, vmax = 2.4, 4.0
    #elif idx_size == 2:
    #  if idx_field == 0:
    #    vmin, vmax = 0, 6
    #  else:
    #    vmin, vmax = 3.6, 4.0
    #elif idx_size == 3:
    #  if idx_field == 0:
    #    vmin, vmax = 2, 7
    #  else:
    #    vmin, vmax = 3.8, 4.0
    #vmin, vmax = 4, 8
    im = ax.imshow(MyImage.img, cmap = pr.utils.get_colormap(field), extent = [0, MyImage.xbins, 0, MyImage.ybins], vmin = vmin, vmax = vmax)

    if idx_size != nsizes - 1:
      dx, dy = float(sizes[idx_size+1]) / float(size) * MyImage.xbins, float(sizes[idx_size+1]) / float(size) * MyImage.ybins
      xi, xf = MyImage.xbins / 2. - dx / 2., MyImage.xbins / 2. + dx / 2.
      yi, yf = MyImage.ybins / 2. - dy / 2., MyImage.ybins / 2. + dy / 2.
      # Square
      square = pl.Rectangle((xi, yi), dx, dy, facecolor ='none', edgecolor='k')
      ax.add_patch(square)
      x1, y1 = np.array([xi, 0]), np.array([yf, 0])
      x2, y2 = np.array([xf, MyImage.xbins]), np.array([yf, 0])
      # Lines connecting edges
      ax.plot(x1, y1, 'k', linewidth=1.0)
      ax.plot(x2, y2, 'k', linewidth=1.0)
      # Fixing limits
      ax.set_xlim(0, MyImage.xbins)
      ax.set_ylim(0, MyImage.ybins)

    ax.set_xticklabels([])
    ax.set_yticklabels([])
  
    # Colorbar
    if idx_field == 0:
      cax = pl.axes([left-0.03, bottom+(nsizes-idx_size-1)*height+0.01, 0.02, height-0.02])
      #cax = pl.axes([left+(idx_ax-1)*width+0.01, bottom+nrows*height+0.04, width-0.02, 0.02])
      cb = pl.colorbar(im, cax = cax, orientation='vertical', ticks=np.arange(vmin, vmax + 1).astype(int))
      cax.yaxis.set_ticks_position('left')
      cax.yaxis.set_label_position('left')
      cax.set_xlabel(pr.utils.get_label(field), fontsize=30, rotation=90, x=-2.5, labelpad=-300)
      cb.ax.tick_params(labelsize=24)
    elif idx_field == 1:
      if idx_size == 0 or idx_size == 1:
        step = 0.4
      else:
        step = 0.1
      cax = pl.axes([left+ncols*width+0.01, bottom+(nsizes-idx_size-1)*height+0.01, 0.02, height-0.02])
      #cax = pl.axes([left+(idx_ax-4)*width+0.01, bottom - 0.07, width-0.02, 0.02])
      cb = pl.colorbar(im, cax = cax, orientation='vertical', ticks=np.arange(vmin, vmax + 0.05, step))
      cax.yaxis.set_ticks_position('right')
      cax.yaxis.set_label_position('right')
      cax.set_xlabel(pr.utils.get_label(field), fontsize=30, rotation=90, x=3.5, labelpad=-270)
      cb.ax.tick_params(labelsize=24)

pl.subplots_adjust(left=left, right=right, bottom=bottom, top=top, wspace=0.0, hspace=0.0)
fig.set_size_inches(nfields*8, nsizes*8)

fig.savefig('/n/home02/fbecerra/pacha_new/images/zoom_view_'+sim+'.pdf', dpi=100)
pl.show()
pl.close()

end_total = time.time()
print 'Took: %f seconds' %(end_total - start_total) 

import matplotlib
matplotlib.use('Agg')
import src as pr
import numpy as np
import pylab as pl
import time
from mpl_toolkits.axes_grid1 import ImageGrid

def get_label(sim):
  if sim == 'nahw1r4sm3' or sim == 'nahw1r4ad3':
    label = r'$n_{\rm th} = 10^{12}\,{\rm cm}^{-3}$'
  elif sim == 'nahw1r4sm2' or sim == 'nahw1r4ad2':
    label = r'$n_{\rm th} = 10^{10}\,{\rm cm}^{-3}$'
  elif sim == 'nahw1r4sm1' or sim == 'nahw1r4ad1':
    label = r'$n_{\rm th} = 10^8\,{\rm cm}^{-3}$'
  return label



sims, snaps, cmaps, colors, shifts = ['nahw1r4ad1', 'nahw1r4ad2', 'nahw1r4ad3'], [31, 23, 32], ['Reds', 'Blues', 'Greens'], ['#e41a1c', '#377eb8', '#4daf4a'], [3.4, 4.6, 5.6]
base = '/n/hernquistfs2/fbecerra/'
field = 'temp'
gray = '#999999'

fig = pl.figure()
ax = pl.subplot(1, 1, 1)
MyImage = pr.image.Image()

xmin, xmax = 2, 16
ymin, ymax = 3, 5

for idx_sim, sim in enumerate(sims):
  print 'Sim: %s' %sim
  snap = snaps[idx_sim]
  cmap = cmaps[idx_sim]
  color = colors[idx_sim]
  shift = shifts[idx_sim]

  path = base+sim+'/snapdir_%03i/' %snap
  file = sim+'_%03i' %snap
  snapbase = path+file
    
  MySnap = pr.snap.Snap()
  MySnap.read_header(snapbase)
  MySnap.read_fields(snapbase)
  MySnap.center_box()
  MySnap.rotate_box()
 
  MyImage.pspace_plot(MySnap, field)
  ax.imshow(MyImage.img, cmap = cmap, extent = [MyImage.xmin, MyImage.xmax, MyImage.ymin, MyImage.ymax], aspect='auto', interpolation='none', alpha = 0.7)
  ax.plot(MyImage.valx, MyImage.valy, color, alpha=0.8, label=get_label(sim))

  pl.plot(np.array([xmin, xmax]), 2./3*np.array([xmin,xmax]) - shift, ':', color=gray)

pl.legend(loc='upper left', fontsize=24)

# y-axis
pl.ylim(ymin, ymax)
pl.ylabel(pr.utils.get_label(field), fontsize=32)
locs, labels = pl.yticks()
locs = np.array(locs)
idx = np.where(locs % 1 == 0)[0]
pl.yticks(locs[idx], fontsize=24)

# x-axis
pl.xlim(xmin, xmax)
pl.xlabel(pr.utils.get_label('nh'), fontsize=32)
locs, labels = pl.xticks()
locs = np.array(locs)
idx = np.where(locs % 1 == 0)[0]
pl.xticks(locs[idx], fontsize=24)

fig.set_size_inches(10, 9)
fig.savefig('/n/home02/fbecerra/pacha_new/images/pspace_all_core_'+field+'.pdf')
pl.show()
pl.close()

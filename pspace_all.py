import matplotlib
matplotlib.use('Agg')
import src as pr
import numpy as np
import pylab as pl
import time
from mpl_toolkits.axes_grid1 import ImageGrid

sims = ['ah1w3f2', 'ah1w3r0', 'ah1w3']
snaps = ['001', '022', '157']
base= '/scratch/00025/tgreif/'
field = 'temp'
f2snaps = np.append(np.arange(36,126,8), np.arange(132,154,7))
#f2snaps = np.arange(36, 47, 10)
min_dens = 15

fig = pl.figure()
MyImage = pr.image.Image()
nsnaps = len(f2snaps)
ncols = int(np.ceil(np.sqrt(nsnaps)))
nrows = nsnaps / ncols

for snap_idx in range(nsnaps):
  f2snap = f2snaps[snap_idx]
  snaps[0] = '%03d' %f2snap

  for index in range(len(sims)):
  
    path = base+sims[index]+'/snapdir_'+snaps[index]+'/'
    file = sims[index]+'_'+snaps[index]
    snapbase = path+file
    
    if index == 0:
      MySnap = pr.snap.Snap()
      MySnap.read_header(snapbase)
      MySnap.read_fields(snapbase)
      MySnap.center_box()
      MySnap.rotate_box()
    else:
      MySnap2 = pr.snap.Snap()
      MySnap2.read_header(snapbase)
      MySnap2.read_fields(snapbase)
      MySnap2.center_box()
      MySnap2.rotate_box()
      MySnap.merge_snaps(MySnap2)
      del MySnap2
 
  MySnap.calculate_radius()
  part_idx = np.where(MySnap.fields['nh'] <= min_dens)
  MySnap.remove_particles(part_idx) 

  MyImage.pspace_plot(MySnap, field)
  ax = pl.subplot(nrows, ncols,snap_idx+1)
  ax.imshow(MyImage.img, cmap = 'jet', extent = [MyImage.xmin, MyImage.xmax, MyImage.ymin, MyImage.ymax], aspect='auto')
  ax.plot(MyImage.valx, MyImage.valy, 'k')

  pl.ylim(3, 6)
  pl.xlim(min_dens, 22)
  if (((snap_idx+1) % ncols == 1) and (snap_idx / ncols == nrows -1)):
    # y-axis
    pl.ylabel(pr.utils.get_label(field), fontsize=32)
    locs, labels = pl.yticks()
    locs = np.array(locs)
    idx = np.where(locs % 1 == 0)[0]
    pl.yticks(locs[idx], fontsize=24)

    # x-axis
    pl.xlabel(pr.utils.get_label('nh'), fontsize=32)
    locs, labels = pl.xticks()
    locs = np.array(locs)
    idx = np.where(locs % 1 == 0)[0]
    pl.xticks(locs[idx], fontsize=24)
  else:
    ax.set_xticklabels([])
    ax.set_yticklabels([])

pl.subplots_adjust(wspace=0.0, hspace=0.0)
fig.set_size_inches(ncols*8, nrows*8)
fig.savefig('/scratch/02563/fbecerra/paha/pspace/timelapse_pspace_all_'+field+'.pdf')
pl.show()
pl.close()

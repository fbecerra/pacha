import matplotlib
matplotlib.use('Agg')
import src as pr
import numpy as np
import pylab as pl
import time

#frame, sim, snap_init, snap_final = 1, 'nahw1', 101, 138
frame, sim, snap_init, snap_final = 500, 'nahw1r2', 1, 30

base= '/scratch/00025/tgreif/'

plot_fields = ['nh', 'enc_mass', 'temp', 'abH2']

start = time.time()

if sim == 'nahw1r2':
  previous_snap = '/n/hernquistfs2/fbecerra/nahw1/snapdir_138/nahw1_138'
  PreviousSnap = pr.snap.Snap()
  PreviousSnap.read_header(previous_snap)
  PreviousSnap.read_fields(previous_snap)
  PreviousSnap.center_box()
  PreviousSnap.rotate_box()
  print 'Read previous snap'

for snap in xrange(snap_init, snap_final + 1):

  print 'Snap: ', snap
  base = '/n/hernquistfs2/fbecerra/'
  path = base+sim+'/snapdir_%03i/' %snap
  file = sim+'_%03i' %snap
  snapbase = path+file

  MySnap = pr.snap.Snap()
  MySnap.read_header(snapbase)
  MySnap.read_fields(snapbase)
  MySnap.center_box()
  MySnap.rotate_box()
  if sim == 'nahw1r2':
    MySnap.merge_snaps(PreviousSnap)
  
  MySnap.calculate_radius()
  
  MyRadial = pr.radial.Radial()
  MyRadial.radial_profile(MySnap, plot_fields)
  del MySnap
  
  fig = pl.figure()

  for idx, field in enumerate(plot_fields):
    nfields = len(plot_fields)
    if nfields > 1:
      ncols = 2
    else:
      ncols = 1
    nrows = nfields / ncols
  
    ax = pl.subplot(nrows, ncols,idx+1)
    xval = 'radius'
    pl.plot(MyRadial.radial[xval], MyRadial.radial[field], color='#377eb8')
    #xmin, xmax = np.min(MyRadial.radial[xval]), np.max(MyRadial.radial[xval])
    #xmin, xmax = -10, 2 #kpc
    xmin, xmax = -7, 5
    
    # y-axis
    ymin, ymax = pl.ylim()
    if field == 'nh':
      pl.ylim(-4, 18)
    if field == 'abH2':
      pl.ylim(-17, -4)
    if field == 'temp':
      pl.ylim(1, 4.5)
    if field == 'enc_mass':
      pl.ylim(-6, 12)
    if (ncols > 1 and (idx+1) % ncols == 0):
      ax.yaxis.tick_right()
      ax.yaxis.set_ticks_position('both')
      ax.yaxis.set_label_position('right')
    pl.ylabel(pr.utils.get_label(field), fontsize=24)
    locs, labels = pl.yticks()
    locs = np.array(locs)
    ticks_idx = np.where(locs % 1 == 0)[0]
    pl.yticks(locs[ticks_idx])
    
    # x-axis
    if (idx / ncols == nrows -1):
      pl.xlabel(pr.utils.get_label(xval), fontsize=24)
#      locs, labels = pl.xticks()
#      locs = np.array(locs)
#      idx = np.where(locs % 1 == 0)[0]
#      pl.xticks(locs[idx])
      locs = np.arange(xmin, xmax, dtype=int)
      if len(locs) > 6:
        locs = locs[np.where(locs % 2 == 0)[0]]
      pl.xticks(locs)
    else:
      ax.set_xticklabels([])
    
    ax.tick_params(axis='both', labelsize=20)
    pl.xlim(xmin, xmax)
  
  del MyRadial

  pl.subplots_adjust(wspace=0.0, hspace=0.0)
  #fig.tight_layout()
  fig.set_size_inches(6*ncols, 6*nrows)
  pl.savefig('/n/home02/fbecerra/pacha_new/radial/movie_radial_'+sim+'_%04d.pdf' %snap, dpi=100)

  pl.show()
  pl.close()

  frame += 1

end = time.time()
print 'Total time: %f' %(end - start)

import matplotlib
matplotlib.use('TkAgg')
import src as pr
import numpy as np
import pylab as pl
import time

#path = '/Volumes/DARK-MATTER/arepo/ah1w3f0/snapdir_238/'
#file = 'ah1w3f0_238'
#snapbase = path+file

start_total = time.time()

plot_fields = ['nh', 'enc_mass', 'temp', 'abH2']
#plot_fields = ['nh', 'angmom', 'taugrav', 'taupres', 'tgrav', 'tpres']
#plot_fields = ['nh', 'angmom', 'taugrav', 'taupres']
#plot_fields = ['sigma_gas', 'cs', 'omega_rot', 'q']
#plot_fields = ['geff', 'dnh', 'tcool_ratio', 'tcs_ratio'] 
#plot_fields = ['tcool', 'omega_rot', 'gammie', 'tcool_ratio']

### Sinks
#sim, snap_init, snap_final = 'nahw1r4sm1', 1, 20
#sim, snap_init, snap_final = 'nahw1r4sm2', 7, 16
sim, snap_init, snap_final = 'nahw1r4sm3', 14, 26

### Adiabatic
#sim, snap_init, snap_final = 'nahw1r4ad1', 1, 31
#sim, snap_init, snap_final = 'nahw1r4ad2', 7, 23
#sim, snap_init, snap_final = 'nahw1r4ad3', 14, 26

### Extras
#sim, snap_init, snap_final = 'nahw1', 137, 137

for snap in xrange(snap_init, snap_final+1):
  print 'Snapshot: ', snap
  #path = '/n/hernquistfs3/pmocz/Output/CosmologyAxion/MiniHaloSims/output/mh2/snapdir_%03i/' %snap
  path = '/n/hernquistfs2/fbecerra/'+sim+'/snapdir_%03i/' %snap
  file = sim+'_%03i' %snap
  #file = 'mh2_%03i' %snap
  snapbase = path+file
  MySnap = pr.snap.Snap()
  MySnap.read_header(snapbase)
  MySnap.read_fields(snapbase)
  MySnap.center_box()
  MySnap.rotate_box()
  MySnap.calculate_radius()
  pc = pr.plot_collection.PlotCollection()
  pc.multi_radial(MySnap, plot_fields)
  pc.savefig('./radial/multi_radial_'+file+'_'+'_'.join(plot_fields)+'.pdf')
print 'Took: %f seconds' %(time.time() - start_total)

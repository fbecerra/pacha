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
#plot_fields = ['sigma_gas', 'cs', 'omega_rot', 'q']
#plot_fields = ['geff', 'dnh', 'tcool_ratio', 'tcs_ratio'] 
#plot_fields = ['tcool', 'omega_rot', 'gammie', 'tcool_ratio']

for snap in xrange(34,40):

  path = '/n/hernquistfs3/pmocz/Output/CosmologyAxion/MiniHaloSims/output/mh2/snapdir_%03i/' %snap
  #path = '/n/hernquistfs2/fbecerra/nahw1tgs1/snapdir_%03i/' %snap
  #file = 'nahw1tgs1_%03i' %snap
  file = 'mh2_%03i' %snap
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

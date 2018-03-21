import matplotlib
matplotlib.use('TkAgg')
import src as pr
import numpy as np
import pylab as pl
import time

#snapbase = path+file

start_total = time.time()

sims = ['LLGas512', 'BLGas512', 'LBGas512', 'BBGas512']
path = '/n/hernquistfs3/pmocz/Output/CosmologyAxion/BECDMCosmologicalSimsWithBaryons/output/'
base = 'snap'
plot_fields = ['nh', 'temp']
#plot_fields = ['nh', 'angmom', 'taugrav', 'taupres', 'tgrav', 'tpres']
#plot_fields = ['sigma_gas', 'cs', 'omega_rot', 'q']
#plot_fields = ['geff', 'dnh', 'tcool_ratio', 'tcs_ratio'] 
#plot_fields = ['tcool', 'omega_rot', 'gammie', 'tcool_ratio']

for sim in sims:
  for snap in np.arange(10, 19):
    print 'Sim: %s, Snapshot: %i' %(sim, snap)
    snapbase = path + sim + '/snapdir_%03i/' %snap + base + '_%03i' %snap
    MySnap = pr.snap.Snap()
    MySnap.read_header(snapbase, read_HDF=True)
    MySnap.read_fields(snapbase, read_HDF=True)
    MySnap.center_box()
    MySnap.calculate_radius()
    pc = pr.plot_collection.PlotCollection()
    pc.multi_radial(MySnap, plot_fields)
    pc.savefig('/n/home02/fbecerra/pacha_new/radial/'+sim+'_'+'_'.join(plot_fields)+'%03i.pdf' %snap)
print 'Took: %f seconds' %(time.time() - start_total)

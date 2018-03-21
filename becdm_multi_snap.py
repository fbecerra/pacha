import matplotlib
matplotlib.use('TkAgg')
import src as pr
import numpy as np
import pylab as pl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1 import ImageGrid

import time
import h5py

start_total = time.time()


fields = ['nh', 'temp']#, 'abH2'] 
#sim = 'nahw1s2tb1'
sim = 'nahw1r3tgs1'
#sim = 'sink_test'
#fields = ['temp']
#snaps = range(36, 126, 8) # Fragmentation before binary
#snaps = range(146,178,2) # Protostar follow-up after binary
#snaps = np.append(np.arange(36,126,8), np.arange(132,154,7))
for snaps in np.arange(4,9):
   print 'Snapshot ', snaps
   pc = pr.plot_collection.PlotCollection()
   pc.multi_snap([snaps], fields, '/n/hernquistfs2/fbecerra/', sim, 1e-3, PlotTime=False, PlotSize=True, CenterProto=0, PlotProto=False, Rotate=True)
   pc.savefig('/n/home02/fbecerra/pacha_new/images/multi_snap_'+sim+'_%03i.pdf' %snaps)
end_total = time.time()
print 'Took: %f seconds' %(end_total - start_total) 

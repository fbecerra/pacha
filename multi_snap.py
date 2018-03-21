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

### Sinks 
#sim, snap_init, snap_final, img_size = 'nahw1r4sm1', 1, 20, 3
#sim, snap_init, snap_final, img_size = 'nahw1r4sm2', 7, 16, 0.3
sim, snap_init, snap_final, img_size = 'nahw1r4sm3', 14, 26, 0.03

### Adiabatic 
#sim, snap_init, snap_final, img_size = 'nahw1r4ad1', 1, 31, 3
#sim, snap_init, snap_final, img_size = 'nahw1r4ad2', 7, 23, 0.3
#sim, snap_init, snap_final, img_size = 'nahw1r4ad3', 14, 32, 0.03

for snaps in np.arange(snap_init, snap_final+1):
   print 'Snapshot ', snaps
   pc = pr.plot_collection.PlotCollection()
   pc.multi_snap([snaps], fields, '/n/hernquistfs2/fbecerra/', sim, img_size*1e-3, PlotTime=True, PlotSize=True, CenterProto=0, PlotProto=False, Rotate=True)
   pc.savefig('/n/home02/fbecerra/pacha_new/images/multi_snap_'+sim+'_'+str(img_size)+'pc_%03i.pdf' %snaps)
end_total = time.time()
print 'Took: %f seconds' %(end_total - start_total) 

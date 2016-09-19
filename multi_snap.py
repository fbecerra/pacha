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

fields = ['nh']#, 'temp']#, 'abH2'] 
#fields = ['temp']
#snaps = range(36, 126, 8) # Fragmentation before binary
#snaps = range(146,178,2) # Protostar follow-up after binary
snaps = np.append(np.arange(36,126,8), np.arange(132,154,7))
pc = pr.plot_collection.PlotCollection()
pc.multi_snap(snaps, fields, '/scratch/00025/tgreif/', 'ah1w3f2', 20, PlotTime=True, PlotSize=True, CenterProto=1, PlotProto=False, Rotate=True)
pc.savefig('/scratch/02563/fbecerra/paha/images/multi_snap_central_cluster_'+'_'.join(fields)+'.pdf')
end_total = time.time()
print 'Took: %f seconds' %(end_total - start_total) 

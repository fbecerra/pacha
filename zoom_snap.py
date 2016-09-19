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

#:fields = ['nh', 'temp']#, 'abH2']
#snap = 25
#sizes = [1000, 100, 10]
#pc = pr.plot_collection.PlotCollection()
#pc.zoom_snap(snap, fields, sizes, '/scratch/00025/tgreif/', 'ah1w3f2')
#pc.savefig('/scratch/02563/fbecerra/paha/images/zoom_snap_'+'_'.join(fields)+'.pdf')
#end_total = time.time()
#print 'Took: %f seconds' %(end_total - start_total) 

field = 'nh'
sims = ['ah1w3', 'ah1w3r0', 'ah1w3r0', 'ah1w3f2', 'ah1w3f2', 'ah1w3f2']
snaps = [157,22,22,25,25,25]
sizes = [2e6, 2e5, 2e4, 1000, 100, 10]
pc = pr.plot_collection.PlotCollection()
pc.zoom_snap2(snaps, field, sizes, '/scratch/00025/tgreif/', sims)
pc.savefig('/scratch/02563/fbecerra/paha/images/zoom_snap2_'+field+'.pdf')
end_total = time.time()
print 'Took: %f seconds' %(end_total - start_total) 

import matplotlib
matplotlib.use('TkAgg')
import src as pr
import numpy as np
import pylab as pl
import time


path = '/Users/fbecerra/Dropbox/ACH/snapdir_137/'
file = 'nahw1_137'

snapbase = path+file

start = time.time()
MySnap = pr.snap.Snap()
MySnap.read_header(snapbase)
MySnap.read_fields(snapbase)
MySnap.center_box()
MySnap.rotate_box()
print MySnap.params['boxsize']
print np.min(MySnap.fields['x']), np.max(MySnap.fields['x'])
print np.min(MySnap.fields['y']), np.max(MySnap.fields['y'])
print np.min(MySnap.fields['z']), np.max(MySnap.fields['z'])
end = time.time()
print 'Total time: %f' %(end - start)

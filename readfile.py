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
print np.min(MySnap.fields['nh']), np.max(MySnap.fields['nh'])
print np.min(MySnap.fields['temp']), np.max(MySnap.fields['temp'])
end = time.time()
print 'Total time: %f' %(end - start)
#MySnap.center_box()
#MySnap.rotate_box()

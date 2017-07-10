import matplotlib
matplotlib.use('TkAgg')
import src as pr
import numpy as np
import pylab as pl
import time


#path = '/Users/fbecerra/Dropbox/ACH/snapdir_137/'
#file = 'nahw1_137'

path = '/Users/fbecerra/Dropbox/ACH/snapdir_138/'
file = 'nahw1_138'

snapbase = path+file

start = time.time()
MySnap = pr.snap.Snap()
MySnap.read_header(snapbase)
#MySnap.read_fields(snapbase)
MySnap.read_fields(snapbase, read_dm=True)
MySnap.center_box()
MySnap.rotate_box()

#print MySnap.params['boxsize']
#print np.min(MySnap.fields['x']), np.max(MySnap.fields['x'])
#print np.min(MySnap.fields['y']), np.max(MySnap.fields['y'])
#print np.min(MySnap.fields['z']), np.max(MySnap.fields['z'])

print MySnap.params['masstable'], MySnap.params['npartall']
print np.min(MySnap.new_fields['dm']['x']), np.max(MySnap.new_fields['dm']['x'])
print np.min(MySnap.new_fields['dm']['y']), np.max(MySnap.new_fields['dm']['y'])
print np.min(MySnap.new_fields['dm']['z']), np.max(MySnap.new_fields['dm']['z'])
print np.min(MySnap.new_fields['dm']['mass']), np.max(MySnap.new_fields['dm']['mass'])
print np.sum(MySnap.new_fields['dm']['mass']), np.sum(MySnap.fields['mass'])
print np.sum(MySnap.new_fields['sinks']['mass']), MySnap.new_fields['sinks']['id']

end = time.time()
print 'Total time: %f' %(end - start)

import matplotlib
matplotlib.use('Agg')
import src as pr
import numpy as np
import pylab as pl
import time

sims = [#'ah1w3'],
        'ah1w3r0']#, 
        #'ah1w3f2']
snaps = [#'035'],
         '001']#,
         #'065']

nsims = len(sims)

base= '/scratch/00025/tgreif/'

start = time.time()

for sim_idx in range(nsims):
  sim = sims[sim_idx]
  snap = snaps[sim_idx]

  path = base+sim+'/snapdir_'+snap+'/'
  file = sim+'_'+snap
  snapbase = path+file
  
  MySnap = pr.snap.Snap()
  MySnap.read_header(snapbase)
  MySnap.read_fields(snapbase)
  MySnap.center_box()
  MySnap.rotate_box()
  MySnap.calculate_radius()
  print MySnap.params['time']*pr.constants.UNIT_TIME/pr.constants.SEC_PER_YEAR
  print MySnap.params['redshift']
#  print np.sum(MySnap.params['npartall'])
#  print np.sort(MySnap.fields['mass'][MySnap.fields['mass'] != 0])[0:4]
#  print np.min(MySnap.fields['hsml'])

end = time.time()
print 'Total time: %f' %(end - start)

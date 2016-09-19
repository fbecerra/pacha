import matplotlib
matplotlib.use('Agg')
import src as pr
import numpy as np
import pylab as pl
import time

base= '/scratch/00025/tgreif/'
sim = 'ah1w3f2'
snaps = range(1,182,30)
fields = ['mbe_ratio']
fig = pl.figure()
for snap in snaps:
  print snap
  field = fields[0]
  path = base+sim+'/snapdir_%03d/' %snap
  file = sim+'_%03d' %snap
  snapbase = path+file

  MySnap = pr.snap.Snap()
  MySnap.read_header(snapbase)
  MySnap.read_fields(snapbase)
  MySnap.center_box()
  MySnap.rotate_box()
  MySnap.calculate_radius()

  MyRadial = pr.radial.Radial()
  MyRadial.radial_profile(MySnap, fields)
  pl.plot(MyRadial.radial['enc_mass'], MyRadial.radial[field], label='%.2f yr' %(MySnap.params['time']*pr.constants.UNIT_TIME/pr.constants.SEC_PER_YEAR))

xmin, xmax = np.min(MyRadial.radial['enc_mass']), np.max(MyRadial.radial['enc_mass'])

# y-axis
pl.ylabel(pr.utils.get_label(field))
locs, labels = pl.yticks()
locs = np.array(locs)
idx = np.where(locs % 1 == 0)[0]
pl.yticks(locs[idx])

# x-axis
pl.xlabel(pr.utils.get_label('enc_mass'))
locs = np.arange(xmin, xmax, dtype=int)
if len(locs) > 6:
  locs = locs[np.where(locs % 2 == 0)[0]]
pl.xticks(locs)
pl.xlim(xmin, xmax)

pl.plot(np.array([xmin, xmax]), np.array([0,0]), 'k--')
pl.legend(loc='lower right', prop={'size':14})
fig.savefig('/scratch/02563/fbecerra/paha/radial/radial_timelapse_'+'_'.join(fields)+'.pdf')
pl.show()
pl.close()

import matplotlib
matplotlib.use('Agg')
import src as pr
import numpy as np
import pylab as pl
from scipy import integrate
import time

base= '/n/hernquistfs2/fbecerra/'

### Sinks
sim, snap_init, snap_final = 'nahw1r4sm1', 1, 20
#sim, snap_init, snap_final = 'nahw1r4sm2', 9, 16
#sim, snap_init, snap_final = 'nahw1r4sm3', 14, 26

### Adiabatic
##sim, snap_init, snap_final, xmin, xmax = 'nahw1r4ad1', 1, 31, -7, -3
##sim, snap_init, snap_final, xmin, xmax = 'nahw1r4ad2', 7, 23, -8, -3
##sim, snap_init, snap_final, xmin, xmax = 'nahw1r4ad3',

snaps = range(snap_init, snap_final+1)
rif = np.zeros(len(snaps))
masses = np.zeros(len(snaps))
field = 'Nion'

fig = pl.figure()

for idx, snap in enumerate(snaps):
  print snap
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
  MyRadial.radial_profile(MySnap, [field])

  # Version 1
  dnion = MyRadial.radial[field]
  if np.isnan(dnion).any():
    nans = np.argwhere(np.isnan(dnion))[-1]
    x = MyRadial.radial['radius'][nans+1:]
    y = dnion[nans+1:]
  else:
    x = MyRadial.radial['radius']
    y = dnion
  Nion1 = integrate.cumtrapz(y, x=10**x * pr.constants.UNIT_LENGTH * 1e3)

  # Version 2
  idx_sink = np.argmax(MySnap.new_fields['sinks']['mass'])
  sink_mass = MySnap.new_fields['sinks']['mass'][idx_sink]
  Nion2 = 2.37e54 * (sink_mass / 1e6)

  f = Nion1 - Nion2
  idx_rif = np.where(f > 0)[0][0] + 1

  rif[idx] = x[idx_rif]
  masses[idx] = sink_mass

  # Plot
  #pl.plot(x[1:], np.log10(Nion1), label='%.2f yr' %(MySnap.params['time']*pr.constants.UNIT_TIME/pr.constants.SEC_PER_YEAR))
  #pl.plot(x[1:], np.log10(Nion1), label='%.2f Msun' %sink_mass)
  #xmin, xmax = np.min(MyRadial.radial['radius']), np.max(MyRadial.radial['radius'])
  #pl.plot(np.array([xmin, xmax]), np.log10(np.array([Nion2, Nion2])), 'k--')

pl.plot(np.log10(masses), rif, label=sim)
print masses
print rif

# y-axis
pl.ylabel('log (rIF/pc)')
#locs, labels = pl.yticks()
#locs = np.array(locs)
#idx = np.where(locs % 1 == 0)[0]
#pl.yticks(locs[idx])

# x-axis
pl.xlabel('log (Mstar/Msun)')
#locs = np.arange(xmin, xmax, dtype=int)
#if len(locs) > 6:
#  locs = locs[np.where(locs % 2 == 0)[0]]
#pl.xticks(locs)
#pl.xlim(xmin, xmax)

pl.legend(loc='lower right', prop={'size':14})
fig.savefig('/n/home02/fbecerra/pacha_new/radial/radial_rif_' + sim + '.pdf')
pl.show()
pl.close()

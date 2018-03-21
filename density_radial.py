import matplotlib
matplotlib.use('Agg')
import src as pr
import numpy as np
import pylab as pl
import time
from itertools import cycle

def n_isothermal(temp, radius):
  return 2.3e3 * (temp / 300) * (radius)**(-2.)

def sound_speed(mu, temp):
    return np.sqrt(pr.constants.BOLTZMANN * temp / mu / pr.constants.PROTONMASS)

fontsize = 20
lines = cycle(['-', '--', '-.', ':'])
colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf']
gray = '#999999'

base= '/n/hernquistfs2/fbecerra/'

sims = ['nahw1r4sm1', 'nahw1r4sm2', 'nahw1r4sm3', 'nahw1']
snaps = [20, 16, 26, 137]

field = 'nh'
fig = pl.figure(figsize=(7,6))
for idx, snap in enumerate(snaps):
  sim = sims[idx]
  print sim, snap
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
  radius = MyRadial.radial['radius']

  pl.plot(MyRadial.radial['radius'], MyRadial.radial[field], label=sim)
  print MyRadial.radial['radius']
  print MyRadial.radial[field]

xmin, xmax = -4, 3
#radius = np.array([xmin, xmax])
# Plot isothermal fit
Tvir = 34515 
pl.plot(radius, np.log10(n_isothermal(Tvir, 10**radius)), '--', color=colors[idx+1], label='Tvir')
T = 8000
pl.plot(radius, np.log10(n_isothermal(T, 10**radius)), '--', color=colors[idx+2], label='T=8000K')

# y-axis
ymin, ymax = -1, 12
pl.ylim(ymin, ymax)
pl.ylabel(pr.utils.get_label(field))
locs, labels = pl.yticks()
locs = np.array(locs)
idx = np.where(locs % 1 == 0)[0]
pl.yticks(locs[idx])

# x-axis
pl.xlabel(pr.utils.get_label('radius'))
locs = np.arange(xmin, xmax, dtype=int)
if len(locs) > 6:
  locs = locs[np.where(locs % 2 == 0)[0]]
pl.xticks(locs)
pl.xlim(xmin, xmax)

pl.legend(loc='lower right', prop={'size':14})
pl.tight_layout()
fig.savefig('/n/home02/fbecerra/pacha_new/radial/density_all.pdf')
pl.show()
pl.close()

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
lines = cycle(['-', '--'])
#colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf']
colors = ['#08519c', '#9ecae1', '#cb181d', '#fc9272']
gray = '#999999'
time0 = 1.27917e-6

base= '/n/hernquistfs2/fbecerra/'

sim = 'nahw1r4rt1'
snaps = [4, 8, 12, 15]

field = 'phodens'
fig = pl.figure(figsize=(7,6))
for idx, snap in enumerate(snaps):

  print  snap
  path = base+sim+'/advection/snapdir_%03i/' %snap
  file = sim+'_%03d' %snap
  snapbase = path+file

  MySnap = pr.snap.Snap()
  MySnap.read_header(snapbase)
  MySnap.read_fields(snapbase)
  MySnap.center_box()
  MySnap.rotate_box()
  MySnap.calculate_radius()

  time = (MySnap.params['time'] - time0) * pr.constants.UNIT_TIME / pr.constants.SEC_PER_YEAR

  MyRadial = pr.radial.Radial()
  MyRadial.radial_profile(MySnap, [field])

  linestyle = next(lines)
  color = colors[idx]
  pl.plot(MyRadial.radial['radius'], MyRadial.radial[field], linestyle=linestyle, label='%i yr' %time, color=color)

xmin, xmax = -1, 1.2
radius = np.array([xmin, xmax])
# Plot isothermal fit
pl.plot(np.array([xmin, xmax]), -2*np.array([xmin,xmax]) + 9, linestyle=':', color=gray)
pl.text(0.4, 8.3, r'$n_{\gamma} \propto r^{-2}$', fontsize=18)

# y-axis
ymin, ymax = 3, 10
pl.ylim(ymin, ymax)
pl.ylabel(pr.utils.get_label(field), fontsize=20)
locs, labels = pl.yticks()
locs = np.array(locs)
idx = np.where(locs % 1 == 0)[0]
pl.yticks(locs[idx])

# x-axis
pl.xlabel(pr.utils.get_label('radius'), fontsize=20)
locs = np.arange(xmin, xmax, dtype=int)
if len(locs) > 6:
  locs = locs[np.where(locs % 2 == 0)[0]]
pl.xticks(locs)
pl.xlim(xmin, xmax)

pl.legend(loc='lower left', prop={'size':14})
pl.tight_layout()
fig.savefig('/n/home02/fbecerra/pacha_new/radial/photon_advection_radial.pdf')
pl.show()
pl.close()

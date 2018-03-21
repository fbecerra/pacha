import matplotlib
matplotlib.use('Agg')
import src as pr
import numpy as np
import pylab as pl
import time
from itertools import cycle

fontsize = 20
lines = cycle(['-', '--', '-.', ':'])

base= '/n/hernquistfs2/fbecerra/'

### Sinks
#sim, snap_init, snap_final, img_size = 'nahw1r4sm1', 1, 20, 3
#sim, snap_init, snap_final, img_size = 'nahw1r4sm2', 7, 16, 0.3
#sim, snap_init, snap_final, img_size = 'nahw1r4sm3', 14, 26, 0.03

### Adiabatic
#sim, snap_init, snap_final, img_size = 'nahw1r4ad1', 1, 31, 3
#sim, snap_init, snap_final, img_size = 'nahw1r4ad2', 7, 23, 0.3
#sim, snap_init, snap_final, img_size = 'nahw1r4ad3',

snaps = range(snap_init, snap_final+1)
fields = ['tff', 'tgrav', 'tpres']
fig = pl.figure()

fig, ax = pl.subplots(1, 2, figsize=(16, 8))

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
  MyRadial.radial_profile(MySnap, fields)

  linestyle = next(lines)
  ax[0].plot(MyRadial.radial['radius'], np.log10(np.abs(MyRadial.radial['tgrav'])/MyRadial.radial['tff']), linestyle=linestyle, 
             label='%.2f yr' %(MySnap.params['time']*pr.constants.UNIT_TIME/pr.constants.SEC_PER_YEAR))
  ax[1].plot(MyRadial.radial['radius'], np.log10(np.abs(MyRadial.radial['tpres'])/MyRadial.radial['tff']), linestyle=linestyle,
             label='%.2f yr' %(MySnap.params['time']*pr.constants.UNIT_TIME/pr.constants.SEC_PER_YEAR))

#xmin, xmax = np.min(MyRadial.radial['enc_mass']), np.max(MyRadial.radial['enc_mass'])

ax[1].legend(loc='lower right', prop={'size':14})

# y-axis
ax[0].set_ylabel(r'$\log \left(t_{\rm grav}/t_{\rm ff}\right)$', fontsize=fontsize)
ax[1].set_ylabel(r'$\log \left(t_{\rm pres}/t_{\rm ff}\right)$', fontsize=fontsize)

# x-axis
ax[0].set_xlabel(r'$\log(r/{\rm kpc})$', fontsize=fontsize)
ax[1].set_xlabel(r'$\log(r/{\rm kpc})$', fontsize=fontsize)

xmin, xmax = -6, -3
ax[0].set_xlim(xmin, xmax)
ax[1].set_xlim(xmin, xmax)

ax[0].set_ylim(-2, 5)
ax[1].set_ylim(-2, 5)

ax[0].plot(np.array([xmin, xmax]), np.array([0,0]), 'k--')
ax[1].plot(np.array([xmin, xmax]), np.array([0,0]), 'k--')

fig.savefig('./radial/radial_timelapse_'+sim+'_'.join(fields)+'.pdf')
pl.show()
pl.close()

import matplotlib
matplotlib.use('Agg')
import src as pr
import numpy as np
import pylab as pl
import time

base= '/n/hernquistfs2/fbecerra/'
sim = 'nahw1r3tgs2'
snaps = range(11,19)
fields = ['taugrav', 'taupres', 'tgrav', 'tpres']
fig = pl.figure()

fig, ax = pl.subplots(1, 2, figsize=(12,5))

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
  ax[0].plot(MyRadial.radial['radius'], np.log10(MyRadial.radial['taugrav']/MyRadial.radial['taupres']), 
             label='%.2f yr' %(MySnap.params['time']*pr.constants.UNIT_TIME/pr.constants.SEC_PER_YEAR))
  ax[1].plot(MyRadial.radial['radius'], MyRadial.radial['tgrav']/MyRadial.radial['tpres'], 
             label='%.2f yr' %(MySnap.params['time']*pr.constants.UNIT_TIME/pr.constants.SEC_PER_YEAR))

#xmin, xmax = np.min(MyRadial.radial['enc_mass']), np.max(MyRadial.radial['enc_mass'])

ax[0].legend(loc='lower right', prop={'size':14})

# y-axis
ax[0].set_ylabel(r'$\log (\tau_{\rm grav}/\tau_{\rm pres})$')
ax[1].set_ylabel(r'$t_{\rm grav}/t_{\rm pres}$')

# x-axis
ax[0].set_xlabel(r'$\log(r/{\rm pc})$')
ax[1].set_xlabel(r'$\log(r/{\rm pc})$')

xmin, xmax = -6, -3
ax[0].set_xlim(xmin, xmax)
ax[1].set_xlim(xmin, xmax)

ax[0].set_ylim(-2, 2)
ax[1].set_ylim(-50, 50)

ax[0].plot(np.array([xmin, xmax]), np.array([0,0]), 'k--')
ax[1].plot(np.array([xmin, xmax]), np.array([1,1]), 'k--')

fig.savefig('./radial/radial_timelapse_'+sim+'_'.join(fields)+'.pdf')
pl.show()
pl.close()

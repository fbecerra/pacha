import matplotlib
matplotlib.use('Agg')
import src as pr
import numpy as np
import pylab as pl
import time

base= '/n/hernquistfs2/fbecerra/'
sim = 'nahw1r3tgs2'
#snaps = range(4, 11)
snaps = range(11, 19)
fields = ['infall'] #, 'cs', 'mass']
fig = pl.figure(figsize=(7,6))
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
  ##cs = 10**MyRadial.radial['cs']
  ##cs3G = (cs * 1e5)**3. / pr.constants.GRAVITY / pr.constants.SOLAR_MASS * pr.constants.SEC_PER_YEAR
  ##shu_rate = 0.975 * cs3G
  ##LP_rate = 46.9 * cs3G
  ##print shu_rate
  ##print LP_rate

  # Previous snapshot
  path = base+sim+'/snapdir_%03d/' %(snap-1)
  file = sim+'_%03d' %(snap-1)
  snapbase = path+file

#  MySnapPrev = pr.snap.Snap()
#  MySnapPrev.read_header(snapbase)
#  MySnapPrev.read_fields(snapbase)
#  MySnapPrev.center_box()
#  MySnapPrev.rotate_box()
#  MySnapPrev.calculate_radius()
#
#  MyRadialPrev = pr.radial.Radial()
#  MyRadialPrev.radial_profile(MySnapPrev, ['mass'])
#
#  diff_rate = - (10**MyRadial.radial['mass'] - 10**MyRadialPrev.radial['mass'])  / ((MySnap.params['time'] - MySnapPrev.params['time']) * pr.constants.UNIT_TIME / pr.constants.SEC_PER_YEAR)

  pl.plot(MyRadial.radial['radius'], np.log10(MyRadial.radial[field]), label='%.2f yr' %(MySnap.params['time']*pr.constants.UNIT_TIME/pr.constants.SEC_PER_YEAR)) #, label='%.2f yr' %(MySnap.params['time']*pr.constants.UNIT_TIME/pr.constants.SEC_PER_YEAR)) #  label=r'$\dot{M}\simeq -4\pi r^2 \rho v_{\rm rad}$'
  ##pl.plot(MyRadial.radial['radius'], np.log10(shu_rate), color='#377eb8', linestyle=':', label=r'$\dot{M}_{\rm Shu}\simeq 0.975c_{\rm s}^3/G$')
  ##pl.plot(MyRadial.radial['radius'], np.log10(LP_rate), color='#4daf4a', linestyle=':', label=r'$\dot{M}_{\rm LP}\simeq 46.9c_{\rm s}^3/G$')
  #pl.plot(MyRadial.radial['radius'], diff_rate, color='#4daf4a', linestyle=':', label=r'$\dot{M} = dM / dt$')

#xmin, xmax = np.min(MyRadial.radial['enc_mass']), np.max(MyRadial.radial['enc_mass'])
#
# y-axis
pl.ylabel(pr.utils.get_label(field), fontsize=20)
#locs, labels = pl.yticks()
#locs = np.array(locs)
#idx = np.where(locs % 1 == 0)[0]
#pl.yticks(locs[idx])
pl.yticks(fontsize=16)
pl.ylim(-5, 2)

# x-axis
pl.xlabel(pr.utils.get_label('radius'), fontsize=20)
#locs = np.arange(xmin, xmax, dtype=int)
#if len(locs) > 6:
#  locs = locs[np.where(locs % 2 == 0)[0]]
#pl.xticks(locs)
pl.xticks(fontsize=16)
#pl.xlim(-1, 2)
#
#pl.plot(np.array([xmin, xmax]), np.array([0,0]), 'k--')
pl.legend(loc=4, prop={'size':12}, frameon=True)
#pl.axvline(x=np.log10(0.18714), linestyle='dashed', color='k')
pl.tight_layout()
fig.savefig('/n/home02/fbecerra/pacha_new/radial/radial_'+sim+'_'+'_'.join(fields)+'.pdf')
pl.show()
pl.close()

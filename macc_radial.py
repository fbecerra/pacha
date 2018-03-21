import matplotlib
matplotlib.use('Agg')
import src as pr
import numpy as np
import pylab as pl
import time
from itertools import cycle

def get_label(sim):
  if sim == 'nahw1r4sm1' or sim == 'nahw1r4sm2' or sim == 'nahw1r4sm3':
    label1 = 'sink'
  elif sim == 'nahw1r4ad1' or sim == 'nahw1r4ad2' or sim == 'nahw1r4ad3':
    label1 = 'core'

  if sim == 'nahw1r4sm3' or sim == 'nahw1r4ad3':
    label2 = r'$n_{\rm th} = 10^{12}\,{\rm cm}^{-3}$'
  elif sim == 'nahw1r4sm2' or sim == 'nahw1r4ad2':
    label2 = r'$n_{\rm th} = 10^{10}\,{\rm cm}^{-3}$'
  elif sim == 'nahw1r4sm1' or sim == 'nahw1r4ad1':
    label2 = r'$n_{\rm th} = 10^8\,{\rm cm}^{-3}$'

  #return label1+', '+label2

  if sim == 'nahw1r4sm1' or sim == 'nahw1r4sm2' or sim == 'nahw1r4sm3':
    return label2
  elif sim == 'nahw1r4sm1' or sim == 'nahw1r4sm2' or sim == 'nahw1r4sm3':
    return ''

def get_color(sim):
  if sim == 'nahw1r4sm3' or sim == 'nahw1r4ad3':
    color = '#4daf4a'
  elif sim == 'nahw1r4sm2' or sim == 'nahw1r4ad2':
    color = '#377eb8'
  elif sim == 'nahw1r4sm1' or sim == 'nahw1r4ad1':
    color = '#e41a1c'
  return color

def get_linestyle(sim):
  if sim == 'nahw1r4sm1' or sim == 'nahw1r4sm2' or sim == 'nahw1r4sm3':
    linestyle = '-'
  elif sim == 'nahw1r4ad1' or sim == 'nahw1r4ad2' or sim == 'nahw1r4ad3':
    linestyle = '--'
  return linestyle

fontsize = 20
lines = cycle(['-', '--', '-.', ':'])

base= '/n/hernquistfs2/fbecerra/'
### Sinks
#sim, snap_init, snap_final, xmin, xmax = 'nahw1r4sm1', 1, 20, -2, 0
#sim, snap_init, snap_final, xmin, xmax = 'nahw1r4sm2', 7, 16, -3, 0
#sim, snap_init, snap_final, xmin, xmax = 'nahw1r4sm3', 14, 26, -4, 0

### Adiabatic
#sim, snap_init, snap_final, xmin, xmax = 'nahw1r4ad1', 1, 31, -4, 0
#sim, snap_init, snap_final, xmin, xmax = 'nahw1r4ad2', 7, 23, -5, 0
#sim, snap_init, snap_final, xmin, xmax = 'nahw1r4ad3', 14, 32, -6, 0

sims = ['nahw1r4sm1', 'nahw1r4sm2', 'nahw1r4sm3', 'nahw1r4ad1', 'nahw1r4ad2', 'nahw1r4ad3']
snaps = [20, 16, 26, 31, 23, 32]

#snaps = range(snap_init, snap_final+1)
fields = ['infall'] #, 'cs', 'mass']
fig = pl.figure(figsize=(7,6))
for idx_sim, sim in enumerate(sims):
  snap = snaps[idx_sim]
  print 'Sim: %s, snap %i' %(sim,snap)

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
  #path = base+sim+'/snapdir_%03d/' %(snap-1)
  #file = sim+'_%03d' %(snap-1)
  #snapbase = path+file

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

  #linestyle = next(lines)
  pl.plot(MyRadial.radial['radius'], np.log10(MyRadial.radial[field]), linestyle=get_linestyle(sim), label=get_label(sim), color=get_color(sim))
  #pl.plot(MyRadial.radial['radius'], np.log10(MyRadial.radial[field]), linestyle=linestyle,
  #         label='%.2f yr' %(MySnap.params['time']*pr.constants.UNIT_TIME/pr.constants.SEC_PER_YEAR)) 
           #, label='%.2f yr' %(MySnap.params['time']*pr.constants.UNIT_TIME/pr.constants.SEC_PER_YEAR)) #  label=r'$\dot{M}\simeq -4\pi r^2 \rho v_{\rm rad}$'
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
pl.ylim(-3, 1.5)

# x-axis
pl.xlabel(pr.utils.get_label('radius'), fontsize=20)
#locs = np.arange(xmin, xmax, dtype=int)
#if len(locs) > 6:
#  locs = locs[np.where(locs % 2 == 0)[0]]
#pl.xticks(locs)
pl.xticks(fontsize=16)
pl.xlim(-5, 0)
#pl.xlim(xmin, xmax)
#
#pl.plot(np.array([xmin, xmax]), np.array([0,0]), 'k--')
pl.legend(loc=4, prop={'size':14})
#pl.axvline(x=np.log10(0.18714), linestyle='dashed', color='k')
pl.tight_layout()
fig.savefig('/n/home02/fbecerra/pacha_new/radial/radial_timelapse_allsims_'+'_'.join(fields)+'.pdf')
pl.show()
pl.close()

import matplotlib
matplotlib.use('Agg')
import src as pr
import numpy as np
import pylab as pl
import time
from itertools import cycle

def get_label(sim, idx):
  if sim == 'nahw1r4sm3' or sim == 'nahw1r4ad3':
    label1 = r'$n_{\rm th} = 10^{12}\,{\rm cm}^{-3}$'
  elif sim == 'nahw1r4sm2' or sim == 'nahw1r4ad2':
    label1 = r'$n_{\rm th} = 10^{10}\,{\rm cm}^{-3}$'
  elif sim == 'nahw1r4sm1' or sim == 'nahw1r4ad1':
    label1 = r'$n_{\rm th} = 10^8\,{\rm cm}^{-3}$'

  if idx == 0:
    label2 = 'initial'
  elif idx == 1:
    label2 = 'final'

  #return label1+', '+label2
  if idx == 0:
    return label1
  else:
    return ''

def get_color(sim):
  if sim == 'nahw1r4sm3' or sim == 'nahw1r4ad3':
    color = '#4daf4a'
  elif sim == 'nahw1r4sm2' or sim == 'nahw1r4ad2':
    color = '#377eb8'
  elif sim == 'nahw1r4sm1' or sim == 'nahw1r4ad1':
    color = '#e41a1c'
  return color

fontsize = 20
lines = cycle(['-', '--', '-.', ':'])
gray = '#999999'
base= '/n/hernquistfs2/fbecerra/'
### Sinks
#sim, snap_init, snap_final, xmin, xmax = 'nahw1r4sm1', 1, 20, 3, 6
#sim, snap_init, snap_final, xmin, xmax = 'nahw1r4sm2', 7, 16, 2, 5
#sim, snap_init, snap_final, xmin, xmax = 'nahw1r4sm3', 14, 26, 1, 5

### Adiabatic
#sim, snap_init, snap_final, xmin, xmax = 'nahw1r4ad1', 1, 31, 1, 6
#sim, snap_init, snap_final, xmin, xmax = 'nahw1r4ad2', 7, 23, 0, 5
#sim, snap_init, snap_final, xmin, xmax = 'nahw1r4ad3', 14, 32

sims = ['nahw1r4ad1', 'nahw1r4ad2', 'nahw1r4ad3']
snaps = [[1, 31], [7, 23], [14, 32]]
#snaps = [[1], [7], [14]]
linestyles = ['-', '--']

#snaps = range(snap_init, snap_final+1)
fields = ['mbe_ratio']
fig = pl.figure(figsize=(7,6))

for idx_sim, sim in enumerate(sims):
  for idx_snap, snap in enumerate(snaps[idx_sim]):
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
  
    linestyle = linestyles[idx_snap]
    pl.plot(MyRadial.radial['enc_mass'], MyRadial.radial[field], linestyle=linestyle, label=get_label(sim, idx_snap), color=get_color(sim))


# y-axis
ymin, ymax = -1.2, 1.2
pl.ylim(ymin, ymax)
pl.ylabel(pr.utils.get_label(field), fontsize=20)
locs, labels = pl.yticks()
locs = np.array(locs)
idx = np.where(locs % 1 == 0)[0]
pl.yticks(locs[idx])

# x-axis
pl.xlabel(pr.utils.get_label('enc_mass'), fontsize=20)
#xmin, xmax = np.min(MyRadial.radial['enc_mass']), np.max(MyRadial.radial['enc_mass'])
xmin, xmax = 0, 5 
locs = np.arange(xmin, xmax, dtype=int)
if len(locs) > 6:
  locs = locs[np.where(locs % 2 == 0)[0]]
pl.xticks(locs)
pl.xlim(xmin, xmax)

pl.plot(np.array([xmin, xmax]), np.array([0,0]), ':', color=gray)
pl.legend(loc='lower right', prop={'size':14})
pl.tight_layout()
fig.savefig('/n/home02/fbecerra/pacha_new/radial/radial_timelapse_allad_'+'_'.join(fields)+'.pdf')
pl.show()
pl.close()

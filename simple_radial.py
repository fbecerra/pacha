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
sim = 'nahw1'
#snaps = range(4, 11)
snaps = range(137, 138)
fields = ['temp'] #, 'cs', 'mass']
#fields = ['nh', 'temp', 'abHM', 'abH2', 'abHII', 'gamma']
fig = pl.figure(figsize=(6.5,6))
for snap in snaps:
  #print snap
  #field = fields[0]
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

  for field in fields:
    print field+': '
    print MyRadial.radial[field]

  print MyRadial.radial['radius']
  pl.plot(MyRadial.radial['radius'], MyRadial.radial[field], linestyle='-', color='red') 

#xmin, xmax = np.min(MyRadial.radial['enc_mass']), np.max(MyRadial.radial['enc_mass'])
xmin, xmax = -2, 5
ymin, ymax = 1, 5

Teff = 6000
Mass = 1e4 * pr.constants.SOLAR_MASS #10^4 solar mass object
Tinf = 1e4 #Tvir
sound_speed = np.sqrt(pr.constants.BOLTZMANN * Tinf / pr.constants.PROTONMASS)
r_a = 2 * pr.constants.GRAVITY * Mass / sound_speed**2.

def temp_var(radius):
    return (Tinf/4.) * (radius/r_a)**(-1)

pl.plot(np.array([xmin, xmax]), np.log10(np.array([Teff, Teff])), '--', color='black')
pl.plot(np.array([xmin, xmax]), np.log10(np.array([temp_var(10**xmin * pr.constants.UNIT_LENGTH / 1e3), temp_var(10**xmax * pr.constants.UNIT_LENGTH / 1e3)])), ':', color='black')

rB5 = 8.68386752702
rB4 = 0.868386752702
rB3 = 0.0868386752702

size = 30
pl.scatter([np.log10(rB5)], [np.log10(Teff)], marker='s', color='k', label=r'$R_{\rm B}(M_\star = 10^5\,{\rm M}_\odot)$', s=size)
pl.scatter([np.log10(rB5)], [np.log10(temp_var(rB5 * pr.constants.UNIT_LENGTH / 1e3))], marker='s', color='k', s=size)

pl.scatter([np.log10(rB4)], [np.log10(Teff)], marker='o', color='k', label=r'$R_{\rm B}(M_\star = 10^4\,{\rm M}_\odot)$', s=size)
pl.scatter([np.log10(rB4)], [np.log10(temp_var(rB4 * pr.constants.UNIT_LENGTH / 1e3))], marker='o', color='k', s=size)

pl.scatter([np.log10(rB3)], [np.log10(Teff)], marker='^', color='k', label=r'$R_{\rm B}(M_\star = 10^3\,{\rm M}_\odot)$', s=size)
pl.scatter([np.log10(rB3)], [np.log10(temp_var(rB3 * pr.constants.UNIT_LENGTH / 1e3))], marker='^', color='k', s=size)

pl.text(2.5, 1.3, r'$\gamma = 5/3$', fontsize=16)
pl.text(3.3, 3.9, r'$T_\star \simeq 6000\,{\rm K}$', fontsize=16)

# y-axis
pl.ylabel(pr.utils.get_label(field), fontsize=20)
#locs, labels = pl.yticks()
#locs = np.array(locs)
#idx = np.where(locs % 1 == 0)[0]
#pl.yticks(locs[idx])
pl.yticks(fontsize=16)
pl.ylim(ymin, ymax)

# x-axis
pl.xlabel(pr.utils.get_label('radius'), fontsize=20)
#locs = np.arange(xmin, xmax, dtype=int)
#if len(locs) > 6:
#  locs = locs[np.where(locs % 2 == 0)[0]]
#pl.xticks(locs)
pl.xticks(fontsize=16)
pl.xlim(xmin, xmax)
#
#pl.plot(np.array([xmin, xmax]), np.array([0,0]), 'k--')
pl.legend(loc='lower left', numpoints=1, prop={'size':14})
#pl.axvline(x=np.log10(0.18714), linestyle='dashed', color='k')
pl.tight_layout()
fig.savefig('/n/home02/fbecerra/pacha_new/radial/radial_'+sim+'_'+'_'.join(fields)+'.pdf')
pl.show()
pl.close()

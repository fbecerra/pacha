import matplotlib
matplotlib.use('Agg')
import src as pr
import numpy as np
import pylab as pl
from scipy import integrate
import time

base= '/n/hernquistfs2/fbecerra/'

### Sinks
#sim, snap_init, snap_final, n_th = 'nahw1r4sm1', 1, 20 
#sim, snap_init, snap_final, n_th = 'nahw1r4sm2', 7, 16
#sim, snap_init, snap_final, n_th = 'nahw1r4sm3', 14, 26

### Adiabatic
#sim, snap_init, snap_final, n_th = 'nahw1r4ad1', 1, 31, 8.
#sim, snap_init, snap_final, n_th = 'nahw1r4ad2', 7, 23, 10.
sim, snap_init, snap_final, n_th = 'nahw1r4ad3', 14, 32, 12.

snaps = range(snap_init, snap_final+1)
fields = ['nh', 'enc_mass']

times = np.array([])
core_radius = np.array([])
core_mass = np.array([])

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

  f = n_th - MyRadial.radial['nh']
  idx_core = np.where(f > 0)[0][0] + 1
  core_radius = np.append(core_radius, MyRadial.radial['radius'][idx_core])
  core_mass = np.append(core_mass, MyRadial.radial['enc_mass'][idx_core])
  times = np.append(times, MySnap.params['time'] * pr.constants.UNIT_TIME / pr.constants.SEC_PER_YEAR)

print times
print core_radius
print core_mass

import matplotlib
matplotlib.use('Agg')
import src as pr
import numpy as np
import pylab as pl
from scipy import integrate
import time

base= '/n/hernquistfs2/fbecerra/'

sim, snap = 'nahw1', 137
Rvir = 1400
fields = ['angmom', 'rho', 'tff']

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

radius = 10**MyRadial.radial['radius']

f = radius - Rvir
idx_rvir = np.where(f > 0)[0][0]

angmom = 10**MyRadial.radial['angmom'][idx_rvir]
dens = 10**MyRadial.radial['rho'][idx_rvir]
tff = 10**MyRadial.radial['tff'][idx_rvir]

print 'Angular momentum: ', angmom
print 'Volume density: ', dens
print 'Free-fall time: ', tff

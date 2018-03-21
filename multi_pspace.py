import matplotlib
matplotlib.use('TkAgg')
import src as pr
import numpy as np
import pylab as pl
import time
from scipy.signal import argrelmax


path = '/n/hernquistfs2/fbecerra/'
base = 'nahw1r4ad2'

plot_fields = ['temp']#, 'abH2'] #, 'tff', 'tcool']
#plot_fields = ['nh', 'temp', 'kappa']
#plot_fields = ['nh', 'enc_mass', 'temp', 'abH2', 'vrad', 'vrot', 'vturb', 'esc_frac', 'tcool', 'geff', 'tff', 'dnh', 'sigma_gas', 'cs', 'omega_rot', 'q', 'entro']
#plot_fields = ['vrot', 'vrad', 'vkep']

for snap in np.arange(7, 12):
  snapbase = path + base + '/snapdir_%03i/' %snap + base + '_%03i' %snap
  
  print snap
  MySnap = pr.snap.Snap()
  MySnap.read_header(snapbase)
  MySnap.read_fields(snapbase)
  MySnap.center_box()
  MySnap.rotate_box()
  
  pc = pr.plot_collection.PlotCollection()
  pc.multi_pspace(MySnap, plot_fields)
  pc.savefig('/n/home02/fbecerra/pacha_new/images/multi_pspace_'+base+'_'+'_'.join(plot_fields)+'_%03i.pdf' %snap)

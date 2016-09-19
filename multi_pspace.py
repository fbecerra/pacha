import matplotlib
matplotlib.use('TkAgg')
import src as pr
import numpy as np
import pylab as pl
import time
from scipy.signal import argrelmax

path = '/Volumes/DARK-MATTER/arepo/ah1w3f0/snapdir_238/'
file = 'ah1w3f0_238'
snapbase = path+file

MySnap = pr.snap.Snap()
MySnap.read_header(snapbase)
MySnap.read_fields(snapbase)
MySnap.center_box()
MySnap.rotate_box()

plot_fields = ['temp', 'abH2'] #, 'tff', 'tcool']
#plot_fields = ['nh', 'temp', 'kappa']
#plot_fields = ['nh', 'enc_mass', 'temp', 'abH2', 'vrad', 'vrot', 'vturb', 'esc_frac', 'tcool', 'geff', 'tff', 'dnh', 'sigma_gas', 'cs', 'omega_rot', 'q', 'entro']
#plot_fields = ['vrot', 'vrad', 'vkep']

pc = pr.plot_collection.PlotCollection()
pc.multi_pspace(MySnap, plot_fields)
pc.savefig('/scratch/02563/fbecerra/paha/pspace/mult_pspace_'+file+'_'+'_'.join(plot_fields)+'.pdf')

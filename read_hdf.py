import matplotlib
matplotlib.use('TkAgg')
from matplotlib import cm
import src as pr
import numpy as np
import pylab as pl
import h5py
import time

path = '/path/to/output/'
file = 'lcdm256'
initial_snap = 15
final_snap = 15
n_snaps = final_snap - initial_snap + 1

for num_fof in range(initial_snap, final_snap + 1):

  groupbase = path + '/' + file + '/groups_%03i/fof_tab_%03i' %(num_fof, num_fof)

  MyGroup = pr.fof.Fof()
  MyGroup.read_groups(groupbase, read_HDF=True)
  print 'Searching at redshift z=', MyGroup.params['redshift']

  halo_masses = MyGroup.fields['mass'] / MyGroup.params['hubbleparam'] * pr.constants.UNIT_MASS / pr.constants.SOLAR_MASS
  min_mass = 5
  max_mass = 13
  nbins = 300

  y, bin_edges = np.histogram(halo_masses, bins=np.logspace(min_mass, max_mass, nbins))
  bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
  pl.plot(bin_centers, y, '-', label='z='+str(MyGroup.params['redshift']))

pl.gca().set_xscale('log')
pl.legend(loc='upper right')
pl.savefig('./images/hmf.pdf', dpi=100)
pl.show()
pl.close()

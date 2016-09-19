import matplotlib
matplotlib.use('TkAgg')
import src as pr
import numpy as np
import pylab as pl
import time
from scipy.signal import argrelmax
import h5py

def find_max(array):
  idx = np.argmax(array)
  value = array[idx]
  return idx, value

#def calculate_radius(x, y, z, center):
#  return 

snaps = [238, 241]
f = h5py.File('protostars_ah1w3f0.hdf5', 'w')

for snap in snaps:

  snapshot = '%03i' %snap
  path = '/Volumes/DARK-MATTER/arepo/ah1w3f0/snapdir_'+snapshot+'/'
  file = 'ah1w3f0_'+snapshot
  snapbase = path+file
  
  #start = time.time()
  MySnap = pr.snap.Snap()
  MySnap.read_header(snapbase)
  MySnap.read_fields(snapbase)
  #end = time.time()
  #print 'Total time: %f' %(end - start)
  #MySnap.center_box()
  #MySnap.rotate_box()
  #MySnap.calculate_radius()
  
  plot_fields = ['kappa']
  #plot_fields = ['sigma_gas', 'geff', 'dnh']
  #plot_fields = ['vrot', 'vrad', 'vkep']
   
  MyRadial = pr.radial.Radial()
  
  count = 1
  max_dens_idx, max_dens_val = find_max(MySnap.fields['nh'])
  
  while max_dens_val > 1e19:
    center = np.array([MySnap.fields['x'][max_dens_idx], MySnap.fields['y'][max_dens_idx], MySnap.fields['z'][max_dens_idx]])
    MySnap.fields['radius'] = np.sqrt((MySnap.fields['x'] - center[0])**2 + (MySnap.fields['y'] - center[1])**2 + (MySnap.fields['z'] - center[2])**2)
    MyRadial.radial_profile(MySnap, ['kappa'])
    diff = 10**MyRadial.radial['kappa'][2:] - 10**MyRadial.radial['kappa'][:-2]
    rad_idx =  min(np.argsort(diff)[:3])
    protostar_radius = 10**MyRadial.radial['radius'][rad_idx]
    protostar_idx = np.where(MySnap.fields['radius'] < protostar_radius)[0]
    f.create_dataset(file+'/Proto%i/Position' %count, data = center)
    f.create_dataset(file+'/Proto%i/Radius' %count, data = protostar_radius)
    #f.create_dataset(file+'/Proto%i/CM' %count, data = pr.utils.center_of_mass(MySnap, protostar_idx))
    #f.create_dataset(file+'/Proto%i/VCM' %count, data = pr.utils.vel_center_of_mass(MySnap, protostar_idx))
    f.create_dataset(file+'/Proto%i/Mass' %count, data = pr.utils.total_mass(MySnap, protostar_idx))
    f.create_dataset(file+'/Proto%i/IDs' %count, data = MySnap.fields['id'][protostar_idx])
    print 'Protostar ', count
    print 'Position ', center
    print 'Radius: ', protostar_radius
    #print 'Center of Mass: ', pr.utils.center_of_mass(MySnap, protostar_idx)
    #print 'Velocity center of mass: ', pr.utils.vel_center_of_mass(MySnap, protostar_idx)
    print 'Total mass: ', pr.utils.total_mass(MySnap, protostar_idx)
    MySnap.remove_particles(protostar_idx)
  #  MySnap.calculate_radius()
    max_dens_idx, max_dens_val = find_max(MySnap.fields['nh'])
    pl.plot(MyRadial.radial['radius'], MyRadial.radial['kappa'])
    ymin, ymax = pl.ylim()
    pl.plot(np.array([np.log10(protostar_radius), np.log10(protostar_radius)]), np.array([ymin, ymax]), '--')
    pl.savefig('kappa_'+file+'_'+str(count)+'.png')
    pl.show()
    pl.close()
    count += 1
f.close()

import matplotlib
matplotlib.use('TkAgg')
import src as pr
import numpy as np
import pylab as pl
import time

#path = '/scratch/00025/tgreif/ah1dm/snapdir_005/'
#file = 'ah1dm_005'

path = '/scratch/00025/tgreif/ah1w3/snapdir_001/'
file = 'ah1w3_001'

snapbase = path+file

start = time.time()
MySnap = pr.snap.Snap()
MySnap.read_header(snapbase)
#MySnap.read_fields(snapbase)
print 'Number of particles: ', MySnap.params['npartall'], np.sum(MySnap.params['npartall'])
#print 'Mass resolution:', np.min(MySnap.fields['mass']), np.max(MySnap.fields['mass'])
print MySnap.params['masstable']*pr.constants.UNIT_MASS/pr.constants.SOLAR_MASS/MySnap.params['hubbleparam']
end = time.time()
print 'Total time: %f' %(end - start)
#MySnap.center_box()
#MySnap.rotate_box()
#
#plot_fields = ['nh', 'temp', 'kappa']
##plot_fields = ['sigma_gas', 'geff', 'dnh']
##plot_fields = ['vrot', 'vrad', 'vkep']
# 
#MyRadial = pr.radial.Radial()
#count = 1
#
#while max(MySnap.fields['nh']) > 1e19:
#  MyRadial.radial_profile(MySnap, plot_fields)
#  for field_plot in plot_fields:
#    pl.plot(MyRadial.radial['radius'], MyRadial.radial[field_plot])
#    pl.savefig(field_plot+'_'+file+'_'+str(count)+'.png')
#    pl.show()
#    pl.close()
#    if field_plot == 'kappa':
#      diff = 10**MyRadial.radial[field_plot][2:] - 10**MyRadial.radial[field_plot][:-2]
#      protostar_radius = 10**MyRadial.radial['radius'][np.argmin(diff)]
#      protostar_idx = np.where(MySnap.fields['radius'] < protostar_radius)[0]
#      print 'Protostar ', count
#      print 'Position ', MySnap.Center
#      print 'Radius: ', protostar_radius
#      print 'Center of Mass: ', pr.utils.center_of_mass(MySnap, protostar_idx)
#      print 'Velocity center of mass: ', pr.utils.vel_center_of_mass(MySnap, protostar_idx)
#      print 'Total mass: ', pr.utils.total_mass(MySnap, protostar_idx)
#      MySnap.remove_particles(protostar_idx)
#      MySnap.center_box()
#      MySnap.rotate_box()
#      MySnap.calculate_radius()
#      count += 1

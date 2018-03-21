import matplotlib
matplotlib.use('TkAgg')
import src as pr
import numpy as np
import pylab as pl
import time

path = '/n/hernquistfs2/fbecerra/'
#base = 'nahw1r4sm3'
#base = 'sink_test'
base = 'nahw1r4ad3'

start = time.time()

def print_minmax(array):
  print np.min(array), np.max(array)

#f = open('./outputs/sink_mass_'+base+'.txt', 'w')

for snap in np.arange(32, 33):
  print 'Snapshot: ', snap
  snapbase = path + base + '/snapdir_%03i/' %snap + base + '_%03i' %snap
#  snapbase = '/n/home00/fmarinacci/mvogelsfs1/SINKS/output5/snapdir_%03i/' %snap + base + '_%03i' %snap
  MySnap = pr.snap.Snap()
  MySnap.read_header(snapbase)
  MySnap.read_fields(snapbase)
  MySnap.calculate_radius()

  print_minmax(MySnap.fields['radius'])

#  zeroids = np.where(MySnap.fields['id'] == 0)[0]
#  try:
#    MySnap.new_fields['sinks']['id']
#  except:
#    print 'No sinks'
#    continue
#  for idx, id in enumerate(MySnap.new_fields['sinks']['id']):
#    f.write(str(MySnap.params['time']) + ' ' + str(id) + ' ' + str(MySnap.new_fields['sinks']['mass'][idx]) + '\n')

  ####### Jeans number ####
  ####everypart = 50
  ####MySnap.calculate_fields('cs')
  ####MySnap.calculate_fields('tff')
  ####jeans_length = MySnap.derived_fields['cs'] * 1e5 * MySnap.derived_fields['tff']
  ####cell_size = MySnap.fields['hsml'] 
  #####idx = np.where((MySnap.fields['mass'] > 0) & (MySnap.fields['allowref'] == 0))[0]
  ####idx = np.where(MySnap.fields['mass'] > 0)[0]
  ####jeans_number = jeans_length[idx] / cell_size[idx]
  ####density = MySnap.fields['nh'][idx]

  ####try:
  ####  xsink, ysink, zsink = MySnap.new_fields['sinks']['x'][0], MySnap.new_fields['sinks']['y'][0], MySnap.new_fields['sinks']['z'][0]
  ####except:
  ####  max_nh = np.argmax(MySnap.fields['nh'])
  ####  xsink, ysink, zsink = MySnap.fields['x'][max_nh],  MySnap.fields['y'][max_nh],  MySnap.fields['z'][max_nh]
  ####x, y, z = MySnap.fields['x'][idx], MySnap.fields['y'][idx], MySnap.fields['z'][idx]
  ####dist_to_sink = np.sqrt((x-xsink)**2. + (y-ysink)**2. + (z-zsink)**2.)

  ####fig = pl.figure()
  #####pl.loglog(jeans_number[0::everypart], density[0::everypart], ',', alpha=0.5)
  #####pl.loglog(dist_to_sink[0::everypart], density[0::everypart], ',', alpha=0.5)
  ####pl.loglog(dist_to_sink[0::everypart], MySnap.fields['mass'][idx][0::everypart], ',', alpha=0.5)
  #####pl.semilogx(dist_to_sink[0::everypart], MySnap.fields['allowref'][idx][0::everypart], ',', alpha=0.5)
  #####fig.savefig('./images/jeans_number_output4_%03i.pdf' %snap, dpi=100)
  #####fig.savefig('./images/dist_to_sink_massallowref_output5_%03i.pdf' %snap, dpi=100)
  ####fig.savefig('./images/dist_to_sink_mass_'+base+'_%03i.pdf' %snap, dpi=100)
  ####pl.xlabel('Distance to sink')
  ####pl.ylabel('Number density')
  ####pl.show()
  ####pl.close()

#f.close()
end = time.time()
print 'Total time: %f' %(end - start)

#####plot_fields = ['gravacc']
#####MyRadial = pr.radial.Radial()
#####MyRadial.radial_profile(MySnap, plot_fields)
#####fig  = pl.figure()
#####gravacc = 10**MyRadial.radial['gravacc'] * pr.constants.UNIT_LENGTH / pr.constants.UNIT_TIME**2 / 1e5 * 1e6 * pr.constants.SEC_PER_YEAR
#####pl.loglog(10**MyRadial.radial['radius'], gravacc)
#####fig.savefig('./gravacc_nahw1.pdf', dpi=100)
#####pl.show()
#####pl.close()
###
###print MyRadial.radial['gravacc']
###print MyRadial.radial['radius']

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

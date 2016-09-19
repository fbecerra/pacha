import matplotlib
matplotlib.use('Agg')
import src as pr
import numpy as np
import matplotlib.pyplot as pl
import time
import h5py

def find_max(array, exclude=np.array([])):
  mask = np.ones(array.shape, dtype=np.bool)
  all_idx = np.arange(len(mask))
  if len(exclude) > 0:
    mask[exclude] = False
  idx_mask = np.argmax(array[mask])
  idx = all_idx[mask][idx_mask]
  value = array[idx]
  return idx, value

init_snap = 181
final_snap = 182
snaps = np.arange(init_snap, final_snap)
if init_snap == 1:
  protocount = 1
else:
  protocount = 17
n_thres = 10**18.5
fraction_thres = 0.5
merger_radius = 0.1
diff_thres = -500
r_merger = 0.1

start = time.time()

for snap in snaps:

  print snap
  snapshot = '%03i' %snap
  path = '/scratch/00025/tgreif/ah1w3f2/snapdir_'+snapshot+'/'
  file = 'ah1w3f2_%03i' %snap
  f = h5py.File('/scratch/02563/fbecerra/paha/protostars/'+file+'.hdf5', 'w')
  snapbase = path+file
  
  MySnap = pr.snap.Snap()
  MySnap.read_header(snapbase)
  MySnap.read_fields(snapbase)
  del MySnap.fields['abHM'], MySnap.fields['abH2'], MySnap.fields['abHII']
#  MySnap.center_box()
#  MySnap.rotate_box()
  x, y, z = MySnap.fields['x'], MySnap.fields['y'], MySnap.fields['z']
  vx, vy, vz = MySnap.fields['vx'], MySnap.fields['vy'], MySnap.fields['vz']
  nh = MySnap.fields['nh']
  excluded_particles = np.array([], dtype=np.int32)
  
  radial_fields = ['kappa', 'rho', 'vrad']
     
  max_dens_idx, max_dens_val = find_max(nh)
  local_count = 0

  while max_dens_val > n_thres:

    MyRadial = pr.radial.Radial()

    mask = np.ones(nh.shape, dtype=np.bool)
    if len(excluded_particles) > 0:
      mask[excluded_particles] = False
    print len(np.where(nh[mask] > n_thres)[0])

    if len(np.where(nh[mask] > n_thres)[0]) < 80000:
      break
    max_mass = 0
    del mask

    # Get protostars properties
    protostar_center = np.array([x[max_dens_idx], y[max_dens_idx], z[max_dens_idx]])
    MySnap.fields['x'] = x - protostar_center[0]
    MySnap.fields['y'] = y - protostar_center[1]
    MySnap.fields['z'] = z - protostar_center[2]
    MySnap.fields['vx'], MySnap.fields['vy'], MySnap.fields['vz'] = vx, vy, vz
    MySnap.calculate_radius()
    MySnap.rotate_box()

    MyRadial.radial_profile(MySnap, radial_fields)
#    del MySnap.fields['temp'], MySnap.fields['u'], MySnap.fields['gamma']
    diff = 10**MyRadial.radial['kappa'][2:] - 10**MyRadial.radial['kappa'][:-2]
#    rad_idx =  min(np.argsort(diff)[:3])
#    diff = np.diff(10**MyRadial.radial['kappa'], n=2)
    index = 0
    rad_idx = np.where(diff < diff_thres)[0][index]
    protostar_radius = 10**MyRadial.radial['radius'][rad_idx]
    while protostar_radius < r_merger:
      index += 1
      rad_idx = np.where(diff < diff_thres)[0][index]
      protostar_radius = 10**MyRadial.radial['radius'][rad_idx]
    del diff
    protostar_idx = np.where(MySnap.fields['radius'] < protostar_radius)[0]
    protostar_mass = pr.utils.total_mass(MySnap, protostar_idx)
    protostar_vel = pr.utils.vel_center_of_mass(MySnap, protostar_idx)
    protostar_ids = MySnap.fields['id'][protostar_idx]
    protostar_id = protocount
    protostar_rho = 10**MyRadial.radial['rho'][rad_idx]
    del MyRadial.radial['rho']
    protostar_vrad = MyRadial.radial['vrad'][rad_idx]
    protostar_progenitors = np.array([])
    excluded_particles = np.append(excluded_particles, protostar_idx)
    del protostar_idx

    # Tracking
    if snap != 1:
      previous_file = 'ah1w3f2_%03i' %(snap-1)
      previous_f = h5py.File('/scratch/02563/fbecerra/paha/protostars/'+previous_file+'.hdf5', 'r')
      for key in previous_f.keys():
        previous_ids = np.array(previous_f[key+'/PartIDs'])#, dtype = np.int64)
        #Match
        part_in_common = np.intersect1d(protostar_ids, previous_ids)
        fraction = float(len(part_in_common)) / float(len(previous_ids))
        fraction2 = float(len(part_in_common)) / float(len(protostar_ids))
        if ((fraction > fraction_thres) | (fraction2 > fraction_thres)):
          # Add progenitor
          progenitor_id = np.int32(key.split('Proto')[1])
          progenitor_mass = np.float64(previous_f[key+'/Mass'])
          if progenitor_mass > max_mass:
            protostar_id = progenitor_id
            max_mass = progenitor_mass
          protostar_progenitors = np.append(protostar_progenitors, progenitor_id)
      del previous_ids, part_in_common, max_mass
      previous_f.close()

    # Save data
    print protostar_id
    try:
      f.create_dataset('Proto%i/Position' %protostar_id, data = protostar_center)
    except:
      continue
    f.create_dataset('Proto%i/Velocity' %protostar_id, data = protostar_vel)
    f.create_dataset('Proto%i/Radius' %protostar_id, data = protostar_radius)
    f.create_dataset('Proto%i/Mass' %protostar_id, data = protostar_mass)
    f.create_dataset('Proto%i/PartIDs' %protostar_id, data = protostar_ids)
    f.create_dataset('Proto%i/Rho' %protostar_id, data = protostar_rho)
    f.create_dataset('Proto%i/Vrad' %protostar_id, data = protostar_vrad)
    f.create_dataset('Proto%i/Progenitors' %protostar_id, data = protostar_progenitors)
    print 'Protostar ', protostar_id
    print 'Position ', protostar_center
    print 'Radius: ', protostar_radius
    print 'Progenitors ', protostar_progenitors
    print 'Total mass: ', protostar_mass
    del protostar_center, protostar_vel, protostar_mass, protostar_ids, protostar_rho, protostar_vrad

    # Find new maximum
    max_dens_idx, max_dens_val = find_max(nh, exclude=excluded_particles)

    # Plotting stuff
    pl.plot(MyRadial.radial['radius'], MyRadial.radial['kappa'])
    del MyRadial.radial['kappa']
    ymin, ymax = pl.ylim()
    pl.plot(np.array([np.log10(protostar_radius), np.log10(protostar_radius)]), np.array([ymin, ymax]), '--')
    pl.savefig('/scratch/02563/fbecerra/paha/protostars/kappa_'+file+'_'+str(protostar_id)+'.pdf')
    pl.show()
    pl.close()

    pl.plot(MyRadial.radial['radius'], MyRadial.radial['vrad'])
    ymin, ymax = pl.ylim()
    pl.plot(np.array([np.log10(protostar_radius), np.log10(protostar_radius)]), np.array([ymin, ymax]), '--')
    pl.savefig('/scratch/02563/fbecerra/paha/protostars/vrad_'+file+'_'+str(protostar_id)+'.pdf')
    pl.show()
    pl.close()
    del MyRadial

    pl.clf()

    if len(protostar_progenitors) == 0:
      protocount += 1

  f.close()

end = time.time()
print 'Total time: %f secs' %(end-start)

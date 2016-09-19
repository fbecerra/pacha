import matplotlib
matplotlib.use('Agg')
import src as pr
import numpy as np
import matplotlib.pyplot as pl
import time
import h5py
import os

def find_max(array, include=np.array([])):
  mask = np.zeros(array.shape, dtype=np.bool)
  all_idx = np.arange(len(mask))
  if len(include) > 0:
    mask[include] = True
  idx_mask = np.argmax(array[mask])
  idx = all_idx[mask][idx_mask]
  value = array[idx]
  return idx, value

def update_data(path, file, id, progenitors):

  f = h5py.File(fpath+file+'.hdf5', 'r')
  fnew = h5py.File(fpath+'fixed_'+file+'.hdf5', 'w')
  proto = 'Proto' + str(np.int32(id))
  for key in f.keys():
    protostar_center = np.array(f[key+'/Position'])
    protostar_vel = np.array(f[key+'/Velocity'])
    protostar_radius = np.array(f[key+'/Radius'])
    protostar_mass = np.array(f[key+'/Mass'])
    protostar_ids = np.array(f[key+'/PartIDs'])
    protostar_rho = np.array(f[key+'/Rho'])
    protostar_vrad = np.array(f[key+'/Vrad'])      
    if key == proto:
      protostar_progenitors = np.array([progenitors]) ###progenitors
    else:
      try:
        protostar_progenitors = np.array(f[key+'/Progenitors'])
      except:
        protostar_progenitors = np.array([])
    fnew.create_dataset(key+'/Position', data = protostar_center)
    fnew.create_dataset(key+'/Velocity', data = protostar_vel)
    fnew.create_dataset(key+'/Radius', data = protostar_radius)
    fnew.create_dataset(key+'/Mass', data = protostar_mass)
    fnew.create_dataset(key+'/PartIDs', data = protostar_ids)
    fnew.create_dataset(key+'/Rho', data = protostar_rho)
    fnew.create_dataset(key+'/Vrad', data = protostar_vrad)
    fnew.create_dataset(key+'/Progenitors', data = protostar_progenitors)
  f.close()
  fnew.close()
  os.system('cp '+fpath+'fixed_'+file+'.hdf5 '+fpath+file+'.hdf5')
  os.system('rm '+fpath+'fixed_'+file+'.hdf5')

def add_data(path, file, id, pos, vel, rad, mass, ids, rho, vrad, progenitors):

  f = h5py.File(fpath+file+'.hdf5', 'r')
  fnew = h5py.File(fpath+'fixed_'+file+'.hdf5', 'w')
  proto = 'Proto' + str(np.int32(id))

  for key in f.keys():
    protostar_center = np.array(f[key+'/Position'])
    protostar_vel = np.array(f[key+'/Velocity'])
    protostar_radius = np.array(f[key+'/Radius'])
    protostar_mass = np.array(f[key+'/Mass'])
    protostar_ids = np.array(f[key+'/PartIDs'])
    protostar_rho = np.array(f[key+'/Rho'])
    protostar_vrad = np.array(f[key+'/Vrad'])
    try:
      protostar_progenitors = np.array(f[key+'/Progenitors'])
    except:
      protostar_progenitors = np.array([])

    fnew.create_dataset(key+'/Position', data = protostar_center)
    fnew.create_dataset(key+'/Velocity', data = protostar_vel)
    fnew.create_dataset(key+'/Radius', data = protostar_radius)
    fnew.create_dataset(key+'/Mass', data = protostar_mass)
    fnew.create_dataset(key+'/PartIDs', data = protostar_ids)
    fnew.create_dataset(key+'/Rho', data = protostar_rho)
    fnew.create_dataset(key+'/Vrad', data = protostar_vrad)
    fnew.create_dataset(key+'/Progenitors', data = protostar_progenitors)

  fnew.create_dataset(proto+'/Position', data = pos)
  fnew.create_dataset(proto+'/Velocity', data = vel)
  fnew.create_dataset(proto+'/Radius', data = rad)
  fnew.create_dataset(proto+'/Mass', data = mass)
  fnew.create_dataset(proto+'/PartIDs', data = ids)
  fnew.create_dataset(proto+'/Rho', data = rho)
  fnew.create_dataset(proto+'/Vrad', data = vrad)
  fnew.create_dataset(proto+'/Progenitors', data = progenitors)

  f.close()
  fnew.close()
  os.system('cp '+fpath+'fixed_'+file+'.hdf5 '+fpath+file+'.hdf5')
  os.system('rm '+fpath+'fixed_'+file+'.hdf5')

initial_snap = 147
final_snap = 148
snaps = np.arange(final_snap, initial_snap, -1)
n_thres = 10**18.3
fraction_thres = 0.5
merger_radius = 0.1
diff_thres = -500
r_merger = 0.1
radial_fields = ['kappa', 'rho', 'vrad']

for snap in snaps:

  print snap
  path = '/scratch/00025/tgreif/ah1w3f2/snapdir_%03i/' %(snap-1)
  previous_file = 'ah1w3f2_%03i' %(snap-1)
  file = 'ah1w3f2_%03i' %snap
  snapbase = path+previous_file

  fpath = '/scratch/02563/fbecerra/paha/protostars/'
  f = h5py.File(fpath+file+'.hdf5', 'r')
#####  fnew = h5py.File(fpath+'fixed_'+file+'.hdf5', 'w')
#####
#####  for key in f.keys():
#####    id = np.float64(key.split('Proto')[1])
#####    protostar_center = np.array(f[key+'/Position'])
#####    protostar_vel = np.array(f[key+'/Velocity'])
#####    protostar_radius = np.array(f[key+'/Radius'])
#####    protostar_mass = np.array(f[key+'/Mass'])
#####    protostar_ids = np.array(f[key+'/PartIDs'])
#####    protostar_rho = np.array(f[key+'/Rho'])
#####    protostar_vrad = np.array(f[key+'/Vrad'])
#####    try:
#####      protostar_progenitors = np.array(f[key+'/Progenitors'])
#####    except:
#####      protostar_progenitors = np.array([])
#####
#####    if key == 'Proto7':
#####      key = 'Proto5'
#####    elif key == 'Proto8':
#####      key = 'Proto6'
#####    elif key == 'Proto9':
#####      key = 'Proto7'
#####    elif key == 'Proto11':
#####      key = 'Proto8'
#####    elif key == 'Proto15':
#####      key = 'Proto13'
#####
#####    idx = np.where(protostar_progenitors == id)[0]
#####    newid = np.float64(key.split('Proto')[1])
#####    if len(idx) != 0:
#####      protostar_progenitors[idx] = newid
#####
#####    fnew.create_dataset(key+'/Position', data = protostar_center)
#####    fnew.create_dataset(key+'/Velocity', data = protostar_vel)
#####    fnew.create_dataset(key+'/Radius', data = protostar_radius)
#####    fnew.create_dataset(key+'/Mass', data = protostar_mass)
#####    fnew.create_dataset(key+'/PartIDs', data = protostar_ids)
#####    fnew.create_dataset(key+'/Rho', data = protostar_rho)
#####    fnew.create_dataset(key+'/Vrad', data = protostar_vrad)
#####    fnew.create_dataset(key+'/Progenitors', data = protostar_progenitors)
#####  f.close()
#####  fnew.close()
#####  os.system('cp '+fpath+'fixed_'+file+'.hdf5 '+fpath+file+'.hdf5')
#####  os.system('rm '+fpath+'fixed_'+file+'.hdf5')


  for key in f.keys():

    print key
    try:
      progenitors = np.array(f[key+'/Progenitors'])
    except:
      progenitors = np.array([])

    if len(progenitors) == 0:
###    if key == 'Proto6':
###
###      protostar_id = np.float64(key.split('Proto')[1])
###      protostar_progenitors = np.array([6., 10.])
###      update_data(fpath, file, protostar_id, protostar_progenitors)

      MySnap = pr.snap.Snap()
      MySnap.read_header(snapbase)
      MySnap.read_fields(snapbase)

      ids = np.array(f[key+'/PartIDs'])
      particles = np.in1d(MySnap.fields['id'], ids)

      x, y, z = MySnap.fields['x'], MySnap.fields['y'], MySnap.fields['z']
      vx, vy, vz = MySnap.fields['vx'], MySnap.fields['vy'], MySnap.fields['vz']
      nh = MySnap.fields['nh']
      include = np.array([], dtype=np.int64)

      cmx = np.sum(nh[particles] * x[particles]) / np.sum(nh[particles])
      cmy = np.sum(nh[particles] * y[particles]) / np.sum(nh[particles])
      cmz = np.sum(nh[particles] * z[particles]) / np.sum(nh[particles])
      particles = np.where((x-cmx)**2 + (y-cmy)**2 + (z-cmz)**2 < 0.5**2)[0]

      xmin, xmax = np.min(x[particles]), np.max(x[particles])
      ymin, ymax = np.min(y[particles]), np.max(y[particles])
      zmin, zmax = np.min(z[particles]), np.max(z[particles])
      particles = np.where((xmin < x) & (x < xmax) & (ymin < y) & (y < ymax) & (zmin < z) & (z < zmax))[0]

      include = np.append(include, particles)

      max_dens_idx, max_dens_val = find_max(nh, include=include)

      if max_dens_val > n_thres:

        print 'Protostar detected'

        # Get protostars properties
        protostar_center = np.array([x[max_dens_idx], y[max_dens_idx], z[max_dens_idx]])
        MySnap.fields['x'] = x - protostar_center[0]
        MySnap.fields['y'] = y - protostar_center[1]
        MySnap.fields['z'] = z - protostar_center[2]
        MySnap.fields['vx'], MySnap.fields['vy'], MySnap.fields['vz'] = vx, vy, vz
        MySnap.calculate_radius()
        MySnap.rotate_box()

        print 'Calculating its properties'

        MyRadial = pr.radial.Radial()
        MyRadial.radial_profile(MySnap, radial_fields)
        diff = 10**MyRadial.radial['kappa'][2:] - 10**MyRadial.radial['kappa'][:-2]
        idx_diff = np.where(diff < diff_thres)[0]
        if idx_diff[len(idx_diff) - 1] < r_merger:
          continue
        index = 0
        rad_idx = idx_diff[index]
        protostar_radius = 10**MyRadial.radial['radius'][rad_idx]
        while protostar_radius < r_merger:
          index += 1
          rad_idx = idx_diff[index]
          protostar_radius = 10**MyRadial.radial['radius'][rad_idx]
        protostar_idx = np.where(MySnap.fields['radius'] < protostar_radius)[0]
        protostar_mass = pr.utils.total_mass(MySnap, protostar_idx)
        protostar_vel = pr.utils.vel_center_of_mass(MySnap, protostar_idx)
        protostar_ids = MySnap.fields['id'][protostar_idx]
        protostar_id = np.float64(key.split('Proto')[1])
        protostar_rho = 10**MyRadial.radial['rho'][rad_idx]
        protostar_vrad = MyRadial.radial['vrad'][rad_idx]
        protostar_progenitors = np.array([])

        f_prev = h5py.File('/scratch/02563/fbecerra/paha/protostars/'+previous_file+'.hdf5', 'r')

        # Did the protostar exist before?
        for key_prev in f_prev.keys():

          print 'Checking if it already existed'

          ids_prev = np.array(f_prev[key_prev+'/PartIDs'])
          part_in_common = np.intersect1d(protostar_ids, ids_prev)
          fraction = float(len(part_in_common)) / float(len(ids_prev))
          fraction2 = float(len(part_in_common)) / float(len(protostar_ids))
          if ((fraction > fraction_thres) | (fraction2 > fraction_thres)):
            # Add progenitor
            protostar_progenitors = np.float64(key_prev.split('Proto')[1])
            f.close()
            f_prev.close()
            print 'Yes, updating progenitors'
            update_data(fpath, file, protostar_id, protostar_progenitors)
            f = h5py.File(fpath+file+'.hdf5', 'r')
            break
        else:
          f_prev.close()
          print 'No, creating new dataset'
          add_data(fpath, previous_file, protostar_id, protostar_center, protostar_vel, protostar_radius, protostar_mass,
                                         protostar_ids, protostar_rho, protostar_vrad, protostar_progenitors)
          update_data(fpath, file, protostar_id, protostar_id)

  f.close()

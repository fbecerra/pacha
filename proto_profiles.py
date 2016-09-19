import matplotlib
matplotlib.use('TkAgg')
import src as pr
import numpy as np
import pylab as pl
import h5py

import time

#pc = pr.plot_collection.PlotCollection()
#pc.proto_profiles('/scratch/00025/tgreif/', 'ah1w3f2')
#pc.savefig('/scratch/02563/fbecerra/paha/protostars/time_profiles.pdf')


colors = ['b', 'r', 'g', 'c', 'm', 'y', '#006666', '#FF6600', '#6699FF', '#993300', '#FF3399']
initial_snap = 1
final_snap = 182
snaps = np.arange(initial_snap, final_snap)
color_index = 0
times = {}
protostar = {}
color = {}
keys = np.array([])

#for snap in snaps:
#
#  print snap
#  path = '/scratch/00025/tgreif/ah1w3f2/snapdir_%03i/' %snap
#  previous_path = '/scratch/00025/tgreif/ah1w3f2/snapdir_%03i/' %(snap-1)
#  file = 'ah1w3f2_%03i' %snap
#  previous_file = 'ah1w3f2_%03i' %(snap-1)
#  previous_snapbase = path + previous_file
#  snapbase = path + file
#
#  fpath = '/scratch/02563/fbecerra/paha/protostars/'
#  f = h5py.File(fpath+file+'.hdf5', 'r')
#
##  MySnap = pr.snap.Snap()
##  MySnap.read_header(snapbase)
##  snaptime = MySnap.params['time'] * pr.constants.UNIT_TIME / pr.constants.SEC_PER_YEAR
#
#  for key in f.keys():
#
#    print key
#    if key not in keys:
#      keys = np.append(keys, key)
#      times[key] = np.array([])
#      protostar[key] = np.array([])
#    number = np.int32(key.split('Proto')[1])
#    try:
#      progenitors = np.array(f[key+'/Progenitors'])
#    except:
#      progenitors = np.array([])
#    nprogenitors = len(progenitors)
#
#    if nprogenitors == 0:
#      color[key] = colors[color_index]
#      color_index += 1
#    elif nprogenitors == 1:
#      progenitor = np.int32(progenitors[0])
#      prog = 'Proto'+str(progenitor)
#      color[key] = color[prog]
#      times[key] = np.append(times[key], snap-1)
#      protostar[key] = np.append(protostar[key], progenitor)
#    else:
#      for progenitor in progenitors:
#        prog = 'Proto'+str(np.int32(progenitor))
#        if prog != key:
#          color[prog] = color[key]
#          times[prog] = np.append(times[prog], snap-1)
#          protostar[prog] = np.append(protostar[prog], progenitor)
#          times[prog] = np.append(times[prog], snap)
#          protostar[prog] = np.append(protostar[prog], number)
#        else:
#          times[key] = np.append(times[key], snap-1)
#          protostar[key] = np.append(protostar[key], progenitor)
#
#for key in keys:
#  pl.plot(times[key], protostar[key], color=color[key])
#
#pl.ylim(0,11)
#pl.savefig('/scratch/02563/fbecerra/paha/protostars/merger_tree.pdf')
#pl.show()
#pl.close()
    

base = '/scratch/00025/tgreif/'
sim = 'ah1w3f2'
distance = np.array([])
time = np.array([])
velocity = np.array([])
mass1 = np.array([])
mass2 = np.array([])
vescape1 = np.array([])
vescape2 = np.array([])
snaps = range(1, 182)
for snap in snaps:

  print snap
  path = base+sim+'/snapdir_%03i/' %snap
  file = sim+'_%03i' %snap
  snapbase = path+file

  MySnap = pr.snap.Snap()
  MySnap.read_header(snapbase)
  snaptime = MySnap.params['time'] * pr.constants.UNIT_TIME / pr.constants.SEC_PER_YEAR
  f = h5py.File('/scratch/02563/fbecerra/paha/protostars/'+file+'.hdf5', 'r')
  key1, key2 = 'Proto1', 'Proto3'
  if key2 in f.keys():
    d1 = np.array(f[key1+'/Position'])
    d2 = np.array(f[key2+'/Position'])
    dist = np.linalg.norm(d2-d1)
    distance = np.append(distance, dist)
    v1 = np.array(f[key1+'/Velocity'])
    v2 = np.array(f[key2+'/Velocity'])
    velocity = np.append(velocity, np.linalg.norm(v2-v1))

    MySnap.read_fields(snapbase)
    x, y, z = MySnap.fields['x'], MySnap.fields['y'], MySnap.fields['z']

    idx1 = np.where((x - d1[0])**2 + (y - d1[1])**2 + (z - d1[2])**2 < dist**2)[0]
    menc1 = np.sum(MySnap.fields['mass'][idx1])
    mass1 = np.append(mass1, menc1)

    idx2 = np.where((x - d2[0])**2 + (y - d2[1])**2 + (z - d2[2])**2 < dist**2)[0]
    menc2 = np.sum(MySnap.fields['mass'][idx2])
    mass2 = np.append(mass2, menc2)

    #radius = np.array(f[key1+'/Radius'])
    time = np.append(time, snaptime)
    f.close()
  else:
    f.close()
    continue

vescape1 = np.sqrt(2 * pr.constants.GRAVITY * mass1 * pr.constants.SOLAR_MASS / distance / pr.constants.ASTRONOMICAL_UNIT) / 1e5
vescape2 = np.sqrt(2 * pr.constants.GRAVITY * mass2 * pr.constants.SOLAR_MASS / distance / pr.constants.ASTRONOMICAL_UNIT) / 1e5

#print distance
#print velocity
#print time
#print vescape1
#print vescape2

# Plot
fig = pl.figure()
ax = pl.subplot(121)
ax.plot(time, distance)
pl.xlim(np.min(time), np.max(time))
pl.ylabel(r'$d$ [au]', fontsize=24)
pl.xlabel(r'$t$ [yr]', fontsize=24)
pl.tick_params(labelsize=18)
#ax.set_xticklabels([])
ax = pl.subplot(122)
ax.plot(time, velocity, label=r'$v_{\rm rel}$')
ax.plot(time, vescape1, '--', label=r'$v_{\rm esc, 1}$')
ax.plot(time, vescape2, ':', label=r'$v_{\rm esc, 2}$')
ax.yaxis.tick_right()
ax.yaxis.set_ticks_position('both')
ax.yaxis.set_label_position('right')
pl.xlim(np.min(time), np.max(time))
pl.legend(loc='upper right', prop={'size':18}, frameon=True)
pl.ylabel(r'$v$ [km/s]', fontsize=24)
pl.xlabel(r'$t$ [yr]', fontsize=24)
pl.tick_params(labelsize=18)
fig.set_size_inches(14,7)
pl.subplots_adjust(wspace=0.0, hspace=0.0)
pl.savefig('/scratch/02563/fbecerra/paha/protostars/distance.pdf')
pl.show()
pl.close()

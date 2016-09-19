import matplotlib
matplotlib.use('TkAgg')
import src as pr
import numpy as np
import pylab as pl
#from scipy.signal import argrelmax
import time
import h5py

plot_protostar = True
#CenterProto = 15
ImgWidth = 30
field = 'nh'
Rotate = True

for snap in xrange(161, 182):

  print snap
  path = '/scratch/00025/tgreif/ah1w3f2/snapdir_%03i/' %snap
  file = 'ah1w3f2_%03i' %snap
  snapbase = path+file

  MySnap = pr.snap.Snap()
  MySnap.read_header(snapbase)
  MySnap.read_fields(snapbase)

  x, y, z = MySnap.fields['x'], MySnap.fields['y'], MySnap.fields['z']
#  vx, vy, vz = MySnap.fields['vx'], MySnap.fields['vy'], MySnap.fields['vz']

  if plot_protostar:
    f = h5py.File('/scratch/02563/fbecerra/paha/protostars/'+file+'.hdf5', 'r')
  
  for CenterProto in xrange(1, 17):

    if plot_protostar:
      try:
        Center = np.array(f['Proto%i/Position' %CenterProto])
      except:
        continue
    
    if plot_protostar:
      MySnap.fields['x'] = x - Center[0]
      MySnap.fields['y'] = y - Center[1]
      MySnap.fields['z'] = z - Center[2]
#      MySnap.fields['vx'], MySnap.fields['vy'], MySnap.fields['vz'] = vx, vy, vz
    else:
      MySnap.center_box()
    if Rotate:
      MySnap.rotate_box()
    
    MyImage = pr.image.Image()
    MyImage.calculate_image(MySnap, field, ImgWidth=ImgWidth)
    img = MyImage.img
    fig = pl.figure()
    ax = fig.add_subplot(111)
    ax.imshow(img, cmap = 'jet', extent = [0, MyImage.width, 0, MyImage.height], vmin = np.min(img), vmax = np.max(img))

    if plot_protostar:
      for key in f.keys():
        protostar_center = np.array(f[key+'/Position'])
        protostar_center -= Center
        if Rotate:
          x_proto = MySnap.params['rot12'][0][0] * protostar_center[0] + MySnap.params['rot12'][0][1] * protostar_center[1] + MySnap.params['rot12'][0][2] * protostar_center[2]
          y_proto = MySnap.params['rot12'][1][0] * protostar_center[0] + MySnap.params['rot12'][1][1] * protostar_center[1] + MySnap.params['rot12'][1][2] * protostar_center[2]
        else:
          x_proto = protostar_center[0]
          y_proto = protostar_center[1]
        protostar_radius = np.float64(f[key+'/Radius'])
        x_proto += MyImage.width / 2
        y_proto += MyImage.height / 2
        if ((x_proto > 0) & (x_proto < MyImage.width) & (y_proto > 0) & (y_proto < MyImage.height)):
          #ax.plot(x_proto, y_proto, 'k+')
          ax.text(x_proto, y_proto, key.split('Proto')[1], fontsize=8)
          circ = pl.Circle((x_proto, y_proto), radius=protostar_radius, edgecolor='k', fill=False)
          ax.add_patch(circ)

    pl.savefig('/scratch/02563/fbecerra/paha/images/proto%i_' %CenterProto +field+'_'+file+'.pdf')
    pl.show()
    pl.close()

  if plot_protostar:
    f.close()

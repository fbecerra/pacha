import matplotlib
matplotlib.use('TkAgg')
import src as pr
import numpy as np
import pylab as pl
#from scipy.signal import argrelmax
import time
import h5py

UNIT_LENGTH = 3.085678e21
UNIT_MASS = 1.989e43
UNIT_VELOCITY = 1e5
UNIT_TIME = UNIT_LENGTH / UNIT_VELOCITY
SEC_PER_YEAR = 3.155e7

frame = 9999
plot_protostar = True
CenterProto = 1
ImgWidth = 20
field = 'nh'
Rotate = True
sim = 'ah1w3f2'


#path = '/scratch/00025/tgreif/'+sim+'/snapdir_126/'
#file = sim+'_126' 
#if plot_protostar:
#  f = h5py.File('/scratch/02563/fbecerra/paha/protostars/'+file+'.hdf5', 'r')
#  try:
#    Center = np.array(f['Proto%i/Position' %CenterProto])
#  except:
#    print 'Protostar not located'
#  f.close()

for snap in xrange(126, 127):

  print snap
  path = '/scratch/00025/tgreif/'+sim+'/snapdir_%03i/' %snap
  file = sim+'_%03i' %snap
  snapbase = path+file

  MySnap = pr.snap.Snap()
  MySnap.read_header(snapbase)
  MySnap.read_fields(snapbase)

  x, y, z = MySnap.fields['x'], MySnap.fields['y'], MySnap.fields['z']
#  vx, vy, vz = MySnap.fields['vx'], MySnap.fields['vy'], MySnap.fields['vz']

  if plot_protostar:
    f = h5py.File('/scratch/02563/fbecerra/paha/protostars/'+file+'.hdf5', 'r')
  
  if plot_protostar:
    try:
      Center = np.array(f['Proto%i/Position' %CenterProto])
    except:
      continue
  
  if plot_protostar:
    MySnap.fields['x'] = x - Center[0]
    MySnap.fields['y'] = y - Center[1]
    MySnap.fields['z'] = z - Center[2]
#    MySnap.fields['vx'], MySnap.fields['vy'], MySnap.fields['vz'] = vx, vy, vz
  else:
    MySnap.center_box()
  if Rotate:
    MySnap.rotate_box()
 
  for ImgWidth in np.arange(20, 30, 10):
 
    MyImage = pr.image.Image()
    MyImage.calculate_image(MySnap, field, ImgWidth=ImgWidth)
    img = MyImage.img
    fig = pl.figure(frameon=False)
    ax = fig.add_subplot(111)
    vmin, vmax = np.min(img), np.max(img)
    im = ax.imshow(img, cmap = 'jet', extent = [0, MyImage.width, 0, MyImage.height], vmin = vmin, vmax = vmax)
  #  ax.text(MyImage.width/20., MyImage.height*(1- 1/20.), 'z = %.2f' %MySnap.params['redshift'], fontsize=24, color='w')
#    ax.text(MyImage.width/20., MyImage.height*(1- 1/20.), 't = %.2f yr' %(MySnap.params['time']*UNIT_TIME/SEC_PER_YEAR), fontsize=24, color='w')
#    ax.text(MyImage.width/20., MyImage.height/20., '%i au' %ImgWidth, fontsize=24, color='w')
    ax.set_xticklabels([])
    ax.set_yticklabels([])
#    cax = fig.add_axes([0.9, 0.2, 0.02, 0.6])
#    cb = pl.colorbar(im, cax = cax, orientation='vertical', ticks=np.arange(vmin, vmax + 1).astype(int))
#    for t in cb.ax.get_yticklabels():
#       t.set_color('w')
#       t.set_fontsize(18)
#    cb.outline.set_color('w')
#    cb.set_label(r'${\rm log}\,n_{\rm H}\,[{\rm cm}^{-3}]$', fontsize=24, color='w', rotation=90, labelpad=-80)
#    cb.ax.tick_params(labelsize=18, color='w')
    fig.tight_layout()
    fig.set_size_inches(10, 10)
    pl.savefig('/scratch/02563/fbecerra/paha/images/movie_%04d.pdf' %frame, bbox_inches='tight', pad_inches=0)
    pl.show()
    pl.close()
  
    frame += 1

  if plot_protostar:
    f.close()

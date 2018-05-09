import matplotlib
matplotlib.use('TkAgg')
import src as pr
import numpy as np
import pylab as pl
#from scipy.signal import argrelmax
import time
import h5py

def get_size(width):
  if width >= 1:
    label = '%i pc' %width
  elif width >= 0.1:
    label = '%.1f pc' %width
  elif width >= 0.01:
    label = '%.2f pc' %width
  elif width >= 0.001:
    label = '%.3f pc' %width
  elif width > 0.00009:
    label = '%.4f pc' %width
  else:
    label = 'Too small'
  return label

UNIT_LENGTH = 3.085678e21
UNIT_MASS = 1.989e43
UNIT_VELOCITY = 1e5
UNIT_TIME = UNIT_LENGTH / UNIT_VELOCITY
SEC_PER_YEAR = 3.155e7

plot_protostar = False
CenterProto = 1
field = 'abH2'
#field = 'temp'
#field = 'nh'
Rotate = False

#path = '/scratch/00025/tgreif/'+sim+'/snapdir_126/'
#file = sim+'_126' 
#if plot_protostar:
#  f = h5py.File('/scratch/02563/fbecerra/paha/protostars/'+file+'.hdf5', 'r')
#  try:
#    Center = np.array(f['Proto%i/Position' %CenterProto])
#  except:
#    print 'Protostar not located'
#  f.close()

#frame, sim, width_init, width_final, delta_width, snap_init, snap_final, id_center =   1,   'nahw1', 5000, 5000,  -100,   1,  60, 1003198624 
frame, sim, width_init, width_final, delta_width, snap_init, snap_final, id_center = 100,   'nahw1', 5000,  600, -200,  60,  60, 0
#frame, sim, width_init, width_final, delta_width, snap_init, snap_final, id_center = 200,   'nahw1',  600,  600, -100,  60, 100, 0 
#frame, sim, width_init, width_final, delta_width, snap_init, snap_final, id_center = 300,   'nahw1',  600,   20,  -20, 100, 100, 0 
#frame, sim, width_init, width_final, delta_width, snap_init, snap_final, id_center = 400,   'nahw1',   20,   20, -100, 100, 138, 1011053801 
#frame, sim, width_init, width_final, delta_width, snap_init, snap_final, id_center = 500,   'nahw1',   20,    1,   -1, 138, 138, 0
#frame, sim, width_init, width_final, delta_width, snap_init, snap_final, id_center = 600, 'nahw1r2',    1,    1, -100,   1,  15, 0
#frame, sim, width_init, width_final, delta_width, snap_init, snap_final, id_center = 700, 'nahw1r2',    1,  0.1, -0.1,   15,  15, 0
#frame, sim, width_init, width_final, delta_width, snap_init, snap_final, id_center = 710, 'nahw1r2',  0.09, 0.01, -0.01,   15,  15, 0
#frame, sim, width_init, width_final, delta_width, snap_init, snap_final, id_center = 800, 'nahw1r2',  0.01, 0.01, -0.01,   15,  30, 0
#frame, sim, width_init, width_final, delta_width, snap_init, snap_final, id_center = 900, 'nahw1r2',  0.01, 0.001, -0.001,   30,  30, 0
#frame, sim, width_init, width_final, delta_width, snap_init, snap_final, id_center = 910, 'nahw1r2',  0.0009, 0.00005, -0.0001,   30,  30, 0


for snap in xrange(snap_init, snap_final + 1):

  print 'Snap: ', snap
  base = '/n/hernquistfs2/fbecerra/'
  path = base+sim+'/snapdir_%03i/' %snap
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
  elif id_center:
    idx_center = np.where(MySnap.fields['id'] == id_center)[0]
    MySnap.fields['x'] = x - MySnap.fields['x'][idx_center] 
    MySnap.fields['y'] = y - MySnap.fields['y'][idx_center] 
    MySnap.fields['z'] = z - MySnap.fields['z'][idx_center]
  else:
    MySnap.center_box() 
    #center = MySnap.Center
    #print center, np.array(center) * (1 + MySnap.params['redshift']) * MySnap.params['hubbleparam']
  if Rotate:
    MySnap.rotate_box()
 
  for ImgWidth in np.arange(width_init, width_final+delta_width, delta_width):

    print 'Width: ', ImgWidth
 
    MyImage = pr.image.Image()
    MyImage.calculate_image(MySnap, field, ImgWidth=ImgWidth*1e-3)
    img = MyImage.img
    fig = pl.figure(frameon=False)
    ax = fig.add_subplot(111)
    vmin, vmax = np.min(img), np.max(img)
    #if field == 'temp':
    #  vmin, vmax = 1.6, 4.4
    im = ax.imshow(img, cmap = pr.utils.get_colormap(field), extent = [0, MyImage.width, 0, MyImage.height], vmin = vmin, vmax = vmax)

    ## Text
    #ax.text(MyImage.width/20., MyImage.height*(1- 1/20.), 'z = %.2f' %MySnap.params['redshift'], fontsize=24, color='w')
    #ax.text(MyImage.width/20., MyImage.height*(1- 1/20.), 't = %.2f yr' %(MySnap.params['time']*UNIT_TIME/SEC_PER_YEAR), fontsize=24, color='w')
    ax.text(MyImage.width/20., MyImage.height/20., get_size(ImgWidth), fontsize=24, color='w')

    ax.set_xticklabels([])
    ax.set_yticklabels([])

    ## Colorbar
    cax = fig.add_axes([0.88, 0.2, 0.02, 0.6])
    if field == 'nh':
      cb = pl.colorbar(im, cax = cax, orientation='vertical', ticks=np.arange(vmin, vmax + 1).astype(int))
    elif field == 'temp' or field == 'abH2':
      delta = vmax - vmin
      if delta > 2:
        cb = pl.colorbar(im, cax = cax, orientation='vertical', ticks=np.arange(vmin, vmax + 1).astype(int))
      elif delta > 0.2:
        cb = pl.colorbar(im, cax = cax, orientation='vertical', ticks=np.arange(np.around(vmin, decimals=1), np.around(vmax, decimals=1) + 0.1, 0.1))
      else:
        cb = pl.colorbar(im, cax = cax, orientation='vertical', ticks=np.arange(np.around(vmin, decimals=2), np.around(vmax, decimals=2) + 0.01, 0.01))
    for t in cb.ax.get_yticklabels():
       t.set_color('w')
       t.set_fontsize(18)
    cb.outline.set_edgecolor('w')
    cb.set_label(pr.utils.get_label(field), fontsize=24, color='w', rotation=90, labelpad=-80)
    cb.ax.tick_params(labelsize=18, color='w')

    fig.tight_layout()
    fig.set_size_inches(10, 10)
    pl.savefig('/n/home02/fbecerra/pacha_new/movies/movie_'+field+'_%04d.pdf' %frame, bbox_inches='tight', pad_inches=0)
    pl.show()
    pl.close()
  
    frame += 1

  if plot_protostar:
    f.close()

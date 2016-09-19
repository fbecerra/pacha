import matplotlib
matplotlib.use('TkAgg')
import src as pr
import numpy as np
import pylab as pl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1 import ImageGrid

import time
import h5py

start_total = time.time()

field = 'nh'
sim = 'ah1w3f2'
snap = 181
base = '/scratch/00025/tgreif/'
path = base+sim+'/snapdir_%03i/' %snap
file = sim+'_%03i' %snap
snapbase = path+file
Rotate = False

MySnap = pr.snap.Snap()
MySnap.read_header(snapbase)
MySnap.read_fields(snapbase)
MySnap.center_box()
x, y, z = MySnap.fields['x'], MySnap.fields['y'], MySnap.fields['z']
if Rotate:
  MySnap.rotate_box()

fig = pl.figure()

ImgWidth1 = 300
MyImage = pr.image.Image()
MyImage.calculate_image(MySnap, field, ImgWidth=ImgWidth1)
ax1 = pl.axes([0.085, 0.25, 0.5, 0.5])
#vmin, vmax = np.min(MyImage.img), np.max(MyImage.img)
vmin, vmax = 11, 21
im = ax1.imshow(MyImage.img, cmap = pr.utils.get_colormap(field), extent = [0, MyImage.width, 0, MyImage.height], vmin=vmin, vmax=vmax)
ax1.text(MyImage.width*(1-2.5/10.), MyImage.height/20., '%i au' %ImgWidth1, fontsize=36, color='w')
ax1.text(MyImage.width/20., MyImage.height*(1-1./10), '%.2f yr' %(MySnap.params['time']*pr.constants.UNIT_TIME/pr.constants.SEC_PER_YEAR), fontsize=36, color='w')
ax1.set_xticklabels([])
ax1.set_yticklabels([])
#Colorbar
cax = pl.axes([0.035, 0.3, 0.02, 0.4])
cb = pl.colorbar(im, cax = cax, orientation='vertical', ticks=np.arange(vmin, vmax + 1).astype(int))
cb.set_label(pr.utils.get_label(field), fontsize=30, labelpad=-85)
cb.ax.tick_params(labelsize=24)

# Primary clump
f = h5py.File('/scratch/02563/fbecerra/paha/protostars/'+file+'.hdf5', 'r')
protostar_center1 = np.array(f['Proto1/Position'])
protostar_center2 = np.array(f['Proto3/Position'])
protostar_center = (protostar_center1 + protostar_center2) / 2.
protostar_center -= MySnap.Center
if Rotate:
  x_proto = MySnap.params['rot12'][0][0] * protostar_center[0] + MySnap.params['rot12'][0][1] * protostar_center[1] + MySnap.params['rot12'][0][2] * protostar_center[2]
  y_proto = MySnap.params['rot12'][1][0] * protostar_center[0] + MySnap.params['rot12'][1][1] * protostar_center[1] + MySnap.params['rot12'][1][2] * protostar_center[2]
  z_proto = MySnap.params['rot12'][2][0] * protostar_center[0] + MySnap.params['rot12'][2][1] * protostar_center[1] + MySnap.params['rot12'][2][2] * protostar_center[2]
else:
  x_proto, y_proto, z_proto = protostar_center[0], protostar_center[1], protostar_center[2]
MySnap.fields['x'] = x - x_proto
MySnap.fields['y'] = y - y_proto
MySnap.fields['z'] = z - z_proto
ImgWidth2 = 60
dx = dy = ImgWidth2
xi, yi = x_proto - dx/2. + ImgWidth1/2., y_proto - dy/2. + ImgWidth1/2.
square = pl.Rectangle((xi, yi), dx, dy, facecolor ='none', edgecolor='k')
ax1.add_patch(square)
MyImage = pr.image.Image()
MyImage.calculate_image(MySnap, field, ImgWidth=ImgWidth2)
ax2 = pl.axes([0.595, 0.05, 0.4, 0.4])
#vmin, vmax = np.min(MyImage.img), np.max(MyImage.img)
im = ax2.imshow(MyImage.img, cmap = pr.utils.get_colormap(field), extent = [0, MyImage.width, 0, MyImage.height], vmin=vmin, vmax=vmax)
ax2.text(MyImage.width*(1-2.5/10.), MyImage.height/20., '%i au' %ImgWidth2, fontsize=36, color='w')
ax2.set_xticklabels([])
ax2.set_yticklabels([])

##Colorbar
#cax = pl.axes([0.595, 0.03, 0.4, 0.02])
#cb = pl.colorbar(im, cax = cax, orientation='horizontal', ticks=np.arange(vmin, vmax + 1).astype(int))
#cb.set_label(pr.utils.get_label(field), fontsize=30, labelpad=-85)
#cb.ax.tick_params(labelsize=24)

transFigure = fig.transFigure.inverted()
coord1 = transFigure.transform(ax1.transData.transform([xi,yi]))
coord2 = transFigure.transform(ax2.transData.transform([0,0]))
line = matplotlib.lines.Line2D((coord1[0],coord2[0]),(coord1[1],coord2[1]), transform=fig.transFigure, color='k')
fig.lines.append(line)

coord1 = transFigure.transform(ax1.transData.transform([xi,yi+dy]))
coord2 = transFigure.transform(ax2.transData.transform([0,ImgWidth2]))
line = matplotlib.lines.Line2D((coord1[0],coord2[0]),(coord1[1],coord2[1]), transform=fig.transFigure, color='k')
fig.lines.append(line)

# Secondary clump
f = h5py.File('/scratch/02563/fbecerra/paha/protostars/'+file+'.hdf5', 'r')
protostar_center = np.array(f['Proto9/Position'])
protostar_center -= MySnap.Center
if Rotate:
  x_proto = MySnap.params['rot12'][0][0] * protostar_center[0] + MySnap.params['rot12'][0][1] * protostar_center[1] + MySnap.params['rot12'][0][2] * protostar_center[2]
  y_proto = MySnap.params['rot12'][1][0] * protostar_center[0] + MySnap.params['rot12'][1][1] * protostar_center[1] + MySnap.params['rot12'][1][2] * protostar_center[2]
  z_proto = MySnap.params['rot12'][2][0] * protostar_center[0] + MySnap.params['rot12'][2][1] * protostar_center[1] + MySnap.params['rot12'][2][2] * protostar_center[2]
else:
  x_proto, y_proto, z_proto = protostar_center[0], protostar_center[1], protostar_center[2]
MySnap.fields['x'] = x - x_proto
MySnap.fields['y'] = y - y_proto
MySnap.fields['z'] = z - z_proto
ImgWidth3 = 20
dx = dy = ImgWidth3
xi, yi = x_proto - dx/2. + ImgWidth1/2., y_proto - dy/2. + ImgWidth1/2.
square = pl.Rectangle((xi, yi), dx, dy, facecolor ='none', edgecolor='k')
ax1.add_patch(square)
MyImage = pr.image.Image()
MyImage.calculate_image(MySnap, field, ImgWidth=ImgWidth3)
#vmin, vmax = np.min(MyImage.img), np.max(MyImage.img)
ax3 = pl.axes([0.595, 0.55, 0.4, 0.4])
ax3.imshow(MyImage.img, cmap = pr.utils.get_colormap(field), extent = [0, MyImage.width, 0, MyImage.height], vmin=vmin, vmax=vmax)
ax3.text(MyImage.width*(1-2.5/10.), MyImage.height/20., '%i au' %ImgWidth3, fontsize=36, color='w')
ax3.set_xticklabels([])
ax3.set_yticklabels([])

##Colorbar
#cax = pl.axes([0.595, 0.93, 0.4, 0.02])
#cb = pl.colorbar(im, cax = cax, orientation='horizontal', ticks=np.arange(vmin, vmax + 1).astype(int))
#cb.set_label(pr.utils.get_label(field), fontsize=30, labelpad=-85)
#cb.ax.tick_params(labelsize=24)

transFigure = fig.transFigure.inverted()
coord1 = transFigure.transform(ax1.transData.transform([xi,yi]))
coord2 = transFigure.transform(ax3.transData.transform([0,0]))
line = matplotlib.lines.Line2D((coord1[0],coord2[0]),(coord1[1],coord2[1]), transform=fig.transFigure, color='k')
fig.lines.append(line)

coord1 = transFigure.transform(ax1.transData.transform([xi,yi+dy]))
coord2 = transFigure.transform(ax3.transData.transform([0,ImgWidth3]))
line = matplotlib.lines.Line2D((coord1[0],coord2[0]),(coord1[1],coord2[1]), transform=fig.transFigure, color='k')
fig.lines.append(line)

f.close()
fig.set_size_inches(16,16)
fig.savefig('/scratch/02563/fbecerra/paha/images/big_zoom.pdf')
pl.show()
pl.close()

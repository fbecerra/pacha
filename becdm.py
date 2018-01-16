import matplotlib
matplotlib.use('TkAgg')
import src as pr
import numpy as np
import pylab as pl
import time

sims = ['LLGas512', 'BLGas512', 'LBGas512', 'BBGas512']
path = '/n/hernquistfs3/pmocz/Output/CosmologyAxion/BECDMCosmologicalSimsWithBaryons/output/'
base = 'snap'
field = 'rho'
ImgWidth = 200

start = time.time()

def print_minmax(array):
  print np.min(array), np.max(array)


for sim in sims: 
  for snap in np.arange(15, 17):
    print 'Sim: %s, Snapshot: %i' %(sim, snap)
    snapbase = path + sim + '/snapdir_%03i/' %snap + base + '_%03i' %snap
    MySnap = pr.snap.Snap()
    MySnap.read_header(snapbase, read_HDF=True)
    MySnap.read_fields(snapbase, read_HDF=True)
    MySnap.center_box()
  
  #  print_minmax(MySnap.fields['id'])
  #  print_minmax(MySnap.fields['mass'])
  #  print_minmax(MySnap.fields['x'])
  #  print_minmax(MySnap.fields['rho'])
  #  print_minmax(MySnap.fields['nh'])
  
    MyImage = pr.image.Image()
    MyImage.calculate_image(MySnap, field, ImgWidth=ImgWidth)
    img = MyImage.img
    fig = pl.figure()
    ax = fig.add_subplot(111)
    ax.imshow(img, cmap = pr.colormaps.viridis, extent = [0, MyImage.width, 0, MyImage.height], vmin = np.min(img), vmax = np.max(img))
  
    pl.savefig('/n/home02/fbecerra/pacha_new/images/'+sim+'_'+field+'_%03i.pdf' %snap)
    pl.show()
    pl.close()

end = time.time()
print 'Total time: %f' %(end - start)


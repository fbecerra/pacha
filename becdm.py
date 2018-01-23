import matplotlib
matplotlib.use('TkAgg')
import src as pr
import numpy as np
import pylab as pl
import time

sims = ['LLGas512', 'LBGas512', 'BLGas512', 'BBGas512']
#sims = ['LLGas512']
rho_min, rho_max = 7e-31, 6e-24
temp_min, temp_max =  2, 3.5e6

path = '/n/hernquistfs3/pmocz/Output/CosmologyAxion/BECDMCosmologicalSimsWithBaryons/output/'
base = 'snap'
field = 'rho'
ImgWidth = 200

start = time.time()

def print_minmax(array):
  print np.min(array), np.max(array)

idx = 1

fig = pl.figure()

for sim in sims: 
  for snap in np.arange(16, 17):
    print 'Sim: %s, Snapshot: %i' %(sim, snap)
    snapbase = path + sim + '/snapdir_%03i/' %snap + base + '_%03i' %snap
    MySnap = pr.snap.Snap()
    MySnap.read_header(snapbase, read_HDF=True)
    MySnap.read_fields(snapbase, read_HDF=True)
    MySnap.center_box()
  
    print_minmax(MySnap.fields['rho'])
    print_minmax(MySnap.fields['temp'])
  
    MyImage = pr.image.Image()
    MyImage.calculate_image(MySnap, field, ImgWidth=ImgWidth)
    img = MyImage.img

    ax = fig.add_subplot(2, 2, idx) #, axisbg='#081d58')
    ax.imshow(img, cmap = pr.utils.get_colormap(field), extent = [0, MyImage.width, 0, MyImage.height], vmin = rho_min, vmax = rho_max) #np.min(img), np.max(img)
    #try:
    #  nsinks = len(MySnap.new_fields['sinks']['id'])
    #  if nsinks > 0:
    #    x_sinks = pr.params.ImgXBins / 2. + pr.params.ImgXBins * MySnap.new_fields['sinks']['x'] / ImgWidth
    #    y_sinks = pr.params.ImgXBins / 2. + pr.params.ImgXBins * MySnap.new_fields['sinks']['y'] / ImgWidth
    #    index = np.where ((x_sinks > 0) & (x_sinks < pr.params.ImgXBins) & (y_sinks > 0) & (y_sinks < pr.params.ImgXBins))[0]
    #    ax.scatter(x_sinks[index], y_sinks[index], c='#ffffd9', marker=',', edgecolors='none')
    #except:
    #  print 'No stars'
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xlim(0, pr.params.ImgXBins)
    ax.set_ylim(0, pr.params.ImgXBins)
    idx += 1

fig.set_size_inches(16, 16)
left, right, bottom, top = 0.1, 0.9, 0.1, 0.9
pl.subplots_adjust(left=left, right=right, bottom=bottom, top=top, wspace=0.0, hspace=0.0)
pl.savefig('/n/home02/fbecerra/pacha_new/images/all_'+field+'_%03i.pdf' %snap)
pl.show()
pl.close()

end = time.time()
print 'Total time: %f' %(end - start)


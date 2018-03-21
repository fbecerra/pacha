import matplotlib
matplotlib.use('TkAgg')
import src as pr
import numpy as np
import pylab as pl
import time
from mpl_toolkits.axes_grid1 import ImageGrid

sims = ['LLGas512', 'LBGas512', 'BLGas512', 'BBGas512']
sims = ['LLGas512']
rho_min, rho_max = 7e-31, 6e-24
temp_min, temp_max =  2, 3.5e6

mins = [7e-31, 2]
maxs = [6e-24, 3.5e6]

path = '/n/hernquistfs3/pmocz/Output/CosmologyAxion/BECDMCosmologicalSimsWithBaryons/output/'
base = 'snap'
fields = ['rho', 'temp']
ImgWidth = 200

start = time.time()

def print_minmax(array):
  print np.min(array), np.max(array)

for sim in sims: 
  for snap in np.arange(10, 11):
    print 'Sim: %s, Snapshot: %i' %(sim, snap)
    fig = pl.figure()
    snapbase = path + sim + '/snapdir_%03i/' %snap + base + '_%03i' %snap
    MySnap = pr.snap.Snap()
    MySnap.read_header(snapbase, read_HDF=True)
    MySnap.read_fields(snapbase, read_HDF=True)
    MySnap.center_box()
  
    print_minmax(MySnap.fields['rho'])
    print_minmax(MySnap.fields['temp'])

    grid = ImageGrid(fig, 111, nrows_ncols = [2, 1], axes_pad = 0.0, label_mode='L', share_all = False,
                     cbar_location = 'right', cbar_mode = 'edge', cbar_pad = '5%', cbar_size = '5%')

    for idx, field in enumerate(fields): 
      MyImage = pr.image.Image()
      MyImage.calculate_image(MySnap, field, ImgWidth=ImgWidth)
    
      vmin, vmax = np.min(MyImage.img), np.max(MyImage.img)
      print vmin, vmax

      im = grid[idx].imshow(MyImage.img, cmap = pr.utils.get_colormap(field), extent = [0, MyImage.xbins, 0, MyImage.ybins], vmin = vmin, vmax = vmax)

      grid[idx].set_xlim(0, MyImage.xbins)
      grid[idx].set_ylim(0, MyImage.ybins)
      grid[idx].set_xticklabels([])
      grid[idx].set_yticklabels([])

      cb = grid.cbar_axes[idx].colorbar(im)
      #cb.set_label_text(pr.utils.get_label(field), fontsize = 24)
      cb.ax.tick_params(labelsize=18)

    fig.set_size_inches(8, 16)

    fig.savefig('/n/home02/fbecerra/pacha_new/images/'+sim+'_'+'_'.join(fields)+'_%03i.pdf' %snap, dpi=100)
    pl.show()
    pl.close()

end = time.time()
print 'Total time: %f' %(end - start)


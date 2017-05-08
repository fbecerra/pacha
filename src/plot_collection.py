from include import *
from mpl_toolkits.axes_grid1 import ImageGrid
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1 import make_axes_locatable
from radial import Radial
from image import Image
from snap import Snap
from utils import *

def smooth(x,window_len=11,window='hanning'):
        if x.ndim != 1:
                raise ValueError, "smooth only accepts 1 dimension arrays."
        if x.size < window_len:
                raise ValueError, "Input vector needs to be bigger than window size."
        if window_len<3:
                return x
        if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
                raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
        s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
        if window == 'flat': #moving average
                w=np.ones(window_len,'d')
        else:
                w=eval('np.'+window+'(window_len)')
        y=np.convolve(w/w.sum(),s,mode='same')
        return y[window_len:-window_len+1]

class PlotCollection:

  def __init__(self):
    self.fig = pl.figure()

  def multi_radial(self, snap, fields):
    MyRadial = Radial()
    MyRadial.radial_profile(snap, fields)
    nfields = len(fields)
    #ncols = int(np.ceil(np.sqrt(nfields)))
    ncols = 2
    nrows = nfields / ncols
    
    for index in range(nfields):
      field = fields[index]
      ax = pl.subplot(nrows, ncols,index+1)
      if field ==  'cs':
        MyRadial.radial[field] = 10**MyRadial.radial[field]
      if field == 'geff':
        smooth_idx = np.where(MyRadial.radial[field] < 0)[0]
        MyRadial.radial[field][smooth_idx] = (MyRadial.radial[field][smooth_idx+1] + MyRadial.radial[field][smooth_idx-1]) / 2.
        MyRadial.radial[field] = smooth(MyRadial.radial[field], window_len=6)
        idx = range(4,20)
        subsmooth = smooth(MyRadial.radial[field][idx], window_len=10)
        MyRadial.radial[field][idx] = subsmooth
      pl.plot(MyRadial.radial['radius'], MyRadial.radial[field], 'k')
      xmin, xmax = np.min(MyRadial.radial['radius']), np.max(MyRadial.radial['radius'])

      if field == 'mbe_ratio' or field == 'q' or field == 'radial_tcool_ratio' or field == 'tcs_ratio' or field == 'gammie':
        pl.plot(np.array([xmin, xmax]), np.array([0,0]), 'k--')
      if field == 'nh':
        ymin, ymax = pl.ylim()
        pl.plot(np.array([xmin, xmax]), -2*np.array([xmin,xmax]) + 20, 'k--')
        pl.text(4,13, r'$n_{\rm H} \propto r^{-2}$', fontsize=18)
        pl.ylim(ymin, ymax)
    
      # y-axis
      if field == 'geff':
        pl.ylim(0.5, 2)
      if field == 'radial_tcool_ratio':
        pl.ylim(-2, 4)
      if field == 'tcool_ratio':
        ymin, ymax = pl.ylim()
        pl.ylim(0, 3)
      if field == 'gammie':
        pl.ylim(-2, 3)
      if (ncols > 1 and (index+1) % ncols == 0):
        ax.yaxis.tick_right()
        ax.yaxis.set_ticks_position('both')
        ax.yaxis.set_label_position('right')
      pl.ylabel(get_label(field), fontsize = 46)
      locs, labels = pl.yticks()
      locs = np.array(locs)
      if field ==  'cs':
        idx = np.where(locs % 10 == 0)[0]
      else:
        idx = np.where(locs % 1 == 0)[0]
      pl.yticks(locs[idx])
    
      # x-axis
      if (index / ncols == nrows -1):
        pl.xlabel(get_label('radius'), fontsize=46)
#        locs, labels = pl.xticks()
#        locs = np.array(locs)
#        idx = np.where(locs % 1 == 0)[0]
#        pl.xticks(locs[idx])
        locs = np.arange(xmin, xmax, dtype=int)
        if len(locs) > 6:
          locs = locs[np.where(locs % 2 == 0)[0]]
        pl.xticks(locs)
      else:
        ax.set_xticklabels([])

      pl.xlim(xmin, xmax)
      ax.tick_params(axis='both', labelsize=38)
    
    pl.subplots_adjust(wspace=0.0, hspace=0.0)
    self.fig.set_size_inches(ncols*7, nrows*7)
    del MyRadial

  def multi_pspace(self, snap, fields):
    MyImage = Image()
    nfields = len(fields)
    ncols = int(np.ceil(np.sqrt(nfields)))
    nrows = nfields / ncols

    for index in range(nfields):
      field = fields[index]
      MyImage.pspace_plot(snap, field)
      ax = pl.subplot(nrows, ncols, index+1)
      im = ax.imshow(MyImage.img, cmap = 'jet', extent = [MyImage.xmin, MyImage.xmax, MyImage.ymin, MyImage.ymax], aspect = 'auto')
      ax.plot(MyImage.valx, MyImage.valy, 'k')

      # y-axis
      if field == 'temp':
        minvaly = np.min(MyImage.valy)
	ymin, ymax = pl.ylim()
        xmin, xmax = pl.xlim()
        pl.plot(np.array([xmin, xmax]), 2./3*np.array([xmin,xmax]) - 28./3, 'k--')
        pl.text(10,2, r'$T \propto n_{\rm H}^{2/3}$', fontsize=20)
        pl.ylim(1, ymax)
      if field == 'abH2':
        ymin, ymax = pl.ylim()
        pl.ylim(-20, ymax)
      if ((index+1) % ncols == 0):
        ax.yaxis.tick_right()
        ax.yaxis.set_ticks_position('both')
        ax.yaxis.set_label_position('right')
      pl.ylabel(get_label(field), fontsize=22)
      locs, labels = pl.yticks()
      locs = np.array(locs)
      idx = np.where(locs % 1 == 0)[0]
      pl.yticks(locs[idx], fontsize=20)

      # x-axis
      if (index / ncols == nrows -1):
        pl.xlabel(get_label('nh'), fontsize=22)
        locs, labels = pl.xticks()
        locs = np.array(locs)
        idx = np.where(locs % 1 == 0)[0]
        pl.xticks(locs[idx], fontsize=20)
      else:
        ax.set_xticklabels([])

    pl.subplots_adjust(top=0.9, wspace=0.0, hspace=0.0)
#    vmin, vmax = np.min(MyImage.img[MyImage.img != -np.inf]), np.max(MyImage.img)
    vmin, vmax = -20.5, -1.0
    cax = pl.axes([0.15, 0.93, 0.7, 0.025])
    locs = np.arange(vmin, vmax + 1).astype(int)
    locs = locs[np.where(locs % 2 == 0)[0]]
    cb = pl.colorbar(im, cax = cax, orientation='horizontal', ticks=locs)
    cax.set_xlabel(r'${\rm log}\left(M_{\rm bin}/M_{\rm tot}\right)$', fontsize=22, labelpad=-70)
    cb.ax.tick_params(labelsize=20)
    self.fig.set_size_inches(ncols*6, nrows*6)
    del MyImage

  def multi_snap(self, snaps, fields, base, sim, ImgWidth, PlotTime=False, PlotSize=False, CenterProto=None, PlotProto=False, Rotate=False):
    nsnaps = len(snaps)
    nfields = len(fields)
    if nfields == 1:
      ncols = np.ceil(np.sqrt(nsnaps))
      nrows = np.ceil(nsnaps / ncols)
      ncols, nrows = int(ncols), int(nrows)
      cbar_mode = 'single'
    else:
      nrows, ncols = nfields, nsnaps
      cbar_mode = 'edge'
    grid = ImageGrid(self.fig, 111, nrows_ncols = [nrows, ncols], axes_pad = 0.0, label_mode='L', share_all = False,
                     cbar_location = 'right', cbar_mode = cbar_mode, cbar_pad = '5%', cbar_size = '5%')
   
    vmin, vmax = 13, 20.5
 
    for idx_snap in range(nsnaps):
    
      snap = snaps[idx_snap]
      path = base+sim+'/snapdir_%03i/' %snap
      file = sim+'_%03i' %snap
      snapbase = path+file
    
      MySnap = Snap()
      MySnap.read_header(snapbase)
      MySnap.read_fields(snapbase)
      if CenterProto:
        if snap < 126:
          f = h5py.File('/scratch/02563/fbecerra/paha/protostars/'+file+'.hdf5', 'r')
          try:
            Center = np.array(f['Proto%i/Position' %CenterProto])
          except:
            f.close()
            continue
          f.close()
          ImgWidth = 20
        else:
          ImgWidth = 50
        MySnap.fields['x'] -= Center[0]
        MySnap.fields['y'] -= Center[1]
        MySnap.fields['z'] -= Center[2]
      else:
        MySnap.center_box()
      if Rotate:
        MySnap.rotate_box()
    
      for idx_field in range(nfields):
        field = fields[idx_field]
        MyImage = Image()
        MyImage.calculate_image(MySnap, field, ImgWidth=ImgWidth)
    
        idx_ax = nsnaps * idx_field + idx_snap
        # vmin, vmax = np.min(MyImage.img), np.max(MyImage.img)
        im = grid[idx_ax].imshow(MyImage.img, cmap = get_colormap(field), extent = [0, MyImage.xbins, 0, MyImage.ybins], vmin = vmin, vmax = vmax)
        if len(MySnap.new_fields['sinks']['id']) > 0:
          x_sinks = MyImage.xbins / 2. + MyImage.xbins * MySnap.new_fields['sinks']['x'] / MyImage.width
          y_sinks = MyImage.ybins / 2. + MyImage.ybins * MySnap.new_fields['sinks']['y'] / MyImage.height
          for cx, cy in zip(x_sinks, y_sinks):
            if ((cx > 0) & (cx < MyImage.xbins) & (cy > 0) & (cy < MyImage.ybins)):
              grid[idx_ax].plot(cx, cy, 'w+', ms=5, mew=3)
        if PlotProto:
          f = h5py.File('/scratch/02563/fbecerra/paha/protostars/'+file+'.hdf5', 'r')
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
              grid[idx_ax].plot(x_proto, y_proto, 'k+')
              circ = pl.Circle((x_proto, y_proto), radius=protostar_radius, edgecolor='k', fill=False)
              grid[idx_ax].add_patch(circ)
          f.close()
        if PlotTime:
          grid[idx_ax].text(MyImage.xbins/10., MyImage.ybins/10., '%.2f yr' %(MySnap.params['time']*UNIT_TIME/SEC_PER_YEAR), fontsize=36, color='w')
        if PlotSize:
          grid[idx_ax].text(MyImage.xbins*(1-2.5/10), MyImage.ybins*(1-1./10), '%i au' %ImgWidth, fontsize=36, color='w')
        grid[idx_ax].set_xticklabels([])
        grid[idx_ax].set_yticklabels([])
        if nfields == 1 and idx_snap == 0:
          cb = grid.cbar_axes[idx_ax].colorbar(im)
          cb.set_label_text(get_label(field), fontsize = nrows*17)
          cb.ax.tick_params(labelsize=nrows*14)
          #cb.set_clim()
        if ((idx_snap == nsnaps - 1) and (nfields != 1)):
          idx_cb = idx_ax - (idx_field+1) * (nsnaps-1)
          cb = grid.cbar_axes[idx_cb].colorbar(im)
          cb.set_label_text(get_label(field), fontsize = 24)
          cb.ax.tick_params(labelsize=18)
          #cb.set_clim()
    
    self.fig.set_size_inches(ncols*8, nrows*8)
    del MySnap, MyImage

  def zoom_snap2(self, snaps, field, sizes, base, sims):
    nsizes = len(sizes)
    ncols = np.ceil(np.sqrt(nsizes))
    nrows = np.ceil(nsizes / ncols)
    ncols, nrows = int(ncols), int(nrows)
    left, right, bottom, top = 0.1, 0.9, 0.1, 0.9
    width = (right - left)/ncols
    height = (top - bottom)/nrows

    for idx_size in range(nsizes):
      sim = sims[idx_size]
      snap = snaps[idx_size]
      size = sizes[idx_size]
      if idx_size == 0:
        read_snap = True
      else:
        if sims[idx_size -1] == sim and snaps[idx_size -1] == snap:
          read_snap = False
        else:
          read_snap = True
      path = base+sim+'/snapdir_%03i/' %snap
      file = sim+'_%03i' %snap
      snapbase = path+file
 
      if read_snap:
        MySnap = Snap()
        MySnap.read_header(snapbase)
        MySnap.read_fields(snapbase)
        MySnap.center_box()
        MySnap.rotate_box()
  
      MyImage = Image()
      MyImage.calculate_image(MySnap, field, ImgWidth=size)
  
      if idx_size == 0 or idx_size == 1 or idx_size == 2:
        idx_subplot = idx_size + 1
      elif idx_size == 3:
        idx_subplot = 6
      elif idx_size == 4:
        idx_subplot = 5
      elif idx_size == 5:
        idx_subplot = 4
  
      ax = self.fig.add_subplot(nrows, ncols, idx_subplot)
      ax.set(aspect=1)
      #vmin, vmax = np.min(MyImage.img), np.max(MyImage.img)
      if idx_size == 0:
        vmin, vmax = 3,7
      elif idx_size == 1:
        vmin, vmax = 6,9.5
      elif idx_size == 2:
        vmin, vmax = 8,12
      elif idx_size == 3:
        vmin, vmax = 10,14
      elif idx_size == 4:
        vmin, vmax = 12,16
      elif idx_size == 5:
        vmin, vmax = 14,19
      im = ax.imshow(MyImage.img, cmap = get_colormap(field), extent = [0, MyImage.xbins, 0, MyImage.ybins], vmin = vmin, vmax = vmax)
      if size > 1e4:
        size_label = str(size/2e5)+' pc'
      else:
        size_label = '%i au' %size
      ax.text(MyImage.xbins/10., MyImage.ybins/10., size_label, fontsize=36, color='w')
      ax.set_xticklabels([])
      ax.set_yticklabels([])
 
      # Zoom effect 
      if idx_size != nsizes - 1:
        dx, dy = float(sizes[idx_size+1]) / float(size) * MyImage.xbins, float(sizes[idx_size+1]) / float(size) * MyImage.ybins
        xi, xf = MyImage.xbins / 2. - dx / 2., MyImage.xbins / 2. + dx / 2.
        yi, yf = MyImage.ybins / 2. - dy / 2., MyImage.ybins / 2. + dy / 2.
        # Square
        square = pl.Rectangle((xi, yi), dx, dy, facecolor ='none', edgecolor='k')
        ax.add_patch(square)
        if idx_size == 0 or idx_size == 1:
          x1, y1 = np.array([xi, MyImage.xbins]), np.array([yi, 0])
          x2, y2 = np.array([xi, MyImage.xbins]), np.array([yf, MyImage.ybins])
        elif idx_size == 2:
          x1, y1 = np.array([xi, 0]), np.array([yf, 0])
          x2, y2 = np.array([xf, MyImage.xbins]), np.array([yf, 0])
        elif idx_size == 3 or idx_size ==4:
          x1, y1 = np.array([xf, 0]), np.array([yi, 0])
          x2, y2 = np.array([xf, 0]), np.array([yf, MyImage.ybins])
        # Lines connecting edges
        ax.plot(x1, y1, 'k')
        ax.plot(x2, y2, 'k')
        # Fixing limits
        ax.set_xlim(0, MyImage.xbins)
        ax.set_ylim(0, MyImage.ybins)

      # Colorbar
      if idx_size == 0 or idx_size == 1 or idx_size == 2:
        cax = pl.axes([left+(idx_subplot-1)*width+0.01, bottom+nrows*height+0.04, width-0.02, 0.02])
        cb = pl.colorbar(im, cax = cax, orientation='horizontal', ticks=np.arange(vmin, vmax + 1).astype(int)) 
        cax.set_xlabel(get_label(field), fontsize=30, labelpad=-85)
        cb.ax.tick_params(labelsize=24) 
      elif idx_size == 3 or idx_size == 4 or idx_size == 5:
        cax = pl.axes([left+(idx_subplot-4)*width+0.01, bottom - 0.07, width-0.02, 0.02])
        cb = pl.colorbar(im, cax = cax, orientation='horizontal', ticks=np.arange(vmin, vmax + 1).astype(int))
        cax.set_xlabel(get_label(field), fontsize=30, labelpad=-85)
        cb.ax.tick_params(labelsize=24)
 
    pl.subplots_adjust(left=left, right=right, bottom=bottom, top=top, wspace=0.0, hspace=0.0) 
    self.fig.set_size_inches(ncols*8, nrows*8)
    del MySnap, MyImage

  def zoom_snap(self, snap, fields, sizes, base, sim):
    path = base+sim+'/snapdir_%03i/' %snap
    file = sim+'_%03i' %snap
    snapbase = path+file

    MySnap = Snap()
    MySnap.read_header(snapbase)
    MySnap.read_fields(snapbase)
    MySnap.center_box()
    #MySnap.rotate_box()

    nsizes = len(sizes)
    nfields = len(fields)
    nrows, ncols = nfields, nsizes
    cbar_mode = 'edge'
    grid = ImageGrid(self.fig, 111, nrows_ncols = [nrows, ncols], axes_pad = 0.0, label_mode='all', share_all = False,
                     cbar_location = 'right', cbar_mode = cbar_mode, cbar_pad = '5%', cbar_size = '5%')

    for idx_size in range(nsizes):

      size = sizes[idx_size]
      for idx_field in range(nfields):
        field = fields[idx_field]
        MyImage = Image()
        MyImage.calculate_image(MySnap, field, ImgWidth=size)

        idx_ax = nsizes * idx_field + idx_size
        im = grid[idx_ax].imshow(MyImage.img, cmap = get_colormap(field), extent = [0, MyImage.xbins, 0, MyImage.ybins], vmin = np.min(MyImage.img), vmax = np.max(MyImage.img))
        grid[idx_ax].text(MyImage.xbins/10., MyImage.ybins/10., '%i au' %size, fontsize=36, color='w')
        grid[idx_ax].set_xticklabels([])
        grid[idx_ax].set_yticklabels([])
        if ((idx_ax + 1) % nsizes != 0 ):
          dx, dy = float(sizes[idx_size+1]) / float(size) * MyImage.xbins, float(sizes[idx_size+1]) / float(size) * MyImage.ybins
          xi, xf = MyImage.xbins / 2. - dx / 2., MyImage.xbins / 2. + dx / 2. 
          yi, yf = MyImage.ybins / 2. - dy / 2., MyImage.ybins / 2. + dy / 2.
          # Square
          square = pl.Rectangle((xi, yi), dx, dy, facecolor ='none', edgecolor='w')
          grid[idx_ax].add_patch(square)
          # Lines connecting edges
          grid[idx_ax].plot(np.array([xi, MyImage.xbins]), np.array([yi, 0]), 'w')
          grid[idx_ax].plot(np.array([xi, MyImage.xbins]), np.array([yf, MyImage.ybins]), 'w')
          # Fixing limits
          grid[idx_ax].set_xlim(0, MyImage.xbins)
          grid[idx_ax].set_ylim(0, MyImage.ybins)

#        if nfields == 1 and idx_size == 0:
#          cb = grid.cbar_axes[idx_ax].colorbar(im)
#          cb.set_label_text(get_label(field), fontsize = 52)
#          cb.ax.tick_params(labelsize=40)
#          #cb.set_clim()
#        if ((idx_size == nsizes - 1) and (nfields != 1)):
#          idx_cb = idx_ax - (idx_field+1) * (nsizes-1)
#          cb = grid.cbar_axes[idx_cb].colorbar(im)
#          cb.set_label_text(get_label(field), fontsize = 24)
#          cb.ax.tick_params(labelsize=18)
#          #cb.set_clim()

    self.fig.set_size_inches(ncols*8, nrows*8)
    del MySnap, MyImage

  def proto_profiles(self, base, sim):
    mass = {}
    time = {}
    radius = {}
    rho = {}
    vrad = {}
    mass_flux = {}
    keys = np.array([])
    snaps = range(1, 182)

    for snap in snaps:

      print snap

      path = base+sim+'/snapdir_%03i/' %snap
      file = sim+'_%03i' %snap
      snapbase = path+file

      MySnap = Snap()
      MySnap.read_header(snapbase)
      MySnap.read_fields(snapbase)
      snaptime = MySnap.params['time']
      f = h5py.File('/scratch/02563/fbecerra/paha/protostars/'+file+'.hdf5', 'r')
      for key in f.keys():
        if key not in keys:
          keys = np.append(keys, key)
          mass[key] = np.array([])
          time[key] = np.array([])
          radius[key] = np.array([])
          rho[key] = np.array([])
          vrad[key] = np.array([])
          mass_flux[key] = np.array([])
        mass[key] = np.append(mass[key], np.array(f[key+'/Mass']))
        radius[key] = np.append(radius[key], np.array(f[key+'/Radius']))
        time[key] = np.append(time[key], snaptime)

        rho[key] = np.append(rho[key], np.array(f[key+'/Rho']))
        vrad[key] = np.append(vrad[key], np.array(f[key+'/Vrad']))
        
      f.close()

    for key in keys:
      # Check units
      mass_flux[key] = -4 * np.pi * (radius[key] * ASTRONOMICAL_UNIT)**2 * rho[key] * vrad[key] * 1e5 * SEC_PER_YEAR / SOLAR_MASS
      radius[key] = radius[key] * ASTRONOMICAL_UNIT / SOLAR_RADIUS
      time[key] = time[key] * UNIT_TIME / SEC_PER_YEAR

      print key, time[key]

      # Plots
      ax = pl.subplot(311)
      ax.plot(time[key], np.log10(mass[key]))
      print key, mass[key]
      ax.set_xticklabels([])
      pl.ylabel(r'${\rm log}\left({\rm M}_\star/{\rm M}_\odot\right)$')
      if key == 'Proto1':
        pl.xlim(np.min(time[key]), np.max(time[key]))

      ax = pl.subplot(312)
      ax.plot(time[key], np.log10(radius[key]))
      print key, radius[key]
      ax.set_xticklabels([])
      pl.ylabel(r'${\rm log}\left({\rm R}_\star/{\rm R}_\odot\right)$')
      if key == 'Proto1':
        pl.xlim(np.min(time[key]), np.max(time[key]))

      ax = pl.subplot(313)
      ax.plot(time[key], mass_flux[key])
      print key, mass_flux[key]
      pl.xlabel(r'Time [yr]')
      pl.ylabel(r'${\rm dM}_\star / {\rm dt} [{\rm M}_\odot\,{\rm yr}^{-1}]$')
      if key == 'Proto1':
        pl.xlim(np.min(time[key]), np.max(time[key]))

    self.fig.set_size_inches(8,24)

  def savefig(self, label):
    self.fig.savefig(label, dpi=100)
    pl.show()
    pl.close()

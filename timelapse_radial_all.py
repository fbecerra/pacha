import matplotlib
matplotlib.use('Agg')
import src as pr
import numpy as np
import pylab as pl
import time

def get_label(field):
  if field == 'vrad':
    label = r'$|v_{\rm rad}|$'
  elif field == 'vrot':
    label = r'$v_{\rm rot}$'
  elif field == 'vturb':
    label = r'$v_{\rm turb}$'
  elif field == 'vkep':
    label = r'$v_{\rm kep}$'
  elif field == 'cs':
    label = r'$c_{\rm s}$'
  if field == 'vrad_ratio':
    label = r'$|v_{\rm rad}|/c_{\rm s}$'
  elif field == 'vrot_ratio':
    label = r'$v_{\rm rot}/c_{\rm s}$'
  elif field == 'vturb_ratio':
    label = r'$v_{\rm turb}/c_{\rm s}$'
  elif field == 'vkep_ratio':
    label = r'$v_{\rm kep}/c_{\rm s}$'
  return label

sims = [['ah1w3'],
        ['ah1w3'],
        ['ah1w3'],
        ['ah1w3r0', 'ah1w3'],
        ['ah1w3f2', 'ah1w3r0', 'ah1w3'],
        ['ah1w3f2', 'ah1w3r0', 'ah1w3'],
        ['ah1w3f2', 'ah1w3r0', 'ah1w3']]
snaps = [['002'],
         ['035'],
         ['153'],
         ['015', '157'],
         ['001', '022', '157'],
         ['065', '022', '157'],
         ['181', '022', '157']]

linestyles = [#'solid',
              [8, 4],
              [8, 4, 2, 4],
              [2, 4],
              [8, 4],
              [8, 4, 2, 4],
              [2, 4],
              #[8, 4, 2, 4, 2, 4],
              #[16, 4],
              #[16, 4, 2, 4]]
              'solid']

nsims = len(sims)

base= '/scratch/00025/tgreif/'
lookback_times = np.zeros(nsims, float)

plot_fields = [['nh', 'enc_mass', 'temp', 'abH2', 'abHM', 'abHII'],
               ['mbe_ratio']]
#               ['geff', 'dnh', 'radial_tcool_ratio', 'tcs_ratio']]
all_fields = [item for sublist in plot_fields for item in sublist]

start = time.time()

for sim_idx in range(nsims):
  sim = sims[sim_idx]
  snap = snaps[sim_idx]

  for index in range(len(sim)):
  
    path = base+sim[index]+'/snapdir_'+snap[index]+'/'
    file = sim[index]+'_'+snap[index]
    snapbase = path+file
    
    if index == 0:
      MySnap = pr.snap.Snap()
      MySnap.read_header(snapbase)
      MySnap.read_fields(snapbase)
      MySnap.center_box()
      print MySnap.Center
      MySnap.rotate_box()
      print MySnap.params['time']
      lookback_times[sim_idx] += MySnap.params['time']
    else:
      MySnap2 = pr.snap.Snap()
      MySnap2.read_header(snapbase)
      MySnap2.read_fields(snapbase)
      MySnap2.center_box()
      print MySnap2.Center
      MySnap2.rotate_box()
      print MySnap.params['time']
      MySnap.merge_snaps(MySnap2)
      lookback_times[sim_idx] += MySnap2.params['time']
      del MySnap2
  
  end = time.time()
  print 'Total time: %f' %(end - start)
  MySnap.calculate_radius()
  
  MyRadial = pr.radial.Radial()
  MyRadial.radial_profile(MySnap, all_fields)
  #MyRadial.radial_profile(MySnap, ['temp', 'rho', 'mu', 'gamma', 'enc_mass'])
  #print MyRadial.radial['temp']
  #print MyRadial.radial['rho']
  #print MyRadial.radial['mu']
  #print MyRadial.radial['gamma']
  #print MyRadial.radial['enc_mass']
  del MySnap
  
  for idx in range(len(plot_fields)):
    fig = pl.figure(idx + 1)
    fields = plot_fields[idx]
    nfields = len(fields)
    #ncols = int(np.ceil(np.sqrt(nfields)))
    if nfields > 1:
      ncols = 2
    else:
      ncols = 1
    nrows = nfields / ncols
  
    for index in range(nfields):
      field = fields[index]
      ax = pl.subplot(nrows, ncols,index+1)
      if field ==  'cs':
        MyRadial.radial[field] = 10**MyRadial.radial[field]
      if field == 'mbe_ratio':
        xval = 'enc_mass'
      else:
        xval = 'radius'
      line, = pl.plot(MyRadial.radial[xval], MyRadial.radial[field])
      line.set_dashes(linestyles[sim_idx])
      xmin, xmax = np.min(MyRadial.radial[xval]), np.max(MyRadial.radial[xval])
 
      if sim_idx == nsims-1: 
        if field == 'tcool_ratio' or field == 'mbe_ratio' or field == 'q' or field == 'radial_tcool_ratio' or field == 'tcs_ratio':
          pl.plot(np.array([xmin, xmax]), np.array([0,0]), 'k--')
        if field == 'nh':
          ymin, ymax = pl.ylim()
          pl.plot(np.array([xmin, xmax]), -2*np.array([xmin,xmax]) + 20, 'k--')
          pl.text(4,13, r'$n_{\rm H} \propto r^{-2}$', fontsize=24)
          pl.ylim(ymin, 23)
    
        # y-axis
        ymin, ymax = pl.ylim()
        if field == 'radial_tcool_ratio' or fields == 'geff':
          pl.ylim(-2, 4)
        if fields == 'temp':
          pl.ylim(1, ymax)
        if fields == 'abH2':
          pl.ylim(ymin, -3)
        if fields == 'abHM':
          pl.ylim(ymin, -2)
        if fields == 'mbe_ratio':
          pl.ylim(-5.8, 6.5)
        if (ncols > 1 and (index+1) % ncols == 0):
          ax.yaxis.tick_right()
          ax.yaxis.set_ticks_position('both')
          ax.yaxis.set_label_position('right')
        pl.ylabel(pr.utils.get_label(field), fontsize=24)
        locs, labels = pl.yticks()
        locs = np.array(locs)
        if field ==  'cs':
          ticks_idx = np.where(locs % 10 == 0)[0]
        else:
          ticks_idx = np.where(locs % 1 == 0)[0]
        pl.yticks(locs[ticks_idx])
    
        # x-axis
        if (index / ncols == nrows -1):
          pl.xlabel(pr.utils.get_label(xval), fontsize=24)
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
    
        ax.tick_params(axis='both', labelsize=20)
        pl.xlim(-2, xmax)
  
    if sim_idx == nsims -1:
      ax = pl.subplot(nrows, ncols, 1)
      if idx == 0:
        lookback_times *= pr.constants.UNIT_TIME/pr.constants.SEC_PER_YEAR
        lookback_times = np.abs(lookback_times - lookback_times[sim_idx])
        legends = ['%.2e yr' %lookback_time for lookback_time in lookback_times]
      if idx == 0:
        pl.legend(legends, loc = 'lower left', prop={'size':14}, frameon=True)
      elif idx == 1:
        pl.legend(legends, loc = 'upper left', prop={'size':14}, frameon=True)
      pl.subplots_adjust(wspace=0.0, hspace=0.0)
      fig.set_size_inches(ncols*7, nrows*7)
  
  del MyRadial


for idx in range(len(plot_fields)):
  fig = pl.figure(idx + 1)
  fig.savefig('/scratch/02563/fbecerra/paha/radial/timelapse_multi_radial_all_'+'_'.join(plot_fields[idx])+'.pdf', dpi=100)

end = time.time()
print 'Total time: %f' %(end - start)

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

sims = ['ah1w3f2', 'ah1w3r0', 'ah1w3']
snaps = ['181', '022', '157']

base= '/scratch/00025/tgreif/'

start = time.time()
for index in range(len(sims)):

  path = base+sims[index]+'/snapdir_'+snaps[index]+'/'
  file = sims[index]+'_'+snaps[index]
  snapbase = path+file
  
  if index == 0:
    MySnap = pr.snap.Snap()
    MySnap.read_header(snapbase)
    MySnap.read_fields(snapbase)
    MySnap.center_box()
    print MySnap.Center
    MySnap.rotate_box()
  else:
    MySnap2 = pr.snap.Snap()
    MySnap2.read_header(snapbase)
    MySnap2.read_fields(snapbase)
    MySnap2.center_box()
    print MySnap2.Center
    MySnap2.rotate_box()
    MySnap.merge_snaps(MySnap2)
    del MySnap2

end = time.time()
print 'Total time: %f' %(end - start)
MySnap.calculate_radius()

#plot_fields = ['nh', 'temp', 'abH2', 'abHII', 'esc_frac', 'cool_rate']
#pc = pr.plot_collection.PlotCollection()
#pc.multi_radial(MySnap, plot_fields)
#pc.savefig('/scratch/02563/fbecerra/paha/radial/multi_radial_all_'+snaps[0]+'_'.join(plot_fields)+'.pdf')

## nH, M, T, yH2, yHI, yHII
#plot_fields = ['nh', 'enc_mass', 'temp', 'abH2', 'abHM', 'abHII']
#pc = pr.plot_collection.PlotCollection()
#pc.multi_radial(MySnap, plot_fields)
#pc.savefig('/scratch/02563/fbecerra/paha/radial/multi_radial_all_'+'_'.join(plot_fields)+'.pdf')

## Phase space
#plot_fields = ['temp', 'abH2', 'abHM', 'abHII']
#pc = pr.plot_collection.PlotCollection()
#pc.multi_pspace(MySnap, plot_fields)
#pc.savefig('/scratch/02563/fbecerra/paha/pspace/multi_pspace_all_'+'_'.join(plot_fields)+'.pdf')

# Velocities
plot_fields = ['vrot', 'vrad', 'vkep', 'vturb', 'cs']
plot_fields_right = ['vrad', 'vrot', 'vkep', 'vturb']
plot_fields_left = ['vrad', 'vrot', 'vkep', 'vturb', 'cs']
linestyles = ['solid',
              [8, 4],
              [8, 4, 2, 4],
              [2, 4],
              [8, 4, 2, 4, 2, 4],
              [16, 4],
              [16, 4, 2, 4]]
MyRadial = pr.radial.Radial()
MyRadial.radial_profile(MySnap, plot_fields)
fig = pl.figure()

# Right
ax = pl.subplot(122)
cs = 10**MyRadial.radial['cs']
for idx in range(len(plot_fields_right)):
  field_plot = plot_fields_left[idx]
  if field_plot == 'vrad':
    y = np.abs(MyRadial.radial[field_plot])
  else:
    y = MyRadial.radial[field_plot]
  line, = pl.plot(MyRadial.radial['radius'], np.log10(y/cs), label= get_label(field_plot+'_ratio'))
  line.set_dashes(linestyles[idx])
ax.yaxis.tick_right()
ax.yaxis.set_ticks_position('both')
ax.yaxis.set_label_position('right')
pl.xlim(np.min(MyRadial.radial['radius']), np.max(MyRadial.radial['radius']))
pl.ylim(-2, 2)
pl.plot(np.array([np.min(MyRadial.radial['radius']), np.max(MyRadial.radial['radius'])]), np.array([0,0]), 'k--')
pl.xlabel(r'${\rm log}\left(r/{\rm au}\right)$', fontsize=24)
pl.ylabel(r'${\rm log}\,\mathcal{M}$', fontsize=24)
pl.legend(loc='lower center', prop={'size':18}, frameon=True)
pl.tick_params(labelsize=18)

# Left
ax = pl.subplot(121)
for idx in range(len(plot_fields_left)):
  field_plot = plot_fields_left[idx]
  if field_plot == 'vrad':
    y = np.abs(MyRadial.radial[field_plot])
  elif field_plot == 'cs':
    y = cs
  else:
    y = MyRadial.radial[field_plot]
  line, = pl.plot(MyRadial.radial['radius'], np.log10(y), label= get_label(field_plot))
  line.set_dashes(linestyles[idx])
pl.xlim(np.min(MyRadial.radial['radius']), np.max(MyRadial.radial['radius']))
ymin, ymax = pl.ylim()
pl.ylim(-1.5, ymax)
pl.xlabel(r'${\rm log}\left(r/{\rm au}\right)$', fontsize=24)
pl.ylabel(r'${\rm log}\left(v/{\rm km\,s}^{-1}\right)$', fontsize=24)
pl.legend(loc='lower center',prop={'size':18}, frameon=True)

pl.subplots_adjust(wspace=0.0, hspace=0.0)
fig.set_size_inches(2*7, 7)
pl.tick_params(labelsize=18)
pl.savefig('/scratch/02563/fbecerra/paha/radial/velocities_all.pdf', dpi=100)
pl.show()
pl.close()

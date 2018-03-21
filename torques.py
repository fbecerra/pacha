import matplotlib
matplotlib.use('Agg')
import src as pr
import numpy as np
import pylab as pl
import time
from itertools import cycle

def get_label(sim):
  if sim == 'nahw1r4sm1' or sim == 'nahw1r4sm2' or sim == 'nahw1r4sm3':
    label1 = 'sink'
  elif sim == 'nahw1r4ad1' or sim == 'nahw1r4ad2' or sim == 'nahw1r4ad3':
    label1 = 'core'

  if sim == 'nahw1r4sm3' or sim == 'nahw1r4ad3':
    label2 = r'$n_{\rm th} = 10^{12}\,{\rm cm}^{-3}$'
  elif sim == 'nahw1r4sm2' or sim == 'nahw1r4ad2':
    label2 = r'$n_{\rm th} = 10^{10}\,{\rm cm}^{-3}$'
  elif sim == 'nahw1r4sm1' or sim == 'nahw1r4ad1':
    label2 = r'$n_{\rm th} = 10^8\,{\rm cm}^{-3}$'

  #return label1+', '+label2

  if sim == 'nahw1r4sm1' or sim == 'nahw1r4sm2' or sim == 'nahw1r4sm3':
    return label2
  elif sim == 'nahw1r4sm1' or sim == 'nahw1r4sm2' or sim == 'nahw1r4sm3':
    return ''

def get_color(sim):
  if sim == 'nahw1r4sm3' or sim == 'nahw1r4ad3':
    color = '#4daf4a'
  elif sim == 'nahw1r4sm2' or sim == 'nahw1r4ad2':
    color = '#377eb8'
  elif sim == 'nahw1r4sm1' or sim == 'nahw1r4ad1':
    color = '#e41a1c'
  return color

def get_linestyle(sim):
  if sim == 'nahw1r4sm1' or sim == 'nahw1r4sm2' or sim == 'nahw1r4sm3':
    linestyle = '-'
  elif sim == 'nahw1r4ad1' or sim == 'nahw1r4ad2' or sim == 'nahw1r4ad3':
    linestyle = '--'
  return linestyle

def smooth(x,window_len=6,window='hanning'):
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

fontsize = 24
colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf']
colors = ['#08306b', '#08519c', '#2171b5', '#4292c6', '#6baed6', '#9ecae1', '#c6dbef']
colors = ['#0868ac', '#2b8cbe', '#4eb3d3', '#7bccc4', '#a8ddb5', '#ccebc5', '#e0f3db']
colors_sm = ['#1f78b4', '#33a02c', '#e31a1c', '#ff7f00', '#6a3d9a', '#b15928']
colors_ad = ['#a6cee3', '#b2df8a', '#fb9a99', '#fdbf6f', '#cab2d6', '#ffff99']

gray = '#999999'
linestyles = ['-', '--']

base= '/n/hernquistfs2/fbecerra/'

#sims, snaps, xmin, xmax = ['nahw1r4sm1', 'nahw1r4ad1'], [[3], [2]], -2.5, 0
#sims, snaps, xmin, xmax, lmin, lmax = ['nahw1r4sm1', 'nahw1r4ad1'], [[3, 6, 7, 11, 16, 20], [2, 13, 31]], -2.5, 0, 20, 25
#sims, snaps, xmin, xmax, lmin, lmax = ['nahw1r4sm2', 'nahw1r4ad2'], [[8, 9, 10, 13, 16], [9, 11, 23]], -3.5, -1, 19, 24
#sims, snaps, xmin, xmax, lmin, lmax = ['nahw1r4sm3', 'nahw1r4ad3'], [[14, 16, 17, 20, 23, 26], [14, 19, 26]], -4.5, -2, 18, 23

#sims = ['nahw1r4ad1', 'nahw1r4ad2', 'nahw1r4ad3']
#snaps = [[1, 31], [7, 23], [14, 32]]

sims = ['nahw1r4sm1', 'nahw1r4sm2', 'nahw1r4sm3', 'nahw1r4ad1', 'nahw1r4ad2', 'nahw1r4ad3']
snaps = [20, 16, 26, 31, 23, 32]
xmin, xmax = -4.5, 1
lmin, lmax = 19, 25

#snaps = range(snap_init-1, snap_init)
fields = ['angmom', 'taugrav', 'taupres', 'tff', 'tgrav', 'tpres']
fig = pl.figure()

fig, ax = pl.subplots(3, 2, figsize=(16, 24))

for idx, sim in enumerate(sims):

#  for idx_snap, snap in enumerate(snaps[idx]):
    snap = snaps[idx]
    print sim, snap
    path = base+sim+'/snapdir_%03d/' %snap
    file = sim+'_%03d' %snap
    snapbase = path+file
  
    MySnap = pr.snap.Snap()
    MySnap.read_header(snapbase)
    MySnap.read_fields(snapbase)
    MySnap.center_box()
    MySnap.rotate_box()
    MySnap.calculate_radius()
  
    MyRadial = pr.radial.Radial()
    MyRadial.radial_profile(MySnap, fields)

    radius = MyRadial.radial['radius']
    time = MySnap.params['time']*pr.constants.UNIT_TIME/pr.constants.SEC_PER_YEAR
    tff = 10**MyRadial.radial['tff']

    #color = colors[idx_snap]
    #linestyle = linestyles[idx]

    ax[0,0].plot(radius, MyRadial.radial['taugrav'], linestyle=get_linestyle(sim), color=get_color(sim), label=get_label(sim))
    ax[0,1].plot(radius, MyRadial.radial['taupres'], linestyle=get_linestyle(sim), color=get_color(sim), label=get_label(sim))
    ax[1,0].plot(radius, np.log10(10**MyRadial.radial['taugrav']/10**MyRadial.radial['taupres']), linestyle=get_linestyle(sim), color=get_color(sim), label=get_label(sim))
    ax[1,1].plot(radius, MyRadial.radial['angmom'], linestyle=get_linestyle(sim), color=get_color(sim), label=get_label(sim))
    ax[2,0].plot(radius, smooth(np.log10(np.abs(MyRadial.radial['tgrav'])/tff)), linestyle=get_linestyle(sim), color=get_color(sim), label=get_label(sim))
    ax[2,1].plot(radius, smooth(np.log10(np.abs(MyRadial.radial['tpres'])/tff)), linestyle=get_linestyle(sim), color=get_color(sim), label=get_label(sim))

#xmin, xmax = np.min(MyRadial.radial['enc_mass']), np.max(MyRadial.radial['enc_mass'])

ax[1,1].legend(loc='lower right', prop={'size':20})

# y-axis
ax[0,0].set_ylabel(r'$\log \left(\tau_{\rm grav}/{\rm cm}^2\,{\rm s}^2\right)$', fontsize=fontsize)
ax[0,1].set_ylabel(r'$\log \left(\tau_{\rm pres}/{\rm cm}^2\,{\rm s}^2\right)$', fontsize=fontsize)
ax[1,0].set_ylabel(r'$\log (\tau_{\rm grav}/\tau_{\rm pres})$', fontsize=fontsize)
ax[1,1].set_ylabel(r'$\log \left(l/{\rm cm}^2\,{\rm s}\right)$', fontsize=fontsize)
#ax[1].set_ylabel(r'$t_{\rm grav}/t_{\rm pres}$')
ax[2,0].set_ylabel(r'$\log \left(t_{\rm grav}/t_{\rm ff}\right)$', fontsize=fontsize)
ax[2,1].set_ylabel(r'$\log \left(t_{\rm pres}/t_{\rm ff}\right)$', fontsize=fontsize)

# x-axis
ax[0,0].set_xlabel(pr.utils.get_label('radius'), fontsize=fontsize)
ax[0,1].set_xlabel(pr.utils.get_label('radius'), fontsize=fontsize)
ax[1,0].set_xlabel(pr.utils.get_label('radius'), fontsize=fontsize)
ax[1,1].set_xlabel(pr.utils.get_label('radius'), fontsize=fontsize)
ax[2,0].set_xlabel(pr.utils.get_label('radius'), fontsize=fontsize)
ax[2,1].set_xlabel(pr.utils.get_label('radius'), fontsize=fontsize)

ax[0,0].set_xlim(xmin, xmax)
ax[0,1].set_xlim(xmin, xmax)
ax[1,0].set_xlim(xmin, xmax)
ax[1,1].set_xlim(xmin, xmax)
ax[2,0].set_xlim(xmin, xmax)
ax[2,1].set_xlim(xmin, xmax)

ax[0,0].set_ylim(9.5, 13.5)
ax[0,1].set_ylim(9.5, 13.5)
ax[1,0].set_ylim(-1.5, 2)
ax[1,1].set_ylim(lmin, lmax)
ax[2,0].set_ylim(-2, 3)
ax[2,1].set_ylim(-2, 3)

ax[1,0].plot(np.array([xmin, xmax]), np.array([0,0]), ':', color=gray)
ax[2,0].plot(np.array([xmin, xmax]), np.array([0,0]), ':', color=gray)
ax[2,1].plot(np.array([xmin, xmax]), np.array([0,0]), ':', color=gray)

fig.savefig('./radial/radial_timelapse_'+sim+'_'+'_'.join(fields)+'.pdf')
pl.show()
pl.close()

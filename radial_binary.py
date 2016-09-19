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

base= '/scratch/00025/tgreif/'
sim = 'ah1w3f2'
snaps = range(1,26)
for snap in snaps:
  print snap
  path = base+sim+'/snapdir_'+snap+'/'
  file = sim+'_'+snap
  snapbase = path+file

  MySnap = pr.snap.Snap()
  MySnap.read_header(snapbase)
  MySnap.read_fields(snapbase)
  MySnap.center_box()
  MySnap.rotate_box()
  MySnap.calculate_radius()

  plot_fields = ['geff', 'dnh', 'radial_tcool_ratio', 'tcs_ratio']
  pc = pr.plot_collection.PlotCollection()
  pc.multi_radial(MySnap, plot_fields)
  pc.savefig('/scratch/02563/fbecerra/paha/radial/radial_'+file+'_'.join(plot_fields)+'.pdf')

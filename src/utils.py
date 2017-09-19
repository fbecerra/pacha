from include import *

def center_of_mass(Snap, idx):
  # Returns center of mass
  rho = Snap.fields['rho'][idx]
  totrho = np.sum(rho)
  cmx = np.sum(rho * Snap.fields['x'][idx]) / totrho
  cmy = np.sum(rho * Snap.fields['y'][idx]) / totrho
  cmz = np.sum(rho * Snap.fields['z'][idx]) / totrho
  del rho, totrho
  return np.array([cmx, cmy, cmz])

def vel_center_of_mass(Snap, idx):
  # Returns ivelocity of center of mass
  rho = Snap.fields['rho'][idx]
  totrho = np.sum(rho)
  vcmx = np.sum(rho * Snap.fields['vx'][idx]) / totrho
  vcmy = np.sum(rho * Snap.fields['vy'][idx]) / totrho
  vcmz = np.sum(rho * Snap.fields['vz'][idx]) / totrho
  del rho, totrho
  return np.array([vcmx, vcmy, vcmz])

def total_mass(Snap, idx):
  # Returns total mass
  return np.sum(Snap.fields['mass'][idx])

def get_colormap(field):
  if field == 'nh':
    return viridis
  elif field == 'temp':
    return matplotlib.cm.gist_heat
  elif field == 'gravacc':
    return matplotlib.cm.jet
  elif field == 'gradp':
    return matplotlib.cm.jet
  elif field == 'abhm':
    return matplotlib.cm.gist_stern
  elif field == 'abH2':
    return matplotlib.cm.gist_stern
  elif field == 'abhii':
    return matplotlib.cm.get_cmap('own1')
  elif field == 'abdii':
    return matplotlib.cm.gist_stern
  elif field == 'abhd':
    return matplotlib.cm.gist_stern
  elif field == 'abheii':
    return matplotlib.cm.gist_stern
  elif field == 'abheiii':
    return matplotlib.cm.gist_stern
  elif field == 'gamma':
    return matplotlib.cm.jet
  elif field == 'pdvrate':
    return matplotlib.cm.gist_stern
  elif field == 'h2rate':
    return matplotlib.cm.gist_stern
  elif field == 'cierate':
    return matplotlib.cm.gist_stern
  elif field == 'chemrate':
    return matplotlib.cm.gist_stern
  elif field == 'esc_frac':
    return matplotlib.cm.gist_stern
  elif field == 'allowref':
    return matplotlib.cm.gist_stern
  elif field == 'divvel':
    return matplotlib.cm.gist_stern
  elif field == 'cool':
    return matplotlib.cm.gist_stern
  elif field == 'collapse':
    return matplotlib.cm.get_cmap('own1')

def get_label(field):
  if field == 'x':
    if LengthUnit == 0:
      label = r'$x\,[{\rm pc}]$'
    elif LengthUnit == 1:
      label = r'$x\,[{\rm kpc}]$'
    else:
      label = r'$x\,[{\rm AU}]$'
  elif field == 'y':
    if LengthUnit == 0:
      label = r'$y\,[{\rm pc}]$'
    elif LengthUnit == 1:
      label = r'$y\,[{\rm kpc}]$'
    else:
      label = r'$y\,[{\rm AU}]$'
  elif field == 'z':
    if LengthUnit == 0:
      label = r'$z\,[{\rm pc}]$'
    elif LengthUnit == 1:
      label = r'$z\,[{\rm kpc}]$'
    else: 
      label = r'$z\,[{\rm AU}]$'
  elif field == 'radius':
    if LengthUnit == 0:
      label = r'${\rm log}\left(r/{\rm pc}\right)$'
    elif LengthUnit == 1:
      label = r'${\rm log}\left(r/{\rm kpc}\right)$'
    else:
      label = r'${\rm log}\left(r/{\rm au}\right)$'
  elif field == 'enc_mass':
    label = r'${\rm log}\left(M_{\rm enc}/{\rm M}_\odot\right)$'
  elif field == 'mbe_ratio':
    label = r'${\rm log}\left(M_{\rm enc}/M_{\rm BE}\right)$'
  elif field == 'acc_ratio':
    label = r'${\rm log}\left(t_{\rm ff}/t_{\rm acc}\right)$'
  elif field == 'vrad':
    label = r'$v_{\rm rad}\,[{\rm km}\,{\rm s}^{-1}]$'
  elif field == 'frac_rad':
    label = r'$v_{\rm rad}/c_{\rm s}$'
  elif field == 'vrot':
    label = r'$v_{\rm rot}\,[{\rm km}\,{\rm s}^{-1}]$'
  elif field == 'frac_rot':
    label = r'$v_{\rm rot}/c_{\rm s}$'
  elif field == 'vturb':
    label = r'$v_{\rm turb}\,[{\rm km}\,{\rm s}^{-1}]$'
  elif field == 'frac_turb':
    label = r'$v_{\rm turb}/\left|v_{\rm rad}\right|$'
  elif field == 'vkep':
    label = r'$v_{\rm kep}\,[{\rm km}\,{\rm s}^{-1}]$'
  elif field == 'kep_ratio':
    label = r'$v_{\rm rot}/v_{\rm kep}$'
  elif field == 'mach':
    label = r'$\mathcal{M}$'
  elif field == 'angmom':
    label == r'$l\,[?]'
  elif field == 'taugrav':
    label == r'$\tau_{\rm grav}\,[?]'
  elif field == 'taupres':
    label == r'$\tau_{\rm pres}\,[?]'
  elif field == 'nh':
    label = r'${\rm log}\left(n_{\rm H}/{\rm cm}^{-3}\right)$'
  elif field == 'dnh':
    label = r'${\rm log}\,\sigma_\delta$'
  elif field == 'temp':
    label = r'${\rm log}\left(T/{\rm K}\right)$'
  elif field == 'gravacc':
    label = r'$a_{\rm grav}$'
  elif field == 'gradp':
    label = r'$\nabla P$'
  elif field == 'abHM':
    label = r'${\rm log}\,{\rm y}_{{\rm H}^-}$'
  elif field == 'abH2':
    label = r'${\rm log}\,{\rm y}_{{\rm H}_2}$'
  elif field == 'abHII':
    label = r'${\rm log}\,{\rm y}_{\rm HII}$'
  elif field == 'abdii':
    label = r'${\rm log}\,{\rm y}_{\rm DII}$'
  elif field == 'abhd':
    label = r'${\rm log}\,{\rm y}_{\rm HD}$'
  elif field == 'abheii':
    label = r'${\rm log}\,{\rm y}_{\rm HeII}$'
  elif field == 'abheiii':
    label = r'${\rm log}\,{\rm y}_{\rm HeIII}$'
  elif field == 'gamma':
    label = r'$\gamma$'
  elif field == 'geff':
    label = r'$\gamma_{\rm eff}$'
  elif field == 'esc_frac' or field == 'radial_esc_frac':
    label = r'${\rm log}\,f_{\rm esc}$'
  elif field == 'tcool_ratio' or field =='radial_tcool_ratio':
    label = r'${\rm log}\left(t_{\rm cool}/t_{\rm ff}\right)$'
  elif field == 'tcs_ratio':
    label = r'${\rm log}\left(t_{\rm ff}/t_{\rm cs}\right)$'
  elif field == 'tcool':
    label = r'${\rm log}\,t_{\rm cool}$'
  elif field == 'tff':
    label = r'${\rm log}\,t_{\rm ff}$'
  elif field == 'collapse':
    label = r'${\rm log}\left(t_{\rm ff}/t_{\rm c_s}\right)$'
  elif field == 'sigma_gas':
    label = r'${\rm log}\left(\Sigma/{\rm g}\,{\rm cm}^{-2}\right)$'
  elif field == 'cs':
    label = r'${\rm c}_{\rm s}\,[{\rm km}\,{\rm s}^{-1}]$'
  elif field == 'omega_rot':
    label = r'${\rm log}\left(\Omega/{\rm yr}^{-1}\right)$'
  elif field == 'q':
    label = r'${\rm log}\,{\rm Q}$'
  elif field == 'gammie':
    label = r'${\rm log}\left(t_{\rm cool}/3\Omega^{-1} \right)$'

  return label

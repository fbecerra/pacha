IO_POS		= 1   
IO_VEL		= 1
IO_ID		= 1
IO_MASS		= 1
IO_U		= 1
IO_RHO		= 1
IO_VOL		= 1
IO_GRADV	= 0
IO_DELAUNAY	= 0
IO_GRAVACC	= 1
IO_GRADP	= 1
IO_PHODENS      = 0
IO_CHEM		= 1
IO_GAMMA	= 1
IO_ALLOWREF     = 1
IO_DIVVEL	= 1

fields_list         = ['x', 'y', 'z', 'radius', 'vx', 'vy', 'vz', 'id', 'mass', 'u', 'temp', 'rho', 'nh', 'vol', 'hsml', 'delaunay', 
                       'gravaccx', 'gravaccy', 'gravaccz', 'gravacc', 'gradpx', 'gradpy', 'gradpz', 'abHM', 'abH2', 'abHII', 'mu', 'gamma', 'phodens']
derived_fields_list = ['tff', 'cool_rate', 'esc_frac', 'tcool_ratio', 'entro', 'kappa', 'cs', 'mbe', 'tcs_ratio', 'tcool']
radial_fields_list  = ['enc_mass', 'sigma_gas', 'geff', 'dnh', 'vrad', 'vrot', 'vturb', 'omega_rot', 'vkep', 'q', 'radial_esc_frac', 'radial_cool_rate', 'radial_tcool_ratio', 'mbe_ratio', 'gammie', 'infall', 'angmom', 'taugrav', 'taupres', 'tgrav', 'tpres', 'Nion'] 

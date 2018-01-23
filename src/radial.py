from include import *
from scipy import interpolate

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

class Radial:

  def __init__(self):
    self.radial = dict()

  def radial_profile(self, Snap, fields):
    radius = Snap.fields['radius']
    #logrmin = np.log10(RadMin)
    logrmin = np.log10(2*np.min(radius[radius != 0]))
    logrmax = np.log10(np.max(radius))
    rfac = RadialBins / (logrmax - logrmin)
    logr = np.log10(radius)
    r_idx = (rfac*(logr - logrmin)).astype(int)
    r_idx[logr == -np.inf] = -1

    self.radial['radius'] = logrmin + np.arange(RadialBins) / rfac

    for field in fields:
      # Define array
      self.radial[field] = np.zeros(RadialBins, float)

      # Field read from file
      if field in fields_list:
        try:
          Snap.fields[field][0]
        except:
          print 'Field '+field+' was not read!'
        else:
          for j in np.unique(r_idx[(r_idx >= 0) & (r_idx < RadialBins)]):
            idx = np.where(r_idx == j)[0]
            if field == 'mass':
              self.radial[field][j] = np.log10(np.sum(Snap.fields['mass'][idx]))
            else:
              self.radial[field][j] = np.sum(Snap.fields['mass'][idx]*np.log10(Snap.fields[field][idx])) / np.sum(Snap.fields['mass'][idx])
      # Field derived from basic fields
      elif field in derived_fields_list:
        try:
          Snap.derived_fields[field][0]
        except:
          Snap.calculate_fields(field)
        for j in np.unique(r_idx[(r_idx >= 0) & (r_idx < RadialBins)]):
          idx = np.where(r_idx == j)[0]
          self.radial[field][j] = np.sum(Snap.fields['mass'][idx]*np.log10(Snap.derived_fields[field][idx])) / np.sum(Snap.fields['mass'][idx])
      # Field needs radial values
      elif field in radial_fields_list:
        # Enclosed mass
        if field == 'enc_mass':
          try:
            self.radial['mass'][0]
          except:
            self.radial_profile(Snap, ['mass'])
          self.radial[field] = np.log10(np.cumsum(10**self.radial['mass']))

        # Gas surface density
        if field == 'sigma_gas':
          # Convert units
          if LengthUnit == 0:
            fac = UNIT_LENGTH
          elif LengthUnit == 1:
            fac = UNIT_LENGTH / 1e3
          elif LengthUnit == 2:
            fac = ASTRONOMICAL_UNIT
          else:
            print 'Length unit not implemented, we will assume pc'
            fac = UNIT_LENGTH
          if not ComovingUnits:
            fac /= (1 + Snap.params['redshift'])
          fac /= Snap.params['hubbleparam']

          try:
            self.radial['mass'][0]
          except:
            self.radial_profile(Snap, ['mass'])

          rad = 10**self.radial['radius']
          self.radial[field][0] = np.log10(10**self.radial['mass'][0] * SOLAR_MASS / np.pi / (fac * rad[0])**2)
          self.radial[field][1:RadialBins] = np.log10(10**self.radial['mass'][1:RadialBins] * SOLAR_MASS / np.pi / \
                                             (fac * fac *(rad[1:RadialBins]**2 - rad[0:RadialBins-1]**2)))
          del rad

        # Effective gamma
        if field == 'geff':
          try:
            self.radial['mu'][0]
          except:
            self.radial_profile(Snap, ['mu'])
          try:
            self.radial['temp'][0]
          except:
            self.radial_profile(Snap, ['temp'])
          try:
            self.radial['nh'][0]
          except:
            self.radial_profile(Snap, ['nh'])
          mu   = 10**self.radial['mu']
          temp = 10**self.radial['temp']
          nh   = 10**self.radial['nh']
          rad  = self.radial['radius']
          #smooth_mu   = smooth(mu, window_len=21)
          #smooth_temp = smooth(temp, window_len=21)
          #smooth_nh   = smooth(nh, window_len=21)
          self.radial[field][1:RadialBins] = np.log10(nh[1:RadialBins] * temp[1:RadialBins] / mu[1:RadialBins] / \
                                                     (nh[0:RadialBins-1] * temp[0:RadialBins-1] / mu[0:RadialBins-1])) / \
                                             np.log10(nh[1:RadialBins] / nh[0:RadialBins-1])
          #self.radial[field][1:RadialBins] = np.log10(smooth_nh[1:RadialBins] * smooth_temp[1:RadialBins] / smooth_mu[1:RadialBins] / \
          #                                           (smooth_nh[0:RadialBins-1] * smooth_temp[0:RadialBins-1] / smooth_mu[0:RadialBins-1])) / \
          #                                   np.log10(smooth_nh[1:RadialBins] / smooth_nh[0:RadialBins-1])
          self.radial[field][0] = self.radial[field][1]
          #self.radial[field] = smooth(self.radial[field], window_len=6)
          del mu, temp, nh
          #del smooth_mu, smooth_temp, smooth_nh
 
        # Density dispersion
        if field == 'dnh':
          try:
            self.radial['nh'][0]
          except:
            self.radial_profile(Snap, ['nh'])
          try:
            self.radial['mass'][0]
          except:
            self.radial_profile(Snap, ['mass'])
          for j in np.unique(r_idx[(r_idx >= 0) & (r_idx < RadialBins)]):
            idx = np.where(r_idx == j)[0]
            nh = 10**self.radial['nh'][j]
#            mass = 10**self.radial['mass'][j]
#            dnh = sum(Snap.fields['mass'][idx] * ((Snap.fields['nh'][idx] - nh) / nh)**2)
#            self.radial[field][j] = np.log10(np.sqrt(dnh / mass))
            self.radial[field][j] = np.log10(np.sqrt(np.sum(Snap.fields['mass'][idx] * ((Snap.fields['nh'][idx] - nh) / nh)**2) / \
                                             10**self.radial['mass'][j]))
                                      #sum(Snap.fields['mass'][idx])))
          del nh

        # Radial and rotational velocities
        if field == 'vrad' or field == 'vrot' or field == 'vturb' or field == 'angmom':
          if IO_VEL:
            for j in np.unique(r_idx[(r_idx >= 0) & (r_idx < RadialBins)]):
              idx = np.where(r_idx == j)[0]
              mass = Snap.fields['mass'][idx]
              totmass = np.sum(mass)
              x, y, z = Snap.fields['x'][idx], Snap.fields['y'][idx], Snap.fields['z'][idx]
              vx, vy, vz = Snap.fields['vx'][idx], Snap.fields['vy'][idx], Snap.fields['vz'][idx]
              radius = 10**self.radial['radius'][j]
              if field == 'vrot' or field == 'vturb':
                try:
                  self.radial['mass'][0]
                except:
                  self.radial_profile(Snap, ['mass'])
                rad_mass = 10**self.radial['mass'][j]
              vel_cm_x = np.sum(mass * vx) / totmass
              vel_cm_y = np.sum(mass * vy) / totmass
              vel_cm_z = np.sum(mass * vz) / totmass

              if field == 'vrad' or field == 'vturb':
                vradx = x * (vx - vel_cm_x)
                vrady = y * (vy - vel_cm_y)
                vradz = z * (vz - vel_cm_z)
              if field == 'vrad':
                self.radial['vrad'][j] = np.sum(mass * (vradx + vrady + vradz) / Snap.fields['radius'][idx]) / totmass
                del vradx, vrady, vradz
              if field == 'vrot' or field == 'angmom':
                angmom_x = np.sum(mass * (y * (vz - vel_cm_z) - z * (vy - vel_cm_y)))
                angmom_y = np.sum(mass * (z * (vx - vel_cm_x) - x * (vz - vel_cm_z)))
                angmom_z = np.sum(mass * (x * (vy - vel_cm_y) - y * (vx - vel_cm_x)))
		if field == 'angmom':
                  self.radial[field][j] = np.log10(np.sqrt(angmom_x**2 + angmom_y**2 + angmom_z**2) / totmass * UNIT_VELOCITY / UNIT_LENGTH)
		if field == 'vrot':
                  self.radial[field][j] = np.sqrt(angmom_x**2 + angmom_y**2 + angmom_z**2) / radius / rad_mass
                del angmom_x, angmom_y, angmom_z
              if field == 'vturb':
                vradx = vradx / Snap.fields['radius'][idx]
                vrady = vrady / Snap.fields['radius'][idx]
                vradz = vradz / Snap.fields['radius'][idx]

                angmom_x = y * (vz - vel_cm_z) - z * (vy - vel_cm_y)
                angmom_y = z * (vx - vel_cm_x) - x * (vz - vel_cm_z)
                angmom_z = x * (vy - vel_cm_y) - y * (vx - vel_cm_x)
                vrotx = (angmom_y * z - angmom_z * y) / Snap.fields['radius'][idx]**2
                vroty = (angmom_z * x - angmom_x * z) / Snap.fields['radius'][idx]**2
                vrotz = (angmom_x * y - angmom_y * x) / Snap.fields['radius'][idx]**2
                del angmom_x, angmom_y, angmom_z

                vturbx = vx - vradx - vrotx
                vturby = vy - vrady - vroty
                vturbz = vz - vradz - vrotz
                del vradx, vrady, vradz
                del vrotx, vroty, vrotz

                self.radial[field][j] = np.sum(mass * np.sqrt(vturbx**2 + vturby**2 + vturbz**2))/ rad_mass
                del vturbx, vturby, vturbz
              del vel_cm_x, vel_cm_y, vel_cm_z
              del x, y, z, radius
              del vx, vy, vz
              del mass, totmass 
          else:
            print 'Velocity was not read!'

        # Angular velocity 1
        if field == 'omega_rot':
          try:
            self.radial['vrot'][0]
          except:
            self.radial_profile(Snap, ['vrot'])
          self.radial[field] = np.log10(self.radial['vrot'] / 10**self.radial['radius'])

#        # Angular velocity 2
#        if field == 'omega_rot':
#          if IO_VEL:
#            for j in np.unique(r_idx[(r_idx >= 0) & (r_idx < RadialBins)]):
#              idx = np.where(r_idx == j)[0]
#              mass = Snap.fields['mass'][idx]
#              totmass = np.sum(mass)
#              x, y, z = Snap.fields['x'][idx], Snap.fields['y'][idx], Snap.fields['z'][idx]
#              vx, vy, vz = Snap.fields['vx'][idx], Snap.fields['vy'][idx], Snap.fields['vz'][idx]
#              radius = Snap.fields['radius'][idx]
#              vel_cm_x = np.sum(mass * vx) / totmass
#              vel_cm_y = np.sum(mass * vy) / totmass
#              vel_cm_z = np.sum(mass * vz) / totmass
#              angmom_x = y * (vz - vel_cm_z) - z * (vy - vel_cm_y)
#              angmom_y = z * (vx - vel_cm_x) - x * (vz - vel_cm_z)
#              angmom_z = x * (vy - vel_cm_y) - y * (vx - vel_cm_x)
#              omegarot = np.sqrt(angmom_x**2 + angmom_y**2 + angmom_z**2) / radius / radius
#              self.radial[field][j] = np.sum(mass * np.log10(omegarot)) / totmass
#              del angmom_x, angmom_y, angmom_z, omegarot
 
        # Keplerian velocity
        if field == 'vkep':
          try:
            self.radial['enc_mass'][0]
          except:
            self.radial_profile(Snap, ['enc_mass'])
          # Convert units
          if LengthUnit == 0:
            fac = UNIT_LENGTH
          elif LengthUnit == 1:
            fac = UNIT_LENGTH / 1e3
          elif LengthUnit == 2:
            fac = ASTRONOMICAL_UNIT
          else:
            print 'Length unit not implemented, we will assume pc'
            fac = UNIT_LENGTH
          if not ComovingUnits:
            fac /= (1 + Snap.params['redshift'])
          fac /= Snap.params['hubbleparam']

          self.radial[field] = np.sqrt(GRAVITY * 10**self.radial['enc_mass'] * SOLAR_MASS / (fac * 10**self.radial['radius'])) / 1.0e5

        # Gravitational torque
        if field == 'taugrav' or field == 'tgrav':
          if IO_GRAVACC:
            for j in np.unique(r_idx[(r_idx >= 0) & (r_idx < RadialBins)]):
              idx = np.where(r_idx == j)[0]
              mass = Snap.fields['mass'][idx]
              totmass = np.sum(mass)
              x, y, z = Snap.fields['x'][idx], Snap.fields['y'][idx], Snap.fields['z'][idx]
              ax, ay, az = Snap.fields['gravaccx'][idx], Snap.fields['gravaccy'][idx], Snap.fields['gravaccz'][idx]
              taugrav_x = np.sum(mass * (y * az - z * ay))
              taugrav_y = np.sum(mass * (z * ax - x * az))
              taugrav_z = np.sum(mass * (x * ay - y * ax))
              if field == 'taugrav':
                self.radial[field][j] = np.log10(np.sqrt(taugrav_x**2 + taugrav_y**2 + taugrav_z**2) / totmass)
              if field == 'tgrav':
                vx, vy, vz = Snap.fields['vx'][idx], Snap.fields['vy'][idx], Snap.fields['vz'][idx]

                vel_cm_x = np.sum(mass * vx) / totmass
                vel_cm_y = np.sum(mass * vy) / totmass
                vel_cm_z = np.sum(mass * vz) / totmass

                angmom_x = np.sum(mass * (y * (vz - vel_cm_z) - z * (vy - vel_cm_y)))
                angmom_y = np.sum(mass * (z * (vx - vel_cm_x) - x * (vz - vel_cm_z)))
                angmom_z = np.sum(mass * (x * (vy - vel_cm_y) - y * (vx - vel_cm_x)))
                angmomsq = angmom_x**2. + angmom_y**2. + angmom_z**2.

                self.radial[field][j] = angmomsq / (angmom_x * taugrav_x + angmom_y * taugrav_y + angmom_z * taugrav_z)
                del vx, vy, vz
                del vel_cm_x, vel_cm_y, vel_cm_z
                del angmom_x, angmom_y, angmom_z, angmomsq
            del x, y, z
            del ax, ay, az
            del taugrav_x, taugrav_y, taugrav_z
            del mass, totmass
          else:
            print 'Gravitational acceleration was not read!'

        # Pressure torque
        if field == 'taupres' or field == 'tpres':
          if IO_GRADP:
            for j in np.unique(r_idx[(r_idx >= 0) & (r_idx < RadialBins)]):
              idx = np.where(r_idx == j)[0]
              mass = Snap.fields['mass'][idx]
              rho = Snap.fields['rho'][idx]
              totmass = np.sum(mass)
              x, y, z = Snap.fields['x'][idx], Snap.fields['y'][idx], Snap.fields['z'][idx]
              gradpx, gradpy, gradpz = Snap.fields['gradpx'][idx], Snap.fields['gradpy'][idx], Snap.fields['gradpz'][idx]
              taupres_x = np.sum(mass * (y * gradpz - z * gradpy) / rho)
              taupres_y = np.sum(mass * (z * gradpx - x * gradpz) / rho)
              taupres_z = np.sum(mass * (x * gradpy - y * gradpx) / rho)
              if field == 'taupres':
                self.radial[field][j] = np.log10(np.sqrt(taupres_x**2 + taupres_y**2 + taupres_z**2) / totmass)
              if field == 'tpres':
                vx, vy, vz = Snap.fields['vx'][idx], Snap.fields['vy'][idx], Snap.fields['vz'][idx]

                vel_cm_x = np.sum(mass * vx) / totmass
                vel_cm_y = np.sum(mass * vy) / totmass
                vel_cm_z = np.sum(mass * vz) / totmass

                angmom_x = np.sum(mass * (y * (vz - vel_cm_z) - z * (vy - vel_cm_y)))
                angmom_y = np.sum(mass * (z * (vx - vel_cm_x) - x * (vz - vel_cm_z)))
                angmom_z = np.sum(mass * (x * (vy - vel_cm_y) - y * (vx - vel_cm_x)))
                angmomsq = angmom_x**2 + angmom_y**2 + angmom_z**2

                self.radial[field][j] = angmomsq / (angmom_x * taupres_x + angmom_y * taupres_y + angmom_z * taupres_z)
                del vx, vy, vz
                del vel_cm_x, vel_cm_y, vel_cm_z
                del angmom_x, angmom_y, angmom_z, angmomsq
            del x, y, z
            del gradpx, gradpy, gradpz
            del taupres_x, taupres_y, taupres_z
            del rho, mass, totmass
          else:
            print 'Pressure gradient was not read!'


        # Toomre's Q parameter
        if field == 'q':
          try:
            self.radial['cs'][0]
          except:
            self.radial_profile(Snap, ['cs'])
          try:
            self.radial['omega_rot'][0]
          except:
            self.radial_profile(Snap, ['omega_rot'])
          try:
            self.radial['sigma_gas'][0]
          except:
            self.radial_profile(Snap, ['sigma_gas'])

          self.radial[field] = np.log10(10**self.radial['cs'] * 1.0e5 * 10**self.radial['omega_rot'] / SEC_PER_YEAR / 
                                        np.pi / GRAVITY / 10**self.radial['sigma_gas'])

#          # Version 2: Disk
#          disk_radius = np.sqrt(Snap.fields['x']**2 + Snap.fields['y']**2)
#          logdiskr = np.log10(disk_radius)
#          r_idx = (rfac*(logdiskr - logrmin)).astype(int)
#          r_idx[logdiskr == -np.inf] = -1
#          scale_height = 0.1 # AU
#          if LengthUnit == 0:
#            fac = UNIT_LENGTH
#          elif LengthUnit == 1:
#            fac = UNIT_LENGTH / 1e3
#          elif LengthUnit == 2:
#            fac = ASTRONOMICAL_UNIT
#          else:
#            print 'Length unit not implemented, we will assume pc'
#            fac = UNIT_LENGTH
#          if not ComovingUnits:
#            fac /= (1 + Snap.params['redshift'])
#          fac /= Snap.params['hubbleparam']
#          scale_height *= fac
#          rad = 10**self.radial['radius']
#          for j in np.unique(r_idx[(r_idx >= 0) & (r_idx < RadialBins)]):
#            idx = np.where((r_idx == j) & (Snap.fields['z'] < scale_height))[0]
#            mass = Snap.fields['mass'][idx]
#            totmass = np.sum(mass)
#            x, y, z = Snap.fields['x'][idx], Snap.fields['y'][idx], Snap.fields['z'][idx]
#            vx, vy, vz = Snap.fields['vx'][idx], Snap.fields['vy'][idx], Snap.fields['vz'][idx]
#            radius = rad[j]
#            vel_cm_x = np.sum(mass * vx) / totmass
#            vel_cm_y = np.sum(mass * vy) / totmass
#            vel_cm_z = np.sum(mass * vz) / totmass
#
#            angmom_x = np.sum(mass * (y * (vz - vel_cm_z) - z * (vy - vel_cm_y)))
#            angmom_y = np.sum(mass * (z * (vx - vel_cm_x) - x * (vz - vel_cm_z)))
#            angmom_z = np.sum(mass * (x * (vy - vel_cm_y) - y * (vx - vel_cm_x)))
#            omega_rot = np.sqrt(angmom_x**2 + angmom_y**2 + angmom_z**2) / radius / totmass / radius
#            del vel_cm_x, vel_cm_y, vel_cm_z
#            del angmom_x, angmom_y, angmom_z
#
#            cs = np.sum(mass * np.log10(np.sqrt(Snap.fields['gamma'][idx] * BOLTZMANN *Snap.fields['temp'][idx] / Snap.fields['mu'][idx] / PROTONMASS) / 1.0e5)) / totmass
#
#            if j == 0:
#              sigma_gas = totmass * SOLAR_MASS / np.pi / (fac * rad[0])**2
#            else:
#              sigma_gas = totmass * SOLAR_MASS / np.pi / (fac * fac *(rad[j]**2 - rad[j-1]**2))
#
#            self.radial[field][j] = np.log10(10**cs * 1.0e5 * omega_rot / SEC_PER_YEAR / np.pi / GRAVITY / sigma_gas)

        # Bonnor-Ebert mass
        if field == 'mbe_ratio':
          if IO_U and IO_RHO and IO_CHEM and IO_GAMMA:
            mbe = np.zeros(RadialBins, float)
            for j in np.unique(r_idx[(r_idx >= 0) & (r_idx < RadialBins)]):
              idx = np.where(r_idx == j)[0]
              mbe[j] = np.sum(Snap.fields['mass'][idx] * 2.6e-11 * Snap.fields['temp'][idx]**(3./2.) * Snap.fields['rho'][idx]**(-1./2.) \
                                                    * Snap.fields['mu'][idx]**(-3./2.) * Snap.fields['gamma'][idx]**2)
            try:
              self.radial['enc_mass'][0]
            except:
              self.radial_profile(Snap, ['enc_mass'])
            self.radial[field] = np.log10(10**self.radial['enc_mass'] / (np.cumsum(mbe) / 10**self.radial['enc_mass']))
            del mbe
          else:
            print 'Energy, Density, Chemistry or Gamma not read!'

        # Radial escape fraction
        if field == 'radial_esc_frac':
          if IO_RHO:
            try:
              self.radial['nh'][0]
            except:
              self.radial_profile(Snap, ['nh'])
            self.radial[field] = np.exp(-10**self.radial['nh'] / 1.5e16)
          else:
            print 'Density was not read!'

        # Radial cooling rate
        if field == 'radial_cool_rate':
          if IO_U and IO_RHO and IO_CHEM:
            try:
              self.radial['radial_esc_frac'][0]
            except:
              self.radial_profile(Snap, ['radial_esc_frac'])
            try:
              self.radial['abHII'][0]
            except:
              self.radial_profile(Snap, ['abHII'])
            try:
              self.radial['abH2'][0]
            except:
              self.radial_profile(Snap, ['abH2'])
            try:
              self.radial['temp'][0]
            except:
              self.radial_profile(Snap, ['temp'])
            abe = 10**self.radial['abHII']
            abHI = np.maximum(1. - 2. * 10**self.radial['abH2'] - 10**self.radial['abHII'], 0.)
            self.radial[field] = 7.5e-19 * np.exp(-1.18348e5 / 10**self.radial['temp']) / (1. + np.sqrt(10**self.radial['temp'] / 1.0e5)) \
                                         * abHI * abe * 10**self.radial['nh'] * 10**self.radial['nh'] * self.radial['radial_esc_frac']
            del abe, abHI
          else:
            print 'Energy, Density or Chemistry was not read!'

        # Radial cooling time
        if field == 'radial_tcool_ratio':
          if IO_U and IO_RHO:
            try:
              self.radial['radial_cool_rate'][0]
            except:
              self.radial_profile(Snap, ['radial_cool_rate'])
            try:
              self.radial['rho'][0]
            except:
              self.radial_profile(Snap, ['rho'])
            try:
              self.radial['u'][0]
            except:
              self.radial_profile(Snap, ['u'])
            tff = np.sqrt(3. * np.pi / 32. / GRAVITY / 10**self.radial['rho'])
            self.radial[field] = np.log10(10**self.radial['u'] * 10**self.radial['rho'] / self.radial['radial_cool_rate'] / tff)
          else:
            print ' Energy or Density was not read!'

        # Radial Gammie criterion
        if field == 'gammie':
          if IO_U and IO_RHO:
            try:
              self.radial['tcool'][0]
            except:
              self.radial_profile(Snap, ['tcool'])
            try:
              self.radial['omega_rot'][0]
            except:
              self.radial_profile(Snap, ['omega_rot'])
            try:
              self.radial['radial_cool_rate'][0]
            except:
              self.radial_profile(Snap, ['radial_cool_rate'])
            try:
              self.radial['rho'][0]
            except:
              self.radial_profile(Snap, ['rho'])
            try:
              self.radial['u'][0]
            except:
              self.radial_profile(Snap, ['u'])
#          self.radial[field] = np.log10(10**self.radial['tcool'] * 10**self.radial['omega_rot'] / 3.)
            self.radial[field] = np.log10(10**self.radial['u'] * 10**self.radial['rho'] / self.radial['radial_cool_rate'] * 10**self.radial['omega_rot'] / SEC_PER_YEAR / 3.)

        # Infall rate
        if field == 'infall':
          if IO_VEL and IO_RHO:
            try:
              self.radial['vrad'][0]
            except:
              self.radial_profile(Snap, ['vrad'])
            try:
              self.radial['rho'][0]
            except:
              self.radial_profile(Snap, ['rho'])

            # Convert units
            if LengthUnit == 0:
              fac = UNIT_LENGTH
            elif LengthUnit == 1:
              fac = UNIT_LENGTH / 1e3
            elif LengthUnit == 2:
              fac = ASTRONOMICAL_UNIT
            else:
              print 'Length unit not implemented, we will assume pc'
              fac = UNIT_LENGTH

            self.radial[field] =  (-4) * np.pi * (fac * 10**self.radial['radius'])**2. * 10**self.radial['rho'] * (self.radial['vrad'] * 1e5) / SOLAR_MASS * SEC_PER_YEAR
          else:
            print ' Energy or Density was not read!'

  def find_max(self, field):
    # Returns index of cell with maximum value of a field
    return np.argmax(self.radial[field])

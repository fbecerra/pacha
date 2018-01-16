from include import *
from rosseland import *

class Snap:

  def __init__(self):
    self.params = dict()
    self.fields = dict()
    self.new_fields = dict()
    self.new_fields['sinks'] = dict()
    self.new_fields['dm'] = dict()
    self.derived_fields = dict()
    self.quantities = dict()
    self.keys = np.array([])

  def read_header(self, snapbase, read_HDF=False):
    # Read header
    self.read_HDF = read_HDF
    self.params['nsnaptypes'] = 6
    self.snapfile = snapbase+'.0'

    if self.read_HDF:

      self.file = h5py.File(self.snapfile+'.hdf5', 'r')

      self.params['masstable'] = np.array(self.file['Header'].attrs['MassTable'], dtype=np.dtype('d'))
      self.params['time'] = np.array(self.file['Header'].attrs['Time'], dtype=np.dtype('d'))
      self.params['redshift'] = np.array(self.file['Header'].attrs['Redshift'], dtype=np.dtype('d'))
      self.params['flag_sfr'] = np.array(self.file['Header'].attrs['Flag_Sfr'], dtype=np.dtype('i'))
      self.params['flag_feedback'] = np.array(self.file['Header'].attrs['Flag_Feedback'], dtype=np.dtype('i'))
      self.params['npartall'] = np.array(self.file['Header'].attrs['NumPart_Total'], dtype=np.dtype('i'))
      self.params['flag_cooling'] = np.array(self.file['Header'].attrs['Flag_Cooling'], dtype=np.dtype('i'))
      self.params['numfiles'] = np.array(self.file['Header'].attrs['NumFilesPerSnapshot'], dtype=np.dtype('i'))
      self.params['boxsize'] = np.array(self.file['Header'].attrs['BoxSize'], dtype=np.dtype('d'))
      self.params['omega_m'] = np.array(self.file['Header'].attrs['Omega0'], dtype=np.dtype('d'))
      self.params['omega_lambda'] = np.array(self.file['Header'].attrs['OmegaLambda'], dtype=np.dtype('d'))
      self.params['hubbleparam'] = np.array(self.file['Header'].attrs['HubbleParam'], dtype=np.dtype('d'))

      self.file.close()

    else:

      self.file = open(self.snapfile, mode = 'rb')
      self.file.seek(28)
      self.params['masstable'] = np.fromfile(self.file, dtype= np.dtype('d'), count=self.params['nsnaptypes'])
      self.params['time'] = np.fromfile(self.file, dtype=np.dtype('d'), count=1)[0]
      self.params['redshift'] = np.fromfile(self.file, dtype=np.dtype('d'), count=1)[0]
      self.params['flag_sfr'] = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
      self.params['flag_feedback'] = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
      self.params['npartall'] = np.fromfile(self.file, dtype= np.dtype('i'), count=self.params['nsnaptypes'])
      self.params['flag_cooling'] = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
      self.params['numfiles'] = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
      self.params['boxsize'] = np.fromfile(self.file, dtype=np.dtype('d'), count=1)[0]
      self.params['omega_m'] = np.fromfile(self.file, dtype=np.dtype('d'), count=1)[0]
      self.params['omega_lambda'] = np.fromfile(self.file, dtype=np.dtype('d'), count=1)[0]
      self.params['hubbleparam'] = np.fromfile(self.file, dtype=np.dtype('d'), count=1)[0]
      self.file.close()

  def read_fields(self, snapbase, read_dm=False, read_HDF=False):
    self.read_dm = read_dm
    self.read_HDF = read_HDF

    # Number of sinks
    nsinks = self.params['npartall'][5]

    # Define arrays to use
    if IO_POS:
      self.fields['x'] = np.array([], dtype = np.float64)
      self.fields['y'] = np.array([], dtype = np.float64)
      self.fields['z'] = np.array([], dtype = np.float64)
      self.fields['radius'] = np.array([], dtype = np.float64)
      if self.read_dm:
        self.new_fields['dm']['x'] = np.array([], dtype = np.float64)
        self.new_fields['dm']['y'] = np.array([], dtype = np.float64)
        self.new_fields['dm']['z'] = np.array([], dtype = np.float64)
      if nsinks:
        self.new_fields['sinks']['x'] = np.array([], dtype = np.float64)
        self.new_fields['sinks']['y'] = np.array([], dtype = np.float64)
        self.new_fields['sinks']['z'] = np.array([], dtype = np.float64)
    if IO_VEL:
      self.fields['vx'] = np.array([], dtype = np.float64)
      self.fields['vy'] = np.array([], dtype = np.float64)
      self.fields['vz'] = np.array([], dtype = np.float64)
      if self.read_dm:
        self.new_fields['dm']['vx'] = np.array([], dtype = np.float64)
        self.new_fields['dm']['vy'] = np.array([], dtype = np.float64)
        self.new_fields['dm']['vz'] = np.array([], dtype = np.float64)
      if nsinks:
        self.new_fields['sinks']['vx'] = np.array([], dtype = np.float64)
        self.new_fields['sinks']['vy'] = np.array([], dtype = np.float64)
        self.new_fields['sinks']['vz'] = np.array([], dtype = np.float64)
    if IO_ID:
      self.fields['id'] = np.array([], dtype = np.int64)
      if self.read_dm:
        self.new_fields['dm']['id'] = np.array([], dtype = np.int64)
      if nsinks:
        self.new_fields['sinks']['id'] = np.array([], dtype = np.int64)
    if IO_MASS:
      self.fields['mass'] = np.array([], dtype = np.float64)
      if self.read_dm:
        self.new_fields['dm']['mass'] = np.array([], dtype = np.float64)
      if nsinks:
        self.new_fields['sinks']['mass'] = np.array([], dtype = np.float64)
    if IO_U:
      self.fields['u'] = np.array([], dtype = np.float64)
      self.fields['temp'] = np.array([], dtype = np.float64)
    if IO_RHO:
      self.fields['rho'] = np.array([], dtype = np.float64)
      self.fields['nh'] = np.array([], dtype = np.float64)
    if IO_VOL:
      self.fields['vol'] = np.array([], dtype = np.float64)
      self.fields['hsml'] = np.array([], dtype = np.float64)
    if IO_DELAUNAY:
      self.fields['delaunay'] = np.array([], dtype = np.float64)
    if IO_GRAVACC:
      self.fields['gravaccx'] = np.array([], dtype = np.float64)
      self.fields['gravaccy'] = np.array([], dtype = np.float64)
      self.fields['gravaccz'] = np.array([], dtype = np.float64)
      self.fields['gravacc'] = np.array([], dtype = np.float64)
    if IO_GRADP:
      self.fields['gradpx'] = np.array([], dtype = np.float64)
      self.fields['gradpy'] = np.array([], dtype = np.float64)
      self.fields['gradpz'] = np.array([], dtype = np.float64)
    if IO_CHEM:
      self.fields['abHM'] = np.array([], dtype = np.float64)
      self.fields['abH2'] = np.array([], dtype = np.float64)
      self.fields['abHII'] = np.array([], dtype = np.float64)
      self.fields['mu'] = np.array([], dtype = np.float64)
    if IO_GAMMA:
      self.fields['gamma'] = np.array([], dtype = np.float64)

    # Check if header has been read
    try:
      self.params['numfiles']
    except:
      self.read_header(snapbase, read_HDF=self.read_HDF)

    # Iterate over files to read all particles
     ########## vvvvvv HDF vvvvvvvv #########
    if self.read_HDF:

      for i in range(self.params['numfiles']):

        self.snapfile = snapbase + '.' + repr(i)
        self.file = h5py.File(self.snapfile+'.hdf5', 'r')
        self.npart = np.array(self.file['Header'].attrs['NumPart_ThisFile'], dtype=np.dtype('d'))

        # Read positions
        if IO_POS:
          # Convert units
          if LengthUnit == 0:
            fac = 1
          elif LengthUnit == 1:
            fac = 1e3
          elif LengthUnit == 2:
            fac = UNIT_LENGTH / ASTRONOMICAL_UNIT
          else:
            print 'Length unit not implemented, we will assume pc'
            fac = 1
          if not ComovingUnits:
            fac /= (1 + self.params['redshift'])
          fac /= self.params['hubbleparam']
          # Read values
          pos = fac * np.array(self.file['PartType0/Coordinates'], dtype=np.dtype('d'))
          self.fields['x'] = np.append(self.fields['x'], pos[:,0])
          self.fields['y'] = np.append(self.fields['y'], pos[:,1])
          self.fields['z'] = np.append(self.fields['z'], pos[:,2])
          if self.read_dm:
            pos = fac * np.array(self.file['PartType1/Coordinates'], dtype=np.dtype('d'))
            self.new_fields['dm']['x'] = np.append(self.new_fields['dm']['x'], pos[:,0])
            self.new_fields['dm']['y'] = np.append(self.new_fields['dm']['y'], pos[:,1])
            self.new_fields['dm']['z'] = np.append(self.new_fields['dm']['z'], pos[:,2])
          if nsinks:
            pos = fac * np.array(self.file['PartType4/Coordinates'], dtype=np.dtype('d'))
            self.new_fields['sinks']['x'] = np.append(self.new_fields['sinks']['x'], pos[:,0])
            self.new_fields['sinks']['y'] = np.append(self.new_fields['sinks']['y'], pos[:,1])
            self.new_fields['sinks']['z'] = np.append(self.new_fields['sinks']['z'], pos[:,2])
          del pos

        # Read velocities
        if IO_VEL:
          # Convert units
          fac = 1 / np.sqrt(1 + self.params['redshift'])
          # Read values
          vel = fac * np.array(self.file['PartType0/Velocities'], dtype=np.dtype('d'))
          self.fields['vx'] = np.append(self.fields['vx'], vel[:,0])
          self.fields['vy'] = np.append(self.fields['vy'], vel[:,1])
          self.fields['vz'] = np.append(self.fields['vz'], vel[:,2])
          if self.read_dm:
            vel = fac * np.array(self.file['PartType1/Velocities'], dtype=np.dtype('d'))
            self.new_fields['dm']['vx'] = np.append(self.new_fields['dm']['vx'], vel[:,0])
            self.new_fields['dm']['vy'] = np.append(self.new_fields['dm']['vy'], vel[:,1])
            self.new_fields['dm']['vz'] = np.append(self.new_fields['dm']['vz'], vel[:,2])
          if nsinks:
            vel = fac * np.array(self.file['PartType4/Velocities'], dtype=np.dtype('d'))
            self.new_fields['sinks']['vx'] = np.append(self.new_fields['sinks']['vx'], vel[:,0])
            self.new_fields['sinks']['vy'] = np.append(self.new_fields['sinks']['vy'], vel[:,1])
            self.new_fields['sinks']['vz'] = np.append(self.new_fields['sinks']['vz'], vel[:,2])
          del vel

        # Read IDs
        if IO_ID:
          id = np.array(self.file['PartType0/ParticleIDs'], dtype=np.dtype('i'))
          self.fields['id'] = np.append(self.fields['id'], id)
          if self.read_dm:
            id = np.array(self.file['PartType1/ParticleIDs'], dtype=np.dtype('i'))
            self.new_fields['dm']['id'] = np.append(self.new_fields['dm']['id'], id)
          if nsinks:
            id = np.array(self.file['PartType4/ParticleIDs'], dtype=np.dtype('i'))
            self.new_fields['sinks']['id'] = np.append(self.new_fields['sinks']['id'], id)
          del id

        # Read mass
        if IO_MASS:
          # Convert units
          fac = UNIT_MASS / SOLAR_MASS / self.params['hubbleparam']
          # Read values
          mass = fac * np.array(self.file['PartType0/Masses'], dtype=np.dtype('d'))
          self.fields['mass'] = np.append(self.fields['mass'], mass)
          if self.read_dm:
            ndmpart = self.npart[1]
            dmpartmass = self.params['masstable'][1]
            mass = fac * partmass * np.ones(ndmpart, dtype = np.float64)
            self.new_fields['dm']['mass'] = np.append(self.new_fields['dm']['mass'], mass)
          if nsinks:
            mass = fac * np.array(self.file['PartType4/Masses'], dtype=np.dtype('d'))
            self.new_fields['sinks']['mass'] = np.append(self.new_fields['sinks']['mass'], mass)
          del mass

        # Read energy
        if IO_U:
          # Convert units
          fac = UNIT_ENERGY / UNIT_MASS
          # Read values
          u = fac * np.array(self.file['PartType0/InternalEnergy'], dtype=np.dtype('d'))
          self.fields['u'] = np.append(self.fields['u'], u)
          del u

        # Read density
        if IO_RHO:
          # Convert units
          fac = UNIT_MASS / UNIT_LENGTH**3 * (1 + self.params['redshift'])**3 * self.params['hubbleparam']**2
          # Read values
          rho = fac * np.array(self.file['PartType0/Density'], dtype=np.dtype('d'))
          self.fields['rho'] = np.append(self.fields['rho'], rho)
          del rho

     ############# ^^^^ HDF ^^^^^ #########
    else:

      for i in range(self.params['numfiles']):
      
        self.snapfile = snapbase + '.' + repr(i)
        self.file = open(self.snapfile, mode = 'rb')
  
        self.nbytes = 256 - 6*4
        self.file.seek(4)
        self.npart = np.fromfile(self.file, dtype=np.dtype('i'), count=self.params['nsnaptypes'])
        self.file.seek(self.nbytes,1)
        self.nbytes = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
  
        [ngas, nhalo, ndisk, nstars, nboundary, nsinks] = self.npart
        ndm = nhalo + ndisk + nstars
  
        # Read positions
        if IO_POS:
          self.nbytes = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
          # Convert units
          if LengthUnit == 0:
            fac = 1
          elif LengthUnit == 1:
            fac = 1e3
          elif LengthUnit == 2:
            fac = UNIT_LENGTH / ASTRONOMICAL_UNIT
          else:
            print 'Length unit not implemented, we will assume pc'
            fac = 1
          if not ComovingUnits:
            fac /= (1 + self.params['redshift'])
          fac /= self.params['hubbleparam']
          # Read values
          pos = fac * np.fromfile(self.file, dtype=np.dtype('d'), count=3*ngas)
          self.fields['x'] = np.append(self.fields['x'], pos[0::3])
          self.fields['y'] = np.append(self.fields['y'], pos[1::3])
          self.fields['z'] = np.append(self.fields['z'], pos[2::3])
          if self.read_dm:
            pos = fac * np.fromfile(self.file, dtype=np.dtype('d'), count=3*ndm)
            self.new_fields['dm']['x'] = np.append(self.new_fields['dm']['x'], pos[0::3])
            self.new_fields['dm']['y'] = np.append(self.new_fields['dm']['y'], pos[1::3])
            self.new_fields['dm']['z'] = np.append(self.new_fields['dm']['z'], pos[2::3])
            self.file.seek(self.nbytes-3*(ngas+ndm+nsinks)*8,1)
          else:
            self.file.seek(self.nbytes-3*(ngas+nsinks)*8,1)
          if nsinks:
            pos = fac * np.fromfile(self.file, dtype=np.dtype('d'), count=3*nsinks)
            self.new_fields['sinks']['x'] = np.append(self.new_fields['sinks']['x'], pos[0::3])
            self.new_fields['sinks']['y'] = np.append(self.new_fields['sinks']['y'], pos[1::3])
            self.new_fields['sinks']['z'] = np.append(self.new_fields['sinks']['z'], pos[2::3])
          del pos
          self.nbytes = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
  
        # Read velocities
        if IO_VEL:
          self.nbytes = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
          # Convert units
          fac = 1 / np.sqrt(1 + self.params['redshift'])
          # Read values
          vel = fac * np.fromfile(self.file, dtype=np.dtype('d'), count=3*ngas)
          self.fields['vx'] = np.append(self.fields['vx'], vel[0::3])
          self.fields['vy'] = np.append(self.fields['vy'], vel[1::3])
          self.fields['vz'] = np.append(self.fields['vz'], vel[2::3])
          if self.read_dm:
            vel = fac * np.fromfile(self.file, dtype=np.dtype('d'), count=3*ndm)
            self.new_fields['dm']['vx'] = np.append(self.new_fields['dm']['vx'], vel[0::3])
            self.new_fields['dm']['vy'] = np.append(self.new_fields['dm']['vy'], vel[1::3])
            self.new_fields['dm']['vz'] = np.append(self.new_fields['dm']['vz'], vel[2::3])
            self.file.seek(self.nbytes-3*(ngas+ndm+nsinks)*8,1)
          else:
            self.file.seek(self.nbytes-3*(ngas+nsinks)*8,1)
          if nsinks:
            vel = fac * np.fromfile(self.file, dtype=np.dtype('d'), count=3*nsinks)
            self.new_fields['sinks']['vx'] = np.append(self.new_fields['sinks']['vx'], vel[0::3])
            self.new_fields['sinks']['vy'] = np.append(self.new_fields['sinks']['vy'], vel[1::3])
            self.new_fields['sinks']['vz'] = np.append(self.new_fields['sinks']['vz'], vel[2::3])
          del vel
          self.nbytes = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
  
        # Read IDs
        if IO_ID:
          self.nbytes = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
          id = np.fromfile(self.file, dtype=np.dtype('i'), count=ngas)
          self.fields['id'] = np.append(self.fields['id'], id)
          if self.read_dm:
            id = np.fromfile(self.file, dtype=np.dtype('i'), count=ndm)
            self.new_fields['dm']['id'] = np.append(self.new_fields['dm']['id'], id)
            self.file.seek(self.nbytes-(ngas+ndm+nsinks)*4,1)
          else:
            self.file.seek(self.nbytes-(ngas+nsinks)*4,1)
          if nsinks:
            id = np.fromfile(self.file, dtype=np.dtype('i'), count=nsinks)
            self.new_fields['sinks']['id'] = np.append(self.new_fields['sinks']['id'], id)
          del id
          self.nbytes = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
  
        # Read mass
        if IO_MASS:
          self.nbytes = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
          # Convert units
          fac = UNIT_MASS / SOLAR_MASS / self.params['hubbleparam']
          # Read values
          mass = fac * np.fromfile(self.file, dtype=np.dtype('d'), count=ngas)
          self.fields['mass'] = np.append(self.fields['mass'], mass)
          if self.read_dm:
            ndmpart_total = 0
            for idx in xrange(1, self.params['nsnaptypes']):
              ndmpart = self.params['npartall'][idx]
              partmass = self.params['masstable'][idx]
              if partmass:
                mass = fac * partmass * np.ones(ndmpart, dtype = np.float64)
                self.new_fields['dm']['mass'] = np.append(self.new_fields['dm']['mass'], mass)
  #            else:
  #              mass = fac * np.fromfile(self.file, dtype=np.dtype('d'), count=ndmpart)
  #              ndmpart_total += ndmpart
  #            self.new_fields['dm']['mass'] = np.append(self.new_fields['dm']['mass'], mass)
            self.file.seek(self.nbytes-(ngas+ndmpart_total+nsinks)*8,1)
          else:
            self.file.seek(self.nbytes-(ngas+nsinks)*8,1)
          if nsinks:
            mass = fac * np.fromfile(self.file, dtype=np.dtype('d'), count=nsinks)
            self.new_fields['sinks']['mass'] = np.append(self.new_fields['sinks']['mass'], mass)
          del mass
          self.nbytes = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
  
        # Read energy
        if IO_U:
          self.nbytes = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
          # Convert units
          fac = UNIT_ENERGY / UNIT_MASS
          # Read values
          u = fac * np.fromfile(self.file, dtype=np.dtype('d'), count=ngas)
          self.fields['u'] = np.append(self.fields['u'], u)
          del u
          self.file.seek(self.nbytes-ngas*8,1)
          self.nbytes = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
  
        # Read density
        if IO_RHO:
          self.nbytes = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
          # Convert units
          fac = UNIT_MASS / UNIT_LENGTH**3 * (1 + self.params['redshift'])**3 * self.params['hubbleparam']**2
          # Read values
          rho = fac * np.fromfile(self.file, dtype=np.dtype('d'), count=ngas)
          self.fields['rho'] = np.append(self.fields['rho'], rho)
          del rho
          self.file.seek(self.nbytes-ngas*8,1)
          self.nbytes = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
  
        # Read delaunay
        if IO_DELAUNAY:
          self.nbytes = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
          # Convert units
          if LengthUnit == 0:
            fac = 1
          elif LengthUnit == 1:
            fac = 1e9
          elif LengthUnit == 2:
            fac = (UNIT_LENGTH / ASTRONOMICAL_UNIT)**3
          else:
            print 'Length unit not implemented, we will assume pc'
            fac = 1
          if not ComovingUnits:
            fac /= (1 + self.params['redshift'])**3
          fac /= self.params['hubbleparam']**3
          # Read values
          delaunay = fac * np.fromfile(self.file, dtype=np.dtype('d'), count=ngas)
          self.fields['delaunay'] = np.append(self.fields['delaunay'], delaunay)
          del delaunay
          self.file.seek(self.nbytes-ngas*8,1)
          self.nbytes = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
  
        # Read gravitational acceleration
        if IO_GRAVACC:
          self.nbytes = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
          # Convert units
          # fac = (1 + self.params['redshift'])**(-2)
          fac = UNIT_LENGTH / UNIT_TIME**2.
          # Read values
          gravacc = fac * np.fromfile(self.file, dtype=np.dtype('d'), count=3*ngas)
          self.fields['gravaccx'] = np.append(self.fields['gravaccx'], gravacc[0::3])
          self.fields['gravaccy'] = np.append(self.fields['gravaccy'], gravacc[1::3])
          self.fields['gravaccz'] = np.append(self.fields['gravaccz'], gravacc[2::3])
          del gravacc
          self.file.seek(self.nbytes-3*ngas*8,1)
          self.nbytes = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
  
        # Read pressure gradient
        if IO_GRADP:
          self.nbytes = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
          # Convert units
          fac = (1 + self.params['redshift'])**4 * self.params['hubbleparam']**3
          fac *= UNIT_MASS / UNIT_LENGTH**2. / UNIT_TIME**2.
          # Read values
          gradp = fac * np.fromfile(self.file, dtype=np.dtype('d'), count=3*ngas)
          self.fields['gradpx'] = np.append(self.fields['gradpx'], gradp[0::3])
          self.fields['gradpy'] = np.append(self.fields['gradpy'], gradp[1::3])
          self.fields['gradpz'] = np.append(self.fields['gradpz'], gradp[2::3])
          del gradp
          self.file.seek(self.nbytes-3*ngas*8,1)
          self.nbytes = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
  
        # Read chemistry
        if IO_CHEM:
          self.nbytes = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
          chem = np.fromfile(self.file, dtype=np.dtype('d'), count=3*ngas)
          self.fields['abHM'] = np.append(self.fields['abHM'], chem[0::3])
          self.fields['abH2'] = np.append(self.fields['abH2'], chem[1::3])
          self.fields['abHII'] = np.append(self.fields['abHII'], chem[2::3])
          del chem
          self.file.seek(self.nbytes-3*ngas*8,1)
          self.nbytes = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
  
        # Read gamma
        if IO_GAMMA:
          self.nbytes = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
          gamma = np.fromfile(self.file, dtype=np.dtype('d'), count=ngas)
          self.fields['gamma'] = np.append(self.fields['gamma'], gamma)
          del gamma
          self.file.seek(self.nbytes-ngas*8,1)
          self.nbytes = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
  
        self.file.close()
  
#    if IO_POS and IO_VEL:
#      self.center_box()
#      self.rotate_box()
#      self.calculate_radius()

    if IO_POS:
      # Convert units
      if LengthUnit == 0:
        fac = 1
      elif LengthUnit == 1:
        fac = 1e3
      elif LengthUnit == 2:
        fac = UNIT_LENGTH / ASTRONOMICAL_UNIT
      else:
        print 'Length unit not implemented, we will assume pc'
        fac = 1
      if not ComovingUnits:
        fac /= (1 + self.params['redshift'])
      fac /= self.params['hubbleparam']
      self.params['boxsize'] = fac * self.params['boxsize']

    if IO_RHO:
      self.fields['nh'] = self.fields['rho'] * HYDROGEN_MASSFRAC / PROTONMASS

    if IO_VOL:
      if LengthUnit == 0:
        fac = UNIT_LENGTH
      elif LengthUnit == 1:
        fac = UNIT_LENGTH / 1e3
      elif LengthUnit == 2:
        fac = ASTRONOMICAL_UNIT
      self.fields['vol'] = self.fields['mass'] * SOLAR_MASS / self.fields['rho']
      self.fields['hsml'] = (3.*self.fields['vol']/4./np.pi)**(1./3.) / fac

    if IO_GRAVACC:
      self.fields['gravacc'] = np.sqrt(self.fields['gravaccx']**2. + self.fields['gravaccy']**2. + self.fields['gravaccz']**2.)

    if IO_CHEM:
      abe = self.fields['abHII']
      self.fields['mu'] = (1. + 4. * HE_ABUND) / (1 + HE_ABUND - self.fields['abH2'] + abe)
      del abe

    if IO_U:
      if IO_CHEM:
        mu = self.fields['mu']
      else:
        mu = MU_NEUTRAL
      if IO_GAMMA:
        gamma = self.fields['gamma']
      else:
        gamma = GAMMA_ADB
      self.fields['temp'] = mu * (gamma - 1) * PROTONMASS / BOLTZMANN * self.fields['u']
      del mu, gamma
    
  def calculate_radius(self):
      self.fields['radius'] = np.sqrt(self.fields['x']**2 + self.fields['y']**2 + self.fields['z']**2)

  def calculate_quantities(self, field):
    if field == 'center_of_mass':
      # Returns center of mass
      cmx = np.sum(self.fields['rho'] * self.fields['x']) / np.sum(self.fields['rho'])
      cmy = np.sum(self.fields['rho'] * self.fields['y']) / np.sum(self.fields['rho'])
      cmz = np.sum(self.fields['rho'] * self.fields['z']) / np.sum(self.fields['rho']) 
      self.quantities[field] = np.array([cmx, cmy, cmz])

  def find_max(self, field):
    # Returns index of cell with maximum value of a field
    if field in fields_list:
      return np.argmax(self.fields[field])
    elif field in derived_fields_list:
      return np.argmax(self.derived_fields[field])

#    if field = 'maximum_value':
#      # Returns index of cell with maximum value of a field
#      try:
#        self.quantities
#      self.quantities[field] = self.fields[field][max_index]

  def show_fields(self):
    for key in self.fields.keys():
      self.keys = np.append(self.keys, key)
    for key in self.derived_fields.keys():
      self.keys = np.append(self.keys, key)
    self.keys = np.unique(self.keys)
    print self.keys

  def remove_particles(self, idx):
    mask = np.ones(self.fields['x'].shape, dtype=np.bool)
    mask[idx] = False
    # Delete elements with indices in idx
    for key in self.fields.keys():
      self.fields[key] = self.fields[key][mask]
    for key in self.derived_fields.keys():
      self.derived_fields[key] = self.derived_fields[key][mask]

  def merge_snaps(self, snap):
    box_fraction = 10
    boxsize = self.params['boxsize']
    # Assume second snapshot bigger
    idx1 = np.where((np.abs(self.fields['x']) < boxsize / 2 / box_fraction) &
                    (np.abs(self.fields['y']) < boxsize / 2 / box_fraction) &
                    (np.abs(self.fields['z']) < boxsize / 2 / box_fraction))[0]
    idx2 = np.where((np.abs(snap.fields['x']) >= boxsize / 2 / box_fraction) &
                    (np.abs(snap.fields['y']) >= boxsize / 2 / box_fraction) &
                    (np.abs(snap.fields['z']) >= boxsize / 2 / box_fraction))[0]
    self.params['boxsize'] = snap.params['boxsize']
    for key in self.fields.keys():
      if len(self.fields[key]) != 0:
        self.fields[key] = np.append(self.fields[key][idx1], snap.fields[key][idx2])
    for key in self.derived_fields.keys():
      if len(self.derived_fields[key]) != 0:
        self.derived_fields[key] = np.append(self.derived_fields[key][idx1], snap.derived_fields[key][idx2])

  def center_box(self):
    if FlagCenter == 0:
      # Highest Density
      self.max_dens_idx = self.find_max('rho')
      self.Center = [self.fields['x'][self.max_dens_idx], self.fields['y'][self.max_dens_idx], self.fields['z'][self.max_dens_idx]]
    elif FlagCenter == 1:
      # Fixed Center
      if LengthUnit == 0:
        fac = 1.
      elif LengthUnit == 1:
        fac = 1e3
      elif LengthUnit == 2:
        fac = UNIT_LENGTH / ASTRONOMICAL_UNIT
      else:
        print 'Length unit not implemented, we will assume pc'
        fac = 1.
      if not ComovingUnits:
        fac /= (1 + self.params['redshift'])
      fac /= self.params['hubbleparam']
      self.Center = fac * np.array(Center)
    elif FlagCenter == 2:
      # Center of mass
      try:
        self.quantities['center_of_mass']
      except:
        self.calculate_quantities('center_of_mass')
      self.Center = self.quantities['center_of_mass']
#      self.Center[0] = sum(self.fields['rho'] * self.fields['x']) / sum(self.fields['rho'])
#      self.Center[1] = sum(self.fields['rho'] * self.fields['y']) / sum(self.fields['rho'])
#      self.Center[2] = sum(self.fields['rho'] * self.fields['z']) / sum(self.fields['rho'])

    # Update positions
    self.fields['x'] -= self.Center[0]
    self.fields['y'] -= self.Center[1]
    self.fields['z'] -= self.Center[2]

    #Update DM positions
    if self.read_dm:
      self.new_fields['dm']['x'] -= self.Center[0]
      self.new_fields['dm']['y'] -= self.Center[1]
      self.new_fields['dm']['z'] -= self.Center[2]

    # Update sink positions
    if len(self.new_fields['sinks']) > 0:
      self.new_fields['sinks']['x'] -= self.Center[0]
      self.new_fields['sinks']['y'] -= self.Center[1]
      self.new_fields['sinks']['z'] -= self.Center[2]

  def rotate_box(self):
    if FlagRotate == 0:
      # Rotation baed on velocity vector
      if LengthUnit == 0:
        fac = 1
      elif LengthUnit == 1:
        fac = 1e3
      elif LengthUnit == 2:
        fac = UNIT_LENGTH / ASTRONOMICAL_UNIT
      else:
        print 'Length unit not implemented, we will assume pc'
        fac = 1
      if not ComovingUnits:
        fac /= (1 + self.params['redshift'])
      fac /= self.params['hubbleparam']

      if Usage == 0:
        r_search = ImgSize / 2
      else:
        r_search = fac * self.params['boxsize'] / 2

      part_idx = np.where((self.fields['x'] < r_search) & (self.fields['y'] < r_search) & (self.fields['z'] < r_search))[0]
 
      vel_cm = np.zeros(3)
      vec = np.zeros(3)
      dvx = np.array([])
      dvy = np.array([])
      dvz = np.array([])

      mass = self.fields['mass'][part_idx]
      totmass = np.sum(mass)
      vx, vy, vz = self.fields['vx'][part_idx], self.fields['vy'][part_idx], self.fields['vz'][part_idx]
 
      vel_cm[0] = np.sum(mass * vx) / totmass 
      vel_cm[1] = np.sum(mass * vy) / totmass 
      vel_cm[2] = np.sum(mass * vz) / totmass
  
      dvx = vx - vel_cm[0] 
      dvy = vy - vel_cm[1] 
      dvz = vz - vel_cm[2]
      del vel_cm
      del vx, vy, vz
  
      vec[0] = np.sum(mass * (self.fields['y'][part_idx] * dvz - self.fields['z'][part_idx] * dvy)) 
      vec[1] = np.sum(mass * (self.fields['z'][part_idx] * dvx - self.fields['x'][part_idx] * dvz)) 
      vec[2] = np.sum(mass * (self.fields['x'][part_idx] * dvy - self.fields['y'][part_idx] * dvx))
      del dvx, dvy, dvz
      del part_idx
      del mass, totmass
  
      vec /= np.linalg.norm(vec)
  
    elif FlagRotate == 1:
      # Fixed angle
      theta = 2 * np.pi * Angle / 180
      vec = [0, np.sin(theta), np.cos(theta)]

    alpha = np.arccos(vec[2] / np.sqrt(vec[1]**2 + vec[2]**2))
    if vec[1] >= 0:
      alpha *= -1
    beta = np.arcsin(vec[0])
    del vec

    rot1 = np.array([[1 , 0,  0],
                     [0, np.cos(alpha), np.sin(alpha)],
                     [0, -np.sin(alpha), np.cos(alpha)]])
    rot2 = np.array([[np.cos(beta), 0, -np.sin(beta)],
                     [0, 1, 0],
                     [np.sin(beta), 0, np.cos(beta)]])
    self.params['rot12'] = np.dot(rot2, rot1)
    del alpha, beta
 
    x, y, z = self.fields['x'], self.fields['y'], self.fields['z']
    vx, vy, vz = self.fields['vx'], self.fields['vy'], self.fields['vz']
 
    x_new = self.params['rot12'][0][0] * x + self.params['rot12'][0][1] * y + self.params['rot12'][0][2] * z
    y_new = self.params['rot12'][1][0] * x + self.params['rot12'][1][1] * y + self.params['rot12'][1][2] * z
    z_new = self.params['rot12'][2][0] * x + self.params['rot12'][2][1] * y + self.params['rot12'][2][2] * z
    vx_new = self.params['rot12'][0][0] * vx + self.params['rot12'][0][1] * vy + self.params['rot12'][0][2] * vz
    vy_new = self.params['rot12'][1][0] * vx + self.params['rot12'][1][1] * vy + self.params['rot12'][1][2] * vz
    vz_new = self.params['rot12'][2][0] * vx + self.params['rot12'][2][1] * vy + self.params['rot12'][2][2] * vz
    del rot1, rot2
    del x, y, z
    del vx, vy, vz
  
    self.fields['x'] = x_new
    self.fields['y'] = y_new
    self.fields['z'] = z_new
    self.fields['vx'] = vx_new
    self.fields['vy'] = vy_new
    self.fields['vz'] = vz_new
    del x_new, y_new, z_new, vx_new, vy_new, vz_new

    if len(self.new_fields['sinks']) > 0:
      x_sink, y_sink, z_sink = self.new_fields['sinks']['x'], self.new_fields['sinks']['y'], self.new_fields['sinks']['z']
      vx_sink, vy_sink, vz_sink = self.new_fields['sinks']['vx'], self.new_fields['sinks']['vy'], self.new_fields['sinks']['vz']

      x_new_sink = self.params['rot12'][0][0] * x_sink + self.params['rot12'][0][1] * y_sink + self.params['rot12'][0][2] * z_sink
      y_new_sink = self.params['rot12'][1][0] * x_sink + self.params['rot12'][1][1] * y_sink + self.params['rot12'][1][2] * z_sink
      z_new_sink = self.params['rot12'][2][0] * x_sink + self.params['rot12'][2][1] * y_sink + self.params['rot12'][2][2] * z_sink
      vx_new_sink = self.params['rot12'][0][0] * vx_sink + self.params['rot12'][0][1] * vy_sink + self.params['rot12'][0][2] * vz_sink
      vy_new_sink = self.params['rot12'][1][0] * vx_sink + self.params['rot12'][1][1] * vy_sink + self.params['rot12'][1][2] * vz_sink
      vz_new_sink = self.params['rot12'][2][0] * vx_sink + self.params['rot12'][2][1] * vy_sink + self.params['rot12'][2][2] * vz_sink
      del x_sink, y_sink, z_sink
      del vx_sink, vy_sink, vz_sink

      self.new_fields['sinks']['x'] = x_new_sink
      self.new_fields['sinks']['y'] = y_new_sink
      self.new_fields['sinks']['z'] = z_new_sink
      self.new_fields['sinks']['vx'] = vx_new_sink
      self.new_fields['sinks']['vy'] = vy_new_sink
      self.new_fields['sinks']['vz'] = vz_new_sink
      del x_new_sink, y_new_sink, z_new_sink, vx_new_sink, vy_new_sink, vz_new_sink

    if self.read_dm:
      x_dm, y_dm, z_dm = self.new_fields['dm']['x'], self.new_fields['dm']['y'], self.new_fields['dm']['z']
      vx_dm, vy_dm, vz_dm = self.new_fields['dm']['vx'], self.new_fields['dm']['vy'], self.new_fields['dm']['vz']

      x_new_dm = self.params['rot12'][0][0] * x_dm + self.params['rot12'][0][1] * y_dm + self.params['rot12'][0][2] * z_dm
      y_new_dm = self.params['rot12'][1][0] * x_dm + self.params['rot12'][1][1] * y_dm + self.params['rot12'][1][2] * z_dm
      z_new_dm = self.params['rot12'][2][0] * x_dm + self.params['rot12'][2][1] * y_dm + self.params['rot12'][2][2] * z_dm
      vx_new_dm = self.params['rot12'][0][0] * vx_dm + self.params['rot12'][0][1] * vy_dm + self.params['rot12'][0][2] * vz_dm
      vy_new_dm = self.params['rot12'][1][0] * vx_dm + self.params['rot12'][1][1] * vy_dm + self.params['rot12'][1][2] * vz_dm
      vz_new_dm = self.params['rot12'][2][0] * vx_dm + self.params['rot12'][2][1] * vy_dm + self.params['rot12'][2][2] * vz_dm
      del x_dm, y_dm, z_dm
      del vx_dm, vy_dm, vz_dm

      self.new_fields['dm']['x'] = x_new_dm
      self.new_fields['dm']['y'] = y_new_dm
      self.new_fields['dm']['z'] = z_new_dm
      self.new_fields['dm']['vx'] = vx_new_dm
      self.new_fields['dm']['vy'] = vy_new_dm
      self.new_fields['dm']['vz'] = vz_new_dm
      del x_new_dm, y_new_dm, z_new_dm, vx_new_dm, vy_new_dm, vz_new_dm

  def calculate_fields(self, field):

    if field == 'cs':
      if IO_U and IO_CHEM and IO_GAMMA:
        self.derived_fields[field] = np.sqrt(self.fields['gamma'] * BOLTZMANN *self.fields['temp'] / self.fields['mu'] / PROTONMASS) / 1.0e5
      else:
        print 'Energy, Gamma or Chemistry was not read!'

    # Sound crossing time
    if field == 'tcs_ratio':
      if IO_U and IO_CHEM and IO_GAMMA:
        try:
            self.derived_fields['cs'][0]
        except:
          self.calculate_fields('cs')
        try:
          self.derived_fields['tff'][0]
        except:
          self.calculate_fields('tff')
        if LengthUnit == 0:
          fac = 1
        elif LengthUnit == 1:
          fac = 1e3
        elif LengthUnit == 2:
          fac = UNIT_LENGTH / ASTRONOMICAL_UNIT
        else:
          print 'Length unit not implemented, we will assume pc'
          fac = 1
        self.derived_fields[field] = self.derived_fields['tff'] / (fac * self.fields['radius'] / self.derived_fields['cs'])

    # Free-fall time
    if field == 'tff':
      if IO_RHO:
        self.derived_fields[field] = np.sqrt(3. * np.pi / 32. / GRAVITY / self.fields['rho'])
      else:
        print 'Density was not read!'

    # Cooling rate
    if field == 'cool_rate':
      if IO_U and IO_RHO and IO_CHEM:
        try:
            self.derived_fields['esc_frac'][0]
        except:
          self.calculate_fields('esc_frac')
        abe = self.fields['abHII']
        abHI = np.maximum(1. - 2. * self.fields['abH2'] - self.fields['abHII'], 0.)
        self.derived_fields[field] = 7.5e-19 * np.exp(-1.18348e5 / self.fields['temp']) / (1. + np.sqrt(self.fields['temp'] / 1.0e5)) \
                                     * abHI * abe * self.fields['nh'] * self.fields['nh'] * self.derived_fields['esc_frac']
        del abe, abHI
      else:
        print 'Energy, Density or Chemistry was not read!'

    # Escape fraction
    if field == 'esc_frac':
      if IO_RHO:
        self.derived_fields[field] = np.exp(-self.fields['nh'] / 1.5e16)
      else:
        print 'Density was not read!'

    # Cooling time
    if field == 'tcool_ratio':
      if IO_U and IO_RHO:
        try:
          self.derived_fields['cool_rate'][0]
        except:
          self.calculate_fields('cool_rate')
        try:
          self.derived_fields['tff'][0]
        except:
          self.calculate_fields('tff')
        self.derived_fields[field] = self.fields['u'] * self.fields['rho'] / self.derived_fields['cool_rate'] / self.derived_fields['tff']
      else:
        print ' Energy or Density was not read!'

    if field == 'tcool':
      if IO_U and IO_RHO:
        try:
          self.derived_fields['cool_rate'][0]
        except:
          self.calculate_fields('cool_rate')
        self.derived_fields[field] = self.fields['u'] * self.fields['rho'] / self.derived_fields['cool_rate']
      else:
        print ' Energy or Density was not read!'

    # Entropy
    if field == 'entro':
      if IO_RHO and IO_U:
        self.derived_fields[field] = self.fields['mass'] * self.fields['rho']**(1. - GAMMA_ADB) * self.fields['temp']
      else:
        print ' Energy or Density was not read!'

    # Rosselanad mean opacity
    if field == 'kappa':
      if IO_RHO and IO_U:
        # Index in density
        drho, rho_idx_left = np.modf((np.log10(self.fields['rho']) - KAPPA_RHO[0]) / DELTA_KAPPA_RHO)
        rho_idx_left = rho_idx_left.astype(int)
        rho_idx_right = rho_idx_left + 1
        idx = np.where(rho_idx_left < 0)[0]
        rho_idx_left[idx] = 0
        rho_idx_right[idx] = 0
        drho[idx] = 0
        idx = np.where(rho_idx_right > len(KAPPA_RHO) - 1)[0]
        rho_idx_left[idx] = len(KAPPA_RHO) - 1
        rho_idx_right[idx] = len(KAPPA_RHO) - 1
        drho[idx] = 0

        # Index in temperature
        dtemp, temp_idx_left = np.modf((np.log10(self.fields['temp']) - KAPPA_TEMP[0]) / DELTA_KAPPA_TEMP)
        temp_idx_left = temp_idx_left.astype(int)
        temp_idx_right = temp_idx_left + 1
        idx = np.where(temp_idx_left < 0)[0]
        temp_idx_left[idx] = 0
        temp_idx_right[idx] = 0
        dtemp[idx] = 0
        idx = np.where(temp_idx_right > len(KAPPA_TEMP) - 1)[0]
        temp_idx_left[idx] = len(KAPPA_TEMP) - 1
        temp_idx_right[idx] = len(KAPPA_TEMP) - 1
        dtemp[idx] = 0
        del idx

        # Calculate values
        kappa_left = KAPPA_P[temp_idx_left, rho_idx_left] + drho * (KAPPA_P[temp_idx_left, rho_idx_right] - KAPPA_P[temp_idx_left, rho_idx_left])
        kappa_right = KAPPA_P[temp_idx_right, rho_idx_left] + drho * (KAPPA_P[temp_idx_right, rho_idx_right] - KAPPA_P[temp_idx_right, rho_idx_left])
        self.derived_fields[field] = 10**(kappa_left + dtemp * (kappa_right - kappa_left))
        del drho, rho_idx_left, rho_idx_right
        del dtemp, temp_idx_left, temp_idx_right
        del kappa_left, kappa_right 
      else:
        print 'Energy or Density was not read!'

from include import *

path = '/Volumes/DARK-MATTER/arepo/ah1w3f0/snapdir_241/'
file = 'ah1w3f0_241'
snapbase = path+file
snapfile = snapbase+'.0'
nsnaptypes = 6

f = open(snapfile, mode = 'rb')

# Read header
f.seek(28)
masstable = np.array(struct.unpack(repr(nsnaptypes) + 'd', f.read(48)))
time = struct.unpack('d', f.read(8))[0]
redshift = struct.unpack('d', f.read(8))[0]
flag_sfr = struct.unpack('i', f.read(4))[0]
flag_feedback = struct.unpack('i', f.read(4))[0]
npartall = np.array(struct.unpack(repr(nsnaptypes) + 'i', f.read(24)))
flag_cooling = struct.unpack('i', f.read(4))[0]
numfiles = struct.unpack('i', f.read(4))[0]
boxsize = struct.unpack('d', f.read(8))[0]
omega_m = struct.unpack('d', f.read(8))[0]
omega_lambda = struct.unpack('d', f.read(8))[0]
hubbleparam = struct.unpack('d', f.read(8))[0]
f.close()

npart = np.zeros((numfiles, nsnaptypes), dtype = 'i')

x = np.array([])
y = np.array([])
z = np.array([])
totmass = np.array([])
density = np.array([])
totmu = np.array([])
totgamma = np.array([])
totu = np.array([])

for i in range(numfiles):

  snapfile = snapbase + '.' + repr(i)
  f = open(snapfile, mode = 'rb')

  bytesleft = 256 - 6*4
  f.seek(4)
  npart[i, :] = np.array(struct.unpack(repr(nsnaptypes) + 'i', f.read(24)))
  f.seek(bytesleft,1)
  nbytes = struct.unpack('i', f.read(4))[0]
  
  ngas = npart[i, 0]

  # Read positions
  nbytes = struct.unpack('i', f.read(4))[0]
  if IO_POS:
    pos = np.array(struct.unpack(repr(3*ngas) + 'd', f.read(3*ngas*8)))
    x = np.append(x, pos[0::3]) 
    y = np.append(y, pos[1::3])
    z = np.append(z, pos[2::3])
    f.seek(nbytes-3*ngas*8,1)
  else:
    f.seek(nbytes,1)
  nbytes =  struct.unpack('i', f.read(4))[0]

  # Read velocities
  nbytes =  struct.unpack('i', f.read(4))[0]
  if IO_VEL:
    vel = np.array(struct.unpack(repr(3*ngas) + 'd', f.read(3*ngas*8)))
    velx = vel[0::3]
    vely = vel[1::3]
    velz = vel[2::3]
    f.seek(nbytes-3*ngas*8,1)
  else:
    f.seek(nbytes,1)
  nbytes =  struct.unpack('i', f.read(4))[0]

  # Read IDs
  nbytes = struct.unpack('i', f.read(4))[0]
  if IO_ID:
    id = np.array(struct.unpack(repr(ngas) + 'i', f.read(ngas*4)))
    f.seek(nbytes-ngas*4,1)
  else:
    f.seek(nbytes,1)
  nbytes = struct.unpack('i', f.read(4))[0]

  # Read mass
  nbytes = struct.unpack('i', f.read(4))[0]
  if IO_MASS:
    mass = np.array(struct.unpack(repr(ngas) + 'd', f.read(ngas*8)))
    totmass = np.append(totmass, mass)
    f.seek(nbytes-ngas*8,1)
  else:
    f.seek(nbytes,1)
  nbytes = struct.unpack('i', f.read(4))[0]

  # Read energy
  nbytes = struct.unpack('i', f.read(4))[0]
  if IO_U:
    fac = UNIT_ENERGY / UNIT_MASS
    u = fac*np.array(struct.unpack(repr(ngas) + 'd', f.read(ngas*8)))
    totu = np.append(totu, u)
    f.seek(nbytes-ngas*8,1)
  else:
    f.seek(nbytes,1)
  nbytes = struct.unpack('i', f.read(4))[0]

  # Read density
  nbytes = struct.unpack('i', f.read(4))[0]
  if IO_RHO:
    fac = UNIT_MASS / UNIT_LENGTH**3 * (1 + redshift)**3 * hubbleparam**2
    rho = fac * np.array(struct.unpack(repr(ngas) + 'd', f.read(ngas*8)))
    density = np.append(density, rho)
    f.seek(nbytes-ngas*8,1)
  else:
    f.seek(nbytes,1)
  nbytes = struct.unpack('i', f.read(4))[0]

  # Read volume
  nbytes = struct.unpack('i', f.read(4))[0]
  if IO_VOL:
    vol = np.array(struct.unpack(repr(ngas) + 'd', f.read(ngas*8)))
    f.seek(nbytes-ngas*8,1)
  else:
    f.seek(nbytes,1)
  nbytes = struct.unpack('i', f.read(4))[0]

  # Read delaunay
  nbytes = struct.unpack('i', f.read(4))[0]
  if IO_DELAUNAY:
    delaunay = np.array(struct.unpack(repr(ngas) + 'd', f.read(ngas*8)))
    f.seek(nbytes-ngas*8,1)
  else:
    f.seek(nbytes,1)
  nbytes = struct.unpack('i', f.read(4))[0]

  # Read gravitational acceleration
  nbytes =  struct.unpack('i', f.read(4))[0]
  if IO_GRAVACC:
    gravacc = np.array(struct.unpack(repr(3*ngas) + 'd', f.read(3*ngas*8)))
    #gravaccx = gravacc[0::3]
    #gravaccy = gravacc[1::3]
    #gravaccz = gravacc[2::3]
    f.seek(nbytes-3*ngas*8,1)
  else:
    f.seek(nbytes,1)
  nbytes =  struct.unpack('i', f.read(4))[0]

  # Read pressure gradient
  nbytes =  struct.unpack('i', f.read(4))[0]
  if IO_GRADP:
    gradp = np.array(struct.unpack(repr(3*ngas) + 'd', f.read(3*ngas*8)))
    #gradpx = gradp[0::3]
    #gradpy = gradp[1::3]
    #gradpz = gradp[2::3]
    f.seek(nbytes-3*ngas*8,1)
  else:
    f.seek(nbytes,1)
  nbytes =  struct.unpack('i', f.read(4))[0]

  # Read chemistry
  nbytes =  struct.unpack('i', f.read(4))[0]
  if IO_CHEM:
    chem = np.array(struct.unpack(repr(3*ngas) + 'd', f.read(3*ngas*8)))
    abHM = chem[0::3]
    abH2 = chem[1::3]
    abHII = chem[2::3]
    abe = abHII
    mu = (1. + 4. * HE_ABUND) / (1 + HE_ABUND - abH2 + abe)
    totmu = np.append(totmu, mu)
    f.seek(nbytes-3*ngas*8,1)
  else:
    f.seek(nbytes,1)
  nbytes =  struct.unpack('i', f.read(4))[0]

  # Read gamma
  nbytes = struct.unpack('i', f.read(4))[0]
  if IO_GAMMA:
    gamma = np.array(struct.unpack(repr(ngas) + 'd', f.read(ngas*8)))
    totgamma = np.append(totgamma, gamma)
    f.seek(nbytes-ngas*8,1)
  else:
    f.seek(nbytes,1)
  nbytes = struct.unpack('i', f.read(4))[0]

  f.close()

max_dens = np.where(density == max(density))[0]
centerx = x[max_dens]
centery = y[max_dens]
centerz = z[max_dens]
x = x - centerx
y = y - centery 
z = z - centerz

fac = UNIT_LENGTH / ASTRONOMICAL_UNIT

r = fac*np.sqrt(x**2 + y**2 + z**2)

r = r[np.nonzero(r)]
#logrmin = np.log10(min(r))
logrmin = np.log10(0.01)
logrmax = np.log10(max(r))

radial_bins = 100
rfac = radial_bins/(logrmax - logrmin)
r_idx = (rfac*(np.log10(r) - logrmin)).astype(int)

rad = logrmin + np.arange(radial_bins) / rfac

dens = np.zeros(radial_bins, float)
temp = np.zeros(radial_bins, float)
temperature = totmu * (totgamma - 1) * PROTONMASS / BOLTZMANN * totu

fac =  HYDROGEN_MASSFRAC / PROTONMASS

for j in np.unique(r_idx[r_idx > 0]):
  idx = np.where(r_idx == j)[0]
  dens[j-1] = sum(totmass[idx]*np.log10(fac*density[idx]))/sum(totmass[idx])
  temp[j-1] = sum(totmass[idx]*np.log10(temperature[idx]))/sum(totmass[idx])

pl.plot(rad, temp)
pl.savefig('test.png')
pl.show()

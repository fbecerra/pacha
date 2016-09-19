import numpy as np
import struct

path = '/Users/fbecerra/Classes/CS205/pavoreal/data/'
file = 'merger_330'
snapfile = path+file

f = open(snapfile, mode='rb')

nsnaptypes = 6
bytesleft = 256 - 6*4 - 6*8 - 8 - 8 - 2*4 - 6*4

nbytes = struct.unpack('i', f.read(4))[0]
npart = np.array(struct.unpack(repr(nsnaptypes) + 'i', f.read(24))) 
massarr = np.array(struct.unpack(repr(nsnaptypes) + 'd', f.read(48)))
time = struct.unpack('d', f.read(8))[0]
redshift = struct.unpack('d', f.read(8))[0]
flag_sfr = struct.unpack('i', f.read(4))[0]
flag_feedback = struct.unpack('i', f.read(4))[0]
npartTotal = np.array(struct.unpack(repr(nsnaptypes) + 'i', f.read(24)))

ntot = sum(npart)
ngas = npart[0]

f.seek(bytesleft,1)
nbytes = struct.unpack('i', f.read(4))[0]

# Read positions
nbytes = struct.unpack('i', f.read(4))[0]
pos = np.array(struct.unpack(repr(3*ngas) + 'f', f.read(3*ngas*4)))
x = pos[0::3]
y = pos[1::3]
z = pos[2::3]

# Normalize positions
xmin , xmax = np.min(x), np.max(x)
ymin , ymax = np.min(y), np.max(y)
zmin , zmax = np.min(z), np.max(z)
x = (x - xmin) / (xmax - xmin)
y = (y - ymin) / (ymax - ymin)
z = (z - zmin) / (zmax - zmin)

f.seek(nbytes-3*ngas*4,1)
nbytes =  struct.unpack('i', f.read(4))[0]

# Skip velocities
nbytes =  struct.unpack('i', f.read(4))[0]
#vel = np.array(struct.unpack(repr(3*ntot) + 'f', f.read(3*ntot*4)))
#f.seek(nbytes2-3*ntot*4,1)
f.seek(nbytes,1)
nbytes =  struct.unpack('i', f.read(4))[0]

# Skip IDs
nbytes = struct.unpack('i', f.read(4))[0]
#id = np.array(struct.unpack(repr(ntot) + 'i', f.read(ntot*4)))
#f.seek(nbytes2-ntot*4,1)
f.seek(nbytes,1)
nbytes = struct.unpack('i', f.read(4))[0]

# Skip energy
nbytes = struct.unpack('i', f.read(4))[0]
#u = np.array(struct.unpack(repr(ngas) + 'f', f.read(ngas*4)))
#f.seek(nbytes2-ngas*4,1)
f.seek(nbytes,1)
nbytes = struct.unpack('i', f.read(4))[0]

# Read density
nbytes = struct.unpack('i', f.read(4))[0]
rho = np.array(struct.unpack(repr(ngas) + 'f', f.read(ngas*4)))
rho = np.log10(rho)
f.seek(nbytes-ngas*4,1)
nbytes = struct.unpack('i', f.read(4))[0]

# Skip smoothing length
nbytes = struct.unpack('i', f.read(4))[0]
hsml = np.array(struct.unpack(repr(ngas) + 'f', f.read(ngas*4)))

print np.min(rho), np.max(rho)

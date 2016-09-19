from include import *
from utils import *

base
path
snapnum
snapbase = path + '/' + base + '/' + 'snapdir' + '_' + snapnum + '/' + base + '_' + snapnum
snapfile=snapbase + '.0'

#### Reading Header ####


f = open(self.snapfile, mode = 'rb')

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

for i in range(numfiles):

    if numfiles > 1:

        snapfile = snapbase + '.' + repr(i)

    f = open(snapfile, mode = 'rb')

    f.seek(4)

    npart[i, :] = np.array(struct.unpack(repr(nsnaptypes) + 'i', f.read(24)))

    f.close()

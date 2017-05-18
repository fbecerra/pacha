from include import *

class Fof:

  def __init__(self):
    self.params = dict()
    self.fields = dict()

  def read_header(self, groupfile):

    self.file = open(groupfile, mode = 'rb')
    
    self.params['nbytes'] = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
    self.params['ngroups'] = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
    self.params['nsubgroups'] = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
    self.params['nids'] = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
    self.params['totngroups'] = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
    self.params['totnsubgroups'] = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]

    dummy = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
    self.params['totnids'] = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
    dummy2 = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
    self.params['numfiles'] = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
    dummy3 = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]

    self.params['time'] = np.fromfile(self.file, dtype=np.dtype('d'), count=1)[0]
    self.params['redshift'] = np.fromfile(self.file, dtype=np.dtype('d'), count=1)[0]
    self.params['hubbleparam'] = np.fromfile(self.file, dtype=np.dtype('d'), count=1)[0]
    self.params['boxsize'] = np.fromfile(self.file, dtype=np.dtype('d'), count=1)[0]
    self.params['omega_m'] = np.fromfile(self.file, dtype=np.dtype('d'), count=1)[0]
    self.params['omega_lambda'] = np.fromfile(self.file, dtype=np.dtype('d'), count=1)[0]
    self.params['flag_doubleprecision'] = np.fromfile(self.file, dtype=np.dtype('d'), count=1)[0]
    
    self.file.close()

  def read_groups(self, groupbase):

    try:
      self.params['flag_doubleprecision']
    except:
      self.read_header(groupbase + '.0')

    if self.params['flag_doubleprecision'] == 0:
      float_type = 'f'
    else:
      float_type = 'd'

    self.fields['len'] = np.zeros(self.params['totngroups'], 'i')
    self.fields['mass'] = np.zeros(self.params['totngroups'], float_type)
    self.fields['pos'] = np.zeros(3 * self.params['totngroups'], float_type)
    self.fields['cm'] = np.zeros(3 * self.params['totngroups'], float_type)
    self.fields['vel'] = np.zeros(3 * self.params['totngroups'], float_type)
    self.fields['lentype'] = np.zeros(6 * self.params['totngroups'], 'i')
    self.fields['masstype'] = np.zeros(6 * self.params['totngroups'], float_type)
    self.fields['ids'] = np.zeros(self.params['totnids'], 'i')

    count = 0
    countids = 0

    # Iterate over files to read all groups
    for i in range(self.params['numfiles']):

      self.groupfile = groupbase + '.' + repr(i)
      self.read_header(self.groupfile)
      self.file = open(self.groupfile, mode = 'rb')

      self.blksize = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
      self.ngroups = np.fromfile(self.file, dtype=np.dtype('i'), count=1)[0]
      self.file.seek(self.blksize + 12)

      if self.ngroups != 0:

        self.fields['len'][count: count + self.ngroups] = np.fromfile(self.file, dtype = 'i', count = self.ngroups)
        self.file.seek(8, 1)

        self.fields['mass'][count: count + self.ngroups] = np.fromfile(self.file, dtype = float_type, count = self.ngroups)
        self.file.seek(8, 1)

        self.fields['pos'][3 * count: 3 * count + 3 * self.ngroups] = np.fromfile(self.file, dtype = float_type, count = 3 * self.ngroups)
        self.file.seek(8, 1)

        self.fields['cm'][3 * count: 3 * count + 3 * self.ngroups] = np.fromfile(self.file, dtype = float_type, count = 3 * self.ngroups)
        self.file.seek(8, 1)

        self.fields['vel'][3 * count: 3 * count + 3 * self.ngroups] = np.fromfile(self.file, dtype = float_type, count = 3 * self.ngroups)
        self.file.seek(8, 1)

        self.fields['lentype'][6 * count: 6 * count + 6 * self.ngroups] = np.fromfile(self.file, dtype = 'i', count = 6 * self.ngroups)
        self.file.seek(8, 1)

        self.fields['masstype'][6 * count: 6 * count + 6 * self.ngroups] = np.fromfile(self.file, dtype = float_type, count = 6 * self.ngroups)
        self.file.seek(8, 1)

      if self.params['nids'] != 0:
        self.fields['ids'][countids: countids + self.params['nids']] = np.fromfile(self.file, dtype = 'i', count = self.params['nids'])

      self.file.close()

      count += self.ngroups
      countids += self.params['nids']

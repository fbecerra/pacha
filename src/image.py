from include import *
from cutils import *

class Image:

#  def __init__(self):

  def calculate_image(self, snap, field, ImgWidth=1):
    self.xbins = ImgXBins
    self.ybins = self.xbins / ImgScreenRatio
    self.width = ImgWidth
    self.height = self.width / ImgScreenRatio 
    self.size = max(self.width, self.height)
    self.img = np.zeros([self.xbins, self.ybins], float)
    sum = np.zeros([self.xbins, self.ybins], float)
    if IO_POS and IO_RHO and IO_VOL:
      ax1 = snap.fields['x']
      ax2 = snap.fields['y']
      ax3 = snap.fields['z']
      idx = np.where((np.abs(ax1) < self.width / 2.) & (np.abs(ax2) < self.height / 2.) & (np.abs(ax3) < self.size / 2.))[0]
      min1 = np.minimum(np.maximum((ax1[idx] + self.width / 2. - snap.fields['hsml'][idx]) / self.width * self.xbins, 0).astype(int), self.xbins - 1).astype(int)
      max1 = np.minimum(np.maximum((ax1[idx] + self.width / 2. + snap.fields['hsml'][idx]) / self.width * self.xbins, 0).astype(int), self.xbins - 1).astype(int)
      min2 = np.minimum(np.maximum((ax2[idx] + self.height / 2. - snap.fields['hsml'][idx]) / self.height * self.ybins, 0).astype(int), self.ybins - 1).astype(int)
      max2 = np.minimum(np.maximum((ax2[idx] + self.height / 2. + snap.fields['hsml'][idx]) / self.height * self.ybins, 0).astype(int), self.ybins - 1).astype(int)
      calculate_image(snap, snap.fields[field], min1, max1, min2, max2, idx, sum, self.img)
      del ax1, ax2, ax3
      del idx
      del min1, max1, min2, max2
      self.img = np.rot90(np.log10(self.img / sum))
      del sum
    else:
      print 'Required fields not read!'

  def pspace_plot(self, snap, field):
    self.bins = PSpaceBins
    if IO_RHO:
      self.xmin = np.log10(np.min(snap.fields['nh']))
      self.xmax = np.log10(np.max(snap.fields['nh']))

      if field in fields_list:
        yvalues = snap.fields[field]
      elif field in derived_fields_list:
        try:
          snap.derived_fields[field][0]
        except:
          snap.calculate_fields(field)
        yvalues = snap.derived_fields[field]

      self.ymin = np.log10(np.min(yvalues))
      self.ymax = np.log10(np.max(yvalues))
      self.valy = np.zeros(self.bins, float)
      sumy = np.zeros(self.bins, float)
      self.img = np.zeros([self.bins, self.bins], float)
      sum = np.zeros([self.bins, self.bins], float)
      idx = np.where((self.xmin <= np.log10(snap.fields['nh'])) & (np.log10(snap.fields['nh']) <= self.xmax) &
                     (self.ymin <= np.log10(yvalues)) & (np.log10(yvalues) <= self.ymax))[0]
      xidx = np.maximum(np.minimum((np.log10(snap.fields['nh'][idx]) - self.xmin) / (self.xmax - self.xmin) * self.bins , self.bins - 1).astype(int), 0).astype(int)
      yidx = np.maximum(np.minimum((np.log10(yvalues[idx]) - self.ymin) / (self.ymax - self.ymin) * self.bins , self.bins - 1).astype(int), 0).astype(int)
      calculate_phase_space(snap, np.log10(yvalues), xidx, yidx, idx, self.valy, sumy, self.img, sum)
      totval = np.sum(self.img)

      self.valx = self.xmin + (np.arange(self.bins) + 0.5) * (self.xmax - self.xmin) / self.bins
      self.valy = self.valy / sumy
      self.img = np.rot90(np.log10(self.img / totval))
      del sumy, sum, idx, totval

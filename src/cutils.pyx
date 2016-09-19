cimport cython
cimport numpy as np
import numpy as np

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

@cython.boundscheck(False)
@cython.wraparound(False)
def calculate_image(snap, np.ndarray[DTYPE_t, ndim=1] field,
                    np.ndarray[np.int64_t, ndim=1] min1, np.ndarray[np.int64_t, ndim=1] max1, 
                    np.ndarray[np.int64_t, ndim=1] min2, np.ndarray[np.int64_t, ndim=1] max2, 
                    np.ndarray[np.int64_t, ndim=1] idx, 
                    np.ndarray[DTYPE_t, ndim=2] sum, np.ndarray[DTYPE_t, ndim=2] img):

  cdef unsigned int i, j ,k
  cdef int length = idx.shape[0]
  cdef DTYPE_t sum_part, img_part

  for i in xrange(0, length):
    sum_part = snap.fields['nh'][idx[i]] * snap.fields['nh'][idx[i]] * snap.fields['mass'][idx[i]]
    img_part = sum_part * field[idx[i]]
    for j in xrange(min1[i], max1[i] + 1):
      for k in xrange(min2[i], max2[i] + 1):
        sum[j, k] += sum_part
        img[j, k] += img_part

@cython.boundscheck(False)
@cython.wraparound(False)
def calculate_phase_space(snap, np.ndarray[DTYPE_t, ndim=1] field,
                          np.ndarray[np.int64_t, ndim=1] xidx, np.ndarray[np.int64_t, ndim=1] yidx,
                          np.ndarray[np.int64_t, ndim=1] idx,
                          np.ndarray[DTYPE_t, ndim=1] valy, np.ndarray[DTYPE_t, ndim=1] sumy,
                          np.ndarray[DTYPE_t, ndim=2] val, np.ndarray[DTYPE_t, ndim=2] sum):

  cdef unsigned int i, j
  cdef int length = idx.shape[0]
  cdef DTYPE_t mass

  for i in xrange(0, length):
    j = idx[i]
    mass = snap.fields['mass'][j]
    valy[xidx[i]] += mass * field[j]
    sumy[xidx[i]] += mass
    val[xidx[i], yidx[i]] += mass
    sum[xidx[i], yidx[i]] += 1

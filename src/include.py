import matplotlib
matplotlib.use('Agg')
from constants import *
from fields import *
from params import *
#from mpi4py import MPI
from colormaps import *
import numpy as np
import h5py
import math
import re
import random
import scipy as sc
import pylab as pl
pl.ion()
from matplotlib import rc
fontsize=14
rc('text', usetex=True)
rc('font', **{'family':'serif','serif':'Computer Modern Roman', 'size':fontsize})
rc('axes', labelsize=fontsize, linewidth=1.5)
rc('legend', fontsize=fontsize, numpoints=1, frameon=False)
rc('xtick', labelsize=fontsize)
rc('ytick', labelsize=fontsize)
rc('lines', lw=2.0, mew=0.3)
rc('grid', linewidth=2.0)

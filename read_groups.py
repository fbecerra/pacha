import matplotlib
matplotlib.use('TkAgg')
import src as pr
import numpy as np
import pylab as pl
import time

path = '/path/to/output/folder/'
file = 'mh1dm'
initial_snap = 25
final_snap = 90
threshold_mass = 5e5

for num_fof in range(initial_snap, final_snap):

  groupbase = path + '/' + file + '/groups_%03i/fof_tab_%03i' %(num_fof, num_fof)

  MyGroup = pr.fof.Fof()
  MyGroup.read_groups(groupbase)
  print 'Searching at redshift z=', MyGroup.params['redshift']

  # No halos yet
  if len(MyGroup.fields['mass']) <= 0:
    continue

  # Any halo above mass threshold?
  if np.max(MyGroup.fields['mass']) / MyGroup.params['hubbleparam'] > threshold_mass * pr.constants.SOLAR_MASS / pr.constants.UNIT_MASS:
    print 'Halo found at time t=', MyGroup.params['time'] * pr.constants.UNIT_TIME / pr.constants.SEC_PER_YEAR, 'yr in snapshot ', num_fof
    break
 
# No halos above threshols
if (num_fof == final_snap) and (np.max(MyGroup.fields['mass']) / MyGroup.params['hubbleparam'] < threshold_mass * pr.constants.SOLAR_MASS / pr.constants.UNIT_MASS):
  print "Halo not found :("
  sys.exit()

# Halo properties 
halo_mass = np.max(MyGroup.fields['mass'])
mmg = np.where(MyGroup.fields['mass'] == halo_mass)[0]
groups_pos = MyGroup.fields['pos'].reshape((len(MyGroup.fields['pos'])/3.,-1))
halo_pos  = groups_pos[mmg,:]/MyGroup.params['hubbleparam'] / (1 + MyGroup.params['redshift'])   #[kpc] and comoving coordinates
groups_vel = MyGroup.fields['vel'].reshape((len(MyGroup.fields['vel'])/3.,-1))
halo_vel   = groups_vel[mmg,:] * (1 + MyGroup.params['redshift'])                #[km s-1] and comoving coordinates

# Output properties
print "\t***Halo properties***"
print "Position\t\t:", groups_pos[mmg,:][0], " in code units and ", halo_pos[0], " in kpc"
print "Velocity\t\t:", groups_vel[mmg,:][0], " in code units and ", halo_vel[0], " in km*s-1"
print "Mass\t\t\t:", halo_mass, " in code units and ", halo_mass * pr.constants.UNIT_MASS / pr.constants.SOLAR_MASS / MyGroup.params['hubbleparam'], " in Msun"

group_offset = np.zeros(MyGroup.params['totngroups'], int)
for i in xrange(1, MyGroup.params['totngroups']):
  group_offset[i] = group_offset[i-1] + MyGroup.fields['len'][i-1]
offset = group_offset[mmg][0]
left = MyGroup.fields['len'][mmg][0]

## Read snap
#snapbase = path + file + '/snapdir_%03i/' % num_fof + file + '_%03i' %num_fof
#print snapbase
#MySnap = pr.snap.Snap()
#MySnap.read_header(snapbase)
#MySnap.read_fields(snapbase)

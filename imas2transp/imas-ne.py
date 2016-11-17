import imas
#import numpy
#import sys
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 18})

# This program gets data from the DB entry, just for practicing the PUAL GET command
# It servers also as a nested of 3 level nested AoS (type 3 at the top,
# type 2 below)

shot = 386; run = 1
treename = 'ids'

imas_obj = imas.ids(shot, run, 0, 0)

imas_obj.open()# Create a new instance of database
#imas_obj.open_env(' ', 'test', '2.0')
# write(*,*) 'Opened MDS pulse file, idx = ', idx

ids = imas_obj.core_profiles

print 'Test of the GET_SLICE function'
interpol = 1 # Interpolation mode = closest neighbour

twant = 45.0
ids.getSlice(twant, interpol)

print 'IDS_Properties homogeneous : ' + str(ids.ids_properties.homogeneous_time)
print 'IDS_Properties comment : ' + ids.ids_properties.comment
print 'Size of ids.profiles_1d : ' + str(len(ids.profiles_1d))
# print 'ids.profiles_1d.time : ",ids.profiles_1d(:).tim
print 'Main IDS time : ' + str(ids.time)
print 'Ip = ' + str(ids.global_quantities.ip)

#print len(ids.profiles_1d); sys.exit(0)

i = 0; print dir(ids.profiles_1d[i])
print 'Time slice [' + str(i) + '] = ' + str(ids.profiles_1d[i].time)
#print 'rho = ' + str(ids.profiles_1d[i].grid.rho_tor_norm)
#print 'ne = ' + str(ids.profiles_1d[i].electrons.density)
#sys.exit(0)

plt.plot(ids.profiles_1d[i].grid.rho_tor_norm, ids.profiles_1d[i].electrons.density, 'b')
plt.xlabel('ids.profiles_1d[i].grid.rho_tor_norm')
plt.ylabel('ids.profiles_1d[i].electrons.density')
#plt.grid(True)
plt.ylim([1.8e19, 4.8e19])
plt.tight_layout()
plt.savefig('imas-ne.png', format = 'png')
plt.show()
plt.close()

imas_obj.close()

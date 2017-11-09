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
imas_obj.open() # Create a new instance of database
ids = imas_obj.equilibrium

print 'Test of the GET_SLICE function'
interpol = 1 # Interpolation mode = closest neighbour
twant = 0.0
ids.getSlice(twant, interpol)

i = 0
plt.plot(ids.time_slice[i].global_quantities.magnetic_axis/b_field_tor)
#plt.xlabel('ids.profiles_1d[i].grid.rho_tor_norm')
#plt.ylabel('ids.profiles_1d[i].electrons.density')
#plt.grid(True)
#plt.ylim([1.8e19, 4.8e19])
plt.tight_layout()
plt.savefig('imas-b_field_tor.png', format = 'png')
plt.show()
plt.close()

imas_obj.close()

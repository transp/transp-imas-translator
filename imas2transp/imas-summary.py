import imas
#import numpy
import sys
#import matplotlib
import matplotlib.pyplot as plt

#matplotlib.use('PS')

plt.rcParams.update({'font.size': 16})

# This program gets data from the DB entry, just for practicing the PUAL GET command
# It servers also as a nested of 3 level nested AoS (type 3 at the top,
# type 2 below)

shot = 385; run = 30
treename = 'ids'

imas_obj = imas.ids(shot, run, 0, 0)
imas_obj.open()
ids = imas_obj.summary

print 'Test of the GET_SLICE function'
interpol = 1 # Interpolation mode = closest neighbour
twant = 1000.0
ids.getSlice(twant, interpol)

print 'IDS_Properties homogeneous : ' + str(ids.ids_properties.homogeneous_time)
print 'IDS_Properties comment : ' + ids.ids_properties.comment

print ids.code.name
print ids.code.version, '\n'

print ids.global_quantities.ip.value
print ids.global_quantities.current_ohm.value
print ids.global_quantities.current_non_inductive.value
print ids.global_quantities.current_bootstrap.value, '\n'

print ids.global_quantities.energy_total.value
print ids.global_quantities.energy_thermal.value, '\n'

print ids.global_quantities.tau_energy.value
#print ids.global_quantities.h_98.value
#print ids.global_quantities.h_mode.value, '\n'

sys.exit(0)

plt.plot(ids.global_quantities.ip.value, 'b')
plt.xlabel('ids.global_quantities.ip.value')
#plt.grid(True)
#plt.ylim([1.8e19, 4.8e19])
plt.tight_layout()
plt.savefig('imas-summary.png', format = 'png')
#plt.show()
plt.close()

imas_obj.close()

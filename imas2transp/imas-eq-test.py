import imas
import matplotlib.pyplot as plt

#shot = 386; run = 1
#shot = 120; run = 28
shot = 385; run = 30
itm = imas.ids(shot, run)
itm.open()

# A'ight, let's see what's in there...
print 'Every available IDS:', itm.__dict__.keys()

# We only want the equilibrium IDS for now
ids = itm.equilibrium

# Tell the IDS what time slice you want data for
twant = 1000.0 # Greater or equal than final time
interpol = 1 # Interpolation mode = nearest neighbor
ids.getSlice(twant, interpol)

print '\nDumping last-time-slice contents of each data object in equilibrium IDS that was written by transp2imas:'

print '\nMembers of container object equilibrium are:\n', dir(ids)
print 'equilibrium%time[-1] =', ids.time

print '\nMembers of container object equilibrium.time_slice[-1]%profiles_1d are:\n', dir(ids.time_slice[-1].profiles_1d)
print 'equilibrium.time_slice[-1]%profiles_1d.psi =', ids.time_slice[-1].profiles_1d.psi

for i in range(6):
 print i, ids.time_slice[-1].profiles_1d.psi[i + 1] - ids.time_slice[-1].profiles_1d.psi[i]

plt.plot(ids.time_slice[-1].profiles_1d.psi)
#plt.xlabel('X')
#plt.ylabel('TE')
#plt.grid(True)
#plt.tight_layout()
#plt.savefig('te-final6-140358X63.png', format = 'png')
plt.show()
#plt.close()

itm.close()

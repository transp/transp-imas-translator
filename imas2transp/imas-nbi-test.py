import imas
import matplotlib.pyplot as plt
import sys

shot = 385; run = 30
itm = imas.ids(shot, run)
itm.open()

# A'ight, let's see what's in there...
print 'Every available IDS:', itm.__dict__.keys()

# We only want the nbi IDS for now
ids = itm.nbi

#sys.exit(0)

# NEVER try to tell the nbi IDS what time slice you want data for!
# I did and paging-palooza ensued...
#twant = 1000.0 # Greater or equal than final time
#interpol = 1 # Interpolation mode = nearest neighbor
#ids.getSlice(twant, interpol)

print 'I don\'t know yet how to load nbi IDS data into python. ids.getSlice() causes severe paging so can\'t be called. But it does some hidden initialization that must be done before accessing data.'
sys.exit(0)

print '\nDumping last-time-slice contents of each data object in nbi IDS that was written by transp2imas:'

#sys.exit(0)

print '\nMembers of container object nbi are:\n', dir(ids)
print 'nbi =', ids

#sys.exit(0)

#print '\nMembers of container object nbi.unit are:\n', dir(ids.unit)
#sys.exit(0)
#print '\nMembers of container object nbi.unit[0].beamlets_group[0] are:\n', dir(ids.unit[0].beamlets_group[0])

print 'nbi.unit[0].beamlets_group[0].beamlets.tangency_radii[0] =', ids.unit[0].beamlets_group[0].beamlets.tangency_radii[0]
print 'nbi.unit[1].beamlets_group[0].beamlets.tangency_radii[0] =', ids.unit[1].beamlets_group[0].beamlets.tangency_radii[0]
print 'nbi.unit[0].beamlets_group[0].tangency_radius =', ids.unit[0].beamlets_group[0].tangency_radius
print 'nbi.unit[1].beamlets_group[0].tangency_radius =', ids.unit[1].beamlets_group[0].tangency_radius

itm.close()

import imas
#import numpy
import sys

shot = 386; run = 1
#shot = 120; run = 28
itm = imas.ids(shot, run)
itm.open()

# A'ight, let's see what's in there...
print 'Every available IDS:', itm.__dict__.keys()

# We only want the core_profiles IDS for now
ids = itm.core_profiles

# I expected this to return an array populated with all scalar time slices, but no....
#print 'ids.time =', ids.time, '(why is this an empty array?!?)'

# Tell the IDS what time slice you want data for
#print 'Test of the getSlice() function'
twant = 1000.0 # Greater or equal than final time
interpol = 1 # Interpolation mode = nearest neighbor
ids.getSlice(twant, interpol)

print '\nDumping last-time-slice contents of all data objects in core_profiles IDS that were written by transp2imas:'
#print dir(ids)
# Check that we got the time slice we wanted
print 'core_profiles%time(-1) =', ids.time

#print dir(ids.global_quantities)
print 'core_profiles%global_quantities%v_loop(-1) =', ids.global_quantities.v_loop

itm.close()

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
print 'core_profiles%global_quantities%li(-1) =', ids.global_quantities.li
print 'core_profiles%global_quantities%beta_pol(-1) =', ids.global_quantities.beta_pol

#print dir(ids.profiles_1d[-1])
print 'core_profiles%profiles_1d(-1)%t_i_average =', ids.profiles_1d[-1].t_i_average
print 'core_profiles%profiles_1d(-1)%n_i_total_over_n_e =', ids.profiles_1d[-1].n_i_total_over_n_e
print 'core_profiles%profiles_1d(-1)%zeff =', ids.profiles_1d[-1].zeff
print 'core_profiles%profiles_1d(-1)%j_tor =', ids.profiles_1d[-1].j_tor
print 'core_profiles%profiles_1d(-1)%j_ohmic =', ids.profiles_1d[-1].j_ohmic
print 'core_profiles%profiles_1d(-1)%j_bootstrap =', ids.profiles_1d[-1].j_bootstrap
print 'core_profiles%profiles_1d(-1)%q =', ids.profiles_1d[-1].q
print 'core_profiles%profiles_1d(-1)%magnetic_shear =', ids.profiles_1d[-1].magnetic_shear
print 'core_profiles%profiles_1d(-1)%time =', ids.profiles_1d[-1].time

itm.close()

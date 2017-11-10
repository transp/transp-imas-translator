import imas

shot = 386; run = 1
#shot = 120; run = 28
#shot = 385; run = 30
itm = imas.ids(shot, run)
itm.open()

# A'ight, let's see what's in there...
print 'Every available IDS:', itm.__dict__.keys()

# We only want the core_profiles IDS for now
ids = itm.core_profiles

# I expected this to return an array populated with all scalar time slices, but no....
#print 'ids.time =', ids.time, '(why is this an empty array?!?)'

# Tell the IDS what time slice you want data for
twant = 1000.0 # Greater or equal than final time
interpol = 1 # Interpolation mode = nearest neighbor
ids.getSlice(twant, interpol)

print '\nDumping last-time-slice contents of each data object in core_profiles IDS that was written by transp2imas:'

print '\nMembers of container object core_profiles are:\n', dir(ids)
print 'core_profiles%time[-1] =', ids.time

print '\nMembers of container object core_profiles.global_quantities are:\n', dir(ids.global_quantities)
print 'core_profiles%global_quantities%v_loop[-1] =', ids.global_quantities.v_loop
print 'core_profiles%global_quantities%li[-1] =', ids.global_quantities.li
print 'core_profiles%global_quantities%beta_pol[-1] =', ids.global_quantities.beta_pol

print '\nMembers of container object core_profiles.profiles_1d[-1] are:\n', dir(ids.profiles_1d[-1])
print 'core_profiles%profiles_1d[-1]%t_i_average =', ids.profiles_1d[-1].t_i_average
print 'core_profiles%profiles_1d[-1]%n_i_total_over_n_e =', ids.profiles_1d[-1].n_i_total_over_n_e
print 'core_profiles%profiles_1d[-1]%zeff =', ids.profiles_1d[-1].zeff
print 'core_profiles%profiles_1d[-1]%j_tor =', ids.profiles_1d[-1].j_tor
print 'core_profiles%profiles_1d[-1]%j_ohmic =', ids.profiles_1d[-1].j_ohmic
print 'core_profiles%profiles_1d[-1]%j_bootstrap =', ids.profiles_1d[-1].j_bootstrap
print 'core_profiles%profiles_1d[-1]%q =', ids.profiles_1d[-1].q
print 'core_profiles%profiles_1d[-1]%magnetic_shear =', ids.profiles_1d[-1].magnetic_shear
print 'core_profiles%profiles_1d[-1]%time =', ids.profiles_1d[-1].time

print '\nMembers of container object core_profiles.profiles_1d[-1].grid are:\n', dir(ids.profiles_1d[-1].grid)
print 'core_profiles%profiles_1d[-1]%grid%rho_tor_norm =', ids.profiles_1d[-1].grid.rho_tor_norm
print 'core_profiles%profiles_1d[-1]%grid%rho_tor =', ids.profiles_1d[-1].grid.rho_tor
print 'core_profiles%profiles_1d[-1]%grid%psi =', ids.profiles_1d[-1].grid.psi
print 'core_profiles%profiles_1d[-1]%grid%volume =', ids.profiles_1d[-1].grid.volume
print 'core_profiles%profiles_1d[-1]%grid%area =', ids.profiles_1d[-1].grid.area

print '\nMembers of container object core_profiles.profiles_1d[-1].electrons are:\n', dir(ids.profiles_1d[-1].electrons)
print 'core_profiles%profiles_1d[-1]%electrons%temperature =', ids.profiles_1d[-1].electrons.temperature
print 'core_profiles%profiles_1d[-1]%electrons%density =', ids.profiles_1d[-1].electrons.density

print '\nMembers of container object core_profiles.profiles_1d[-1].ion[1] are:\n', dir(ids.profiles_1d[-1].ion[1])
print 'core_profiles%profiles_1d[-1]%ion[1]%z_ion =', ids.profiles_1d[-1].ion[1].z_ion
print 'core_profiles%profiles_1d[-1]%ion[1]%label =', ids.profiles_1d[-1].ion[1].label
print 'core_profiles%profiles_1d[-1]%ion[1]%temperature =', ids.profiles_1d[-1].ion[1].temperature
print 'core_profiles%profiles_1d[-1]%ion[1]%density =', ids.profiles_1d[-1].ion[1].density
print 'core_profiles%profiles_1d[-1]%ion[2]%density_fast =', ids.profiles_1d[-1].ion[2].density_fast
print 'core_profiles%profiles_1d[-1]%ion[2]%pressure_fast_perpendicular =', ids.profiles_1d[-1].ion[2].pressure_fast_perpendicular
print 'core_profiles%profiles_1d[-1]%ion[2]%pressure_fast_parallel =', ids.profiles_1d[-1].ion[2].pressure_fast_parallel

print '\nMembers of container object core_profiles.profiles_1d[-1].ion[1].element[0] are:\n', dir(ids.profiles_1d[-1].ion[1].element[0])
print 'core_profiles%profiles_1d[-1]%ion[1]%element[0]%a =', ids.profiles_1d[-1].ion[1].element[0].a
print 'core_profiles%profiles_1d[-1]%ion[1]%element[0]%z_n =', ids.profiles_1d[-1].ion[1].element[0].z_n

print '\nMembers of container object core_profiles.profiles_1d[-1].ion[4].state[0] are:\n', dir(ids.profiles_1d[-1].ion[4].state[0])
print 'core_profiles%profiles_1d[-1]%ion[4]%state[0]%label =', ids.profiles_1d[-1].ion[4].state[0].label
print 'core_profiles%profiles_1d[-1]%ion[4]%state[0]%temperature =', ids.profiles_1d[-1].ion[4].state[0].temperature
print 'core_profiles%profiles_1d[-1]%ion[4]%state[0]%density =', ids.profiles_1d[-1].ion[4].state[0].density

print '\nIon species:'
for i in range(len(ids.profiles_1d[-1].ion)):
  print str(i + 1) + ': ' + ids.profiles_1d[-1].ion[i].label + ', a = ' + str(ids.profiles_1d[-1].ion[i].element[0].a) + ', z_n = ' + str(ids.profiles_1d[-1].ion[i].element[0].z_n)

itm.close()

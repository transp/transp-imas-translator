module transp2imas_module

	use ids_schemas
	use ids_routines

	type (ids_equilibrium) :: eq
	type (ids_core_profiles) :: cp
	type (ids_core_transport) :: ct
	type (ids_edge_sources) :: es
	type (ids_nbi) :: nbi
	type (ids_summary) :: sum

	!type (ids_pf_active) :: pf
	!type (ids_actuator) :: actuator
	!type (ids_atomic_data) :: atomic
	!type (ids_controllers) :: controllers
	!type (ids_core_sources) :: cs
	!type (ids_em_coupling) :: em
	!type (ids_magnetics) :: mag
	!type (ids_pf_passive) :: pfpassive
	!type (ids_schedule) :: schedule
	!type (ids_sdn) :: sdn
	!type (ids_simulation) :: sim
	!type (ids_temporary) :: tmp
	!type (ids_tf) :: tf

	integer :: whichtimeslice
	integer :: whichprofile
	real*8 :: timeofinterest

end module transp2imas_module

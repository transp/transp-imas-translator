module transp2ids_module
   use ids_schemas
   use ids_routines

   type (ids_core_profiles) :: cp  ! Declaration of the empty ids cp to be filled
   type (ids_equilibrium) :: eq
   type (ids_pf_active) :: pf  ! Declaration of the empty ids pf to be filled

   type (ids_actuator) :: actuator  ! Declaration of the empty ids pf to be filled
   type (ids_atomic_data) :: atomic  ! Declaration of the empty ids pf to be filled
   type (ids_controllers) :: controllers  ! Declaration of the empty ids pf to be filled
   type (ids_core_sources) :: cs  ! Declaration of the empty ids pf to be filled
   type (ids_core_transport) :: ct  ! Declaration of the empty ids pf to be filled
   type (ids_em_coupling) :: em  ! Declaration of the empty ids pf to be filled
   type (ids_magnetics) :: mag  ! Declaration of the empty ids pf to be filled
   type (ids_pf_passive) :: pfpassive  ! Declaration of the empty ids pf to be filled
   type (ids_schedule) :: schedule  ! Declaration of the empty ids pf to be filled
   type (ids_sdn) :: sdn  ! Declaration of the empty ids pf to be filled
!   type (ids_simulation) :: sim  ! Declaration of the empty ids pf to be filled
   type (ids_temporary) :: tmp  ! Declaration of the empty ids pf to be filled
   type (ids_tf) :: tf  ! Declaration of the empty ids pf to be filled

      integer whichtimeslice
      integer whichprofile
      real*8 timeofinterest
end module transp2ids_module

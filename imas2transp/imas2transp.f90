program imas2transp

   use ids_schemas              !! These are the Fortran type definitions for the Physics Data Model
   use ids_routines             !! These are the Access Layer routines + management of IDS structures

   use transp_ufiles_0d         !! transp ufile scalar type
   use transp_ufiles_1d         !! transp ufile 1D array
   use transp_ufiles_2d         !! transp ufile 2D array
   use transp_ufiles_3d         !! transp ufile 3D array

   implicit none

   !internal variables
   type (ids_core_profiles) :: cp ! Declaration of the empty ids variables to be filled
   type (ids_equilibrium)   :: eq
   type (ids_nbi)           :: nbi
   type (ids_pf_active)     :: pf

   type(transp_ufiles_0d_data) :: t_uf0d
   type(transp_ufiles_1d_data) :: t_uf1d
   type(transp_ufiles_2d_data) :: t_uf2d
   type(transp_ufiles_3d_data) :: t_uf3d

   character(len=132) :: longstring
   character(len=5)   :: treename

   integer :: idx, shot, run
   integer :: i,j,k,ii
   integer :: N_Psi, N_Phi, N_Rho_Tor_Norm, N_Rho_Tor, N_Ions, N_q, N_r, N_z
   integer :: N_X_Points,N_S_Points, N_Outline_Points
   integer :: N_Profiles_2D, N_Dim1, N_Dim2
   integer :: N_C_Dim1, N_C_Dim2
   integer :: N_PF_Coils

   !prepare to write transp ufile
   integer :: ilun
   CHARACTER(len=16) :: prefix, suffix, disk, directory
   integer :: ishot
   CHARACTER(len=10) :: shdate
   CHARACTER(len=4) :: tdev !(tokamak id)
   integer :: ndim
   CHARACTER(len=160) :: comment
   integer :: tlen, xlen, ylen, zlen, tmplen
   integer :: is,jj,kk,ll,Ntmp,ierr
   real :: tmp1, tmp2, tmp3, eps = 1.e-6
   real :: r0 !vacuum major radius
   REAL :: twopi
   integer :: iarg,nargs
   character*3 :: args(2),cshot,crun
   integer :: nsctime, nprtime

   twopi=2.*4.*atan(1.)

   treename = 'ids'
   !shot = 20
   !run = 107
   cshot = ' '
   crun = ' '
   call get_arg_count(nargs)
   if (nargs.lt.2) then
      write(*,*) 'imas2transp should be given shot and run number'
      write(*,*) 'syntax: ./imas2transp shotnumber runnumber'
      write(*,*) 'shotnumber and runnumber are <3 digits integers'
   endif
   iarg=1
   call get_arg(iarg,args(iarg))
   cshot=args(iarg)
   iarg=2
   call get_arg(iarg,args(iarg))
   crun=args(iarg)
   write(*,*) 'imas2transp runid=',trim(cshot), trim(crun)
   !convert string into integer
   read(cshot,*,iostat=ierr) shot
   read(crun,*,iostat=ierr) run
   write(*,'(a,2i3,a)') 'Open shot', shot, run, 'in MDS'

   call imas_open(treename,shot,run,idx)
   call ids_get(idx,"core_profiles", cp)
   call ids_get(idx,"equilibrium", eq)
   call ids_get(idx,"nbi", nbi)
   call ids_get(idx,"pf_active", pf)

   !example test
   !ii       = 2
   !N_Rho_Tor_Norm   = size(cp%profiles_1d(ii)%grid%rho_tor_norm)
   !N_Psi            = size(eq%time_slice(ii)%profiles_1d%psi)
   !N_Ions           = size(cp%profiles_1d(ii)%ion)
   !N_Outline_Points = size(eq%time_slice(ii)%boundary%lcfs%r)
   !N_X_Points       = size(eq%time_slice(ii)%boundary%x_point)
   !N_S_Points       = size(eq%time_slice(ii)%boundary%strike_point)
   !N_Profiles_2D    = size(eq%time_slice(ii)%profiles_2d)

   !set up transp ufile writing
   ilun=33
   disk=''
   directory='./'
   ! These subroutines are hardcoded for ITER. Commenting out. JAC 12/6/16
   ! call write_shotnumber(shot,run,ishot)
   ! tdev='ITER'
   ! call write_shdate(shdate)
   ! call write_omment(comment)
   ! write(*,*) "shot number=",ishot,"shot date=",shdate, "comment=",trim(comment)
   write(*,*)
   write(*,*)
   write(*,*) "========START UFILE WRITING======================"
   write(*,*)
   write(*,*)

   !good 1. A38601.CUR:plasma_current_Amperes(time_seconds)
   ! plasma current:
   ! *CUR(TIM)       { 'plasma current',     'amps' }
   nsctime = size(eq%time)
   tmplen=size(eq%time_slice)

   ndim=1
   prefix='A'
   suffix='CUR'
   t_uf1d%nsc=0
   write(comment,'(a)') 'imported from equilibrium IDS time_slice global_quantities ip'

   t_uf1d%nx=nsctime
   if (t_uf1d%nx > 0) then
      t_uf1d%labelx='time:'
      t_uf1d%unitsx='second'
      t_uf1d%labelf='plasma current:'
      t_uf1d%unitsf='amps'
      t_uf1d%iproc=0
      allocate(t_uf1d%X(t_uf1d%nx))
      allocate(t_uf1d%F(t_uf1d%nx))
      t_uf1d%X=eq%time
      t_uf1d%F(:)=eq%time_slice(:)%global_quantities%ip
      call put_data_to_ufiles(ilun, &
         prefix,suffix,disk,directory, &
         ishot, &
         tdev,ndim,shdate, &
         comment, &
         t_uf0d,t_uf1d,t_uf2d,t_uf3d, &
         ierr)
      deallocate(t_uf1d%X, t_uf1d%F)
   endif
   write(*,*) '  --->',trim(prefix),'.',trim(suffix),' ufile write done.'
   write(*,*)

!not needed 2. A38601.DF4:diffusivity cm**2/sec (time seconds, (sqrt(Phi/Phi(a))))
   !*DF4(TIM,RAD)   { 'He4++ Diffusivity',  'cm**2/sec' }
   !core_transport
   !Effective diffusivity {dynamic} [m^2.s^-1], it=1,51, iphi=1,101

!todo: sign and name 3. A38601.PFC : Free boundary poloidal field coil currents
   !  *PFC(TIM,-) { 'PFC currents', 'A' }
   ! pfc Ampere(time second, ccindex) it=1,358, ind=1,12
   ! pf_active coil(1...N)/current
   !           coil(1...N)/current/data (coil(:)/current/time)
   !           coil(1...N)/current/time
!30. minus sign M38601.PFC
!31. plus  sign P38601.PFC
   !???ask SunHeeKim the sign and name
#if 0
   ylen=size(pf%coil)
   write(*,*) "size pf coil", ylen
   if (ylen <=0) then
      write(*,*) "WARNING : no coil current"
      !go to 101
   else
      do ll=1,ylen
         write(*,*) "  pf coil name ", pf%coil(ll)%name
      enddo
   endif
   !2 additional coils are for internal coils.
   !You have to handle these coils just like other CS/PF coils.
   !These coils are not virtual, they are almost real.
   !Therefore, if possible please use 14 coils instead of 12 coils,
   !with assigning zero currents on the last 2 coils as done in the IDS data .
   !However if you have a difficulty for doing this,
   !you can ignore the last 2 coils as these coils would not affect the scenario simulations
   !which would be main interests for TRANSP.
   if (ylen .gt. 12) ylen=12
   !cjdo kk=2,ylen
   !cj   if (tlen .ne. size(pf%coil(kk)%current%time)) then
   !cj      write(*,*) "WARNING : coil current len not match"
   !cj      go to 101
   !cj   endif
   !cjenddo
   tlen=size(pf%time)    !!coil(1)%current%time) not filled
   write(*,*)"time points=",tlen,"pf coil index",ylen

   ndim=2
   prefix='A'
   suffix='PFC'
   t_uf2d%nsc=0
   write(comment,'(a,i3)') 'imported from ids pf coil current, total',ylen

   t_uf2d%nx=tlen  !time
   t_uf2d%ny=ylen  !ccindex
   t_uf2d%nf1=t_uf2d%nx
   if (t_uf2d%nx > 0 .and. t_uf2d%ny > 0) then
      t_uf2d%labelx='time:'
      t_uf2d%unitsx='second'
      t_uf2d%labely='ccindex:'
      t_uf2d%unitsy='-'
      t_uf2d%labelf='PFC currents:'
      t_uf2d%unitsf='A'
      t_uf2d%iproc=0
      allocate(t_uf2d%X(t_uf2d%nx))
      allocate(t_uf2d%Y(t_uf2d%ny))
      allocate(t_uf2d%F(t_uf2d%nx,t_uf2d%ny))
      t_uf2d%X=pf%time   !!coil(1)%current%time not filled
      do kk=1,t_uf2d%ny
         t_uf2d%Y(kk)=kk
         do jj=1,t_uf2d%nx
            t_uf2d%F(jj,kk)=pf%coil(kk)%current%data(jj)
         enddo
      enddo
      call put_data_to_ufiles(ilun, &
         prefix,suffix,disk,directory, &
         ishot, &
         tdev,ndim,shdate, &
         comment, &
         t_uf0d,t_uf1d,t_uf2d,t_uf3d, &
         ierr)
      deallocate(t_uf2d%X, t_uf2d%Y, t_uf2d%F)
   endif
   write(*,*) '  --->',trim(prefix),'.',trim(suffix),' ufile write done.'
   write(*,*)
#endif
!good 4. A38601.RBZ R*(external toroidal field)
   ! B*R T*m (time seconds) it= 358
   ! core_profiles : cp vacuum_toroidal_field
   !    vacuum_toroidal_field/r0  meter
   !    vacuum_toroidal_field/b0 (time)  Tesla
   !    cp b0 t1 tn=  0.  0.
   ! equilibrium : eq vacuum_toroidal_field
   !    vacuum_toroidal_field/r0  meter
   !    vacuum_toroidal_field/b0 (time)  Tesla

   tlen = size(eq%time)
   tmplen=size(eq%vacuum_toroidal_field%b0)
   if (tlen .ne. size(cp%time)) then
      write(*,*) "WARNING : time points not match"
      go to 101
   endif
   write(*,*) "cp time t1=", cp%time(1), "cp time tn=", cp%time(tlen)
   !write(*,*) "cp b0 t1 tn=", cp%vacuum_toroidal_field%b0(1), cp%vacuum_toroidal_field%b0(tlen)
   write(*,*) "eq time t1=", eq%time(1), "eq time tn=", eq%time(tlen), "size of b0=", tmplen
   write(*,*) "eq r0=",eq%vacuum_toroidal_field%r0
   write(*,*) "eq b0 t1 tn=", eq%vacuum_toroidal_field%b0(1), eq%vacuum_toroidal_field%b0(tlen)

   !use eq%vacuum_toroidal_field%b0 since cp is zero
   ndim=1
   prefix='A'
   suffix='RBZ'
   t_uf1d%nsc=0
   write(comment,'(a)') 'imported from equilibrium IDS vacuum_toroidal_field r0 * b0'

   t_uf1d%nx=size(eq%time)
   if (t_uf1d%nx > 0) then
      t_uf1d%labelx='time:'
      t_uf1d%unitsx='seconds'
      t_uf1d%labelf='B*R:'
      t_uf1d%unitsf='T*cm'
      t_uf1d%iproc=0
      allocate(t_uf1d%X(t_uf1d%nx))
      allocate(t_uf1d%F(t_uf1d%nx))
      t_uf1d%X=eq%time
      !t_uf1d%F=100. * eq%vacuum_toroidal_field%r0 * eq%vacuum_toroidal_field%b0
      r0=6.20
      t_uf1d%F=100. * r0 * eq%vacuum_toroidal_field%b0
      call put_data_to_ufiles(ilun, &
         prefix,suffix,disk,directory, &
         ishot, &
         tdev,ndim,shdate, &
         comment, &
         t_uf0d,t_uf1d,t_uf2d,t_uf3d, &
         ierr)
      deallocate(t_uf1d%X, t_uf1d%F)
   endif
   write(*,*) '  --->',trim(prefix),'.',trim(suffix),' ufile write done.'
   write(*,*)

!not needed 5. A38601.TP4
   ! confinement (time seconds, (sqrt(Phi/Phi(a)))) it=51, iphi=101

!not needed 6. A38601.VC4
   ! *VC4(TIM,RAD)   { 'He4++ v(convective)',        'cm/sec' }
   !  velocity cm/sec  (time seconds, (sqrt(Phi/Phi(a)))) it=51, iphi=101

!good 7. A38601.VSF : surface loop voltage  Volt (time seconds) it=358
   ! *VSF(TIM)       { 'surface voltage',    'volts' }
   ! core_profiles global_quantities/v_loop  (time FLT_1D) LCFS loop voltage {dynamic} [V]
   tlen = size(cp%time)
   tmplen=size(cp%global_quantities%v_loop)
   write(*,*) "cp v_loop t1 tn=", cp%global_quantities%v_loop(1), cp%global_quantities%v_loop(tmplen)

   ndim=1
   prefix='A'
   suffix='VSF'
   t_uf1d%nsc=0
   write(comment,'(a)') 'imported from ids cp global_quantities v_loop'

   t_uf1d%nx=tlen
   if (t_uf1d%nx > 0) then
      t_uf1d%labelx='time:'
      t_uf1d%unitsx='seconds'
      t_uf1d%labelf='surface voltage:'
      t_uf1d%unitsf='volts'
      t_uf1d%iproc=0
      allocate(t_uf1d%X(t_uf1d%nx))
      allocate(t_uf1d%F(t_uf1d%nx))
      t_uf1d%X=cp%time
      t_uf1d%F=cp%global_quantities%v_loop
      call put_data_to_ufiles(ilun, &
         prefix,suffix,disk,directory, &
         ishot, &
         tdev,ndim,shdate, &
         comment, &
         t_uf0d,t_uf1d,t_uf2d,t_uf3d, &
         ierr)
      deallocate(t_uf1d%X, t_uf1d%F)
   endif
   write(*,*) '  --->',trim(prefix),'.',trim(suffix),' ufile write done.'
   write(*,*)

!good set NLZFI2=.T & NRIxxx=8 8. A38601.ZEF : charge into Zeff at core
   ! ! Zeff vs. time
   ! *ZEF(TIM)       { 'Zeff',               ' ' }
   !  zeff (time  seconds) it=358
!good 27. B38601.ZE2
!good 27. B38601.ZFA : change into Zeff at edge
   ! Zeff (time seconds, ((Phi/Phi(a)))
   !*ZFA(TIM)       { 'Zeff at edge',       ' ' }
!good 28. B38601.ZF2 : charges (time, rad)
   ! Zeff (time seconds, ! ((Phi/Phi(a)))
   !*ZF2(TIM,RAD)   { 'Zeff profile',       ' ' }
   ! 144
   ! 81
   ! core_profiles/profiles_1d(time)/zeff(profiles_1d(:)/grid/rho_tor_norm)
   ! Effective charge {dynamic} [-]  FLT_1D  1-

   tlen=size(cp%profiles_1d)
   ylen=size(cp%profiles_1d(1)%grid%rho_tor_norm)
   write(*,*) "cp zeff t1=",cp%profiles_1d(1)%zeff(1),cp%profiles_1d(1)%zeff(ylen)
   write(*,*) "cp zeff tn=",cp%profiles_1d(tlen)%zeff(1),cp%profiles_1d(tlen)%zeff(ylen)

   !write .ZEF
   ndim=1
   prefix='A'
   suffix='ZEF'
   t_uf1d%nsc=0
   write(comment,'(a)') 'imported from ids cp profiles_1d(time) zeff(1:1)'

   t_uf1d%nx=tlen
   if (t_uf1d%nx > 0) then
      t_uf1d%labelx='time:'
      t_uf1d%unitsx='seconds'
      t_uf1d%labelf='Zeff at core:'
      t_uf1d%unitsf='-'
      t_uf1d%iproc=0
      allocate(t_uf1d%X(t_uf1d%nx))
      allocate(t_uf1d%F(t_uf1d%nx))
      t_uf1d%X=cp%time
      do kk=1,tlen
         t_uf1d%F(kk)=cp%profiles_1d(kk)%zeff(1)
      enddo
      call put_data_to_ufiles(ilun, &
         prefix,suffix,disk,directory, &
         ishot, &
         tdev,ndim,shdate, &
         comment, &
         t_uf0d,t_uf1d,t_uf2d,t_uf3d, &
         ierr)
      deallocate(t_uf1d%X, t_uf1d%F)
   endif
   write(*,*) '  --->',trim(prefix),'.',trim(suffix),' ufile write done.'
   write(*,*)

   !write .ZFA
   ndim=1
   prefix='A'
   suffix='ZFA'
   t_uf1d%nsc=0
   write(comment,'(a)') 'imported from ids cp profiles_1d(time) zeff(tlen:tlen)'

   t_uf1d%nx=tlen
   if (t_uf1d%nx > 0) then
      t_uf1d%labelx='time:'
      t_uf1d%unitsx='seconds'
      t_uf1d%labelf='Zeff at edge:'
      t_uf1d%unitsf='-'
      t_uf1d%iproc=0
      allocate(t_uf1d%X(t_uf1d%nx))
      allocate(t_uf1d%F(t_uf1d%nx))
      t_uf1d%X=cp%time
      do kk=1,tlen
         t_uf1d%F(kk)=cp%profiles_1d(kk)%zeff(ylen)
      enddo
      call put_data_to_ufiles(ilun, &
         prefix,suffix,disk,directory, &
         ishot, &
         tdev,ndim,shdate, &
         comment, &
         t_uf0d,t_uf1d,t_uf2d,t_uf3d, &
         ierr)
      deallocate(t_uf1d%X, t_uf1d%F)
   endif
   write(*,*) '  --->',trim(prefix),'.',trim(suffix),' ufile write done.'
   write(*,*)

   !write .ZF2
   ndim=2
   prefix='A'
   suffix='ZF2'
   t_uf2d%nsc=0
   write(comment,'(a)') 'imported from ids cp profiles_1d(time) zeff(1:tlen)'

   t_uf2d%nx=tlen  !time
   t_uf2d%ny=ylen  !normalized toroidal flux
   t_uf2d%nf1=t_uf2d%nx
   if (t_uf2d%nx > 0 .and. t_uf2d%ny > 0) then
      t_uf2d%labelx='time:'
      t_uf2d%unitsx='seconds'
      t_uf2d%labely='normalized toroidal field:'
      t_uf2d%unitsy='-'
      t_uf2d%labelf='Zeff profile:'
      t_uf2d%unitsf='-'
      t_uf2d%iproc=0
      allocate(t_uf2d%X(t_uf2d%nx))
      allocate(t_uf2d%Y(t_uf2d%ny))
      allocate(t_uf2d%F(t_uf2d%nx,t_uf2d%ny))
      t_uf2d%X=cp%time
      t_uf2d%Y=cp%profiles_1d(1)%grid%rho_tor_norm
      do kk=1, t_uf2d%ny
         do jj=1, t_uf2d%nx
            t_uf2d%F(jj,kk)=cp%profiles_1d(jj)%zeff(kk)
         enddo
      enddo
      call put_data_to_ufiles(ilun, &
         prefix,suffix,disk,directory, &
         ishot, &
         tdev,ndim,shdate, &
         comment, &
         t_uf0d,t_uf1d,t_uf2d,t_uf3d, &
         ierr)
      deallocate(t_uf2d%X, t_uf2d%Y, t_uf2d%F)
   endif
   write(*,*) '  --->',trim(prefix),'.',trim(suffix),' ufile write done.'
   write(*,*)

!not needed 9. AA38601.GRB : 38601ITER 2 0 6 ;-SHOT #- F(X) DATA -UF2DWR- 10Jul2014
   ! *GRB(TIM,RAD)   { 'MHD g (R*Bt) profile', 'Tesla*cm' }
   !-1.0000E+35   ;-SCALAR, LABEL FOLLOWS: TMIN:     Tmin      seconds
   ! 1.0000E+35   ;-SCALAR, LABEL FOLLOWS: TMAX:     Tmax      seconds
   ! 0.0000E+00   ;-SCALAR, LABEL FOLLOWS: TSMOOTH:  delta(t)  seconds
   ! 1.0000E+00   ;-SCALAR, LABEL FOLLOWS: FBDY:     init.bdy
   ! 8.0000E-02   ;-SCALAR, LABEL FOLLOWS: CRAT:     min.curv
   ! 3.2000E-01   ;-SCALAR, LABEL FOLLOWS: GS_ERRMAX:allow.err.
   ! 5.0000E-01   ;-SCALAR, LABEL FOLLOWS: JAC_VARMX:allow.var.
   !-1.0000E+00   ;-SCALAR, LABEL FOLLOWS: BTOR_CCW: Dir. of B
   ! 1.0000E+00   ;-SCALAR, LABEL FOLLOWS: ITOR_CCW: Dir. of I
   ! x=(sqrt(phi/philim))          ;-INDEPENDENT VARIABLE LABEL: X-
   ! Time                Seconds   ;-INDEPENDENT VARIABLE LABEL: Y-
   ! EFIT g (R*Bt)       Tesla*m   ;-DEPENDENT VARIABLE LABEL-
   ! 0                             ;-PROC CODE- 0:RAW 1:AVG 2:SM 3:AVG+SM
   !  41                    ;-# OF X PTS-
   ! 144                    ;-# OF Y PTS- X,Y,F(X,Y) DATA FOLLOW:
   ! eq%time_slice(:)%profiles_1d%f

!todo add the scalars, all else good 10. AA38601.LIM
   ! AA38601.LIM:  38601ITER 1 0 6  ;-SHOT #- F(X) DATA -UF1DWR- 10Jul2014
   ! ! /LIM
   ! ! Numerical limiter
   ! *LIM(TIM,-)     { 'Axisymmetric Limiter',  'cm' }
   !-1.0000E+35  ;-SCALAR, LABEL FOLLOWS: TMIN:     Tmin      seconds
   ! 1.0000E+35  ;-SCALAR, LABEL FOLLOWS: TMAX:     Tmax      seconds
   ! 0.0000E+00  ;-SCALAR, LABEL FOLLOWS: TSMOOTH:  delta(t)  seconds
   ! 1.0000E+00  ;-SCALAR, LABEL FOLLOWS: FBDY:     init.bdy
   ! 8.0000E-02  ;-SCALAR, LABEL FOLLOWS: CRAT:     min.curv
   ! 3.2000E-01  ;-SCALAR, LABEL FOLLOWS: GS_ERRMAX:allow.err.
   ! 5.0000E-01  ;-SCALAR, LABEL FOLLOWS: JAC_VARMX:allow.var.
   !-1.0000E+00  ;-SCALAR, LABEL FOLLOWS: BTOR_CCW: Dir. of B
   ! 1.0000E+00  ;-SCALAR, LABEL FOLLOWS: ITOR_CCW: Dir. of I
   ! R of limiter contourm         ;-INDEPENDENT VARIABLE LABEL-
   ! Z of limiter contourm         ;-DEPENDENT VARIABLE LABEL-
   ! 57                    ;-# OF PTS-  X, F(X) DATA FOLLOW:
   ! equilibrium : eq time_slice(:)/boundary/
   ! type ! 0 (limiter) or 1 (diverted) {dynamic} ! INT_0D
   ! active_limiter_point RZ position of the active limiter point (point of the plasma boundary in contact with the limiter)  structure
   ! active_limiter_point/r Major radius {dynamic} [m]  FLT_0D
   ! active_limiter_point/z Height {dynamic} [m]  FLT_0D
   ! ???? check to see if it changes with time
#if 0
   tlen=size(eq%time)

   ndim=1
   prefix='A'
   suffix='LIM'
   t_uf1d%nsc=2
   write(comment,'(a)') 'imported from equilibrium IDS time_slice(time) boundary active_limiter_point r z'
   t_uf1d%nsc=2

   if (t_uf1d%nsc>0) then
      allocate(t_uf1d%SCVAL(t_uf1d%nsc))
      allocate(t_uf1d%labels(t_uf1d%nsc))
      allocate(t_uf1d%unitss(t_uf1d%nsc))
      t_uf1d%SCVAL(1) =  1.
      t_uf1d%SCVAL(2) =  1.
      write(t_uf1d%labels(1),'(a)') 'BTOR_CCW: Dir. of B'
      write(t_uf1d%labels(2),'(a)') 'ITOR_CCW: Dir. of I'
      write(t_uf1d%unitss(1),'(a)') ''
      write(t_uf1d%unitss(2),'(a)') ''
   endif

   t_uf1d%nx=tlen
   if (t_uf1d%nx > 0) then
      t_uf1d%labelx='R of limiter contourm:'
      t_uf1d%unitsx='cm'
      t_uf1d%labelf='Z of limiter contourm:'
      t_uf1d%unitsf='cm'
      t_uf1d%iproc=0
      allocate(t_uf1d%X(t_uf1d%nx))
      allocate(t_uf1d%F(t_uf1d%nx))
      do jj=1,tlen
         t_uf1d%X(jj)=100. * eq%time_slice(jj)%boundary%active_limiter_point%r
         t_uf1d%F(jj)=100. * eq%time_slice(jj)%boundary%active_limiter_point%z
      enddo
      call put_data_to_ufiles(ilun, &
         prefix,suffix,disk,directory, &
         ishot, &
         tdev,ndim,shdate, &
         comment, &
         t_uf0d,t_uf1d,t_uf2d,t_uf3d, &
         ierr)
      deallocate(t_uf1d%X, t_uf1d%F)
   endif
   if (t_uf1d%nsc>0) then
      deallocate(t_uf1d%SCVAL)
      deallocate(t_uf1d%labels)
      deallocate(t_uf1d%unitss)
   endif
   write(*,*) '  --->',trim(prefix),'.',trim(suffix),' ufile write done.'
   write(*,*)

!not needed 11. AA38601.MMX : Equilibrium "Moments" -- Fourier representation
   ! AA38601.MMX:  38601ITER 3 0 6 ;-SHOT #- F(X) DATA -UF3DWR- 10Jul2014
   ! /MMX
   ! *MMX(TIM,-)     {'Equilibrium moments', 'cm' }
   !-1.0000E+35 ;-SCALAR, LABEL FOLLOWS: TMIN:     Tmin      seconds
   ! 1.0000E+35 ;-SCALAR, LABEL FOLLOWS: TMAX:     Tmax      seconds
   ! 0.0000E+00 ;-SCALAR, LABEL FOLLOWS: TSMOOTH:  delta(t)  seconds
   ! 1.0000E+00 ;-SCALAR, LABEL FOLLOWS: FBDY:     init.bdy
   ! 8.0000E-02 ;-SCALAR, LABEL FOLLOWS: CRAT:     min.curv
   ! 3.2000E-01 ;-SCALAR, LABEL FOLLOWS: GS_ERRMAX:allow.err.
   ! 5.0000E-01 ;-SCALAR, LABEL FOLLOWS: JAC_VARMX:allow.var.
   !-1.0000E+00 ;-SCALAR, LABEL FOLLOWS: BTOR_CCW: Dir. of B
   ! 1.0000E+00 ;-SCALAR, LABEL FOLLOWS: ITOR_CCW: Dir. of I
   ! 2.0000E+00 ;-SCALAR, LABEL FOLLOWS: MMX_MAXE: exponent  maximum
   ! Time                Seconds   ;-INDEPENDENT VARIABLE LABEL: X-
   ! x=(sqrt(phi/philim))          ;-INDEPENDENT VARIABLE LABEL: Y-
   ! moment index                  ;-INDEPENDENT VARIABLE LABEL: Z-
   ! equilibrium moments m         ;-DEPENDENT VARIABLE LABEL-
   ! 144                    ;-# OF X PTS-
   !  41                    ;-# OF Y PTS-
   ! 68                    ;-# OF Z PTS- X,Y,F(X,Y)
   ! ???ids doesn't have it ---> change to writing RFS, ZFS
!good. set LEVGEO=8. Does I need to wrie the above scalars? 33. boundary RFS, ZFS
   ! *RFS(TIM,-)     { 'R(theta,xi)', 'cm' }
   ! *ZFS(TIM,-)     { 'Z(theta,xi)', 'cm' }
   ! eq%time_slice(ii)%boundary%lcfs%z(i), i=1,N_Outline_Points)
   ! ???confim with KIM
   tlen=size(eq%time)
   ylen=size(eq%time_slice(1)%boundary%lcfs%r)

   ndim=2
   prefix='A'
   suffix='RFS'
   t_uf2d%nsc=2
   write(comment,'(a)') 'imported from equilibrium IDS time_slice(time) boundary lcfs rz(1:ylen)'

   if (t_uf2d%nsc>0) then
      allocate(t_uf2d%SCVAL(t_uf2d%nsc))
      allocate(t_uf2d%labels(t_uf2d%nsc))
      allocate(t_uf2d%unitss(t_uf2d%nsc))
      t_uf2d%SCVAL(1) =  1.
      t_uf2d%SCVAL(2) =  1.
      write(t_uf2d%labels(1),'(a)') 'BTOR_CCW: Dir. of B'
      write(t_uf2d%labels(2),'(a)') 'ITOR_CCW: Dir. of I'
      write(t_uf2d%unitss(1),'(a)') ''
      write(t_uf2d%unitss(2),'(a)') ''
   endif

   t_uf2d%nx=tlen
   t_uf2d%ny=ylen
   if (t_uf2d%nx > 0 .and. t_uf2d%ny > 0) then
      t_uf2d%labelx='Time:'
      t_uf2d%unitsx='seconds'
      t_uf2d%labely='# of bdy R points:'
      t_uf2d%unitsy='-'
      t_uf2d%labelf='R:'
      t_uf2d%unitsf='cm'
      t_uf2d%iproc=0
      allocate(t_uf2d%X(t_uf2d%nx))
      allocate(t_uf2d%Y(t_uf2d%ny))
      allocate(t_uf2d%F(t_uf2d%nx,t_uf2d%ny))
      t_uf2d%X=eq%time
      do kk=1,t_uf2d%ny
         t_uf2d%Y(kk)=kk
         do jj=1,t_uf2d%nx
            t_uf2d%F(jj,kk)=100.0*eq%time_slice(jj)%boundary%lcfs%r(kk)
         enddo
      enddo

      call put_data_to_ufiles(ilun, &
         prefix,suffix,disk,directory, &
         ishot, &
         tdev,ndim,shdate, &
         comment, &
         t_uf0d,t_uf1d,t_uf2d,t_uf3d, &
         ierr)
      deallocate(t_uf2d%X, t_uf2d%Y, t_uf2d%F)
   endif
   if (t_uf2d%nsc>0) then
      deallocate(t_uf2d%SCVAL)
      deallocate(t_uf2d%labels)
      deallocate(t_uf2d%unitss)
   endif
   write(*,*) '  --->',trim(prefix),'.',trim(suffix),' ufile write done.'
   write(*,*)

   ndim=2
   prefix='A'
   suffix='ZFS'
   t_uf2d%nsc=0
   write(comment,'(a)') 'imported from equilibrium IDS time_slice(time) boundary lcfs z(1:ylen)'

   t_uf2d%nx=tlen
   t_uf2d%ny=ylen
   if (t_uf2d%nx > 0 .and. t_uf2d%ny > 0) then
      t_uf2d%labelx='Time:'
      t_uf2d%unitsx='seconds'
      t_uf2d%labely='# of bdy Z points:'
      t_uf2d%unitsy='-'
      t_uf2d%labelf='Z(theta,xi):'
      t_uf2d%unitsf='meter'
      t_uf2d%iproc=0
      allocate(t_uf2d%X(t_uf2d%nx))
      allocate(t_uf2d%Y(t_uf2d%ny))
      allocate(t_uf2d%F(t_uf2d%nx,t_uf2d%ny))
      t_uf2d%X=eq%time
      do kk=1,t_uf2d%ny
         t_uf2d%Y(kk)=kk
         do jj=1,t_uf2d%nx
            !t_uf2d%F(jj,kk)=100.0*eq%time_slice(jj)%boundary%lcfs%z(kk)
            t_uf2d%F(jj,kk)=100.0*eq%time_slice(1)%boundary%lcfs%z(kk)
         enddo
      enddo
      call put_data_to_ufiles(ilun, &
         prefix,suffix,disk,directory, &
         ishot, &
         tdev,ndim,shdate, &
         comment, &
         t_uf0d,t_uf1d,t_uf2d,t_uf3d, &
         ierr)
      deallocate(t_uf2d%X, t_uf2d%Y, t_uf2d%F)
   endif
   write(*,*) '  --->',trim(prefix),'.',trim(suffix),' ufile write done.'
   write(*,*)
#endif
!good. write scalars? 12. AA38601.PF0 : 38601ITER 1 0 6 ;-SHOT #- F(X) DATA -UF1DWR- 10Jul2014
   ! Poloidal flux at magnetic axis -- e.g. as supplied by EFIT
   !   Internally in TRANSP Psi(0)=0 is presumed; this is the offset to
   !   get to Psi(0) according to a free boundary equilibrium calculation or
   !   analysis:  Psi(free_bdy) = Psi(TRANSP) + PF0
   !   CAUTION: poloidal flux in TRANSP is generally Wb/rad = Wb/(2pi)
   !   so that Bpol = grad(Psi)/R.  Input data is assumed to be "really"
   !   in these units, even if units are labeled as "Wb" or "Webers".  This
   !   is necessary for backwards compatibility reasons.
   ! *PF0(TIM)       { 'Psi(axis)', 'Wb/rad' }
   !-1.0E+35   ;-SCALAR, LABEL FOLLOWS: TMIN:     Tmin      seconds
   ! 1.0E+35   ;-SCALAR, LABEL FOLLOWS: TMAX:     Tmax      seconds
   ! 0.0E+00   ;-SCALAR, LABEL FOLLOWS: TSMOOTH:  delta(t)  seconds
   ! 1.0E+00   ;-SCALAR, LABEL FOLLOWS: FBDY:     init.bdy
   ! 8.0E-02   ;-SCALAR, LABEL FOLLOWS: CRAT:     min.curv
   ! 3.2E-01   ;-SCALAR, LABEL FOLLOWS: GS_ERRMAX:allow.err.
   ! 5.0E-01   ;-SCALAR, LABEL FOLLOWS: JAC_VARMX:allow.var.
   !-1.0E+00   ;-SCALAR, LABEL FOLLOWS: BTOR_CCW: Dir. of B
   ! 1.0E+00   ;-SCALAR, LABEL FOLLOWS: ITOR_CCW: Dir. of I
   ! Time      Seconds   ;-INDEPENDENT VARIABLE LABEL-
   ! EFIT delta(Psi0)    Wb/rad    ;-DEPENDENT VARIABLE LABEL-
   ! 144       ;-# OF PTS-  X, F(X) DATA FOLLOW:

   ! core_profiles : cp profiles_1d(time)/grid/psi
   ! ( FLT_1D profiles_1d(time)/grid/rho_tor_norm) Poloidal magnetic flux {dynamic} [Wb]

   ! equilibrium : eq time_slice(time)/global_quantities/
   ! psi_axis      Poloidal flux at the magnetic axis {dynamic} [Wb] FLT_0D
   ! psi_boundary  Poloidal flux at the selected plasma boundary {dynamic} [Wb]  FLT_0D
   ! psi           Poloidal flux {dynamic} [Wb]  FLT_1D  1- 1...N(=81)
   !time_slice(:)/profiles_2d(:)/psi  Values of the poloidal flux at the grid in the poloidal plane {dynamic} [Wb]  FLT_2D
   !1- time_slice(:)/profiles_2d(:)/grid/dim1 (=65)
   !2- time_slice(:)/profiles_2d(:)/grid/dim2 (=129)

   tlen=size(eq%time)

   ndim=1
   prefix='A'
   suffix='PF0'
   t_uf1d%nsc=2
   write(comment,'(a)') 'imported from equilibrium IDS time_slice(time) global_quantities psi_axis'

   if (t_uf1d%nsc>0) then
      allocate(t_uf1d%SCVAL(t_uf1d%nsc))
      allocate(t_uf1d%labels(t_uf1d%nsc))
      allocate(t_uf1d%unitss(t_uf1d%nsc))
      t_uf1d%SCVAL(1) =  1.
      t_uf1d%SCVAL(2) =  1.
      write(t_uf1d%labels(1),'(a)') 'BTOR_CCW: Dir. of B'
      write(t_uf1d%labels(2),'(a)') 'ITOR_CCW: Dir. of I'
      write(t_uf1d%unitss(1),'(a)') ''
      write(t_uf1d%unitss(2),'(a)') ''
   endif

   t_uf1d%nx=tlen
   if (t_uf1d%nx > 0) then
      t_uf1d%labelx='TIME:'
      t_uf1d%unitsx='seconds'
      t_uf1d%labelf='Psi(axis):'
      t_uf1d%unitsf='Wb/rad'
      t_uf1d%iproc=0
      allocate(t_uf1d%X(t_uf1d%nx))
      allocate(t_uf1d%F(t_uf1d%nx))
      t_uf1d%X=eq%time
      t_uf1d%F(:)=eq%time_slice(:)%global_quantities%psi_axis / twopi
      call put_data_to_ufiles(ilun, &
         prefix,suffix,disk,directory, &
         ishot, &
         tdev,ndim,shdate, &
         comment, &
         t_uf0d,t_uf1d,t_uf2d,t_uf3d, &
         ierr)
      deallocate(t_uf1d%X, t_uf1d%F)
   endif
   if (t_uf1d%nsc>0) then
      deallocate(t_uf1d%SCVAL)
      deallocate(t_uf1d%labels)
      deallocate(t_uf1d%unitss)
   endif
   write(*,*) '  --->',trim(prefix),'.',trim(suffix),' ufile write done.'
   write(*,*)

!good. write scalars? 13. AA38601.PLF :  38601ITER 1 0 6 ;-SHOT #- F(X) DATA -UF1DWR- 10Jul2014
   ! Enclosed poloidal flux from eq. code: meaning the flux difference from
   ! the plasma center (mag. axis) to the plasma boundary
   ! (as used by TRANSP, near the last closed flux surface).
   ! CAUTION: poloidal flux in TRANSP is generally Wb/rad = Wb/(2pi)
   ! so that Bpol = grad(Psi)/R.  Input data is assumed to be "really"
   ! in these units, even if units are labeled as "Wb" or "Webers".  This
   ! is necessary for backwards compatibility reasons.
   ! *PLF(TIM)       { 'Poloidal flux', 'Wb/rad' }
   !-1.0000E+35   ;-SCALAR, LABEL FOLLOWS: TMIN:     Tmin      seconds
   ! 1.0000E+35   ;-SCALAR, LABEL FOLLOWS: TMAX:     Tmax      seconds
   ! 0.0000E+00   ;-SCALAR, LABEL FOLLOWS: TSMOOTH:  delta(t)  seconds
   ! 1.0000E+00   ;-SCALAR, LABEL FOLLOWS: FBDY:     init.bdy
   ! 8.0000E-02   ;-SCALAR, LABEL FOLLOWS: CRAT:     min.curv
   ! 3.2000E-01   ;-SCALAR, LABEL FOLLOWS: GS_ERRMAX:allow.err.
   ! 5.0000E-01   ;-SCALAR, LABEL FOLLOWS: JAC_VARMX:allow.var.
   !-1.0000E+00   ;-SCALAR, LABEL FOLLOWS: BTOR_CCW: Dir. of B
   ! 1.0000E+00   ;-SCALAR, LABEL FOLLOWS: ITOR_CCW: Dir. of I
   ! Time                Seconds   ;-INDEPENDENT VARIABLE LABEL-
   ! EFIT psi(poloidal)  Wb/rad    ;-DEPENDENT VARIABLE LABEL-
   ! 144                    ;-# OF PTS-  X, F(X) DATA FOLLOW:
   tlen=size(eq%time)

   ndim=1
   prefix='A'
   suffix='PLF'
   t_uf1d%nsc=2
   write(comment,'(a)') 'imported from equilibrium IDS time_slice(time) global_quantities psi_boundary-psi_axis'

   if (t_uf1d%nsc>0) then
      allocate(t_uf1d%SCVAL(t_uf1d%nsc))
      allocate(t_uf1d%labels(t_uf1d%nsc))
      allocate(t_uf1d%unitss(t_uf1d%nsc))
      t_uf1d%SCVAL(1) =  1.
      t_uf1d%SCVAL(2) =  1.
      write(t_uf1d%labels(1),'(a)') 'BTOR_CCW: Dir. of B'
      write(t_uf1d%labels(2),'(a)') 'ITOR_CCW: Dir. of I'
      write(t_uf1d%unitss(1),'(a)') ''
      write(t_uf1d%unitss(2),'(a)') ''
   endif

   t_uf1d%nx=tlen
   if (t_uf1d%nx > 0) then
      t_uf1d%labelx='TIME:'
      t_uf1d%unitsx='seconds'
      t_uf1d%labelf='Poloidal flux:'
      t_uf1d%unitsf='Wb/rad'
      t_uf1d%iproc=0
      allocate(t_uf1d%X(t_uf1d%nx))
      allocate(t_uf1d%F(t_uf1d%nx))
      t_uf1d%X=eq%time
      t_uf1d%F(:)= (eq%time_slice(:)%global_quantities%psi_boundary - &
         eq%time_slice(:)%global_quantities%psi_axis) / twopi
      call put_data_to_ufiles(ilun, &
         prefix,suffix,disk,directory, &
         ishot, &
         tdev,ndim,shdate, &
         comment, &
         t_uf0d,t_uf1d,t_uf2d,t_uf3d, &
         ierr)
      deallocate(t_uf1d%X, t_uf1d%F)
   endif
   if (t_uf1d%nsc>0) then
      deallocate(t_uf1d%SCVAL)
      deallocate(t_uf1d%labels)
      deallocate(t_uf1d%unitss)
   endif
   write(*,*) '  --->',trim(prefix),'.',trim(suffix),' ufile write done.'
   write(*,*)

!good. write scalars? 13.1 AA38601.PSI :  38601ITER 1 0 6 ;-SHOT #- F(X) DATA -UF1DWR- 10Jul2014
   ! poloidal flux on the (R,Z) grid.
   ! CAUTION: poloidal flux in TRANSP is generally Wb/rad = Wb/(2pi)
   ! so that Bpol = grad(Psi)/R.  Input data is assumed to be "really"
   ! in these units, even if units are labeled as "Wb" or "Webers".  This
   ! is necessary for backwards compatibility reasons.
   ! *PLF(TIM)       { 'Poloidal flux', 'Wb/rad' }
   !-1.0000E+35   ;-SCALAR, LABEL FOLLOWS: TMIN:     Tmin      seconds
   ! 1.0000E+35   ;-SCALAR, LABEL FOLLOWS: TMAX:     Tmax      seconds
   ! 0.0000E+00   ;-SCALAR, LABEL FOLLOWS: TSMOOTH:  delta(t)  seconds
   ! 1.0000E+00   ;-SCALAR, LABEL FOLLOWS: FBDY:     init.bdy
   ! 8.0000E-02   ;-SCALAR, LABEL FOLLOWS: CRAT:     min.curv
   ! 3.2000E-01   ;-SCALAR, LABEL FOLLOWS: GS_ERRMAX:allow.err.
   ! 5.0000E-01   ;-SCALAR, LABEL FOLLOWS: JAC_VARMX:allow.var.
   !-1.0000E+00   ;-SCALAR, LABEL FOLLOWS: BTOR_CCW: Dir. of B
   ! 1.0000E+00   ;-SCALAR, LABEL FOLLOWS: ITOR_CCW: Dir. of I
   ! Time                Seconds   ;-INDEPENDENT VARIABLE LABEL-
   ! EFIT psi(poloidal)  Wb/rad    ;-DEPENDENT VARIABLE LABEL-
   ! 144                    ;-# OF PTS-  X, F(X) DATA FOLLOW:

   ndim=3
   prefix='A'
   suffix='PSI'
   t_uf3d%nsc=2
   if (t_uf3d%nsc>0) then
      allocate(t_uf3d%SCVAL(t_uf3d%nsc))
      write(*,*) "number of scalars to be written", t_uf3d%nsc
      allocate(t_uf3d%labels(t_uf3d%nsc))
      allocate(t_uf3d%unitss(t_uf3d%nsc))
      is=1
      t_uf3d%SCVAL(is)=1.
      write(t_uf3d%labels(is),'(a)') 'BTOR_CCW: Dir. of B'
      t_uf3d%unitss(is)='-'
      is=2
      t_uf3d%SCVAL(is)=1.
      write(t_uf3d%labels(is),'(a)') 'ITOR_CCW: Dir. of I'
      t_uf3d%unitss(is)='-'
   else
      write(*,*) "no scalars to be written"
   endif

   t_uf3d%nx=size(eq%time)
   t_uf3d%ny=size(eq%time_slice(1)%profiles_2d(1)%grid%dim1)
   t_uf3d%nz=size(eq%time_slice(1)%profiles_2d(1)%grid%dim2)
   t_uf3d%nf1=t_uf3d%nx
   t_uf3d%nf2=t_uf3d%ny
   write(*,*) "3d nx=", t_uf3d%nx, 'ny=', t_uf3d%ny, 'nz=', t_uf3d%nz
   if (t_uf3d%nx > 0 .and. t_uf3d%ny > 0 .and. t_uf3d%nz > 0) then
      t_uf3d%labelx='timex:'
      t_uf3d%unitsx='second'
      t_uf3d%labely='grid_dim1:'
      t_uf3d%unitsy='meter'
      t_uf3d%labelz='grid_dim2:'
      t_uf3d%unitsz='meter'
      t_uf3d%labelf='jinchen'
      t_uf3d%unitsf='f'
      t_uf3d%iproc=0
      allocate(t_uf3d%X(t_uf3d%nx))
      allocate(t_uf3d%Y(t_uf3d%ny))
      allocate(t_uf3d%Z(t_uf3d%nz))
      allocate(t_uf3d%F(t_uf3d%nx,t_uf3d%ny,t_uf3d%nz))
      t_uf3d%X=eq%time
      t_uf3d%Y=eq%time_slice(1)%profiles_2d(1)%grid%dim1
      t_uf3d%Z=eq%time_slice(1)%profiles_2d(1)%grid%dim2
      do ll=1,t_uf3d%nz
         do kk=1,t_uf3d%ny
            do jj=1,t_uf3d%nx
               t_uf3d%F(jj,kk,ll)=eq%time_slice(jj)%profiles_2d(1)%psi(kk,ll)
!            if (jj.eq.1) write(*,'(a,i3.3,x,i3.3,6(x,e10.3))') &
!            'r z psi_rz psi_rz/2pi psi_axis psi_bdy: ', kk,ll,&
!            eq%time_slice(jj)%profiles_2d(1)%grid%dim1(kk), &
!            eq%time_slice(jj)%profiles_2d(1)%grid%dim2(ll), &
!            eq%time_slice(jj)%profiles_2d(1)%psi(kk,ll), &
!            eq%time_slice(jj)%profiles_2d(1)%psi(kk,ll)/twopi, &
!            eq%time_slice(jj)%global_quantities%psi_axis , &
!            eq%time_slice(jj)%global_quantities%psi_boundary
            enddo
         enddo
      enddo
      call put_data_to_ufiles(ilun, &
         prefix,suffix,disk,directory, &
         ishot, &
         tdev,ndim,shdate, &
         comment, &
         t_uf0d,t_uf1d,t_uf2d,t_uf3d, &
         ierr)
      deallocate(t_uf3d%X, t_uf3d%Y, t_uf3d%Z, t_uf3d%F)
   endif
   if (t_uf3d%nsc>0) then
      deallocate(t_uf3d%SCVAL)
      deallocate(t_uf3d%labels)
      deallocate(t_uf3d%unitss)
   endif


!good. write scalars? set NRIxxx=6. 14. AA38601.PRS 38601ITER 2 0 6 ;-SHOT #- F(X) DATA -UF2DWR- 10Jul2014
   ! *PRS(TIM,RAD)   { 'MHD Pressure profile', 'Pascals' }
   !-1.0000E+35  ;-SCALAR, LABEL FOLLOWS: TMIN:     Tmin      seconds
   ! 1.0000E+35  ;-SCALAR, LABEL FOLLOWS: TMAX:     Tmax      seconds
   ! 0.0000E+00  ;-SCALAR, LABEL FOLLOWS: TSMOOTH:  delta(t)  seconds
   ! 1.0000E+00  ;-SCALAR, LABEL FOLLOWS: FBDY:     init.bdy
   ! 8.0000E-02  ;-SCALAR, LABEL FOLLOWS: CRAT:     min.curv
   ! 3.2000E-01  ;-SCALAR, LABEL FOLLOWS: GS_ERRMAX:allow.err.
   ! 5.0000E-01  ;-SCALAR, LABEL FOLLOWS: JAC_VARMX:allow.var.
   !-1.0000E+00  ;-SCALAR, LABEL FOLLOWS: BTOR_CCW: Dir. of B
   ! 1.0000E+00  ;-SCALAR, LABEL FOLLOWS: ITOR_CCW: Dir. of I
   ! x=(sqrt(phi/philim))          ;-INDEPENDENT VARIABLE LABEL: X-
   ! Time                Seconds   ;-INDEPENDENT VARIABLE LABEL: Y-
   ! EFIT Pressure       Pascals   ;-DEPENDENT VARIABLE LABEL-
   !  41                    ;-# OF X PTS-
   ! 144                    ;-# OF Y PTS- X,Y,F(X,Y) DATA FOLLOW:
   !
   ! equilibrium  time_slice(:)/profiles_1d/
   ! psi       Poloidal flux {dynamic} [Wb]  FLT_1D  1- 1...N
   ! phi       Toroidal flux {dynamic} [Wb]  FLT_1D  1- time_slice(:)/profiles_1d/psi
   ! pressure  Pressure {dynamic} [Pa] FLT_1D  1- time_slice(:)/profiles_1d/psi
   ! ???ask which one is mhd pressure. Is it the last line above???
   tlen=size(eq%time)
   ylen=size(eq%time_slice(1)%profiles_1d%psi)

   ndim=2
   prefix='A'
   suffix='PRS'
   t_uf2d%nsc=2
   write(comment,'(a)') 'imported from equilibrium IDS time_slice(time) profiles_1d pressure'

   if (t_uf2d%nsc>0) then
      allocate(t_uf2d%SCVAL(t_uf2d%nsc))
      allocate(t_uf2d%labels(t_uf2d%nsc))
      allocate(t_uf2d%unitss(t_uf2d%nsc))
      t_uf2d%SCVAL(1) =  1.
      t_uf2d%SCVAL(2) =  1.
      write(t_uf2d%labels(1),'(a)') 'BTOR_CCW: Dir. of B'
      write(t_uf2d%labels(2),'(a)') 'ITOR_CCW: Dir. of I'
      write(t_uf2d%unitss(1),'(a)') ''
      write(t_uf2d%unitss(2),'(a)') ''
   endif

   t_uf2d%nx=tlen
   t_uf2d%ny=ylen
   if (t_uf2d%nx > 0 .and. t_uf2d%ny > 0) then
      t_uf2d%labelx='Time:'
      t_uf2d%unitsx='seconds'
      t_uf2d%labely='normalized poloidal flux:'
      t_uf2d%unitsy='Wb'
      t_uf2d%labelf='MHD Pressure profile:'
      t_uf2d%unitsf='Pascals'
      t_uf2d%iproc=0
      allocate(t_uf2d%X(t_uf2d%nx))
      allocate(t_uf2d%Y(t_uf2d%ny))
      allocate(t_uf2d%F(t_uf2d%nx,t_uf2d%ny))
      t_uf2d%X=eq%time
      do kk=1,t_uf2d%ny
         !Be careful that the correct definition of normalized poloidal flux is:
         !(psi-psi(axis))/(psi(boundary)-psi(axis))
         !and not  poloidal flux/flux_at_boundary
         !This is different from the normalized toroidal flux definition.
         !The toroidal flux is always 0 on axis, but the poloidal flux is not.
         !This is the reason why the normalization is different
         !(it is basically the same normalization,
         !but for the toroidal flux reduces to the /boundary) because of its intrinsic 0 value on axis.
         !also check how the poloidal flux is constructed, by plotting it.
         !In Europe it is usually used the convention of a poloidal flux surface that is upside down.
         !If this is the case, the normalization has to invert the values on axis and at the boundary,
         !because in the US we use a convention with the minimum of psi on axis and the maximum at the boundary.
         !Also
         !recently ITER has started to use a definition of poloidal flux that is normalized by a factor 2*pi.
         !We need to make sure that things are consistent.
         t_uf2d%Y(kk)=(eq%time_slice(1)%profiles_1d%psi(kk)   - eq%time_slice(1)%profiles_1d%psi(1)) / &
            (eq%time_slice(1)%profiles_1d%psi(ylen) - eq%time_slice(1)%profiles_1d%psi(1))
         !cj write(*,'(a,i3,x,1e17.10)') "eq time_slices profiles_1d psi ", kk,t_uf2d%Y(kk)
         do jj=1,t_uf2d%nx
            t_uf2d%F(jj,kk)=eq%time_slice(jj)%profiles_1d%pressure(kk)
         enddo
      enddo
      call put_data_to_ufiles(ilun, &
         prefix,suffix,disk,directory, &
         ishot, &
         tdev,ndim,shdate, &
         comment, &
         t_uf0d,t_uf1d,t_uf2d,t_uf3d, &
         ierr)
      deallocate(t_uf2d%X, t_uf2d%Y, t_uf2d%F)
   endif
   if (t_uf2d%nsc>0) then
      deallocate(t_uf2d%SCVAL)
      deallocate(t_uf2d%labels)
      deallocate(t_uf2d%unitss)
   endif
   write(*,*) '  --->',trim(prefix),'.',trim(suffix),' ufile write done.'
   write(*,*)

   ! core_profiles profiles_1d(:)/
   ! ion(:)/p_i Ion pressure for that ion species {dynamic} [Pa]  FLT_1D  1- profiles_1d(:)/grid/rho_tor_norm
   ! /ion(:)/charge_state(:)/p_z Pressure of the charge state considered {dynamic} [Pa]  FLT_1D  1- profiles_1d(:)/grid/rho_tor_norm
   ! /p_i_total  Total ion pressure (sum over the ion species) {dynamic} [Pa]  FLT_1D  1- profiles_1d(:)/grid/rho_tor_norm
   ! /p_i_total_fast_perpendicular Fast (non-thermal) total ion (sum over the ion species) perpendicular pressure {dynamic} [Pa] FLT_1D  1- profiles_1d(:)/grid/rho_tor_norm
!    do jj=1,size(cp%time)
!       do kk=1,size(cp%profiles_1d(jj)%ion)
!          do ll=1, size(cp%profiles_1d(jj)%grid%rho_tor_norm)
!       !Thermal pressure (electrons+ions) {dynamic} [Pa]  FLT_1D  1- profiles_1d(:)/grid/rho_tor_norm
!       !Ion pressure for that ion species {dynamic} [Pa]  FLT_1D  1- profiles_1d(:)/grid/rho_tor_norm
!       !Electron pressure {dynamic} [Pa]  FLT_1D  1- profiles_1d(:)/grid/rho_tor_norm
!       write(85,'(a,3i3,4e13.5)') 'thermal pressure', ll,kk,jj, &
!          cp%profiles_1d(jj)%pressure_thermal(ll), &
!          cp%profiles_1d(jj)%ion(kk)%p_i(ll) , &
!          cp%profiles_1d(jj)%p_e(ll) , &
!          cp%profiles_1d(jj)%pressure_thermal(ll) - &
!        (cp%profiles_1d(jj)%ion(kk)%p_i(ll) + cp%profiles_1d(jj)%p_e(ll))
!          enddo
!       enddo
!    enddo

!good. write scalars? set NRIxxx=8. 15. AA38601.QPR :  38601ITER 2 0 6 ;-SHOT #- F(X) DATA -UF2DWR- 10Jul2014
   ! B field tilt data
   ! *QPR(TIM,RAD)   { 'q profile',          ' ' }
   !-1.0000E+35;-SCALAR, LABEL FOLLOWS: TMIN:     Tmin      seconds
   ! 1.0000E+35 ;-SCALAR, LABEL FOLLOWS: TMAX:     Tmax      seconds
   ! 0.0000E+00 ;-SCALAR, LABEL FOLLOWS: TSMOOTH:  delta(t)  seconds
   ! 1.0000E+00 ;-SCALAR, LABEL FOLLOWS: FBDY:     init.bdy
   ! 8.0000E-02 ;-SCALAR, LABEL FOLLOWS: CRAT:     min.curv
   ! 3.2000E-01 ;-SCALAR, LABEL FOLLOWS: GS_ERRMAX:allow.err.
   ! 5.0000E-01 ;-SCALAR, LABEL FOLLOWS: JAC_VARMX:allow.var.
   !-1.0000E+00 ;-SCALAR, LABEL FOLLOWS: BTOR_CCW: Dir. of B
   ! 1.0000E+00 ;-SCALAR, LABEL FOLLOWS: ITOR_CCW: Dir. of I
   ! x=(sqrt(phi/philim))          ;-INDEPENDENT VARIABLE LABEL: X-
   ! Time                Seconds   ;-INDEPENDENT VARIABLE LABEL: Y-
   ! EFIT q profile                ;-DEPENDENT VARIABLE LABEL-
   !  41                    ;-# OF X PTS-
   ! 144                    ;-# OF Y PTS- X,Y,F(X,Y) DATA FOLLOW:
   ! core_profiles  profiles_1d(time)/q
   !    Safety factor {dynamic} [-] FLT_1D  1- profiles_1d(:)/grid/rho_tor_norm
   ! equilibrium  time_slice(time)/profiles_1d/q
   !    Safety factor {dynamic} [-] FLT_1D  1- time_slice(:)/profiles_1d/psi
   tlen=size(cp%time)
   ylen=size(cp%profiles_1d(1)%grid%rho_tor_norm)

   ndim=2
   prefix='A'
   suffix='QPR'
   t_uf2d%nsc=2
   write(comment,'(a)') 'imported from ids_cp profiles_1d q'

   if (t_uf2d%nsc>0) then
      allocate(t_uf2d%SCVAL(t_uf2d%nsc))
      allocate(t_uf2d%labels(t_uf2d%nsc))
      allocate(t_uf2d%unitss(t_uf2d%nsc))
      t_uf2d%SCVAL(1) =  1.
      t_uf2d%SCVAL(2) =  1.
      write(t_uf2d%labels(1),'(a)') 'BTOR_CCW: Dir. of B'
      write(t_uf2d%labels(2),'(a)') 'ITOR_CCW: Dir. of I'
      write(t_uf2d%unitss(1),'(a)') ''
      write(t_uf2d%unitss(2),'(a)') ''
   endif

   t_uf2d%nx=tlen
   t_uf2d%ny=ylen
   if (t_uf2d%nx > 0 .and. t_uf2d%ny > 0) then
      t_uf2d%labelx='Time:'
      t_uf2d%unitsx='seconds'
      t_uf2d%labely='normalized toroidal field:'
      t_uf2d%unitsy='-'
      t_uf2d%labelf='q profile:'
      t_uf2d%unitsf='-'
      t_uf2d%iproc=0
      allocate(t_uf2d%X(t_uf2d%nx))
      allocate(t_uf2d%Y(t_uf2d%ny))
      allocate(t_uf2d%F(t_uf2d%nx,t_uf2d%ny))
      t_uf2d%X=cp%time
      do kk=1,t_uf2d%ny
         t_uf2d%Y(kk)=cp%profiles_1d(1)%grid%rho_tor_norm(kk)
         do jj=1,t_uf2d%nx
            t_uf2d%F(jj,kk)=cp%profiles_1d(jj)%q(kk)
         enddo
      enddo
      call put_data_to_ufiles(ilun, &
         prefix,suffix,disk,directory, &
         ishot, &
         tdev,ndim,shdate, &
         comment, &
         t_uf0d,t_uf1d,t_uf2d,t_uf3d, &
         ierr)
      deallocate(t_uf2d%X, t_uf2d%Y, t_uf2d%F)
   endif
   if (t_uf2d%nsc>0) then
      deallocate(t_uf2d%SCVAL)
      deallocate(t_uf2d%labels)
      deallocate(t_uf2d%unitss)
   endif
   write(*,*) '  --->',trim(prefix),'.',trim(suffix),' ufile write done.'
   write(*,*)

!good. write scalars? 16. AA38601.TRF : Enclosed toroidal flux from eq. code
   ! AA38601.TRF:  38601ITER 1 0 6 ;-SHOT #- F(X) DATA -UF1DWR- 10Jul2014
   ! *TRF(TIM)       { 'Toroidal flux', 'Webers' }
   !-1.0000E+35 ;-SCALAR, LABEL FOLLOWS: TMIN:     Tmin      seconds
   ! 1.0000E+35 ;-SCALAR, LABEL FOLLOWS: TMAX:     Tmax      seconds
   ! 0.0000E+00 ;-SCALAR, LABEL FOLLOWS: TSMOOTH:  delta(t)  seconds
   ! 1.0000E+00 ;-SCALAR, LABEL FOLLOWS: FBDY:     init.bdy
   ! 8.0000E-02 ;-SCALAR, LABEL FOLLOWS: CRAT:     min.curv
   ! 3.2000E-01 ;-SCALAR, LABEL FOLLOWS: GS_ERRMAX:allow.err.
   ! 5.0000E-01 ;-SCALAR, LABEL FOLLOWS: JAC_VARMX:allow.var.
   !-1.0000E+00 ;-SCALAR, LABEL FOLLOWS: BTOR_CCW: Dir. of B
   ! 1.0000E+00 ;-SCALAR, LABEL FOLLOWS: ITOR_CCW: Dir. of I
   ! Time                Seconds   ;-INDEPENDENT VARIABLE LABEL-
   ! EFIT psi(toroidal)  Webers    ;-DEPENDENT VARIABLE LABEL-
   ! 144                    ;-# OF PTS-  X, F(X) DATA FOLLOW:
   ! equilibrium time_slice(:)/profiles_1d/
   ! phi           Toroidal flux {dynamic} [Wb]  FLT_1D  1- time_slice(:)/profiles_1d/psi
   tlen=size(eq%time)
   tmplen=size(eq%time_slice(1)%profiles_1d%psi)

   ndim=1
   prefix='A'
   suffix='TRF'
   t_uf1d%nsc=2
   write(comment,'(a)') 'imported from ids_eq time_slice profiles_1d phi'

   if (t_uf1d%nsc>0) then
      allocate(t_uf1d%SCVAL(t_uf1d%nsc))
      allocate(t_uf1d%labels(t_uf1d%nsc))
      allocate(t_uf1d%unitss(t_uf1d%nsc))
      t_uf1d%SCVAL(1) =  1.
      t_uf1d%SCVAL(2) =  1.
      write(t_uf1d%labels(1),'(a)') 'BTOR_CCW: Dir. of B'
      write(t_uf1d%labels(2),'(a)') 'ITOR_CCW: Dir. of I'
      write(t_uf1d%unitss(1),'(a)') ''
      write(t_uf1d%unitss(2),'(a)') ''
   endif

   t_uf1d%nx=tlen
   if (t_uf1d%nx > 0) then
      t_uf1d%labelx='TIME:'
      t_uf1d%unitsx='seconds'
      t_uf1d%labelf='Toroidal flux:'
      t_uf1d%unitsf='Wb'
      t_uf1d%iproc=0
      allocate(t_uf1d%X(t_uf1d%nx))
      allocate(t_uf1d%F(t_uf1d%nx))
      t_uf1d%X=eq%time
      do jj=1,tlen
         t_uf1d%F(jj)= eq%time_slice(jj)%profiles_1d%phi(tmplen) - &
            eq%time_slice(jj)%profiles_1d%phi(1)
      enddo
      call put_data_to_ufiles(ilun, &
         prefix,suffix,disk,directory, &
         ishot, &
         tdev,ndim,shdate, &
         comment, &
         t_uf0d,t_uf1d,t_uf2d,t_uf3d, &
         ierr)
      deallocate(t_uf1d%X, t_uf1d%F)
   endif
   write(*,*) '  --->',trim(prefix),'.',trim(suffix),' ufile write done.'
   write(*,*)

!good. set NRIxxx=8. 17. B38601.NER : electron density
   ! *NER(TIM,RAD)   { 'electron density',   'n/cm**3' }
   ! time                   seconds
   ! ((Phi/Phi(a))
   ! density             /cm**3
   ! 144
   ! 81
   ! core_profiles profiles_1d(:)/n_e Electron density (thermal+non-thermal)
   !{dynamic} [m^-3] FLT_1D  1- profiles_1d(:)/grid/rho_tor_norm
   tlen=size(cp%time)
   ylen=size(cp%profiles_1d(1)%grid%rho_tor_norm)

   ndim=2
   prefix='A'
   suffix='NER'
   t_uf2d%nsc=0
   write(comment,'(a)') 'imported from ids cp profiles_1d(time) n_e(1:ylen)'

   t_uf2d%nx=tlen
   t_uf2d%ny=ylen
   if (t_uf2d%nx > 0 .and. t_uf2d%ny > 0) then
      t_uf2d%labelx='Time:'
      t_uf2d%unitsx='seconds'
      t_uf2d%labely='normalized toroidal field:'
      t_uf2d%unitsy='-'
      t_uf2d%labelf='electron density:'
      t_uf2d%unitsf='/cm**3'
      t_uf2d%iproc=0
      allocate(t_uf2d%X(t_uf2d%nx))
      allocate(t_uf2d%Y(t_uf2d%ny))
      allocate(t_uf2d%F(t_uf2d%nx,t_uf2d%ny))
      t_uf2d%X=cp%time
      do kk=1,t_uf2d%ny
         t_uf2d%Y(kk)=cp%profiles_1d(1)%grid%rho_tor_norm(kk)
         do jj=1,t_uf2d%nx
            t_uf2d%F(jj,kk)=1.0e-6 * cp%profiles_1d(jj)%electrons%density(kk)
         enddo
      enddo
      call put_data_to_ufiles(ilun, &
         prefix,suffix,disk,directory, &
         ishot, &
         tdev,ndim,shdate, &
         comment, &
         t_uf0d,t_uf1d,t_uf2d,t_uf3d, &
         ierr)
      deallocate(t_uf2d%X, t_uf2d%Y, t_uf2d%F)
   endif
   write(*,*) '  --->',trim(prefix),'.',trim(suffix),' ufile write done.'
   write(*,*)

!good. set NRIxxx=8. 25. B38601.TER
   ! ! electron temperature
   ! *TER(TIM,RAD)   { 'electron temperature',       'eV' }
   ! time                   seconds
   ! ((Phi/Phi(a))
   ! temperature         eV
   ! 0
   ! 144
   ! 81
   ! core_profiles
   ! profiles_1d(:)/t_e  Electron temperature {dynamic} [eV] FLT_1D  1- profiles_1d(:)/grid/rho_tor_norm
   tlen=size(cp%time)
   ylen=size(cp%profiles_1d(1)%grid%rho_tor_norm)

   ndim=2
   prefix='A'
   suffix='TER'
   t_uf2d%nsc=0
   write(comment,'(a)') 'imported from ids cp profiles_1d(time) t_e(1:ylen)'

   t_uf2d%nx=tlen
   t_uf2d%ny=ylen
   if (t_uf2d%nx > 0 .and. t_uf2d%ny > 0) then
      t_uf2d%labelx='Time:'
      t_uf2d%unitsx='seconds'
      t_uf2d%labely='normalized toroidal field:'
      t_uf2d%unitsy='-'
      t_uf2d%labelf='electron teperature:'
      t_uf2d%unitsf='eV'
      t_uf2d%iproc=0
      allocate(t_uf2d%X(t_uf2d%nx))
      allocate(t_uf2d%Y(t_uf2d%ny))
      allocate(t_uf2d%F(t_uf2d%nx,t_uf2d%ny))
      t_uf2d%X=cp%time
      do kk=1,t_uf2d%ny
         t_uf2d%Y(kk)=cp%profiles_1d(1)%grid%rho_tor_norm(kk)
         do jj=1,t_uf2d%nx
            t_uf2d%F(jj,kk)=cp%profiles_1d(jj)%electrons%temperature(kk)
         enddo
      enddo
      call put_data_to_ufiles(ilun, &
         prefix,suffix,disk,directory, &
         ishot, &
         tdev,ndim,shdate, &
         comment, &
         t_uf0d,t_uf1d,t_uf2d,t_uf3d, &
         ierr)
      deallocate(t_uf2d%X, t_uf2d%Y, t_uf2d%F)
   endif
   write(*,*) '  --->',trim(prefix),'.',trim(suffix),' ufile write done.'
   write(*,*)

!18. B38601.NIR
!19. B38601.NZ1
!20. B38601.NZ2
!21. B38601.NZ3
!22. B38601.NZ4
!23. B38601.Nz5
!24. B38601.NZ5
   ! *NIH(TIM,RAD) { 'H+ ion density', 'n/cm**3' }
   ! *NID(TIM,RAD) { 'D+ ion density', 'n/cm**3' }
   ! *NIT(TIM,RAD) { 'T+ ion density', 'n/cm**3' }
   ! *NI3(TIM,RAD) { 'He3++ ion density',  'n/cm**3' }
   ! *NI4(TIM,RAD) { 'He4++ ion density',  'n/cm**3' }
   ! *NI6(TIM,RAD) { 'Li+++ ion density',  'n/cm**3' }
   ! density /cm**3 (144=time seconds, 81=((Phi/Phi(a)))

   ! core_profiles: cp
   ! profiles_1d(:)/ion(:) Quantities related to the different ion species struct_array [max_size=unbounded]
   ! profiles_1d(:)/ion(:)/label String identifying ion (e.g. H+, D+, T+, He+2, C+, ...) {dynamic} STR_0D
   ! profiles_1d(:)/ion(:)/n_i Ion density (thermal+non-thermal) {dynamic} [m^-3]  FLT_1D  1- profiles_1d(:)/grid/rho_tor_norm
   ! 1- 1...N_Ions
   ! ???which one is NIR
   tlen=size(cp%time)
   ylen=size(cp%profiles_1d(1)%grid%rho_tor_norm)
   N_Ions=size(cp%profiles_1d(1)%ion)
   write(*,*)"# of ion species=", N_Ions

   do ll=1,N_Ions  !=5
      tmp1=cp%profiles_1d(1)%ion(ll)%element(1)%a
      tmp2=cp%profiles_1d(1)%ion(ll)%z_ion
      tmp3=cp%profiles_1d(1)%ion(ll)%element(1)%z_n
      write(*,'(2i3,a,e10.3,a,e10.3,a,e10.3)') jj,ll," atom mass unit ",tmp1," ion charge unit ",&
         tmp2," nuclear charge unit ",tmp3
      do jj=2,t_uf2d%nx
         if (tmp1 .ne. cp%profiles_1d(1)%ion(ll)%element(1)%a .or. &
            tmp2 .ne. cp%profiles_1d(1)%ion(ll)%element(1)%z_n) &
            write(*,*) " WARNING : atom mass not equal"
      enddo
   enddo

   do ll=1,N_Ions  !=5
      write(*,*)"ion",ll,"label:",cp%profiles_1d(1)%ion(ll)%label

      ndim=2
      prefix='A'
      write(suffix,'(a,i1)') 'NZ',ll
      t_uf2d%nsc=3
      write(comment,'(a,i2,a)') 'imported from ids cp profiles_1d(time) ion',ll,&
         'ion cunit: ion charge unit; nuc cunit: nuclear charge unit.'

      if (t_uf2d%nsc>0) then
         allocate(t_uf2d%SCVAL(t_uf2d%nsc))
         allocate(t_uf2d%labels(t_uf2d%nsc))
         allocate(t_uf2d%unitss(t_uf2d%nsc))
         t_uf2d%SCVAL(1) =  cp%profiles_1d(1)%ion(ll)%element(1)%a
         t_uf2d%SCVAL(2) =  cp%profiles_1d(1)%ion(ll)%z_ion
         t_uf2d%SCVAL(3) =  cp%profiles_1d(1)%ion(ll)%element(1)%z_n
         write(t_uf2d%labels(1),'(a)') 'mass unit:'
         write(t_uf2d%labels(2),'(a)') 'ion cunit:'
         write(t_uf2d%labels(3),'(a)') 'nuc cunit:'
         write(t_uf2d%unitss(1),'(a)') '-'
         write(t_uf2d%unitss(2),'(a)') '-'
         write(t_uf2d%unitss(3),'(a)') '-'
      endif

      t_uf2d%nx=tlen
      t_uf2d%ny=ylen
      if (t_uf2d%nx > 0 .and. t_uf2d%ny > 0) then
         t_uf2d%labelx='Time:'
         t_uf2d%unitsx='seconds'
         t_uf2d%labely='normalized toroidal field:'
         t_uf2d%unitsy='-'
         t_uf2d%labelf='ion density:'
         t_uf2d%unitsf='n/cm**3'
         t_uf2d%iproc=0
         allocate(t_uf2d%X(t_uf2d%nx))
         allocate(t_uf2d%Y(t_uf2d%ny))
         allocate(t_uf2d%F(t_uf2d%nx,t_uf2d%ny))
         t_uf2d%X=cp%time
         do kk=1,t_uf2d%ny
            t_uf2d%Y(kk)=cp%profiles_1d(1)%grid%rho_tor_norm(kk)
            do jj=1,t_uf2d%nx
               t_uf2d%F(jj,kk)= 1.0e-6 * cp%profiles_1d(jj)%ion(ll)%density(kk)
            enddo
         enddo
         call put_data_to_ufiles(ilun, &
            prefix,suffix,disk,directory, &
            ishot, &
            tdev,ndim,shdate, &
            comment, &
            t_uf0d,t_uf1d,t_uf2d,t_uf3d, &
            ierr)
         deallocate(t_uf2d%X, t_uf2d%Y, t_uf2d%F)
      endif
      if (t_uf2d%nsc>0) then
         deallocate(t_uf2d%SCVAL)
         deallocate(t_uf2d%labels)
         deallocate(t_uf2d%unitss)
      endif
      write(*,*) '  --->',trim(prefix),'.',trim(suffix),' ufile write done.'
      write(*,*)
   enddo

   !.NIR
   ndim=2
   prefix='A'
   write(suffix,'(a)') 'NIR'
   t_uf2d%nsc=0
   write(comment,'(a,i2,a)') 'imported from ids cp profiles_1d(time) ion, sum over', N_Ions, ' n_i'

   t_uf2d%nx=tlen
   t_uf2d%ny=ylen
   if (t_uf2d%nx > 0 .and. t_uf2d%ny > 0) then
      t_uf2d%labelx='Time:'
      t_uf2d%unitsx='seconds'
      t_uf2d%labely='normalized toroidal field:'
      t_uf2d%unitsy='-'
      t_uf2d%labelf='ion density:'
      t_uf2d%unitsf='n/cm**3'
      t_uf2d%iproc=0
      allocate(t_uf2d%X(t_uf2d%nx))
      allocate(t_uf2d%Y(t_uf2d%ny))
      allocate(t_uf2d%F(t_uf2d%nx,t_uf2d%ny))
      t_uf2d%X=cp%time
      do kk=1,t_uf2d%ny
         t_uf2d%Y(kk)=cp%profiles_1d(1)%grid%rho_tor_norm(kk)
         do jj=1,t_uf2d%nx
            t_uf2d%F(jj,kk)=cp%profiles_1d(jj)%n_i_total_over_n_e(kk) * &
               cp%profiles_1d(jj)%electrons%density(kk) * 1.0e-6
         enddo
      enddo
      call put_data_to_ufiles(ilun, &
         prefix,suffix,disk,directory, &
         ishot, &
         tdev,ndim,shdate, &
         comment, &
         t_uf0d,t_uf1d,t_uf2d,t_uf3d, &
         ierr)
      deallocate(t_uf2d%X, t_uf2d%Y, t_uf2d%F)
   endif
   write(*,*) '  --->',trim(prefix),'.',trim(suffix),' ufile write done.'
   write(*,*)

!26. B38601.TI2
   ! ! ion temperature
   ! *TI2(TIM,RAD)   { 'ion temperature',    'eV' }
   ! time                   seconds
   ! ((Phi/Phi(a))
   ! temperature         eV
   ! 0
   ! 144
   ! 81
   ! profiles_1d(:)/t_i_average  Ion temperature (averaged on charge states and ion species) {dynamic} [eV]  FLT_1D  1- profiles_1d(:)/grid/rho_tor_norm
   ! profiles_1d(:)/ion(:)/t_i Ion temperature of that ion species {dynamic} [eV]  FLT_1D  1- profiles_1d(:)/grid/rho_tor_norm
   ! profiles_1d(:)/ion(:)/charge_state(:)/t_z Temperature of the charge state considered {dynamic} [m^-3] FLT_1D  1- profiles_1d(:)/grid/rho_tor_norm
   tlen=size(cp%time)
   ylen=size(cp%profiles_1d(1)%grid%rho_tor_norm)
   N_Ions=size(cp%profiles_1d(1)%ion)
   write(*,*)"# of ion species=", N_Ions
   do ll=1,N_Ions  !=5
      write(*,*)"ion",ll,"label:",cp%profiles_1d(1)%ion(ll)%label

      ndim=2
      prefix='A'
      write(suffix,'(a,i1)') 'TI',ll
      t_uf2d%nsc=3
      write(comment,'(a,i2)') 'imported from ids cp profiles_1d(time) ion #', ll

      if (t_uf2d%nsc>0) then
         allocate(t_uf2d%SCVAL(t_uf2d%nsc))
         allocate(t_uf2d%labels(t_uf2d%nsc))
         allocate(t_uf2d%unitss(t_uf2d%nsc))
         t_uf2d%SCVAL(1) =  cp%profiles_1d(1)%ion(ll)%element(1)%a
         t_uf2d%SCVAL(2) =  cp%profiles_1d(1)%ion(ll)%z_ion
         t_uf2d%SCVAL(3) =  cp%profiles_1d(1)%ion(ll)%element(1)%z_n
         write(t_uf2d%labels(1),'(a)') 'mass unit:'
         write(t_uf2d%labels(2),'(a)') 'ion cunit:'
         write(t_uf2d%labels(3),'(a)') 'nuc cunit:'
         write(t_uf2d%unitss(1),'(a)') ''
         write(t_uf2d%unitss(2),'(a)') ''
         write(t_uf2d%unitss(3),'(a)') ''
      endif

      t_uf2d%nx=tlen
      t_uf2d%ny=ylen
      if (t_uf2d%nx > 0 .and. t_uf2d%ny > 0) then
         t_uf2d%labelx='Time:'
         t_uf2d%unitsx='seconds'
         t_uf2d%labely='normalized toroidal field:'
         t_uf2d%unitsy='-'
         t_uf2d%labelf='ion temperature:'
         t_uf2d%unitsf='eV'
         t_uf2d%iproc=0
         allocate(t_uf2d%X(t_uf2d%nx))
         allocate(t_uf2d%Y(t_uf2d%ny))
         allocate(t_uf2d%F(t_uf2d%nx,t_uf2d%ny))
         t_uf2d%X=cp%time
         do kk=1,t_uf2d%ny
            t_uf2d%Y(kk)=cp%profiles_1d(1)%grid%rho_tor_norm(kk)
            do jj=1,t_uf2d%nx
               t_uf2d%F(jj,kk)=cp%profiles_1d(jj)%ion(ll)%temperature(kk)
            enddo
         enddo
         call put_data_to_ufiles(ilun, &
            prefix,suffix,disk,directory, &
            ishot, &
            tdev,ndim,shdate, &
            comment, &
            t_uf0d,t_uf1d,t_uf2d,t_uf3d, &
            ierr)
         deallocate(t_uf2d%X, t_uf2d%Y, t_uf2d%F)
      endif
      if (t_uf2d%nsc>0) then
         deallocate(t_uf2d%SCVAL)
         deallocate(t_uf2d%labels)
         deallocate(t_uf2d%unitss)
      endif
      write(*,*) '  --->',trim(prefix),'.',trim(suffix),' ufile write done.'
      write(*,*)
   enddo

!    !.TI2
!    ndim=2
!    prefix='A'
!    write(suffix,'(a)') 'TI2'
!    t_uf2d%nsc=0
!    write(comment,'(a,i2,a)') 'imported from ids cp profiles_1d(time) ion, sum over', N_Ions, ' t_i'
!
!    t_uf2d%nx=tlen
!    t_uf2d%ny=ylen
!    if (t_uf2d%nx > 0 .and. t_uf2d%ny > 0) then
!       t_uf2d%labelx='Time:'
!       t_uf2d%unitsx='seconds'
!       t_uf2d%labely='normalized toroidal field:'
!       t_uf2d%unitsy='-'
!       t_uf2d%labelf='ion temperature:'
!       t_uf2d%unitsf='eV'
!       t_uf2d%iproc=0
!       allocate(t_uf2d%X(t_uf2d%nx))
!       allocate(t_uf2d%Y(t_uf2d%ny))
!       allocate(t_uf2d%F(t_uf2d%nx,t_uf2d%ny))
!       t_uf2d%X=cp%time
!       do kk=1,t_uf2d%ny
!          t_uf2d%Y(kk)=cp%profiles_1d(1)%grid%rho_tor_norm(kk)
!          do jj=1,t_uf2d%nx
!             t_uf2d%F(jj,kk)=cp%profiles_1d(jj)%ion(1)%t_i(kk)
!             do ll=2,N_Ions  !=5
!                t_uf2d%F(jj,kk)=t_uf2d%F(jj,kk)+cp%profiles_1d(jj)%ion(ll)%t_i(kk)
!             enddo
!          enddo
!       enddo
!       call put_data_to_ufiles(ilun, &
!                               prefix,suffix,disk,directory, &
!                               ishot, &
!                               tdev,ndim,shdate, &
!                               comment, &
!                               t_uf0d,t_uf1d,t_uf2d,t_uf3d, &
!                               ierr)
!       deallocate(t_uf2d%X, t_uf2d%Y, t_uf2d%F)
!    endif
!    write(*,*) '  --->',trim(prefix),'.',trim(suffix),' ufile write done.'
!    write(*,*)

!not needed. 32. X38601.NHE ! impurity ion density profile  /cm**3
   ! *NIM(TIM,RAD)   { 'impurity density',   'n/cm**3' }
   ! *ZXR(TIM,RAD)   { 'impurity <Z>', ' ' }
   ! *AXR(TIM,RAD)   { 'impurity <A>', 'AMU' }
   ! /SIM
   ! impurity ion density profile by species
   ! Note:  NEVER decrease the sizes associated with IMPS & IMPT -- it will
   !    cause an incompatibility with old TRDAT output files (dmc).
   !    ALSO: for same reason: do not delete this specification, ever!
   ! *SIM(TIM,RAD,-) { 'impurity species',   'n/cm**3' }
   ! (time seconds, ((Phi/Phi(a)))
   ! 0
   ! 144
   ! 81

   call wrgeqdsk(ishot, eq)

   write(*,*)
   write(*,*)
   write(*,*) "========END UFILE WRITING======================"
   write(*,*)
   write(*,*)

101 continue

   call ids_deallocate(eq)
   call ids_deallocate(cp)
   call ids_deallocate(pf)
   call imas_close(idx)

   stop

end program imas2transp

!write transp ufile comment
subroutine write_omment(comment)
   character(len=*) :: comment
   comment="jin chen transp ids translator from march 27 2015"
end subroutine write_omment

!get shot date 'shdate'
subroutine write_shdate(shdate)
   character(len=*) :: shdate
   integer,dimension(8) :: xvalues
   call date_and_time(VALUES=xvalues)
   write(shdate, '(i4,a,i2.2,a,i2.2)'), xvalues(1),'-',xvalues(2),'-',xvalues(3)
end subroutine write_shdate

!get shot number 'ishot'
subroutine write_shotnumber(shot,run,ishot)
   integer,intent(in) :: shot, run
   integer,intent(out) :: ishot
   if (shot.le.99 .and. shot.ge.10) then
      ishot = shot * 10**(1 + (ceiling(log10(real(run))))) + run
   else if (shot.le.999 .and. shot.ge.1000) then
      ishot = shot * 10**(ceiling(log10(real(run)))) + run
   else
      write(*,*) "... err: wrong shot number!"
   endif
end subroutine write_shotnumber

!write geqdsk file per timeslice
subroutine wrgeqdsk(ishot, eq)
   use ids_routines             !! These are the Access Layer routines + management of IDS structures
   USE EZspline_obj
   USE EZspline
   IMPLICIT NONE
!
!external EZspline_setup, EZlinear_init, EZspline_interp, EZspline_error, EZspline_free
!============
   INTEGER ishot
   type (ids_equilibrium)   :: eq  ! Declaration of the empty ids variables to be filled

!============
   character*10 case(6)
   INTEGER i, j, jj
   INTEGER :: ngeq, ngeqindex, nlim
   INTEGER :: istat = 0
   integer,dimension(8) :: xvalues
   REAL twopi
   character*50 ctmp
   REAL rtmp,ztmp

!============
   REAL rdim,zdim,zmid,bcentr,xdum

!============
   REAL, ALLOCATABLE, DIMENSION(:) :: fpol
   REAL, ALLOCATABLE, DIMENSION(:) :: pres
   REAL, ALLOCATABLE, DIMENSION(:) :: ffprim
   REAL, ALLOCATABLE, DIMENSION(:) :: pprime
   REAL, ALLOCATABLE, DIMENSION(:) :: qpsi
   REAL, ALLOCATABLE, DIMENSION(:,:) :: psirz
   REAL, ALLOCATABLE, DIMENSION(:) :: rbbbs,zbbbs,rlim,zlim

!============
   TYPE (EZspline1_r8) :: spln
   integer :: bcs1(2), nt, n1, n1_new, ierr
!REAL, ALLOCATABLE, DIMENSION(:) :: x1, f1
!REAL, ALLOCATABLE, DIMENSION(:) :: x1_new, f1_new
   real(ezspline_r8), allocatable, dimension(:) :: x1, f1
   real(ezspline_r8), allocatable, dimension(:) :: x1_new, f1_new
   REAL :: deltax, deltax_new

!============
! number of horizontal R /vertical Z grid points (nw)
   INTEGER npsi, idum, nw, nh, nbbbs, nlimitr, izmid
   REAL :: rcentr, rleft, rmaxis, zmaxis, simag, sibry, current
   character*128 fnamegeq, fnamegeqindex



!============
   nt = size(eq%time)
   twopi=2.*4.*atan(1.)
   ngeqindex = 51
   fnamegeqindex='gindex.dat'
   write(*,*) 'opening:', trim(fnamegeqindex), ' to write for ', nt, ' timeslices'
   open(ngeqindex,file=trim(fnamegeqindex),form="formatted",status="unknown",iostat=istat)
   write(ngeqindex,'(a,i3)') 'ntimes=',nt

   do jj=1, nt
      ngeq = 52
      write(fnamegeq,'(a,i6.6,a,i6.6,a)') 'ids_', ishot, '_', jj, '.geq'
      write(*,'(a,i3.3,a,a)') 'at timeslice ', jj, ' write ', trim(fnamegeq)
      write(ngeqindex,'(a,f8.3,2x,a,a)') 'time=',eq%time(jj),'filename=',trim(fnamegeq)

      npsi = size(eq%time_slice(jj)%profiles_1d%psi)
      nw = size(eq%time_slice(jj)%profiles_2d(1)%grid%dim1)
      nh = size(eq%time_slice(jj)%profiles_2d(1)%grid%dim2)
      !nbbbs=size(eq%time_slice(jj)%boundary%lcfs%r)
      nlimitr=56

      if (.not.ALLOCATED(fpol)) ALLOCATE(fpol(nw), STAT=istat)
      if (.not.ALLOCATED(pres)) ALLOCATE(pres(nw), STAT=istat)
      if (.not.ALLOCATED(ffprim)) ALLOCATE(ffprim(nw), STAT=istat)
      if (.not.ALLOCATED(pprime)) ALLOCATE(pprime(nw), STAT=istat)
      if (.not.ALLOCATED(qpsi)) ALLOCATE(qpsi(nw), STAT=istat)
      if (.not.ALLOCATED(psirz)) ALLOCATE(psirz(nw,nh), STAT=istat)
      !if (.not.ALLOCATED(rbbbs)) ALLOCATE(rbbbs(nbbbs), STAT=istat)
      !if (.not.ALLOCATED(zbbbs)) ALLOCATE(zbbbs(nbbbs), STAT=istat)
      !if (.not.ALLOCATED(rlim)) ALLOCATE(rlim(nlimitr), STAT=istat)
      !if (.not.ALLOCATED(zlim)) ALLOCATE(zlim(nlimitr), STAT=istat)

! horizontal/vertical dimension in meter of computational box
      rdim = eq%time_slice(jj)%profiles_2d(1)%grid%dim1(nw) - &
         eq%time_slice(jj)%profiles_2d(1)%grid%dim1(1)
      zdim = eq%time_slice(jj)%profiles_2d(1)%grid%dim2(nh) - &
         eq%time_slice(jj)%profiles_2d(1)%grid%dim2(1)
! r in meter of vacuum toroidal magentic field bcentr
! vacuum toroidal magnetic field at rcentr
      rcentr=eq%vacuum_toroidal_field%r0
      bcentr=eq%vacuum_toroidal_field%b0(jj)
! minimum r in meter of rectangular computational box
      rleft = eq%time_slice(jj)%profiles_2d(1)%grid%dim1(1)
! z center of computational box in meter
      izmid=(nh+1)/2
      zmid= eq%time_slice(jj)%profiles_2d(1)%grid%dim2(izmid)

! r,z of magnetic axis in meter
      rmaxis=eq%time_slice(jj)%global_quantities%magnetic_axis%r
      zmaxis=eq%time_slice(jj)%global_quantities%magnetic_axis%z
      ! 0.590266728E+01 0.552840213E+00
! poloidal flux at magnetic axis / plasma boundary in weber/rad
! since psi_boundary < psi_axis, so inverse it
      simag=-eq%time_slice(jj)%global_quantities%psi_axis / twopi
      sibry=-eq%time_slice(jj)%global_quantities%psi_boundary / twopi
      ! 0.163003605E+03 0.145223672E+03 0.530000000E+01
! plasma current in ampere
      current=eq%time_slice(jj)%global_quantities%ip
      xdum=0.

! poloidal current function in m-T f=rB_t on flux grid
!      fpol(:)=eq%time_slice(jj)%profiles_1d%f(:)
! plasma pressure in nt/m^2 on uniform flux grid
!      pres(:)=eq%time_slice(jj)%profiles_1d%pressure(:)
!      ffprim(:)=eq%time_slice(jj)%profiles_1d%f(:) * &
!             eq%time_slice(jj)%profiles_1d%f_df_dpsi(:)
!      pprime(:)=eq%time_slice(jj)%profiles_1d%dpressure_dpsi(:)
!      qpsi(:)=eq%time_slice(jj)%profiles_1d%q(:)
! poloidal flux in weber/rad on the rectangular grid points
! since psi_boundary < psi_axis, so inverse it by psi_axis-psi-psi_boundary
      psirz(:,:)= - eq%time_slice(jj)%profiles_2d(1)%psi(:,:)
      psirz(:,:)= psirz(:,:) / twopi



!......start........................................................................................
      bcs1 = (/ 0, 0 /)
      n1 = npsi
      allocate(x1(n1), f1(n1))
      x1(:) = (eq%time_slice(jj)%profiles_1d%psi(:)  - eq%time_slice(jj)%profiles_1d%psi(1)) / &
         (eq%time_slice(jj)%profiles_1d%psi(n1) - eq%time_slice(jj)%profiles_1d%psi(1))
      do i = 1, n1
         x1(i) = x1(i) + 1.0e-32 * i ! pspline demands a STRICTLY monotonic profile
         !if (i.lt.3 .or. i.gt.(n1-3)) write(*,*) 'debug: x1 ',jj,i,x1(i)
      enddo

      n1_new=nw
      deltax_new=(x1(n1)-x1(1))/real(n1_new-1)
      allocate(x1_new(n1_new), f1_new(n1_new))
      x1_new = (/ (real(i-1) *deltax_new, i=1, n1_new) /)
      !do i = 1, n1_new
      !   if (i.lt.3 .or. i.gt.(n1_new-3)) write(*,*) 'debug: x1_new ', jj,i,x1_new(i)
      !end do

!old grid : n1(number of grid), x1(grid), f1(data to be interpolated)
!new grid : n1_new(number of new grid), x1_new(new grid), f1_new(new data on new grid)
! 1.fpol
      CALL EZlinear_init(spln,n1,ierr)
      CALL EZspline_error(ierr)

      f1 = eq%time_slice(jj)%profiles_1d%f
      spln%x1 = x1
      CALL EZspline_setup(spln,f1,ierr)
      CALL EZspline_error(ierr)

      CALL EZspline_interp(spln,n1_new,x1_new,f1_new,ierr)
      CALL EZspline_error(ierr)
      fpol=-f1_new

      CALL EZspline_free(spln,ierr)
      CALL EZspline_error(ierr)
!     write(*,'(a,x,i3.3)') 'ezspline interpolate fpol from eq%f',jj

!old grid : n1(number of grid), x1(grid), f1(data to be interpolated)
!new grid : n1_new(number of new grid), x1_new(new grid), f1_new(new data on new grid)
! 2.pres
      CALL EZlinear_init(spln,n1,ierr)
      CALL EZspline_error(ierr)

      f1 = eq%time_slice(jj)%profiles_1d%pressure
      spln%x1 = x1
      CALL EZspline_setup(spln,f1,ierr)
      CALL EZspline_error(ierr)

      CALL EZspline_interp(spln,n1_new,x1_new,f1_new,ierr)
      CALL EZspline_error(ierr)
      pres=f1_new

      CALL EZspline_free(spln,ierr)
      CALL EZspline_error(ierr)
!     write(*,'(a,x,i3.3)') 'ezspline interpolate pres from eq%pressure',jj

!old grid : n1(number of grid), x1(grid), f1(data to be interpolated)
!new grid : n1_new(number of new grid), x1_new(new grid), f1_new(new data on new grid)
! 3.ffprim
      CALL EZlinear_init(spln,n1,ierr)
      CALL EZspline_error(ierr)

      f1 = eq%time_slice(jj)%profiles_1d%f_df_dpsi
      spln%x1 = x1
      CALL EZspline_setup(spln,f1,ierr)
      CALL EZspline_error(ierr)

      CALL EZspline_interp(spln,n1_new,x1_new,f1_new,ierr)
      CALL EZspline_error(ierr)
      ffprim=f1_new * twopi

      CALL EZspline_free(spln,ierr)
      CALL EZspline_error(ierr)
!     write(*,'(a,x,i3.3)') 'ezspline interpolate ffprim from eq%f_df_dpsi',jj

!old grid : n1(number of grid), x1(grid), f1(data to be interpolated)
!new grid : n1_new(number of new grid), x1_new(new grid), f1_new(new data on new grid)
! 4.pprime
      CALL EZlinear_init(spln,n1,ierr)
      CALL EZspline_error(ierr)

      f1 = eq%time_slice(jj)%profiles_1d%dpressure_dpsi
      spln%x1 = x1
      CALL EZspline_setup(spln,f1,ierr)
      CALL EZspline_error(ierr)

      CALL EZspline_interp(spln,n1_new,x1_new,f1_new,ierr)
      CALL EZspline_error(ierr)
      pprime=f1_new * twopi

      CALL EZspline_free(spln,ierr)
      CALL EZspline_error(ierr)
!     write(*,'(a,x,i3.3)') 'ezspline interpolate pprime from eq%dpressure_dpsi',jj

!old grid : n1(number of grid), x1(grid), f1(data to be interpolated)
!new grid : n1_new(number of new grid), x1_new(new grid), f1_new(new data on new grid)
! 5.qpsi
      CALL EZlinear_init(spln,n1,ierr)
      CALL EZspline_error(ierr)

      f1 = eq%time_slice(jj)%profiles_1d%q
      spln%x1 = x1
      CALL EZspline_setup(spln,f1,ierr)
      CALL EZspline_error(ierr)

      CALL EZspline_interp(spln,n1_new,x1_new,f1_new,ierr)
      CALL EZspline_error(ierr)
      qpsi=f1_new

      CALL EZspline_free(spln,ierr)
      CALL EZspline_error(ierr)
!     write(*,'(a,x,i3.3)') 'ezspline interpolate qpsi from eq%q',jj
!.......end.........................................................................................


! r of boundary points in meter
! z of boundary points in meter
      !rbbbs(:)=eq%time_slice(jj)%boundary%lcfs%r(:)
      !zbbbs(:)=eq%time_slice(jj)%boundary%lcfs%z(:)
!     write(*,'(a,x,i3.3)') 'get bdy rz from eq%lcfs',jj
! r of surrounding limiter contour in meter
! z of surrounding limiter contour in meter
      ! not this one rlim=eq%time_slice(jj)%boundary%active_limiter_point%r
      ! not this one zlim=eq%time_slice(jj)%boundary%active_limiter_point%z
      ! should read from excel file write(*,'(a,x,i3.3)') 'get limiter from eq%active_limier',jj
#if 0
      nlim=60
      open(nlim,file='limiter_v3.4.txt',status="unknown",iostat=istat)
      j=0
      do i=1,nlim
         if (i.le.2 .or. (i.ge.22 .and.i.le.23)) then
            read(nlim,*) ctmp
!         if (jj.eq.1) write(*,*) 'limiter_txt:',i, trim(ctmp)
         else
            !read(nlim,*) rlim(j), zlim(j)
            j=j+1
            read(nlim,*) rtmp, ztmp
            if (j.le.19) then
               rlim(19-j+1)=rtmp
               zlim(19-j+1) = ztmp
            else
               rlim(j)=rtmp
               zlim(j)=ztmp
            endif
         endif
      enddo
      if (j.ne.nlimitr) then
         write(*,*) "limiter_txt: err in reading"
      else
         do j=1,nlimitr
            if (jj.eq.1) write(*,*) 'limiter_txt:',j,rlim(j), zlim(j)
         enddo
      endif
      close(nlim)
#endif
!============
!write GEQDSK ASCII file in the format used by EFIT (added 6/04/03)
!============

      call date_and_time(VALUES=xvalues)
      case(1)='ids'
      case(2)='g'
      !case(3)='05/07'
      write(case(3), '(i2.2,a,i2.2)'), xvalues(2),'/',xvalues(3)
      !case(4)='2015'
      write(case(4), '(i4.4)'), xvalues(1)
      case(5)='iter'
      write(case(6), '(f8.3)') eq%time(jj)
      !case(6)='ms'
      open(ngeq,file=trim(fnamegeq),form="formatted",status="unknown",iostat=istat)
      write (ngeq,2000) (case(i),i=1,6),idum,nw,nh
      write (ngeq,2020) rdim,zdim,rcentr,rleft,zmid
      write (ngeq,2020) rmaxis,zmaxis,simag,sibry,bcentr
      write (ngeq,2020) current,simag,xdum,rmaxis,xdum
      write (ngeq,2020) zmaxis,xdum,sibry,xdum,xdum
      write (ngeq,2020) (fpol(i),i=1,nw)
      write (ngeq,2020) (pres(i),i=1,nw)
      write (ngeq,2020) (ffprim(i),i=1,nw)
      write (ngeq,2020) (pprime(i),i=1,nw)
      write (ngeq,2020) ((psirz(i,j),i=1,nw),j=1,nh)
      write (ngeq,2020) (qpsi(i),i=1,nw)
      write (ngeq,2022) nbbbs,nlimitr
      !write (ngeq,2020) (rbbbs(i),zbbbs(i),i=1,nbbbs)
      !write (ngeq,2020) (rlim(i),zlim(i),i=1,nlimitr)
      close(ngeq)

102   continue
      deallocate(x1, f1, x1_new, f1_new)
      !DEALLOCATE(fpol, pres, ffprim, pprime, qpsi, psirz, rbbbs, zbbbs, rlim, zlim)
      DEALLOCATE(fpol, pres, ffprim, pprime, qpsi, psirz)

   enddo  !finish jj timeslice

   close(ngeqindex)

2000 format (6a8,3i4)
2020 format (5e16.9)
2022 format (2i5)
end subroutine wrgeqdsk



!-----------------------------------------------------------------!

!module toray_ga_global
!
!  integer nprof
!
!end module toray_ga_global

subroutine t1mhdeq_geq(runid, time)
  !
  ! test driver for toray-ga
  ! 
  implicit none
  integer, parameter :: r8 = selected_real_kind(12,100)

  character(*) :: runid
  real(r8) :: time, deltat
  !  direction of Btoroidal and Itoroidal
  integer ibccw,ipccw
  ! grid size
  integer nR, nZ, nt1, ns
  integer ier

  !write(*,*)' $Id: toray_ga.f90,v 1.13 2007-01-25 23:28:51 Doug_McCune Exp $ '

  !ier = 0

!  call write_echin_file(runid, time, deltat)

!  call write_geqdsk_file(runid, time, deltat, ibccw, ipccw, nR, nZ, nt1, ns, ier)
!  if(ier/=0) then
!    print *, ' error after mhdcomputgeq'
!    call bad_exit
!  end if

!  call write_toray_namelist(runid, time, ier)
!  if(ier/=0) then
!    print *, ' error after write_toray_namelist'
!    call bad_exit
!  end if

!  ier = 0
!  call torGa_torays   ! dmc -- at present, no status code is returned.
!  ! (ier) argument removed.

!  if(ier/=0) then
!    print *, ' error after torGa_torays'
!    call bad_exit
!  end if
!end subroutine t1mhdeq_geq


!......................................................................................


!subroutine write_toray_namelist(runid, time, ier)
!  use key_access
!  use toray_ga_global
!  implicit none
!  integer, parameter :: r8 = selected_real_kind(12,100)
!  character(*), intent(in) :: runid
!  real(r8), intent(in) :: time
!  integer, intent(out) :: ier
!
!  ier = 0
!  
!  call ka_init
!
!  ! numerics
!  call ka_set('igafit', 1) ! must be 1 for geqdsk 
!  call ka_set('smax', 150.0_r8) ! max arclength in cm
!  call ka_set('ds', 0.5_r8) ! initial step size in cm
!  call ka_set('relerr', 5.0e-5_r8) ! rel ODE error
!  call ka_set('abserr', 5.0e-5_r8) ! absolute error
!  call ka_set('mn0', 32) ! no of points in Gauss integration
!
!  ! physics
!  call ka_set('ezeff', 2.0_r8) ! Zrmeff
!  call ka_set('nharm', 2) ! harmonic number
!  call ka_set('bsratio', 1.0_r8) ! antenna elongation
!  call ka_set('modelc', 4) ! model for ECCD calculation
!  
!  ! rays
!  call ka_set('gauszone', 4) ! no of zones 
!  call ka_set('mray', (/1, 5, 12, 12/)) ! no of rays per zone
!  call ka_set('cr', (/0.0_r8, 0.1_r8, 0.05_r8, -0.05_r8/)) ! azimutal phase for each zone
!
!  ! netcdf output
!  call ka_set('netcdfdat', .TRUE.)
!
!  print *,'  toray namelist'
!  call ka_print
!
!  open(1, file='toray.in', form='formatted')
!  write(1,'(a)') runid
!  close(1)
!  call ka_write_namelist('toray.in', 'EDATA', 'a') ! append mode
!  
!  call ka_free
!
!  call ka_init
!  
!  call ka_set('npts', nprof) ! no of nodes in transport code ! should match dim in echin????
!  call ka_set('ipsi', 1)  ! 1 for geqdsk input
!  call ka_set('gafsep',0.01_r8)  ! ?
!  call ka_set('dsrat', 0.02_r8)  ! ?
!  call ka_set('tolmap',1.e-10_r8) ! tolerance
!  call ka_set('percenflux', 1.0_r8) !0.994_r8) ! to separatrix
!  call ka_set('newbdry', .FALSE.)
!  
!  print *,'gafit namelist'
!  call ka_print
!
!  call ka_write_namelist('gafit.in', 'FITDAT', 'w')
!
!  call ka_free
!
!end subroutine write_toray_namelist

!......................................................................................

!subroutine write_echin_file(runid, time, deltat)
!
!  use toray_ga_global
!  implicit none
!  integer, parameter :: r8 = selected_real_kind(12,100)
!  character(*), intent(in) :: runid
!  real(r8) :: time, deltat
!  
!  character(50) :: rlabel
!  integer nsctime,nprtime,nxmax,nmax,ier,istype,ngot
!  character(100) :: zlabel, zunits
!
!  ! scalars
!  real BZXR, RAXIS
!
!  ! mhd
!  integer nsmax,ntheta,nsgot,i
!  real, allocatable :: theta(:), Rarr(:,:), Zarr(:,:)
!  real, allocatable :: rho(:), psi(:), pmhd(:), qmhd(:), gmhd(:)
!  real tflux, pcur, pi
!
!  ! arrays
!  real, allocatable :: XB(:), TE(:), NE(:), ZEFFI(:)
!
!  ! toray echin stuff
!  integer idamp,j12,nray,nbfld
!  real(r8) :: fmu0,rfmod,x00,z00,thet, phai,bhalf,bsratio,rmajs,b0,rmins
!  real(r8) :: r0
!  integer :: j
!  real(r8),allocatable, dimension(:) :: psinrmr, zef, ene, ete
!
!
!  CALL KCONNECT(runid,rlabel,nsctime,nprtime,nxmax,nmax,ier)
!
!  print *, rlabel
!  print *, nsctime, ' timepoints in scalar functions'
!  print *, nprtime, ' timepoints in profile functions'
!  print *, nxmax, ' x axis points'
!
!  ! read scalars
!
!  call t1scalar('BZXR', zlabel, zunits,real(time),real(deltat),BZXR,ier)
!  call t1scalar('RAXIS', zlabel, zunits,real(time),real(deltat),RAXIS,ier)
!
!  print *, 'Bz R0 = ', BZXR
!  print *, 'Rm    = ', RAXIS
!
!  ! read profiles
!
!  allocate(XB(nxmax))
!  allocate(NE(nxmax))
!  allocate(TE(nxmax))
!  allocate(ZEFFI(nxmax))
!
!  call t1profil('XB',zlabel,zunits,real(time),real(deltat),istype,XB(2),nmax,ngot,ier)
!  print *,'ngot=', ngot, zlabel, zunits
!  XB(1) = 0._r8
!  print *,'XB=', XB(1:ngot)
!
!  call t1profil('NE',zlabel,zunits,real(time),real(deltat),istype,NE(2),nmax,ngot,ier)
!  print *,'ngot=', ngot, zlabel, zunits
!  NE(1) = NE(2)
!  print *,'NE=', NE(1:ngot)
!
!  call t1profil('TE',zlabel,zunits,real(time),real(deltat),istype,TE(2),nmax,ngot,ier)
!  print *,'ngot=', ngot, zlabel, zunits
!  TE(1) = TE(2)
!  print *,'TE=', TE(1:ngot)
!
!  call t1profil('ZEFFI',zlabel,zunits,real(time),real(deltat),istype,ZEFFI(2),nmax,ngot,ier)
!  print *,'ngot=', ngot, zlabel, zunits
!  ZEFFI(1) = ZEFFI(2)
!  print *,'ZEFFI=', ZEFFI(1:ngot)
!
!  ! read r, z grid
!
!  ntheta = 33
!  nsmax = nxmax
!  allocate(theta(ntheta))
!  allocate(Rarr(ntheta, nsmax), Zarr(ntheta, nsmax))
!  allocate(rho(nsmax), psi(nsmax), pmhd(nsmax), qmhd(nsmax), gmhd(nsmax))
!  pi = acos(-1.0)
!  theta = 2.*pi*(/ (real(i-1)/real(ntheta-1), i=1, ntheta) /)
!  call t1mhdeq(real(time),real(deltat),nsmax,ntheta,nsgot,'TRANSP',theta,Rarr,Zarr, &
!       & rho,psi,pmhd,qmhd,gmhd,tflux, pcur, ier)
!  print *, 'nsgot = ', nsgot
!  print *, 'minor radius = ', Rarr(1, nsgot) - Rarr(1,1)
!
!  ! write echin file
!  open(2, file='echin', form='formatted', iostat=ier)
!
!  j12 = nsgot
!  nprof = nsgot
!
!  write(*,*)' Damping model (idamp=2): '
!  write(*,*)' 0: none'
!  write(*,*)' 1: Matsuda-Hu near 1st harmonic, weakly relativistic'
!  write(*,*)' 2: Matsuda-Hu near 2nd harmonic, weakly relativistic'
!  write(*,*)' 3: Matsuda-Hu near 3rd harmonic, weakly relativistic'
!  write(*,*)' 7: Matsuda-Hu between 2nd and 3rd harmonic'
!  write(*,*)' 8: Mazzucato near 2nd harmonic, relativistic'
!  write(*,*)' enter idamp (-1 for default value) '
!  read(5,*) idamp; if(idamp<0 .or. idamp>8 .or. (idamp>=4 .and. idamp<=6)) idamp=2
!
!  write(*,*)' Enter number of rays (nray=30) or -1 for default '
!  read(5,*) nray; if(nray<0) nray=30
!
!  nbfld=3 ! efit model
!
!  write(*,*)' Enter frequency [Hz] (fmu0=110.e9) or -1 for default '
!  read(5,*) fmu0; if(fmu0<0._r8) fmu0 = 110.e9_r8
!
!  write(*,*)' Enter fraction of power launched as O mode (rfmod=1.0) or -1 for default '
!  read(5,*) rfmod; if(rfmod<0._r8 .or. rfmod>1._r8) rfmod = 1.0;
!
!  write(*,*)' Enter position of source [cm] (x00=239.73 z00=69.42) for -1,* for default '
!  read(5,*) x00, z00 
!  if(x00<0._r8) then
!     x00 = 239.73; z00 = 69.42
!  endif
!
!  write(*,*)' Enter poloidal (thet=115.118) and toroidal (phai=204.459) angles [deg] or 360., 360. for default '
!  read(5,*) thet, phai
!  if(thet>=360._r8) thet=115.118_r8; if(phai>=360._r8) phai=204.459_r8;
!  
!  write(*,*)' Enter beam divergence angle bhalf [deg] (1.7) (half angle at half power) or -1 for default '
!  read(5,*) bhalf
!  if(bhalf<0._r8 .or. bhalf>180._r8) bhalf = 1.7_r8
!
!  write(*,*)' Enter beam aspect ratio (bsratio=1.0) or -1. for default '
!  read(5,*) bsratio; if(bsratio<0._r8 .or. bsratio>1._r8)  bsratio=1.0
!  
!  ! keep cm
!  rmajs = RAXIS
!  ! it shouldn't matter whether r0 is raxis ot the geometric center
!  r0 = RAXIS !(Rarr(1, nsgot) + Rarr((ntheta-1)/2, nsgot))/2._r8
!  b0 = BZXR/RAXIS !BZXR/r0
!  rmins = 100._r8*sqrt(tflux/(pi*b0)) ! in cm
!  ! -> Gauss
!  b0 = 1.e4_r8 * b0
!
!  allocate(psinrmr(j12), zef(j12), ene(j12), ete(j12))
!
!  if(size(rho)<j12) print *, ' size error for rho'
!  if(rho(j12)/=1.) print *, ' rho not equal to 1 at edge'
!  if(rho(1  )/=0.) print *, ' rho not equal to 0 on axis'
!  ! rho is sqrt(phi/phi_a)
!  !psinrmr(1:j12) = rho(1:j12)**2  ! original (Alex's input)
!   psinrmr(1:j12) = rho(1:j12)
!
!  if(size(ZEFFI)<j12) print *, ' size error for ZEFFI'
!  zef(1:j12) = ZEFFI(1:j12)
!
!  if(size(NE)<j12) print *, ' size error for NE'
!  ene(1:j12) = NE(1:j12)
!
!  if(size(TE)<j12) print *, ' size error for TE'
!  ! -> keV
!  ete(1:j12) = TE(1:j12)/1.e3_r8
!  
!  write (2, 1001) time
!  write (2, 1000) idamp,j12,nray,nbfld
!  write (2, 1001) fmu0,rfmod,x00,z00,thet, phai,bhalf,bsratio,rmajs, &
!       &   b0,rmins
!
!  write (2, 1001) (psinrmr(j), j=1,j12)
!  write (2, 1001) (zef(j),j=1,j12)
!  write (2, 1001) (ene(j),j=1,j12)
!  write (2, 1001) (ete(j),j=1,j12)
!1000 format (20i4)
!1001 format (5e16.9)
!  
!
!  close(2)
!
!  deallocate(psinrmr, zef, ene, ete)
!
!  deallocate(theta)
!  deallocate(Rarr, Zarr)
!  deallocate(rho, psi, pmhd, qmhd, gmhd)
!  
!  deallocate(XB)
!  deallocate(NE)
!  deallocate(TE)
!  deallocate(ZEFFI)
!
!end subroutine write_echin_file


!......................................................................................


!subroutine write_geqdsk_file(runid, time, deltat, ibccw, ipccw, nR, nZ, nt1, ns, ier)

  ! extract the data from mdsplus and write geqdsk file (=eqdskin)

!  implicit none
!  integer, parameter :: r8 = selected_real_kind(12,100)

!  character(*), intent(in) :: runid
!  real(r8), intent(inout) :: time
!  real(r8), intent(in) :: deltat
!  integer, intent(in) :: ibccw, ipccw, nR, nZ, nt1, ns
!  integer, intent(out) :: ier

  ! locals
  real(r8) :: tmax, tmin
  integer iwarn

  !  GEQDSK id label & other GEQDSK stuff
  !
  character(120) geqdsk_lbl
  !
  integer nh_geqdsk,nv_geqdsk,nb_geqdsk
  real(r8) zcur
  integer :: id_psi=0,id_pmhd=0,id_q=0,id_g=0

  !  (R,Z) grid limits
  real(r8) Rmin,Rmax,Zmin,Zmax
  
  ! misc
  integer irzflag
  character(10) ctime, cdate, cnow


  write(*,*)'entered TRANSP or MDSPlus runid->', trim(runid)
  !write(*,*)"(eg 'MDS+:TRANSPGRID.PPPL.GOV:TRANSP(NSTX.00,11115P07)'"
  !read(5,*) runid
  write(*,*)'entered time slice->', time
  !read(5,*) time
  
  deltat = 0.01_r8
  ibccw  = 1
  ipccw  = 1
  nR = 65    ! this must be an ODD number
  nZ = 65    ! this must be an ODD number
  nt1 = 65
  ns = 65

  ier = 1
  call trx_connect(runid, ier)
  if(ier/=0) return
  ier = ier + 1
  call trx_tlims(tmin, tmax, ier)
  if(ier/=0) return
  write(*,*)'time limits: ', tmin, tmax
  ier = ier + 1
  call trx_time(time, 0.01_r8, iwarn, ier)
  write(*,*)'time of interest: ', time
  if(ier/=0) return
  if(iwarn/=0) then 
     print *,' oops no time slice @ ', time, ' +- 0.01 '
     time = max(min(time, tmax), tmin)
  endif
  ier = ier + 1
  call trx_init_xplasma(ier)
  if(ier/=0) return 
  ier = ier + 1
  call trx_mhd(nt1,ier)
  if(ier/=0) return

  Rmin = 0._r8; Rmax = 0._r8;
  Zmin = 0._r8; Zmax = 0._r8;
  ier = ier + 1
  call trx_wall_RZ(99,Rmin,Rmax,nR,Zmin,Zmax,nZ,ier)
  if(ier/=0) return

  !
  !  get B field components on flux grid;
  !  extrapolate B field components onto the (R,Z) grid
  !
  !  after this call we have Bmod(rho,theta)
  !                          BR(rho,theta) or BR(R,Z)
  !                          BZ(rho,theta) or BZ(R,Z)
  !                                           Bphi(R,Z)
  !
  irzflag=1
  ier = ier + 1
  call trx_bxtr(ibccw,ipccw,irzflag,ier)
  if(ier/=0) return

  write(ctime, '(f10.4)') time
  call date_and_time(date=cdate, time=cnow)
  geqdsk_lbl = runid//'T='//trim(ctime)// &
       & '::'//cdate//'_'//cnow
  !write(*,*) 'Geqdsk file label: ', geqdsk_lbl

  call eq_gfnum('pmhd',id_pmhd)
  call eq_gfnum('q',id_q)

  ! temprorarily set to zero
  zcur = 0._r8
  nh_geqdsk=nZ
  nv_geqdsk=nR
  nb_geqdsk=ns
  open(99, file='eqdskin', form='formatted')
  ier = ier + 1
  call t1mhdeq_compgeq(99,geqdsk_lbl,Rmin,Rmax,Zmin,Zmax,zcur, &
       &        id_pmhd,id_q,nh_geqdsk,nv_geqdsk,nb_geqdsk,ier)
  close(99)
  if(ier/=0) return

  ! ok
  ier = 0
!end subroutine write_geqdsk_file
end subroutine t1mhdeq_geq


!......................................................................................


subroutine t1mhdeq_compgeq(lun_geqdsk,geqdsk_lbl, &
     Rmin,Rmax,Zmin,Zmax,zcur, &
     id_p,id_q,nh,nv,nb,ierr)
  
  !  write GEQDSK file from xplasma module

  !  **correction** dmc 9 Nov 2001:
  !  follow Lang Lao sign convention for G-EQDSK files:
  !  -- psi always increasing
  !  -- direction of current specified in pcur scalar *only*

  use xplasma_obj_instance
  !use eq_module
  implicit NONE

  INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)

  !  input:

  integer lun_geqdsk                ! logical unit where to write file

  character*48 geqdsk_lbl           ! 48 character label for GEQDSK file
  real*8 Rmin,Rmax                  ! (Rmin,Rmax) of Psi(R,Z)
  real*8 Zmin,Zmax                  ! (Zmin,Zmax) of Psi(R,Z)
  real*8 zcur                       ! plasma current (amps)

  !  [Rmin,Rmax]x[Zmin,Zmax] must contain the core plasma but not exceed the
  !  available (R,Z) grids.

  integer id_p                      ! xplasma id:  Pressure profile
  integer id_q                      ! xplasma id:  q profile

  integer nh                        ! #of GEQDSK horizontal pts
  integer nv                        ! #of GEQDSK vertical pts

  !  note nh also controls the number of pts in the 1d profiles

  integer nb                        ! #of pts in bdy contour

  !  output:

  integer ierr                      ! completion code (0=OK, file written).

  !------------------------------------------------------
  ! (old f77 xplasma eq_geqdsk code moved into core modules; most arguments
  ! made optional).

  call xplasma_wr_geqdsk(s,ierr, &
       lun_geqdsk=lun_geqdsk, label=geqdsk_lbl, &
       Rmin=Rmin, Rmax=Rmax, Zmin=Zmin, Zmax=Zmax, &
       cur=zcur, id_pprof=id_p, id_qprof=id_q, &
       nh=nh, nv=nv, nbdy=nb)

  if(ierr.ne.0) then
     write(6,*) ' %eq_geqdsk: error in xplasma_wr_geqdsk:'
     call xplasma_error(s,ierr,6)
  endif

end subroutine t1mhdeq_compgeq


!......................................................................................


module xplasma_profs

  ! module with routines that help to create various types of xplasma profiles

  use xplasma_obj
  use xplasma_ctran
  use xplasma_sol
  use xplasma_rzgeo
  use xplasma_flxint
  use eqi_rzbox_module

  implicit NONE

  private

  public :: xplasma_rzprof,xplasma_rzprof_fun
  public :: xplasma_brz,xplasma_brz_extrap
  public :: xplasma_irhofun
  public :: xplasma_geqdsk_rewrite
  public :: xplasma_wr_geqdsk,xplasma_rhopsi_gen,xplasma_rhopsi_find

  contains

    subroutine xplasma_brz_extrap(s,ier, ispline,sm_edge)

      ! create extrapolated B(R,Z) and Psi(R,Z) by standard method, when
      ! starting only with a prescribed boundary core (inverse representation)
      ! equilibrium.

      ! do not use this if a free boundary equilibrium is available...

      type(xplasma), pointer :: s
      integer, intent(out) :: ier ! exit status code 0=OK

      integer, intent(in), optional :: ispline   ! spline fit control
      !  DEFAULT =1, C1 Akima Hermite
      !   ...use 0 for piecewise linear; use 2 for C2 bicubic splines

      real*8, intent(in),optional :: sm_edge   ! edge smoothing
      !  smoothin in vicinity +/- sm_edge (meters) from plasma boundary
      !  i.e. option to smooth transition btw field internal to plasma and
      !  the field outside; default = no smoothing.

      !------------------------
      external eqm_brz_adhoc
      !------------------------

      call xplasma_brz(s,eqm_brz_adhoc,ier, ispline,sm_edge)

    end subroutine xplasma_brz_extrap

    subroutine xplasma_brz(s,userbvec,ier, ispline,sm_edge)

      ! create B(R,Z) and Psi(R,Z), extending fields already defined
      ! inside the plasma to the region outside-- using user provided
      ! subroutine...

      type(xplasma), pointer :: s

      external userbvec
      integer, intent(out) :: ier ! exit status code 0=OK

      integer, intent(in), optional :: ispline   ! spline fit control
      !  DEFAULT =1, C1 Akima Hermite
      !   ...use 0 for piecewise linear; use 2 for C2 bicubic splines

      real*8, intent(in),optional :: sm_edge   ! edge smoothing
      !  smoothin in vicinity +/- sm_edge (meters) from plasma boundary
      !  i.e. option to smooth transition btw field internal to plasma and
      !  the field outside; default = no smoothing.

      !-------------------------
      !  set up Bphi(R,Z), BR(R,Z), and BZ(R,Z) -- Akima-Hermite interpolation
      !  for axisymmetric case, setup psi(pol) also.

      !    data values from user supplied function  of form:

      !        subroutine userbvec(iv,zR,zZ,zphi,init,BR,BZ,BPHI,Psi,kpsi,ierr)
      !          iv -- vector dimension
      !          zR(iv),zZ(iv),zphi(iv)   ! coordinates at which to evaluate
      !          init -- action code
      !          BR(iv),BZ(iv),BPHI(iv)   ! field components returned
      !          Psi(iv)                  ! Psi values returned
      !          kpsi                     ! =1 if Psi values were set
      !          ierr                     ! completion code, 0=OK

      !  the subroutine must be capable of providing BR,BZ,BPHI, and Psi on
      !  points (R(i),Z(j)) which are on the extrapolated (R,Z) grids
      !  i.e. __RGRID and __ZGRID.  Only the points which are beyond the 
      !  plasma boundary need be filled; others can be zero.
      !

      !     normal calls to userbvec use init=0
      !         returns BR,BZ,BPHI, and Psi
      !     an initialization call uses init=1, returns BR=BZ=BPHI=czero
      !         also returns kpsi=1 if a poloidal flux function is to be
      !         returned. 

      !     a cleanup call uses init=2, returns BR=BZ=BPHI=czero

      !------------------------------
      !  This routine is not required for setting up B and Psi on an extended
      !  (R,Z) grid; the alternative is to call xplasma_rzprof separately for
      !  Psi and for each component of B-- as this routine itself does at the
      !  end...
      !------------------------------

      integer :: id_Rg,id_Zg
      integer :: nR,nZ,iR,iZ
      real*8, dimension(:), allocatable :: Rg,Zg,Rtmp,Ztmp,Phidum
      real*8, dimension(:,:), allocatable :: BRx,BZx,Bphix,Bmodx,Psix
      real*8 :: zdum1(1),zdum2(1),zdum3(1),zdum4(1)
      real*8 :: zsm_edge
      integer :: idum,jvec,ii,iersum,kpsi

      integer :: id_psi,id_BR,id_BZ,id_Bmod,id_Bphi,id_out,jspline,iertmp
      integer :: id_g,bphi_ccw,id_map,lrhomap
      logical :: axisymm,scrapeoff,outside_only

      real*8, dimension(:), pointer :: eqbuf
      integer, dimension(:), pointer :: idata

      real*8, parameter :: ZERO = 0.0d0
      !------------------------------

      sp => s
               
      call xplasma_global_info(s,ier, axisymm=axisymm,scrapeoff=scrapeoff, &
           bphi_ccw=bphi_ccw)
      if(ier.ne.0) return

      if(.not.axisymm) then
         ier=107
         call xplasma_errmsg_append(s,'  xplasma_brz still requires axisymmetry!')
         return
      endif

      if(.not.scrapeoff) then
         ier=109
         call xplasma_errmsg_append(s,'  xplasma_brz needs scrapeoff region defined.')
         return
      endif

      !  find needed resources...

      call xplasma_common_ids(s,ier, &
           id_psi=id_psi,id_BR=id_BR,id_BZ=id_BZ,id_Bmod=id_Bmod,id_g=id_g)
      if(ier.ne.0) return

      call xplasma_find_item(s,'__RGRID',id_Rg,ier)
      if(ier.ne.0) return

      call xplasma_find_item(s,'__ZGRID',id_Zg,ier)
      if(ier.ne.0) return

      call xplasma_find_item(s,'Bphi',id_Bphi,ier)
      if(ier.ne.0) return

      !  fetch grids

      call xplasma_grid_size(s,id_Rg,nR,ier)
      if(ier.ne.0) return

      call xplasma_grid_size(s,id_Zg,nZ,ier)
      if(ier.ne.0) return

      !  initialize userbvec

      zdum1=ZERO; zdum2=ZERO; zdum3=ZERO
      call userbvec(1,zdum1,zdum2,zdum3,1,zdum1,zdum2,zdum3,zdum4,kpsi,ier)

      if(ier.ne.0) then
         call xplasma_errmsg_append(s, &
              ' ?xplasma_brz: external call to userbvec failed.')

         return
      endif

      if(kpsi.eq.1) then
         zsm_edge=ZERO  ! if (smooth Psi) is constructed, no further smoothing
      else
         ! If Psi(ext) imposed, smoothing is an option.
         zsm_edge=ZERO
         if(present(sm_edge)) then
            zsm_edge=sm_edge
         endif
      endif

      call xplasma_find_item(s,'__FASTMAP',id_map,iertmp,nf_noerr=.TRUE.)
      if(id_map.eq.0) then
         call eqi_fastinv_gen(id_Rg,id_Zg,nR,nZ,id_map,ier)
         if(ier.ne.0) then
            call xplasma_errmsg_append(s, &
                 ' ?xplasma_brz: fast map setup failed.')
            return
         endif
      endif

      call xplasma_blackbox_retrieve(s, id_map, iertmp, &
           ia_ptr=idata, r8a_ptr=eqbuf)
      lrhomap=idata(1)

      jspline=1
      if(present(ispline)) jspline=ispline

      !  OK...

      iersum=0

      allocate(Rg(nR),Zg(nZ))
      allocate(Rtmp(nR*nZ),Ztmp(nR*nZ),Phidum(nR*nZ))
      allocate(BRx(nR,nZ),BZx(nR,nZ),Bphix(nR,nZ),Bmodx(nR,nZ),Psix(nR,nZ))

      call xplasma_grid(s,id_Rg,Rg,iertmp)
      call xplasma_grid(s,id_Zg,Zg,iertmp)

      Phidum=ZERO

      ii=0
      do iZ=1,nZ
         Rtmp(ii+1:ii+nR)=Rg     ! (1:nR)
         Ztmp(ii+1:ii+nR)=Zg(iZ)
         ii=ii+nR
      enddo

      call userbvec(nR*nZ,Rtmp,Ztmp,Phidum,0, &
              BRx,BZx,Bphix,Psix,idum,ier)
      iersum=iersum+ier

      !  generally profiles may have been evaluated at points beyond plasma
      !  boundary only.  If so Bphix(iR,iZ) will have zeroes; it is 
      !  advantageous to fill in Bphi now using g(rho)/R formula...

      jvec=0
      do iZ=1,nZ
         do iR=1,nR
            if(Bphix(iR,iZ).eq.ZERO) then
               jvec=jvec+1
               exit
            endif
         enddo
         if(jvec.gt.0) exit
      enddo

      outside_only = (jvec.gt.0)

      !  cleanup call for userbvec...

      zdum1=ZERO; zdum2=ZERO; zdum3=ZERO
      call userbvec(1,zdum1,zdum2,zdum3,2,zdum1,zdum2,zdum3,zdum4,idum,ier)
      iersum=iersum+ier

      if(iersum.ne.0) then
         call xplasma_errmsg_append(s,' ?xplasma_brz: extrapolated field setup failed.')
      else

         call xplasma_author_set(s,xplasma_xmhd,iertmp)

         if(outside_only) then
            ! code "id_fun_in=-1" indicates g(rho)/R
            call chk_auth('Bphi_RZ')
            if(id_Bphi.gt.0) then
               call xplasma_rzprof(s,'Bphi_RZ',id_out,ier, &
                    ispline=jspline, sm_edge=zsm_edge, id_fun_in=id_Bphi, &
                    data=Bphix,label='B_phi on (R,Z) grid',units='T')
            else
               call xplasma_rzprof(s,'Bphi_RZ',id_out,ier, &
                    ispline=jspline, sm_edge=zsm_edge, id_fun_in=-1, &
                    data=Bphix,label='B_phi on (R,Z) grid',units='T')
            endif
            iersum=iersum+ier

            call chk_auth('BR_RZ')
            call xplasma_rzprof(s,'BR_RZ',id_out,ier, &
                 ispline=jspline, sm_edge=zsm_edge, id_fun_in=id_BR, &
                 data=BRx,label='B_R on (R,Z) grid',units='T')
            iersum=iersum+ier

            call chk_auth('BZ_RZ')
            call xplasma_rzprof(s,'BZ_RZ',id_out,ier, &
                 ispline=jspline, sm_edge=zsm_edge, id_fun_in=id_BZ, &
                 data=BZx,label='B_Z on (R,Z) grid',units='T')
            iersum=iersum+ier

            if(kpsi.eq.1) then
               call chk_auth('Psi_RZ')
               call xplasma_rzprof(s,'Psi_RZ',id_out,ier, &
                    ispline=jspline, sm_edge=zsm_edge, id_fun_in=id_Psi, &
                    data=Psix, &
                    label='Poloidal flux on (R,Z) grid',units='Wb/rad')
               iersum=iersum+ier
            endif

         endif

         if(zsm_edge.gt.ZERO) then
            if(.not.outside_only) then
               call eqi_rzsmedg(Rg,nR,Zg,nZ,Bphix,zsm_edge,iertmp)
               call eqi_rzsmedg(Rg,nR,Zg,nZ,BRx,zsm_edge,iertmp)
               call eqi_rzsmedg(Rg,nR,Zg,nZ,BZx,zsm_edge,iertmp)
            endif
         endif

         Bmodx=sqrt(BRx**2+BZx**2+Bphix**2)

         if(.not.outside_only) then
            call chk_auth('BR_RZ')
            call xplasma_create_prof(s,'BR_RZ',id_Rg,id_Zg,BRx,id_out,ier, &
                 ispline=jspline,assoc_id=id_BR, &
                 label='B_R on (R,Z) grid',units='T')
            iersum=iersum+ier

            call chk_auth('BZ_RZ')
            call xplasma_create_prof(s,'BZ_RZ',id_Rg,id_Zg,BZx,id_out,ier, &
                 ispline=jspline,assoc_id=id_BZ, &
                 label='B_Z on (R,Z) grid',units='T')
            iersum=iersum+ier

            call chk_auth('BPhi_RZ')
            call xplasma_create_prof(s,'Bphi_RZ',id_Rg,id_Zg,Bphix,id_out,ier,&
                 ispline=jspline,assoc_id=id_Bphi, &
                 label='Btoroidal on (R,Z) grid',units='T')
            iersum=iersum+ier
         endif

         call chk_auth('BMod_RZ')
         call xplasma_create_prof(s,'Bmod_RZ',id_Rg,id_Zg,Bmodx,id_out,ier, &
              ispline=jspline,assoc_id=id_Bmod, &
              label='mod(B) on (R,Z) grid',units='T')
         iersum=iersum+ier

         call xplasma_author_clear(s,xplasma_xmhd,iertmp)
      endif

      deallocate(Rg,Zg,Rtmp,Ztmp,Phidum)
      deallocate(BRx,BZx,Bphix,Bmodx,Psix)

      if(iersum.gt.0) ier=9999

      nullify(eqbuf,idata)

    CONTAINS
      subroutine chk_auth(zname)
        character*(*), intent(in) :: zname

        !  acquire ownership of existing profile, if necessary...

        !----------------------
        integer :: idp,iertmp
        character*32 :: zauth
        integer :: ildbg = 6
        !----------------------

        call xplasma_profId(s,zname,idp)
        if(idp.eq.0) return

        call xplasma_prof_info(s,idp,iertmp, author=zauth)
        if(zauth.ne.xplasma_xmhd) then
#ifdef __DEBUG
           write(ildbg,*) &
                ' %xplasma_brz(chk_auth): resetting author/owner of '// &
                trim(zname)
           write(ildbg,*) '  from "'//trim(zauth)//'" to "'// &
                trim(xplasma_xmhd)//'".'
#endif
           call xplasma_reset_author(s,idp,zauth,xplasma_xmhd,iertmp)
        endif
      end subroutine chk_auth

    end subroutine xplasma_brz

    !-------------------------------------------------------
    subroutine xplasma_rzprof(s,fname,id_out,ier, &
         id_Rgrid,id_Zgrid, &
         ispline,sm_edge, &
         id_fun_in,lamda, &
         label,units,data)

      ! create a profile f(R,Z) from existing data with extrapolation,
      ! or with array data provided.

      ! either an existing function (id_fun_in) or input data (data) or
      ! both must be provided.

      ! modes of use:
      !   (id_fun_in omitted, data(:,:) provided) -- just use the data given
      !   (id_fun_in provided, data(:,:) omitted) -- use existing profile
      !       at id_fun_in to give variation inside plasma; use extrapolation
      !       based on distance map: f_outside = f_bdy*exp(-d/lamda) where
      !       d is the distance from the plasma
      !   (both id_fun_in and data(:,:) provided) -- use existing profile
      !       at id_fun_in to give variation inside plasma: use data to 
      !       specify the variation beyond the plasma.

      !   ispline -- 0 = Bilinear, 1 = C1 Hermite, 2 = C2 Spline

      !   sm_edge -- Meters -- gives smoothing in vicinity of boundary
      !       a hat function convolution of half width sm_edge is applied
      !       in the vicinity of the boundary

      type(xplasma), pointer :: s

      character*(*), intent(in) :: fname  ! name of profile to create...

      integer, intent(out) :: id_out   ! ID of function just created
      integer, intent(out) :: ier      ! completion code 0=OK

      !------

      integer, intent(in), optional :: id_Rgrid,id_Zgrid  ! R & Z grids to use
      !  (default: __RGRID & __ZGRID)

      integer, intent(in), optional :: ispline  ! interpolation order
      ! 0 (default): piecewise linear; 1: C1 Hermite; 2: C2 Spline

      real*8, intent(in), optional :: sm_edge   ! edge smoothing control
      ! default: no smoothing

      integer,intent(in), optional :: id_fun_in ! f(rho) from which to 
      ! form f(R,Z) (default: 0, none)
      ! note: if id_fun_in = -1, use the formula bphi_ccw*g(rho)/R(rho,theta)

      real*8, intent(in), optional :: lamda     ! scrape off distance
      ! default: huge, making for flat extrapolation

      character*(*), intent(in), optional :: label,units  ! labeling info

      real*8, intent(inout), dimension(:,:), optional :: data  ! f(R,Z) data
      ! default: NONE
      ! data must be defined in extrapolated region-- internal region
      ! will be filled in...

      !--------------------------------
      integer :: id_Rg,id_Zg,nR,nZ,i,iZ,jvec
      integer :: jspline,isource,idf,id_dmap,maptype,iertmp,bphi_ccw
      real*8 :: zsm_edge,zlamda
      logical :: standard_RZ,outside_only,need_dmap,axisymm,scrapeoff,iRflag

      real*8, dimension(:,:), allocatable :: zdata
      real*8, dimension(:), allocatable :: zR,zZ,zZtmp,zwk1,zwk2,zwk3,zwk4
      real*8, dimension(:), allocatable :: zdist,zrho,zth
      logical, dimension(:), allocatable :: inside

      integer :: lrhomap,lchimap,ii,jj,kk,id_map
      real*8, dimension(:), pointer :: eqbuf
      integer, dimension(:), pointer :: idata

      real*8, parameter :: lamda_min=0.0001d0
      real*8, parameter :: zlarge = 1.0d35
      real*8, parameter :: ZERO = 0.0d0
      real*8, parameter :: ONE = 1.0d0
      !--------------------------------

      id_out=0

      id_Rg=0; id_Zg=0
      if(present(id_Rgrid)) id_Rg=id_Rgrid
      if(present(id_Zgrid)) id_Zg=id_Zgrid

      call xplasma_ck_rzgrid(s,id_Rg,nR,id_Zg,nZ,standard_RZ,ier)
      if(ier.ne.0) return
               
      call xplasma_global_info(s,ier, axisymm=axisymm,scrapeoff=scrapeoff, &
           bphi_ccw=bphi_ccw)
      if(ier.ne.0) return

      jspline=0
      if(present(ispline)) jspline=ispline

      zsm_edge=ZERO
      if(present(sm_edge)) zsm_edge=sm_edge

      idf=0
      iRflag=.FALSE.
      if(present(id_fun_in)) then
         if(id_fun_in.gt.0) then
            call xplasma_ck_fun(s,id_fun_in,ier)
            if(ier.ne.0) return
            idf=id_fun_in
         else if(id_fun_in.eq.-1) then
            iRflag=.TRUE.  ! sign*g/R
            call xplasma_common_ids(s,ier, id_g=idf)
         endif
      else
         if(present(lamda)) then

            ier=9999
            call xplasma_errmsg_append(s,'  cannot construct f(R,Z): '//trim(fname))
            call xplasma_errmsg_append(s, &
                 ' ?xplasma_rzprof: cannot use scrapeoff length for function')
            call xplasma_errmsg_append(s, &
                 '  extrapolation: no base function specified.')
            return
         endif
      endif

      zlamda=zlarge
      isource=0
      if(present(lamda)) then
         if(zlamda.le.ZERO) then
            zlamda=zlarge
         else
            zlamda=max(lamda_min,lamda)
         endif
         isource=isource+1
      endif

      if(present(data)) then
         need_dmap=.FALSE.
         isource=isource+1
         if((size(data,1).ne.nR).or.(size(data,2).ne.nZ)) then
            ier=510
            call xplasma_errmsg_append(s,'  cannot construct f(R,Z): '//trim(fname))
            call xplasma_errmsg_append(s, &
                 '  ?xplasma_rzprof: input data array dimensions do not match')
            call xplasma_errmsg_append(s, &
                 '   grid sizes of implicitly or explicitly chosen R&Z grids.')
         endif
      else
         need_dmap=.TRUE.
      endif

      if(isource.gt.1) then
         ier=9999
         call xplasma_errmsg_append(s,'  cannot construct f(R,Z): '//trim(fname))
         call xplasma_errmsg_append(s, &
              '  ?xplasma_rzprof: input data array and scrape-off extrapolation (lamda)')
         call xplasma_errmsg_append(s, &
              '   cannot both be present.')
         call xplasma_errmsg_append(s, &
              '   choose one or the other to define region beyond plasma boundary')
         return

      endif

      !----------------------------------------------
      outside_only = need_dmap.and.standard_RZ.and.(zsm_edge.eq.ZERO)

      allocate(zdata(nR,nZ),zR(nR),zZ(nZ),zZtmp(nR),zth(nR*nZ))
      allocate(zwk1(nR*nZ),zwk2(nR*nZ),zwk3(nR*nZ),zwk4(nR*nZ),zrho(nR*nZ))
      allocate(zdist(nR*nZ),inside(nR*nZ))

      call xplasma_grid(s,id_Rg,zR,ier); if (ier.ne.0) return
      call xplasma_grid(s,id_Zg,zZ,ier); if (ier.ne.0) return

      if(idf.eq.0) then

         ! just use the data provided to create a profile...

         zdata = data

      else

         if(standard_RZ) then
            maptype=3

            !  make sure rhomap is available

            call xplasma_ctrans(sp,iertmp, &
                 R_in=zR(nR/2), Z_in=zZ(nZ/2), &
                 rho_out=zrho(1), theta_out=zth(1), maptype=3)

            call xplasma_find_item(s,'__FASTMAP',id_map,iertmp)
            if(iertmp.ne.0) then
               call xplasma_errmsg_append(s, &
                ' %xplasma_rzprof: __FASTMAP not found: no scrapeoff region?')
               maptype=2
            else
               call xplasma_blackbox_retrieve(s, id_map, iertmp, &
                    ia_ptr=idata, r8a_ptr=eqbuf)
               lrhomap=idata(1)
               lchimap=idata(2)
               jj=lrhomap-nR-1
               kk=lchimap-nR-1
            endif
         else
            maptype=2
         endif
      
         jvec=0
         ii=-nR
         do iZ=1,nZ
            ii=ii+nR

            if(maptype.eq.2) then
               zZtmp=zZ(iZ)
               call xplasma_ctrans(s,.TRUE.,ier, R_in=zR,Z_in=zZtmp, &
                    rho_out=zrho(ii+1:ii+nR),theta_out=zth(ii+1:ii+nR), &
                    maptype=maptype)
               if(ier.ne.0) exit
            else
               jj=jj+nR
               kk=kk+nR
               zrho(ii+1:ii+nR)=eqbuf(jj+1:jj+nR)
               zth(ii+1:ii+nR)=eqbuf(kk+1:kk+nR)
            endif

            if(need_dmap) then
               ! using lamda extrapolation for pts beyond plasma bdy

               do i=1,nR
                  if(zrho(ii+i).gt.ONE) then
                     jvec=jvec+1
                     zwk1(jvec)=zZ(iZ)
                     zwk2(jvec)=zR(i)
                     inside(ii+i)=.FALSE.
                     zrho(ii+i)=ONE
                  else
                     inside(ii+i)=.TRUE.
                     zdist(ii+i)=ZERO
                  endif
               enddo
            else
               do i=1,nR
                  if(zrho(ii+i).le.ONE) then
                     jvec=jvec+1
                     zwk1(jvec)=zR(i)
                     zwk2(jvec)=zrho(ii+i)
                     zwk3(jvec)=zth(ii+i)
                     inside(ii+i)=.TRUE.
                  else
                     inside(ii+i)=.FALSE.
                  endif
               enddo

            endif
         enddo

         if(need_dmap) then
            call xplasma_bdfind(s,jvec,zwk2(1:jvec),zwk1(1:jvec),ier, &
                 maptype=maptype,outside_only=outside_only, &
                 theta_out=zwk3(1:jvec), dist=zwk4(1:jvec))
         else

            ! use idf & iRflag

            call xplasma_eval_prof(s,idf, &
                 xplasma_rho_coord,zwk2(1:jvec), &
                 xplasma_theta_coord,zwk3(1:jvec), &
                 zwk4(1:jvec),ier)

            if(iRflag) then
               zwk4(1:jvec)=bphi_ccw*zwk4(1:jvec)/zwk1(1:jvec)
            endif

         endif

         if(ier.eq.0) then
            if(need_dmap) then
               jvec=0
               ii=-nR
               do iZ=1,nZ
                  ii=ii+nR
                  do i=1,nR
                     if(.not.inside(ii+i)) then
                        jvec=jvec+1
                        zth(ii+i)=zwk3(jvec)
                        zdist(ii+i)=zwk4(jvec)
                     endif
                  enddo
               enddo

               call xplasma_eval_prof(s,idf, &
                    xplasma_rho_coord,zrho, xplasma_theta_coord,zth, &
                    zwk2,ier)

               if(ier.eq.0) then
                  ii=-nR
                  do iZ=1,nZ
                     ii=ii+nR
                     do i=1,nR
                        zdata(i,iZ)=zwk2(ii+i)
                        if((zlamda.lt.zlarge).and.(zdist(ii+i).gt.ZERO)) then
                           zdata(i,iZ)=zdata(i,iZ)*exp(-zdist(ii+i)/zlamda)
                        endif
                     enddo
                  enddo
               endif

            else
               ! using data for pts beyond plasma bdy

               jvec=0
               ii=-nR
               do iZ=1,nZ
                  ii=ii+nR
                  do i=1,nR
                     if(inside(ii+i)) then
                        jvec=jvec+1
                        zdata(i,iZ)=zwk4(jvec)  ! interior-- use function eval
                        data(i,iZ)=zwk4(jvec)   ! pass back also...
                     else
                        zdata(i,iZ)=data(i,iZ)  ! exterior-- use passed data
                     endif
                  enddo
               enddo

            endif

         endif
      endif

      if(ier.eq.0) then

         if(zsm_edge.gt.ZERO) then
            if(idf.eq.0) then
               sp => s
               call eqi_rzsmedg(zR,nR,zZ,nZ,zdata,zsm_edge,iertmp)
            else
               ! use module-internal blending routine; mesh with core data
               call xpblend(s,zR,nR,zZ,nZ,zdata,zsm_edge,idf,iRflag,bphi_ccw, &
                    maptype, iertmp)
            endif
         endif

         call xplasma_create_2dprof(s,fname,id_Rg,id_Zg,zdata,id_out,ier, &
              ispline=jspline,assoc_id=idf,label=label,units=units)

      endif

      deallocate(inside,zdata,zR,zZ,zZtmp,zwk1,zwk2,zwk3,zwk4,zrho,zth,zdist)

      nullify(eqbuf,idata)

    end subroutine xplasma_rzprof

    !-------------------------------------------------------
    subroutine xplasma_rzprof_fun(s,fname,user_fcn,user_iarg,id_out,ier, &
         id_Rgrid,id_Zgrid,ispline,sm_edge,id_fun_in,no_extrap, &
         label,units)

      ! create a profile f(R,Z) from a user provided callable function
      ! possibly combined with existing data.

      ! the callable function covers the entire region inside the limiters
      ! if id_fun_in is omitted; it covers only the region outside the plasma
      ! and inside the limiters if id_fun_in is included.

      ! modes of use:
      !   (id_fun_in omitted) -- just use the given callable function.
      !   (id_fun_in provided)-- use existing profile
      !       at id_fun_in to give variation inside plasma; use the given
      !       function to define the variation beyond the plasma boundary.

      ! the region beyond the limiters is defined by a simple extrapolation.
      ! unless no_extrap is set, in which case the function is used even
      ! beyond the limiters.

      !   sm_edge -- Meters -- gives smoothing in vicinity of boundary
      !       a hat function convolution of half width sm_edge is applied
      !       in the vicinity of the boundary

      type(xplasma), pointer :: s

      character*(*), intent(in) :: fname  ! name of profile to create...

      real*8, external :: user_fcn        ! user provided function...
      integer, intent(in) :: user_iarg    ! argument for user_fcn:

      !  real*8 function user_fcn(user_iarg,R,Z,Phi,ier)

      integer, intent(out) :: id_out   ! ID of function just created
      integer, intent(out) :: ier      ! completion code 0=OK

      !------

      integer, intent(in), optional :: id_Rgrid,id_Zgrid  ! R & Z grids to use
      !  (default: __RGRID & __ZGRID)

      integer, intent(in), optional :: ispline  ! interpolation order
      ! 0 (default): piecewise linear; 1: C1 Hermite; 2: C2 Spline

      real*8, intent(in), optional :: sm_edge   ! edge smoothing control
      ! default: no smoothing

      integer,intent(in), optional :: id_fun_in ! f(rho) from which to 
      ! form f(R,Z) (default: 0, none)
      ! note: if id_fun_in = -1, use the formula bphi_ccw*g(rho)/R(rho,theta)

      logical, intent(in), optional :: no_extrap ! T to suppress extrapolation
      ! beyond limiter.  In this case user_fcn calls are used for (R,Z) even
      ! for points outside the limiter.

      character*(*), intent(in), optional :: label,units  ! labeling info

      !--------------------------------
      integer :: id_Rg,id_Zg,nR,nZ,i,iZ,jvec
      integer :: jspline,idf,maptype,iertmp,bphi_ccw
      real*8 :: zsm_edge,zphi0
      real*8, dimension(:,:), allocatable :: zdata
      real*8, dimension(:), allocatable :: zR,zZ,zwkZ,zrho,zth
      real*8, dimension(:), allocatable :: zdist_lim,zfeval,zwkR
      logical, dimension(:), allocatable :: inside
      logical, dimension(:,:), allocatable :: inside_lim
      logical :: standard_RZ,extrap,axisymm,scrapeoff,iRflag

      real*8, parameter :: ZERO = 0.0d0
      real*8, parameter :: ONE = 1.0d0
      !--------------------------------
      !  dmc: this routine might be sped up by making calls covering
      !   all RxZ instead of just one row at a time.  Other routines
      !   have been optimized in this way... but not yet this one.

      id_out=0

      id_Rg=0; id_Zg=0
      if(present(id_Rgrid)) id_Rg=id_Rgrid
      if(present(id_Zgrid)) id_Zg=id_Zgrid

      call xplasma_ck_rzgrid(s,id_Rg,nR,id_Zg,nZ,standard_RZ,ier)
      if(ier.ne.0) return
               
      call xplasma_global_info(s,ier, axisymm=axisymm,scrapeoff=scrapeoff, &
           bphi_ccw=bphi_ccw)
      if(ier.ne.0) return

      jspline=0
      if(present(ispline)) jspline=ispline

      zsm_edge=0
      if(present(sm_edge)) zsm_edge=sm_edge

      idf=0
      iRflag=.FALSE.
      if(present(id_fun_in)) then
         if(id_fun_in.gt.0) then
            call xplasma_ck_fun(s,id_fun_in,ier)
            if(ier.ne.0) return
            idf=id_fun_in
         else if(id_fun_in.eq.-1) then
            iRflag=.TRUE.  ! sign*g/R
            call xplasma_common_ids(s,ier, id_g=idf)
         endif
      endif

      extrap=.TRUE.
      if(present(no_extrap)) extrap = .not.no_extrap

      !------------------------------------------
      allocate(zdata(nR,nZ),zR(nR),zZ(nZ),zwkZ(nR))
      allocate(zrho(nR),zth(nR),zdist_lim(nR),zfeval(nR),zwkR(nR))
      allocate(inside(nR),inside_lim(nR,nZ))

      call xplasma_grid(s,id_Rg,zR,ier); if (ier.ne.0) return
      call xplasma_grid(s,id_Zg,zZ,ier); if (ier.ne.0) return

      if(standard_RZ) then
         maptype=3
      else
         maptype=2
      endif

      zphi0 = ZERO

      do iZ=1,nZ
         zwkZ=zZ(iZ)

         if(idf.eq.0) then
            if(extrap) then
               call xplasma_lim_distance(s,zR,zwkZ,zdist_lim,ier, &
                    maptype=maptype)
               if(ier.ne.0) exit
            else
               zdist_lim=ZERO
            endif

            do i=1,nR
               if(zdist_lim(i).gt.ZERO) then
                  zdata(i,iZ)=ZERO
                  inside_lim(i,iZ)=.FALSE.
               else
                  zdata(i,iZ) = user_fcn(user_iarg,zR(i),zwkZ(1),ZERO,ier)
                  if(ier.ne.0) exit
                  inside_lim(i,iZ)=.TRUE.
               endif
            enddo
            if(ier.ne.0) exit

         else

            call xplasma_ctrans(s,.TRUE.,ier, R_in=zR,Z_in=zwkZ, &
                 rho_out=zrho,theta_out=zth, maptype=maptype)
            if(ier.ne.0) exit

            if(extrap) then
               call xplasma_lim_distance(s,zR,zwkZ,zdist_lim,ier, &
                    maptype=maptype)
               if(ier.ne.0) exit
            else
               zdist_lim=ZERO
            endif

            jvec=0
            do i=1,nR
               inside_lim(i,iZ)=.TRUE.
               if(zdist_lim(i).gt.ZERO) then
                  zdata(i,iZ)=ZERO
                  inside(i)=.FALSE.
                  inside_lim(i,iZ)=.FALSE.
               else if(zrho(i).gt.ONE) then
                  zdata(i,iZ) = user_fcn(user_iarg,zR(i),zwkZ(1),ZERO,ier)
                  if(ier.ne.0) exit
                  inside(i)=.FALSE.
               else
                  jvec=jvec+1
                  inside(i)=.TRUE.
                  zwkR(jvec)=zR(i)
                  zrho(jvec)=zrho(i)
                  zth(jvec)=zth(i)
               endif
            enddo
            if(ier.ne.0) exit

            if(jvec.gt.0) then
               call xplasma_eval_prof(s,idf, &
                    xplasma_rho_coord,zrho(1:jvec), &
                    xplasma_theta_coord,zth(1:jvec), &
                    zfeval(1:jvec),ier)
               if(ier.ne.0) exit

               if(iRflag) then
                  zfeval(1:jvec) = bphi_ccw*zfeval(1:jvec)/zwkR(1:jvec)
               endif

               jvec=0
               do i=1,nR
                  if(inside(i)) then
                     jvec=jvec+1
                     zdata(i,iZ)=zfeval(jvec)
                  endif
               enddo
            endif
         endif

      enddo

      if(ier.eq.0) then

         sp => s

         if(extrap) then
            ! compute simple extrapolation for beyond-limiter data points
            call eqi_rzfx(zR,nR,zZ,nZ,inside_lim,zdata)
         endif

         if(zsm_edge.gt.ZERO) then
            ! optional smoothing...
            if(idf.eq.0) then
               call eqi_rzsmedg(zR,nR,zZ,nZ,zdata,zsm_edge,iertmp)
            else
               ! use module-internal blending routine; mesh with core data
               call xpblend(s,zR,nR,zZ,nZ,zdata,zsm_edge,idf,iRflag,bphi_ccw, &
                    maptype, iertmp)
            endif
         endif

         call xplasma_create_2dprof(s,fname,id_Rg,id_Zg,zdata,id_out,ier, &
              ispline=jspline,assoc_id=idf,label=label,units=units)

      endif

      deallocate(zdata,zR,zZ,zwkZ,zrho,zth,zdist_lim,zfeval,zwkR)
      deallocate(inside,inside_lim)

    end subroutine xplasma_rzprof_fun

    !-------------------------------------------------------
    subroutine xplasma_ck_fun(s,idf,ier)

      !  ** private **
      !  verify function id -- must be fcn of (rho,theta)

      type(xplasma), pointer :: s

      integer, intent(in) :: idf
      integer, intent(out) :: ier

      !------------------
      integer :: id_grid1,id_grid2,icoord1,icoord2
      character*32 zname1,zname2
      !------------------

      do
         call xplasma_prof_info(s,idf,ier, gridId1=id_grid1, gridId2=id_grid2)
         if(ier.ne.0) exit

         call xplasma_grid_info(s,id_grid1,ier, coord=icoord1, name=zname1)
         if(ier.ne.0) exit

         if(id_grid2.gt.0) then
            call xplasma_grid_info(s,id_grid2,ier, coord=icoord2, name=zname2)
            if(ier.ne.0) exit

         else
            icoord2=0
         endif

         if((icoord1.eq.xplasma_rho_coord).or. &
              (icoord1.eq.xplasma_theta_coord)) then

            continue

         else

            ier=9999
            call xplasma_errmsg_append(s, &
                 '   passed profile needs to be defined vs. flux coordinates.')
            call xplasma_errmsg_append(s, &
                 '   instead it is a function of: '//trim(zname1))
            exit
         endif

         if((icoord2.eq.0).or.(icoord2.eq.xplasma_rho_coord).or. &
              (icoord2.eq.xplasma_theta_coord)) then

            continue

         else

            ier=9999
            call xplasma_errmsg_append(s, &
                 '   passed profile needs to be defined vs. flux coordinates.')
            call xplasma_errmsg_append(s, &
                 '   instead it is a function of: '//trim(zname2))
            exit
         endif

         exit
      enddo

      if(ier.ne.0) then
         call xplasma_errmsg_append(s,' ...error in xplasma_profs module.')
      endif

    end subroutine xplasma_ck_fun

    !-------------------------------------------------------
    subroutine xplasma_ck_rzgrid(s,id_Rgrid,nR,id_Zgrid,nZ,standard_RZ,ier)

      !  ** private ** 
      !  summary info on R & Z grids

      type(xplasma), pointer :: s

      integer, intent(inout) :: id_Rgrid,id_Zgrid ! grid IDs: use default if 0
      integer, intent(out) :: nR,nZ               ! grid sizes
      logical, intent(out) :: standard_RZ         ! T if standard grid is used
      integer, intent(out) :: ier                 ! error code, 0=OK

      !-------------------------
      integer :: icoord,isize,id_Rgrid0,id_Zgrid0
      character*32 :: zname
      !-------------------------

      standard_RZ=.FALSE.
      nR=0
      nZ=0

      call xplasma_find_item(s,'__Rgrid',id_Rgrid0,ier)
      if(ier.ne.0) return

      call xplasma_find_item(s,'__Zgrid',id_Zgrid0,ier)
      if(ier.ne.0) return

      if(id_Rgrid.eq.0) then
         id_Rgrid=id_Rgrid0
      endif

      if(id_Zgrid.eq.0) then
         id_Zgrid=id_Zgrid0
      endif

      if((id_Rgrid.eq.id_Rgrid0).and.(id_Zgrid.eq.id_Zgrid0)) then
         standard_RZ=.TRUE.
      endif

      call xplasma_grid_info(s,id_Rgrid,ier, coord=icoord, size=nR, name=zname)
      if(ier.eq.0) then
         if(icoord.ne.xplasma_R_coord) then
            ier=701
            call xplasma_errmsg_append(s,' passed id grid name: '//trim(zname))
         endif
      endif
      if(ier.ne.0) return

      call xplasma_grid_info(s,id_Zgrid,ier, coord=icoord, size=nZ, name=zname)
      if(ier.eq.0) then
         if(icoord.ne.xplasma_Z_coord) then
            ier=701
            call xplasma_errmsg_append(s,' passed id grid name: '//trim(zname))
         endif
      endif
      if(ier.ne.0) return

    end subroutine xplasma_ck_rzgrid

    !=====================================

    subroutine xpblend(s,zR,nR,zZ,nZ,zdata,zsm,idf,iRflag,ibccw,maptype,ier)

      !  ** PRIVATE **

      !  smooth zdata(R,Z) beyond edge region by matching point and slope
      !  of interior function (idf), in extrapolated (rho,theta) space and
      !  then mapping back to (R,Z).

      !  if iRflag is set, expect idf=id_g and use bphi_ccw*f(rho)/R(rho,theta)
      !      (i.e. the toroidal field) for the interior function.

      type (xplasma), pointer :: s

      integer :: nR,nZ        ! grid sizes
      real*8 :: zR(nR),zZ(nZ) ! grids
      real*8 :: zdata(nR,nZ)  ! data

      real*8 :: zsm           ! smoothing distance

      integer :: idf          ! interior matching function ID

      logical :: iRflag       ! T to use ibccw*f(rho)/R; idf points to f(rho)
      integer :: ibccw        ! +/- 1 for use if (iRflag) true

      integer :: maptype      ! mapping selector for (R,Z) -> (rho,theta)

      integer, intent(out) :: ier    ! status code; 0=OK

      !-----------------------------------------------------------
      integer :: ii,jvec,jvec2,iR,iZ
      real*8 :: del,del1,del2,dR1,dZ1
      real*8 :: dat0,dat1,dat2,datx,b,denom,dxtrap,xx,ff

      real*8, parameter :: ZERO = 0.0d0
      real*8, parameter :: HALF = 0.5d0
      real*8, parameter :: ONE  = 1.0d0
      real*8, parameter :: TWO  = 2.0d0
      real*8, parameter :: THREE  = 3.0d0

      ! data for 1 row @ fixed Z:
      real*8, dimension(:), allocatable :: zzwk,zrho,zth,zRb,zZb,z1vec

      ! collected data for edge blending:
      integer :: isizx,isizf
      real*8, dimension(:), allocatable :: delx,rhox,thx,zRx,datbx, &
           Rinx,Zinx,rho_inx,th_inx,datinx
      integer, dimension(:), allocatable :: iRsave,iZsave
      !-----------------------------------------------------------

      del2 = min( (zR(nR)-zR(1))/max(40,nR), (zZ(nZ)-zZ(1))/max(40,nZ) )
      del1 = HALF*del2

      ier=0

      allocate(zzwk(nR),zrho(nR),zth(nR),zRb(nR),zZb(nR),z1vec(nR))
      z1vec = ONE      ! vector of length nR of 1s...

      isizf=20
10    continue
      if(isizf.eq.20) then
         isizf=10
      else if(isizf.eq.10) then
         isizf=5
      else 
         isizf=isizf-1
         if(isizf.eq.0) then
            call xplasma_errmsg_append(s,' ?xpblend: too many points in range for smoothed blending.')
            ier=9999
            return
         endif
      endif

      isizx = (nR*nZ)/isizf
      if(allocated(delx)) then
         deallocate(delx,rhox,thx,datbx,iRsave,iZsave)
         deallocate(Rinx,Zinx,rho_inx,th_inx,datinx)
      endif
      allocate(delx(isizx),rhox(isizx),thx(isizx),zRx(isizx))
      allocate(datbx(isizx),iRsave(isizx),iZsave(isizx))
      allocate(Rinx(2*isizx),Zinx(2*isizx),rho_inx(2*isizx),th_inx(2*isizx))
      allocate(datinx(2*isizx))

      jvec = 0
      jvec2 = 0

      do iZ=1,nZ
         zzwk = zZ(iZ) ! expand this Z to vector of length nR...

         call xplasma_ctrans(s,.TRUE.,ier, R_in=zR, Z_in=zzwk, &
              rho_out=zrho, theta_out=zth, maptype=maptype)
         if(ier.ne.0) exit

         call xplasma_ctrans(s,.TRUE.,ier, rho_in=z1vec, theta_in=zth, &
              R_out=zRb, Z_out=zZb)
         if(ier.ne.0) exit

         do iR=1,nR
            if(zrho(iR).gt.ONE) then
               del = sqrt((zR(iR)-zRb(iR))**2 + (zZ(iZ)-zZb(iR))**2)
               dR1 = (zRb(iR)-zR(iR))/del
               dZ1 = (zZb(iR)-zZ(iZ))/del
               if((ZERO.lt.del).AND.(del.lt.zsm)) then
                  ! point in range for smoothing

                  jvec = jvec + 1
                  if(jvec.gt.isizx) go to 10  ! check bound

                  delx(jvec)=del
                  iRsave(jvec)=iR
                  iZsave(jvec)=iZ
                  rhox(jvec)=ONE
                  thx(jvec)=zth(iR)
                  zRx(jvec)=zRb(iR)

                  Rinx(jvec2+1)=zRb(iR) + del1*dR1
                  Rinx(jvec2+2)=zRb(iR) + del2*dR1

                  Zinx(jvec2+1)=zZb(iR) + del1*dZ1
                  Zinx(jvec2+2)=zZb(iR) + del2*dZ1

                  jvec2 = jvec2 + 2
               endif
            endif
         enddo

      enddo
      if(ier.ne.0) then
         call xplasma_errmsg_append(s,' ?error 1 in xpblend internal routine')
         return
      endif

      ! evaluate data at boundary

      call xplasma_eval_prof(s,idf, &
                 xplasma_rho_coord,rhox(1:jvec), &
                 xplasma_theta_coord,thx(1:jvec), &
                 datbx(1:jvec),ier)
      if(ier.ne.0) then
         call xplasma_errmsg_append(s,' ?error 2 in xpblend internal routine')
         return
      endif

      if(iRflag) then
         datbx(1:jvec) = ibccw*datbx(1:jvec)/zRx(1:jvec)
      endif

      ! evaluate interior points to allow matching assessment
      ! 1st convert (R,Z) to (rho,theta); use maptype=2 for reasonable accuracy

      call xplasma_ctrans(s,.TRUE.,ier, &
           R_in=Rinx(1:jvec2),Z_in=Zinx(1:jvec2), &
           rho_out=rho_inx(1:jvec2),theta_out=th_inx(1:jvec2), &
           maptype=2)
      if(ier.ne.0) then
         call xplasma_errmsg_append(s,' ?error 3 in xpblend internal routine')
         return
      endif

      ! evaluate interior data 

      call xplasma_eval_prof(s,idf, &
                 xplasma_rho_coord,rho_inx(1:jvec2), &
                 xplasma_theta_coord,th_inx(1:jvec2), &
                 datinx(1:jvec2),ier)
      if(ier.ne.0) then
         call xplasma_errmsg_append(s,' ?error 2 in xpblend internal routine')
         return
      endif

      if(iRflag) then
         datinx(1:jvec2) = ibccw*datinx(1:jvec2)/Rinx(1:jvec2)
      endif

      ! find linear extrapolation along fixed exterior theta lines
      ! i.e. a profile with continuous value and 1st derivative; 2nd
      ! derivative definitely not continuous.

      denom = del1*del2*(del2-del1)

      jvec2=0
      do ii=1,jvec
         dat0=datbx(ii)
         dat1=datinx(jvec2+1)
         dat2=datinx(jvec2+2)
         jvec2 = jvec2 + 2

         !  3 pts: (0,dat0), (del1,dat1), (del2,dat2) -- parabolic fit
         !          bdy       1st int.     2nd int.
         !  d(x) = A*x**2 + B*x + C    d(0)=C=dat0; d'(0)=B
         !    denom = del1*del2*(del2-del1)
         !        A = -[(dat1-dat0)*del2 + (dat2-dat0)*del1]/denom
         !        B =  [(dat1-dat0)*del2**2 + (dat2-dat0)*del1**2]/denom

         B = ((dat1-dat0)*del2**2 - (dat2-dat0)*del1**2)/denom

         del=delx(ii)
         xx=del/zsm             ! normalized extrapolation distance

         dxtrap = dat0 - B*del  ! linear extrapolation of interior data

         iR = iRsave(ii)
         iZ = iZsave(ii)

         ! blend using Hermite cubic f*data + (1-f)*dxtrap
         ! f(0)=0, f'(0)=0, f(1)=1, f'(1)=0
         !     => f(x) = x**2*(-2*x + 3)

         ff = xx*xx*(THREE-TWO*xx)

         zdata(iR,iZ) = ff*zdata(iR,iZ) + (ONE-ff)*dxtrap
      enddo

    end subroutine xpblend

    !------------------------------------------
    subroutine xplasma_irhofun(s,id_axis,zlbl,inprof,zprof,iflag,id,ierr, &
         smdelx)

      !  make XPLASMA profile -- integrated quantity; smooth by 1/2 zone width
      !  to insure smooth derivative from spline

      type (xplasma), pointer :: s

      integer, intent(in) :: id_axis    ! axis id -- must be rho or akin to rho

      character*(*), intent(in) :: zlbl ! name for profile function to create

      integer, intent(in) :: inprof     ! size of the integrated data profile

      real*8, intent(in) :: zprof(inprof) ! the integrated data provided
              ! if inprof = size(x_axis) zprof(1)=0 is expected
              ! if inprof = size(x_axis)-1 the axial data point
              !    is presumed to be omitted.

      integer, intent(in) :: iflag      ! =1 -- volume normalization;
              ! derivative evaluations -> W/m^3, #/sec/m^3, etc.
              ! =2 -- area normalization;
              ! derivative evaluations -> A/m^2 (current density).

              ! >2 -- ID of normalization profile, should be akin to vol(rho)
              !       or area(rho)

      integer, intent(out) :: id        ! id of stored profile (if successful)
      integer, intent(out) :: ierr      ! completion code, 0=OK.

      real*8, intent(in), optional :: smdelx   ! smoothing width (option)

      !  a single zone width hat function smooth (effectively, half a zone
      !  width) is the minimum used by eqi_irhofun.  For more smoothing
      !  than this, give smdelx = [a number greater than x(2)-x(1)] where
      !  x(1:2) are the 1st two grid points identified via id_axis.

      !  Note, the actual smoothing width is proportional to the grid
      !  resolution and so could vary across the grid; a fixed multiplier
      !  of smdelx/(x(2)-x(1)) would be applied in this case.  In the
      !  evaluation routine (eqi_irhofun) there is an upper limit imposed
      !  which corresponds a multiplier that yields (1/4) the profile width
      !  at the grid spacing at rho=0.

      !  if an error occurs and ierr is set, id=0 will be returned.
      !-----------------------------------------------------------------
      integer :: ierr_save
      real*8 :: zsm
      !-----------------------------------------------------------------

      sp => s

      zsm = 0.0d0
      if(present(smdelx)) zsm = smdelx

      call eqi_irhofun(id_axis,zlbl,inprof,zprof,iflag,zsm,id,ierr)

    end subroutine xplasma_irhofun

    !-------------------------------------------------------
    subroutine xplasma_geqdsk_rewrite(filename,ier)

      !  rewrite g-eqdsk file e.g. under new filename, after reading
      !  (e.g. inside an eqi_fromgeqdsk call).

      !  WARNING: the information written is based on the most recent
      !  read of g-eqdsk information.  Information to reconstruct the
      !  g-eqdsk from the data as originally read is NOT saved with each
      !  xplasma object.  Hence: no s pointer argument.

      character*(*), intent(in) :: filename
      integer, intent(out) :: ier

      call eqm_wr_geqdsk(filename,ier)
      if(ier.ne.0) ier=9999

    end subroutine xplasma_geqdsk_rewrite

    !-------------------------------------------------------
    subroutine xplasma_wr_geqdsk(s,ierr, &
         lun_geqdsk, filename, label, Rmin, Rmax, Zmin, Zmax, &
         cur, id_qprof, id_pprof, id_psi_in, &
         nh, nv, nbdy)

      !  Build and write a G-eqdsk file (disk file) from current xplasma
      !  contents.  This differs from "xplasma_geqdsk_rewrite" as the latter
      !  just echoes data read in from another G-eqdsk file or MDSplus data
      !  structure.  This actual constructs the information from current
      !  xplasma contents -- interpolation is involved.

      !  The xplasma object must contain a complete free boundary or 
      !  extrapolated equilibrium, so that Psi(R,Z) covering a rectangle
      !  enclosing the plasma is defined.

      !-----------------
      !  required:

      type (xplasma), pointer :: s  ! XPLASMA object containing equilibrium

      integer, intent(out) :: ierr  ! status code on exit: 0=OK

      !-----------------
      !  optional:

      integer, intent(in), optional :: lun_geqdsk  
      ! FORTRAN LUN on which to write file.  If omitted, the LUN used for
      ! reading G-eqdsk files, in the geqdsk_mds library, is used.
      !   (call geq_getlun(ilun)) (geqdsk_mod default value: 77 as of 7/2006).

      character*(*), intent(in), optional :: filename
      ! default " "; if non-blank, it is the
      ! name of the file to write.  If blank or omitted, the code simply 
      ! writes to lun_geqdsk; it would be up to the caller to OPEN a file.

      character*(*), intent(in), optional :: label
      ! default " "; default means: use xplasma global label; 
      ! if non-blank, the 1st 48 characters are used as a label string in the
      ! G-eqdsk file being written.

      real*8, intent(in), optional :: Rmin,Rmax, Zmin,Zmax
      ! [R,Z] domain over which Psi(R,Z) is written in the G-eqdsk file.
      ! default: use the [R,Z] grid limits already stored in xplasma.  If
      ! explicit limits are provided, overriding the defaults, these must
      ! not extend beyond the grid limits.

      real*8, intent(in), optional :: cur
      ! total plasma current.  Default: use value implied by equilibrium data.

      integer, intent(in), optional :: id_qprof
      ! ID of profile f(rho) defining q(rho).  Default: use value derived from
      ! equilibrium-- calculated here, it will be named Q_EQDSK

      integer, intent(in), optional :: id_pprof
      ! ID of profile f(rho) defining equilibrium pressure.  Default: use value
      ! previously tagged as equilibrium pressure.  If this is not available,
      ! the code will attempt to compute a pressure using the surface averaged
      ! GS equation JxB=grad(P).  This could lead to a non-physical result if
      ! there are errors in the equilibrium data already provided, or if the
      ! assumption of a scalar P is inappropriate.  If a pressure profile is
      ! calculated, it will be named P_EQDSK-- done using an xplasma_gen_p
      ! call.

      integer, intent(in), optional :: id_psi_in
      ! ID of associated pair of Psi(rho) or Psi(R,Z) profiles
      ! If defaulted, the code looks for standard names w/in xplasma object.

      integer, intent(in), optional :: nh,nv
      ! number of horizontal and vertical grid points, respectively.  If
      ! defaulted, the sizes of the [R,Z] grids are used.  NOTE that nh also
      ! controls the size of the 1d profiles f,ff',P,P',and q written in the
      ! G-eqdsk profiles.  These 1d profiles are written over an implied 
      ! evenly spaced Psi grid going from Psi(min) at the axis to Psi(max)
      ! at the boundary.
      
      integer, intent(in), optional :: nbdy
      ! number of points to use to described plasma boundary and limiter.
      ! default: 200.

      !-------------------------------------
      integer :: ilun,istat,ifile,ilen,id_Rc,id_Zc,idp,idq,ii,idpsi
      integer :: inh,inv,inbdy,inum
      integer :: idwk,iertmp

      character*48 glabel
      real*8 :: zRmin,zRmax,zZmin,zZmax,zcur

      !  for plasma current estimate:
      integer, parameter :: inxi=21
      real*8 :: zxi(inxi),zii(inxi)

      integer :: inbdy_dflt=201
      !-------------------------------------

      if(present(lun_geqdsk)) then
         ilun = lun_geqdsk
      else
         call eqi_geq_getlun(ilun)
      endif

      ierr = 0
      ifile= 0
      if(present(filename)) then
         open(unit=ilun,file=trim(filename),status='unknown',iostat=istat)
         if(istat.ne.0) then
            ierr = 9999
            call xplasma_errmsg_append(s, &
                 ' ?xplasma_profs(xplasma_wr_geqdsk): file open failure: '// &
                 trim(filename))
         else
            ifile=1
         endif
      endif
      if(ierr.ne.0) return

      do

         ! label...

         if(present(label)) then
            glabel=label
         else
            glabel='xplasma_wr_geqdsk:'
            ilen=len(trim(glabel))
            call xplasma_global_info(s,ierr,initLabel=glabel(ilen+1:))
         endif
         if(ierr.ne.0) exit
            
         ! Rmin&max, Zmin&max...

         call xplasma_find_item(s,'__R_coord',id_Rc,ierr); if(ierr.ne.0) exit
         call xplasma_find_item(s,'__Z_coord',id_Zc,ierr); if(ierr.ne.0) exit

         call xplasma_coord_info(s,id_Rc,ierr, ngrids=inum); if(ierr.ne.0) exit
         if(inum.eq.0) then
            ierr=615  ! need R & Z rectangle; Psi(R,Z)
            exit
         endif

         call xplasma_coord_info(s,id_Zc,ierr, ngrids=inum); if(ierr.ne.0) exit
         if(inum.eq.0) then
            ierr=615  ! need R & Z rectangle; Psi(R,Z)
            exit
         endif

         call xplasma_coord_info(s,id_Rc,ierr, xmin=zRmin, xmax=zRmax)
         if(ierr.ne.0) exit

         call xplasma_coord_gridId(s,id_Rc,1,idwk,ierr)
         if(ierr.ne.0) exit

         call xplasma_grid_size(s,idwk,inh,ierr)
         if(ierr.ne.0) exit

         call xplasma_coord_info(s,id_Zc,ierr, xmin=zZmin, xmax=zZmax)
         if(ierr.ne.0) exit

         call xplasma_coord_gridId(s,id_Zc,1,idwk,ierr)
         if(ierr.ne.0) exit

         call xplasma_grid_size(s,idwk,inv,ierr)
         if(ierr.ne.0) exit

         if(present(Rmin)) then
            zRmin=max(zRmin,Rmin)
         endif

         if(present(Rmax)) then
            zRmax=min(zRmax,Rmax)
         endif

         if(present(Zmin)) then
            zZmin=max(zZmin,Zmin)
         endif

         if(present(Zmax)) then
            zZmax=min(zZmax,Zmax)
         endif

         !  total current 

         if(present(cur)) then
            zcur=cur
         else
            !  estimate the current

            call xplasma_author_set(s,'xplasma_profs',iertmp)
            do ii=1,inxi
               zxi(ii)=(ii-1)*1.0d0/(inxi-1)
            enddo

            call xplasma_create_integ(s,'__xpasma_profs_tmp_integ', &
                 zxi,idwk,iertmp, cache_enable=.FALSE.)
                 
            call xplasma_rho_zonint(s,idwk,'Itor',zii,ierr)
            if(ierr.ne.0) exit

            zcur=zii(inxi)

            call xplasma_remove_item(s,idwk,iertmp)

            call xplasma_author_clear(s,'xplasma_profs',iertmp)
         endif

         !  P & q

         if(present(id_qprof)) then
            idq=id_qprof
         else
            call xplasma_gen_q(s,'Q_EQDSK',2,idq,ierr)
            if(ierr.ne.0) exit
         endif

         if(present(id_pprof)) then
            idp=id_pprof
         else
            call xplasma_common_ids(s,ierr,id_p=idp)
            if(ierr.ne.0) exit
            if(idp.eq.0) then
               call xplasma_gen_p(s,'P_EQDSK',2,idp,ierr)
               if(ierr.ne.0) exit
            endif
         endif

         if((idp.eq.0).or.(idq.eq.0)) then
            call xplasma_errmsg_append(s, &
                 ' xplasma_profs: zero ids for Pressure or q profiles')
            ierr=9999
            exit
         endif

         !  Psi (id can be zero)

         idpsi=0
         if(present(id_psi_in)) then
            idpsi = id_psi_in
         endif

         !  output grid sizes (default is from R & Z coordinates' 1st grids).

         if(present(nh)) then
            inh=nh
         endif

         if(present(nv)) then
            inv=nv
         endif

         if(present(nbdy)) then
            inbdy=max(61,nbdy)
         else
            inbdy=inbdy_dflt
         endif

         !===============================
         sp => s
         call t1mhdeqi_geqdsk(ilun,glabel,zRmin,zRmax,zZmin,zZmax,zcur, &
              idp,idq,idpsi,inh,inv,inbdy,ierr)
         !===============================
         exit
      enddo

      if(ierr.ne.0) call xplasma_errmsg_append(s,' error in xplasma_wr_geqdsk')
      if(ifile.eq.1) then
         if(ierr.ne.0) then
            close(unit=ilun,status='delete')
         else
            close(unit=ilun)
         endif
      endif

    end subroutine xplasma_wr_geqdsk

    !-------------------------------------------------------
    subroutine xplasma_rhopsi_gen(s,npsi,ierr, tol, psivals, rhovals)

      !  generate an evenly spaced Psi vector (poloidal flux, Wb/rad) that
      !  spans the plasma from axis to edge.  Psi=0 corresponds to the
      !  magnetic axis; Psibdy is taken to be a positive number; the sign
      !  is omitted.  (The sign is available separately-- it corresponds 
      !  to the direction of the toroidal plasma current, jphi_ccw,
      !  available as an optional argument in xplasma_global_info).

      !  also, find the corresponding rho = sqrt(Phi_tor/Phi_tor_bdy)
      !  values that map to the Psi values w/in tolerance (tol).

      !---------------
      !  required arguments:

      type (xplasma), pointer :: s
      integer, intent(in) :: npsi   ! number of Psi values wanted (min. 2)
      integer, intent(out) :: ierr  ! completion status, 0=OK

      !---------------
      !  optional arguments:
      real*8, intent(in), optional :: tol  ! mapping tolerance
      !  default: something close to real*8 machine precision

      real*8, dimension(:), intent(out), optional :: Psivals
      !  The evenly spaced Psi vector generated, Wb/rad, (1:npsi)

      real*8, dimension(:), intent(out), optional :: Rhovals
      !  Rho values satisfying 
      !     abs(Psi(Rhovals(i))-Psivals(i)) <= tol*[Psi at bdy]
      !  (1:npsi) -- sqrt(toroidal_flux/toroidal_flux_at_bdy)
      !  0 on axis, 1 at the edge.

      !-------------------------------------
      real*8 :: psimin,psimax
      real*8, dimension(:), allocatable :: zpsis,zrhos
      integer :: ipsi
      !-------------------------------------

      ierr=0
      if(npsi.lt.2) then
         ierr=9999
         call xplasma_errmsg_append(s, &
              ' ?xplasma_rhopsi_gen: argument error: npsi >= 2 required.')
         return
      endif

      if(present(psivals)) then
         psivals=0
         if(size(psivals).lt.npsi) then
            ierr=612
            call xplasma_errmsg_append(s, &
                ' ?xplasma_rhopsi_gen: "psivals" array size >= npsi required.')
         endif
      endif

      if(present(rhovals)) then
         rhovals=0
         if(size(rhovals).lt.npsi) then
            ierr=612
            call xplasma_errmsg_append(s, &
                ' ?xplasma_rhopsi_gen: "rhovals" array size >= npsi required.')
         endif
      endif

      if(ierr.ne.0) return

      !-----------------------
      !  error checks passed

      call xplasma_psi_range(s,psimin,psimax)

      allocate(zpsis(npsi),zrhos(npsi))

      zpsis(1)=psimin
      zpsis(npsi)=psimax
      do ipsi=2,npsi-1
         zpsis(ipsi)=((npsi-ipsi)*psimin + (ipsi-1)*psimax)/(npsi-1)
      enddo

      if(present(rhovals)) then
         call xplasma_rhopsi_find(s,zpsis,zrhos,ierr, tol=tol)
      endif

      if(ierr.eq.0) then
         if(present(psivals)) psivals(1:npsi)=zpsis(1:npsi)
         if(present(rhovals)) rhovals(1:npsi)=zrhos(1:npsi)
      endif

      deallocate(zpsis,zrhos)

    end subroutine xplasma_rhopsi_gen

    !-------------------------------------------------------
    subroutine xplasma_rhopsi_find(s,psivals,rhovals,ierr, tol, iwarn)

      !  find rho values corresponding to a specified set of Psi values
      !    Psi -- Poloidal flux, Wb/rad
      !    rho -- sqrt(Tor_flux/Tor_flux_at_bdy)

      !  all input Psi values should be in the range [Psimin,Psimax] where
      !  Psimin corresponds to the magnetic axis and Psimax-Psimin
      !  corresponds to the (unsigned) poloidal flux, Wb/rad, enclosed 
      !  within the core plasma.  Usually Psimin=0 is set.

      type (xplasma), pointer :: s
      real*8, dimension(:), intent(in) :: Psivals  ! Psi values in any order
      real*8, dimension(:), intent(out) :: rhovals ! rho values output
      !  sizes of Psivals and rhovals must match

      integer, intent(out) :: ierr   ! status code, =0 on normal exit
      !  error occurs if xplasma is unitialized or contains no MHD equilibrium;
      !  Psi-out-of-range is handled (see iwarn, below).

      real*8, intent(in), optional :: tol     ! accuracy tolerance
      !  on output, rho values satisfy
      !     abs(Psi(rhovals(i))-Psivals(i)) <= tol*[Psi at bdy]
      !  (1:npsi) -- sqrt(toroidal_flux/toroidal_flux_at_bdy)
      !  0 on axis, 1 at the edge.

      integer, intent(out), optional :: iwarn  ! #of Psi values out of range

      !  if Psi <= Psimin, rho=0 is returned; if Psi >= Psimax rho=1 is
      !  returned.

      !------------------------------------
      real*8 :: psimin,psimax,rhomin,rhomax,ztol,psitol
      real*8, dimension(:), allocatable :: zpsi,zrho,zdum,zmina,zmaxa
      integer :: ipsi,inpsi,jwarn,insize
      logical, dimension(:), allocatable :: iok

      external eqi_xpsi_fun
      !------------------------------------

      ierr=0
      if(size(psivals).ne.size(rhovals)) then
         ierr=612
         call xplasma_errmsg_append(s, &
              ' ?xplasma_rhopsi_find: "psivals" and "rhovals" array argument sizes inconsistent.')
         return
      else
         insize=size(psivals)
      endif

      rhovals = 0
      if(present(iwarn)) iwarn=insize

      rhomin = 0
      rhomax = 1

      !-----------------
      !  get bdy & search error tolerance

      if(present(tol)) then
         ztol=tol
      else
         call xplasma_global_info(s,ierr, bdytol=ztol)
         if(ierr.ne.0) return
      endif

      !-----------------
      !  get available range of Psi values
      
      call xplasma_psi_range(s,psimin,psimax)

      psitol = ztol*max(abs(psimin),abs(psimax))

      !-----------------
      !  get the list of interior values for which search is required

      inpsi=0
      allocate(zpsi(insize),zrho(insize),zdum(insize),iok(insize))
      iok = .FALSE.

      jwarn=0

      do ipsi=1,insize
         if(psivals(ipsi).lt.(psimin+psitol)) then
            rhovals(ipsi)=rhomin
            if(psivals(ipsi).lt.(psimin-psitol)) jwarn = jwarn + 1

         else if(psivals(ipsi).gt.(psimax-psitol)) then
            rhovals(ipsi)=rhomax
            if(psivals(ipsi).gt.(psimax+psitol)) jwarn = jwarn + 1

         else
            inpsi=inpsi+1
            zpsi(inpsi)=psivals(ipsi)
            zrho(inpsi)=(psivals(ipsi)-psimin)/(psimax-psimin) ! crude guess
         endif
      enddo

      !  root finder...

      if(inpsi.gt.0) then

         sp => s

         allocate(zmina(inpsi),zmaxa(inpsi)); zmina=rhomin; zmaxa=rhomax

         call zridderx(inpsi, iok, zmina, zmaxa, ztol, psitol, &
              eqi_xpsi_fun, zrho, ierr, inpsi, zpsi, 1, zdum, 1)

         deallocate(zmina,zmaxa)

         if(ierr.ne.0) then

            ierr=9999
            call xplasma_errmsg_append(s, &
                 ' ?xplasma_rhopsi_find: root finder failure.')
            jwarn=insize

         else

            inpsi = 0
            do ipsi=1,insize
               if(psivals(ipsi).lt.(psimin+psitol)) then
                  continue

               else if(psivals(ipsi).gt.(psimax-psitol)) then
                  continue

               else
                  inpsi = inpsi + 1
                  rhovals(ipsi) = zrho(inpsi)
               endif
            enddo

         endif

      endif
      if(present(iwarn)) iwarn = jwarn

      deallocate(zpsi,zrho,zdum,iok)

    end subroutine xplasma_rhopsi_find

end module xplasma_profs


!......................................................................................


subroutine t1mhdeqi_geqdsk(lun_geqdsk,geqdsk_lbl, &
     Rmin,Rmax,Zmin,Zmax,zcur, &
     id_p,id_q,id_psi_in,nh,nv,nb, &
     ierr)
  
  !  write GEQDSK file from xplasma module

  !  **correction** dmc 9 Nov 2001:
  !  follow Lang Lao sign convention for G-EQDSK files:
  !  -- psi always increasing
  !  -- direction of current specified in pcur scalar *only*

  !-----------------------------------
  !  DMC Sep 2010: add code to set magnetic axis to match Psi(R,Z) min
  !
  !  Note: surely no "thread safety" here.  A shared xplasma2 pointer "sp"
  !  from the module "eqi_rzbox_module" is used.  The call to eqi_geq_axis
  !  could involve an excursion into a root finder which needs a module with
  !  data elements that use the SAVE attribute.
  !
  !  At present in the TRANSP/xplasma2 world this routine is called only
  !  from the xplasma_profs module, and this routine makes the only existing
  !  call to eqi_geq_axis.  So... arguments could be changed to pass down
  !  xplasma2 object references and establish thread safety, but... it has
  !  not been attempted as of today.  (DMC noted Jan 2011).
  !-----------------------------------

  use xplasma_definitions
  use eqi_rzbox_module
  USE EZspline_obj
  USE EZspline
  use transp2ids_module

  implicit NONE

  INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)

  !  input:

  integer lun_geqdsk                ! logical unit where to write file

  character*48 geqdsk_lbl           ! 48 character label for GEQDSK file
  real*8 Rmin,Rmax                  ! (Rmin,Rmax) of Psi(R,Z)
  real*8 Zmin,Zmax                  ! (Zmin,Zmax) of Psi(R,Z)
  real*8 zcur                       ! plasma current (amps)

  !  [Rmin,Rmax]x[Zmin,Zmax] must contain the core plasma but not exceed the
  !  available (R,Z) grids.

  integer id_p                      ! xplasma id:  Pressure profile
  integer id_q                      ! xplasma id:  q profile
  integer id_psi_in                 ! xplasma id:  Psi (or 0 to use default)

  integer nh                        ! #of GEQDSK horizontal pts
  integer nv                        ! #of GEQDSK vertical pts

  !  note nh also controls the number of pts in the 1d profiles

  integer nb                        ! #of pts in bdy contour

  !  output:

  integer ierr                      ! completion code (0=OK, file written).

  !  local data arrays...

  real*8 psirz(nh,nv)               ! Psi(R,Z) to be written
  real*8 psi(nh)                    ! local psi grid
  real*8 fpol(nh)                   ! f = R*Bt vs. (psi)
  real*8 ffprime(nh)                ! f*f', ' means d/dpsi, vs. psi
  real*8 pres(nh)                   ! p vs. psi
  real*8 pprime(nh)                 ! p' vs. psi
  real*8 qpsi(nh)                   ! q vs. psi

  real*8 Rcoor(nh,nv)               ! R to feed ids
  real*8 Zcoor(nh,nv)               ! Z to feed ids
  real*8 BR(nh,nv)               ! B_r(R,Z) to be written
  real*8 BZ(nh,nv)               ! B_z(R,Z) to be written
  real*8 Bphi(nh,nv)               ! B_phi(R,Z) to be written
  real*8 Bmod(nh,nv)               ! |B|(R,Z) to be written
  real*8 psirz_rprime(nh,nv)               ! d Psi(R,Z) / d R 
  real*8 psirz_zprime(nh,nv)               ! d Psi(R,Z) / d Z 
  real*8 grz(nh,nv)               ! g(R , Z)
  real*8 tmp2darryr(nh,nv), tmp2darryz(nh,nv)
  real*8 tmp2darryrg(nh,nv), tmp2darryzg(nh,nv)
  real*8 tmp2darryr_prime(nh,nv), tmp2darryz_prime(nh,nv)
  real*8 tmp2darryrg_prime(nh,nv), tmp2darryzg_prime(nh,nv)
  real*8 J_tor(nh,nv), J_pll(nh,nv)

  real*8 xpsi(nh)                   ! rho grid (for equal steps in psi)
  real*8 rgrid(nh)                  ! R grid for psi(R,Z)
  real*8 zgrid(nv)                  ! Z grid for psi(R,Z)
  real*8 zztmp(nh)                  ! array of Z(k) for vector eval
  real*8 rztmp(nv)                  ! array of R(k) for vector eval

  real*8 thgrid(nb)
  real*8 xtmp(nb)
  real*8 rbdy(nb),zbdy(nb)          ! contour boundary
  real*8 rlim(nb),zlim(nb)          ! limiter boundary
  integer nblim                     ! #of limiter points retained

  real*8 :: psi0a                   ! Psi0 offset to machine axis (if avail.)

  !----------------------
  logical, parameter :: intrp_flag = .TRUE.
  !----------------------

  !  scalars

  real*8 rdim,zdim,rleft,zmid,rmaxis,zmaxis
  real*8 rmaxis_save,zmaxis_save ! for debug testing
  real*8 zpsimag,zpsibdy
  real*8 rcentr,bcentr
  real*8 pcur

  real*8 zrmin,zrmax,zzmin,zzmax

  !  dummies

  integer idum
  real*8 xdum


  !------------------------------------------------------
  !  work stuff

  integer :: nsnccwi,nsnccwb,id_g,id_psi,id_R,id_Z

  integer i,j,igot
  integer ifcns(4)
  real*8 zbuf(4)

  real*8 zbufa(nh,4)
  real*8 zbufd(nh,4)

  real*8 zbufb(nb,2)

  real*8 xsmall
  real*8 zdpsidx
  real*8 :: zpsi0,zpsia

  real*8, parameter :: ZERO = 0.0_R8
  real*8, parameter :: C2PI = 6.2831853071795862_R8

  integer id_Bphi,id_BR,id_BZ,id_Bmod
  type(ezspline1_r8) :: spln1
  type(ezspline2_r8) :: spln2
  integer it, iprofile
  real*8 mu0 
  REAL twopi

  !------------------------------------------------------

  idum=0
  xdum=ZERO
  twopi=2.*4.*atan(1.)
  mu0 = 4. * 4. * atan(1.) * 1.e-7

  !  computational region, as per GEQDSK specification

  rdim=Rmax-Rmin
  zdim=Zmax-Zmin
  rleft=Rmin
  zmid=0.5_R8*(Zmin+Zmax)

  call xplasma_common_ids(sp,ierr,id_g=id_g,id_psi=id_psi,id_R=id_R,id_Z=id_Z)
  if(ierr.ne.0) return

  if(id_g.eq.0) then
     ierr=9999
     call xplasma_errmsg_append(sp,' ?eqi_geqdsk: g(rho) not found.')
  endif

  if(id_psi_in.gt.0) then
     id_psi = id_psi_in
  endif

  if(id_psi.le.0) then
     ierr=9999
     call xplasma_errmsg_append(sp,' ?eqi_geqdsk: psi(rho) not found.')
  endif

  if(min(id_R,id_Z).eq.0) then
     ierr=9999
     call xplasma_errmsg_append(sp,' ?eqi_geqdsk: [R,Z](rho,theta) not found.')
  endif
  if(ierr.ne.0) return

  !  magnetic axis

  call xplasma_mag_axis(sp,ierr, raxis=rmaxis,zaxis=zmaxis)
  if(ierr.ne.0) return

  rmaxis_save = rmaxis ! for debug testing
  zmaxis_save = zmaxis

  !  pol. flux & axis & @ bdy -- zpsibdy > zpsimag as per xplasma & G-eqdsk
  !  convention

  call xplasma_psi_range(sp,zpsimag,zpsibdy)

  !  signed vacuum field-- at geometric center (in R) of LCFS

  call xplasma_global_info(sp, ierr, bphi_ccw=nsnccwb, jphi_ccw=nsnccwi)

  call xplasma_RZminmax_plasma(sp, zrmin,zrmax,zzmin,zzmax, ierr)
  if(ierr.ne.0) return

  rcentr=0.5*(zrmin+zrmax)
  bcentr=nsnccwb/rcentr   ! will multiply in bdy R*Bphi later

  !  signed total plasma current

  pcur=nsnccwi*abs(zcur)

  !  ok form profiles:  equispaced psi grid

  psi(1)=zpsimag
  psi(nh)=zpsibdy
  do i=2,nh-1
     psi(i)=abs(zpsimag)+(i-1)*(abs(zpsibdy)-abs(zpsimag))/(nh-1)
  enddo

  !  corresponding x grid

  call xplasma_rhopsi_find(sp,psi,xpsi,ierr)
  if(ierr.ne.0) return

  !  and for finite difference evals near the axis:

  xsmall=xpsi(1)+1.0e-4_R8*(xpsi(2)-xpsi(1))

  !  evaluate profiles & derivatives

  ifcns(1)=id_g
  ifcns(2)=id_p
  ifcns(3)=id_psi
  ifcns(4)=id_q

  !  fcn values

  call xplasma_eval_prof(sp,ifcns(1:4),xpsi,zbufa,ierr)
  if(ierr.ne.0) go to 990

  bcentr = zbufa(nh,1)*bcentr  ! multiply in R*(vacuum field) magnitude

  do i=1,nh
     fpol(i)=zbufa(i,1)*nsnccwb
     pres(i)=zbufa(i,2)
     qpsi(i)=zbufa(i,4)
  enddo

  !  derivatives -- off axis

  call xplasma_eval_prof(sp,ifcns(1:4),xpsi,zbufd,ierr, ideriv1=1)
  if(ierr.ne.0) go to 990

  !  df/dpsi = (df/dx)/(dpsi/dx)
  !   ...ok except @mag. axis where dpsi/dx -> 0

  do i=2,nh
     ffprime(i)=fpol(i)*zbufd(i,1)*nsnccwb/zbufd(i,3)
     pprime(i)=zbufd(i,2)/zbufd(i,3)
  enddo

  !  do finite difference estimate for the axis

  call xplasma_eval_prof(sp,ifcns,xsmall,zbuf,ierr)
  if(ierr.ne.0) go to 990

  zdpsidx=max(1.0e-8_R8*(psi(2)-psi(1))/(xpsi(2)-xpsi(1)), &
       (zbuf(3)-psi(1))/(xsmall-xpsi(1)))
  ffprime(1)=fpol(1)*nsnccwb* &
       ((zbuf(1)-zbufa(1,1))/(xsmall-xpsi(1)))/zdpsidx
  pprime(1)=((zbuf(2)-zbufa(1,2))/(xsmall-xpsi(1)))/zdpsidx

  !  OK... now get psi(R,Z)

  do i=1,nh
     rgrid(i)=Rmin+(i-1)*(Rmax-Rmin)/(nh-1)
  enddo

  do j=1,nv
     zgrid(j)=Zmin+(j-1)*(Zmax-Zmin)/(nv-1)
     zztmp=zgrid(j)
     call xplasma_eval_prof(sp,id_psi, &
          xplasma_R_coord,rgrid,xplasma_Z_coord,zztmp, &
          psirz(1:nh,j),ierr)
     if(ierr.ne.0) go to 990
     
     !  leave psi sign unchanged -- G-EQDSK convention is same as XPLASMA
     !  convention
  enddo

  call xplasma_get_psi0(sp,psi0a,igot,ierr)
  if(ierr.ne.0) go to 990

  !  adjust Psi(R,Z); restoring original EFIT offset; reset extrema as well

  psirz = psirz + psi0a
  zpsimag = zpsimag + psi0a
  zpsibdy = zpsibdy + psi0a

  !-------------------------------------------------------
  !  and now get the plasma boundary specification...

  ifcns(1)=id_R
  ifcns(2)=id_Z
  do i=1,nb
     xtmp(i)=1.0_R8
     thgrid(i)=(i-1)*C2PI/(nb-1)
  enddo

  !  evaluate all but last point; force exact periodicity

  call xplasma_eval_prof(sp,ifcns(1:2), &
       xplasma_theta_coord,thgrid(1:nb-1), xplasma_rho_coord,xtmp(1:nb-1), &
       zbufb(1:nb-1,1:2),ierr)
  if(ierr.ne.0) go to 990
  zbufb(nb,1)=zbufb(1,1)
  zbufb(nb,2)=zbufb(1,2)

  !  copy into output arrays

  rbdy(1:nb)=zbufb(1:nb,1)
  zbdy(1:nb)=zbufb(1:nb,2)

  ! R Z coord to feed ids

  do j=1,nv
  do i=1,nh
     Rcoor(i,j)=rgrid(i)
     Zcoor(i,j)=zgrid(j)
  enddo
  enddo

  ! BR BZ Bmod BPHI

#ifdef METHOS1
  ! method 1, read from sp, got the same number as from s
  call xplasma_common_ids(sp,ierr, &
                          id_BR=id_BR,id_BZ=id_BZ,id_Bmod=id_Bmod)
  call xplasma_find_item(sp,'Bphi',id_Bphi,ierr)
#else
  ! method 2, read from s, got the same number as from sp
  !test_xplasma/test_xplasma.f90
  call eq_gfnum('BR_RZ',id_BR)
  call eq_gfnum('BZ_RZ',id_BZ)
  call eq_gfnum('Bphi_RZ',id_Bphi)
  call eq_gfnum('Bmod_RZ',id_Bmod)
#endif

  write(*,*) " --- id_BR ...=",id_BR,id_BZ,id_Bmod
  write(*,*) " --- id_Bphi=",id_Bphi
  do j=1,nv
     zztmp=zgrid(j)
     call xplasma_eval_prof(sp,id_BR, &
          xplasma_R_coord,rgrid,xplasma_Z_coord,zztmp, BR(1:nh,j),ierr)
          !write(71,'(a,i3,e12.4)') 'BR',j,BR(nh/2,j)
          if(ierr.ne.0) go to 990
     call xplasma_eval_prof(sp,id_BZ, &
          xplasma_R_coord,rgrid,xplasma_Z_coord,zztmp, BZ(1:nh,j),ierr)
          !write(72,'(a,i3,e12.4)') 'BZ',j,BZ(nh/2,j)
          if(ierr.ne.0) go to 990
     call xplasma_eval_prof(sp,id_Bphi, &
          xplasma_R_coord,rgrid,xplasma_Z_coord,zztmp, Bphi(1:nh,j),ierr)
          !write(73,'(a,i3,e12.4)') 'Bphi',j,Bphi(nh/2,j)
          if(ierr.ne.0) go to 990
     call xplasma_eval_prof(sp,id_Bmod, &
          xplasma_R_coord,rgrid,xplasma_Z_coord,zztmp, Bmod(1:nh,j),ierr)
          !write(74,'(a,i3,e12.4)') 'Bmod',j,Bmod(nh/2,j)
          if(ierr.ne.0) go to 990
  enddo

  ! g(R,Z)=RB_phi(R,Z)
  grz=Rcoor*Bphi

  ! J_parallel J_tor

  ! pspline]$/ezspline_derivative.f90 
  ! EZspline_derivative1_array_r8
  ! subroutine EZspline_derivative1_array_r8(spline_o, i1, k1, p1, f, ier)
  ! output : f
  do j=1,nv
     zztmp=zgrid(j)
     call xplasma_eval_prof(sp,id_psi, &
          xplasma_R_coord,rgrid,xplasma_Z_coord,zztmp, &
          psirz_rprime(1:nh,j),ierr, ideriv1=1)
     if(ierr.ne.0) go to 990
  enddo
  do i=1,nh
     rztmp=rgrid(i)
     call xplasma_eval_prof(sp,id_psi, &
          xplasma_R_coord,rztmp,xplasma_Z_coord,zgrid, &
          psirz_zprime(i,1:nv),ierr, ideriv1=1)
     if(ierr.ne.0) go to 990
  enddo

! check geqdsk_mds/geqdsk_mod.f90
! calculate R^2del(1/R^2 grad Psi)
  tmp2darryr=psirz_rprime/Rcoor/Rcoor
  tmp2darryz=psirz_zprime/Rcoor/Rcoor
      CALL EZspline_init(spln2,nh,nv,(/0,0/),(/0,0/), ierr)
      CALL EZspline_error(ierr)
      spln2%x1 = rgrid
      spln2%x2 = zgrid
      CALL EZspline_setup(spln2,tmp2darryr,ierr)
      CALL EZspline_error(ierr)
      CALL ezspline_derivative(spln2,1,0,nh,nv,rgrid,zgrid,tmp2darryr_prime,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_free(spln2,ierr)
      CALL EZspline_error(ierr)

      CALL EZspline_init(spln2,nh,nv,(/0,0/),(/0,0/), ierr)
      CALL EZspline_error(ierr)
      spln2%x1 = rgrid
      spln2%x2 = zgrid
      CALL EZspline_setup(spln2,tmp2darryz,ierr)
      CALL EZspline_error(ierr)
      CALL ezspline_derivative(spln2,0,1,nh,nv,rgrid,zgrid,tmp2darryz_prime,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_free(spln2,ierr)
      CALL EZspline_error(ierr)
  J_tor=Rcoor*Rcoor*(tmp2darryr_prime+tmp2darryz_prime) / mu0 / Rcoor

! calculate g^2del(1/g 1/R^2 grad Psi)
  tmp2darryrg=tmp2darryr/grz
  tmp2darryzg=tmp2darryz/grz
      CALL EZspline_init(spln2,nh,nv,(/0,0/),(/0,0/), ierr)
      CALL EZspline_error(ierr)
      spln2%x1 = rgrid
      spln2%x2 = zgrid
      CALL EZspline_setup(spln2,tmp2darryrg,ierr)
      CALL EZspline_error(ierr)
      CALL ezspline_derivative(spln2,1,0,nh,nv,rgrid,zgrid,tmp2darryrg_prime,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_free(spln2,ierr)
      CALL EZspline_error(ierr)

      CALL EZspline_init(spln2,nh,nv,(/0,0/),(/0,0/), ierr)
      CALL EZspline_error(ierr)
      spln2%x1 = rgrid
      spln2%x2 = zgrid
      CALL EZspline_setup(spln2,tmp2darryzg,ierr)
      CALL EZspline_error(ierr)
      CALL ezspline_derivative(spln2,0,1,nh,nv,rgrid,zgrid,tmp2darryzg_prime,ierr)
      CALL EZspline_error(ierr)
      CALL EZspline_free(spln2,ierr)
      CALL EZspline_error(ierr)
  J_pll=grz*grz*(tmp2darryrg_prime+tmp2darryzg_prime)
  J_pll=J_pll/Bmod / mu0

  !-------------------------------------------------------
  !  and now, finally, the limiter...

  call xplasma_limcon(sp,rlim,zlim,igot,ierr, tol=ZERO)
  if(ierr.ne.0) go to 990

  nblim=nb

  !-------------------------------------------------------
  if(intrp_flag) then

     ! write out and read back from ascii file -- to get the
     ! effect of loss of precision due to ascii I/O in EQDSK format

     rmaxis = rmaxis_save
     zmaxis = zmaxis_save

     call eqi_geq_axis(nh,rgrid,nv,zgrid,psirz,nb,rbdy,zbdy, &
          rmaxis,zmaxis,zPsi0,zPsia,ierr)
     if(ierr.ne.0) then
        call xplasma_errmsg_append(sp, &
             ' ?Newton loop error in eqi_geq_axis, eqi_geqdsk (xplasma).')
        go to 990
     endif

     ! ** adjust Psi(R,Z) so that the true magnetic axis value is zpsimag

     psirz = psirz + (zpsimag - zPsi0)

     ! ** adjust boundary value per result of eqi_geq_axis evaluation
     !    (include the general adjustment to all of psirz(:,:) ...)

     zpsibdy = zPsia + (zpsimag - zPsi0)

  endif

  !  write the GEQDSK file...
     write(*,*) " --- g-eqdsk is written for your reference"

2000 format(a48,3i4)
2001 format(5e16.9)
2002 format(2i5)

  !  from notes (General Atomic, Lang Lao 02/21/00)

  idum=0
  xdum=ZERO

  write(lun_geqdsk,2000) geqdsk_lbl,idum,nh,nv
  write(lun_geqdsk,2001) rdim,zdim,rcentr,rleft,zmid
  write(lun_geqdsk,2001) rmaxis,zmaxis,zpsimag,zpsibdy,bcentr
  write(lun_geqdsk,2001) pcur,zpsimag,xdum,rmaxis,xdum
  write(lun_geqdsk,2001) zmaxis,xdum,zpsibdy,xdum,xdum
  write(lun_geqdsk,2001) (fpol(i),i=1,nh)
  write(lun_geqdsk,2001) (pres(i),i=1,nh)
  write(lun_geqdsk,2001) (ffprime(i),i=1,nh)
  write(lun_geqdsk,2001) (pprime(i),i=1,nh)
  write(lun_geqdsk,2001) ((psirz(i,j),i=1,nh),j=1,nv)
  write(lun_geqdsk,2001) (qpsi(i),i=1,nh)
  write(lun_geqdsk,2002) nb,nblim
  write(lun_geqdsk,2001) (rbdy(i),zbdy(i),i=1,nb)
  write(lun_geqdsk,2001) (rlim(i),zlim(i),i=1,nblim)

  !-------------------------------------------------------
     ! write out in ids eq%timeslice%profiles_2d
     write(*,*) " --- feed ids%eq at time"!, whichtimeslice, whichprofile, timeofinterest

     it=whichtimeslice
     iprofile=whichprofile
        !eq%time_slice(it)%global_quantities%psi_axis= zpsimag * twopi
        !eq%time_slice(it)%global_quantities%psi_boundary = zpsibdy * twopi

        if(associated(eq%time_slice(it)%profiles_2d(iprofile)%grid%dim1)) &
           deallocate(eq%time_slice(it)%profiles_2d(1)%grid%dim1)
        allocate(eq%time_slice(it)%profiles_2d(iprofile)%grid%dim1(nh))
        eq%time_slice(it)%profiles_2d(iprofile)%grid%dim1(:)= rgrid(:)

        if(associated(eq%time_slice(it)%profiles_2d(iprofile)%grid%dim2)) &
           deallocate(eq%time_slice(it)%profiles_2d(1)%grid%dim2)
        allocate(eq%time_slice(it)%profiles_2d(iprofile)%grid%dim2(nv))
        eq%time_slice(it)%profiles_2d(iprofile)%grid%dim2(:)= zgrid(:)

        if(associated(eq%time_slice(it)%profiles_2d(iprofile)%r)) &
           deallocate(eq%time_slice(it)%profiles_2d(iprofile)%r)
        allocate(eq%time_slice(it)%profiles_2d(iprofile)%r(nh,nv))
        eq%time_slice(it)%profiles_2d(iprofile)%r=Rcoor

        if(associated(eq%time_slice(it)%profiles_2d(iprofile)%z)) &
           deallocate(eq%time_slice(it)%profiles_2d(iprofile)%z)
        allocate(eq%time_slice(it)%profiles_2d(iprofile)%z(nh,nv))
        eq%time_slice(it)%profiles_2d(iprofile)%z=Zcoor

        if(associated(eq%time_slice(it)%profiles_2d(iprofile)%psi)) &
           deallocate(eq%time_slice(it)%profiles_2d(iprofile)%psi)
        allocate(eq%time_slice(it)%profiles_2d(iprofile)%psi(nh,nv))
        eq%time_slice(it)%profiles_2d(iprofile)%psi=psirz*twopi

        if(associated(eq%time_slice(it)%profiles_2d(iprofile)%j_parallel)) &
           deallocate(eq%time_slice(it)%profiles_2d(iprofile)%j_parallel)
        allocate(eq%time_slice(it)%profiles_2d(iprofile)%j_parallel(nh,nv))
        eq%time_slice(it)%profiles_2d(iprofile)%j_parallel=J_pll

        if(associated(eq%time_slice(it)%profiles_2d(iprofile)%j_tor)) &
           deallocate(eq%time_slice(it)%profiles_2d(iprofile)%j_tor)
        allocate(eq%time_slice(it)%profiles_2d(iprofile)%j_tor(nh,nv))
        eq%time_slice(it)%profiles_2d(iprofile)%j_tor=J_tor

        if(associated(eq%time_slice(it)%profiles_2d(iprofile)%b_r)) &
           deallocate(eq%time_slice(it)%profiles_2d(iprofile)%b_r)
        allocate(eq%time_slice(it)%profiles_2d(iprofile)%b_r(nh,nv))
        eq%time_slice(it)%profiles_2d(iprofile)%b_r=BR

        if(associated(eq%time_slice(it)%profiles_2d(iprofile)%b_z)) &
           deallocate(eq%time_slice(it)%profiles_2d(iprofile)%b_z)
        allocate(eq%time_slice(it)%profiles_2d(iprofile)%b_z(nh,nv))
        eq%time_slice(it)%profiles_2d(iprofile)%b_z=BZ

        if(associated(eq%time_slice(it)%profiles_2d(iprofile)%b_tor)) &
           deallocate(eq%time_slice(it)%profiles_2d(iprofile)%b_tor)
        allocate(eq%time_slice(it)%profiles_2d(iprofile)%b_tor(nh,nv))
        eq%time_slice(it)%profiles_2d(iprofile)%b_tor=Bphi

        ! ids needs poloidal flux decreases when move from R0 to R0+a
        if(associated(cp%profiles_1d(it)%grid%psi)) &
        cp%profiles_1d(it)%grid%psi(:)= cp%profiles_1d(it)%grid%psi(:) - zpsimag * twopi 

        ! todo use ezspline to interpolate
        if(associated(eq%time_slice(it)%profiles_1d%f_df_dpsi)) &
        eq%time_slice(it)%profiles_1d%f_df_dpsi(:)= ffprime(:)
        if(associated(eq%time_slice(it)%profiles_1d%dpressure_dpsi)) &
        eq%time_slice(it)%profiles_1d%dpressure_dpsi(:)= pprime(:)

  !  that is all

  ierr=0

  go to 1000

  !-----------------------------------
  !  error exit...

990 continue

  ierr=9999
  call xplasma_errmsg_append(sp,' ?data lookup error in eqi_geqdsk (xplasma).')

1000 continue
  return

end subroutine t1mhdeqi_geqdsk

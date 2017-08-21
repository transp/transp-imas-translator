program transp2imas

!---------------------------------------------------------------------
!
!  to run:
!     ./transp2imas <runid>
!
!---------------------------------------------------------------------

   use transp2imas_module
   use ezspline
   use ezspline_obj

   implicit NONE

   integer, parameter :: naxmgf = 126 ! should be same as CPLOTR NAXMGF

   integer imsg,iout                 ! i/o unit numbers

   common/msgs/ imsg,iout

   character*80 arg,argu             ! command line argument
   character*80 rpfile               ! string to establish runid connection

   integer ier,iwarn,istat           ! subroutine completion codes

   character*32 rlabel               ! run label
   character*32 zlabel
   character*16 zunits
!
!  sizes:  scalar & profile databases, max x axis, max item over
!
   integer nsctime,nprtime,nxmax,nmax,iret
!
!  data buffers
!
   real, dimension(:), allocatable :: sctime,prtime,scdata,prdata
   real, dimension(:,:), allocatable :: mgdata,mgslice
   real ztime,dt_avg,zanswer

   ! real*8 is needed to pass time of interest to t1mhdgeq
   real*8 r8ztime
!
!  array for x axis names
!
   integer naxes,nlist               ! number of x axes, no. in list
   integer ngot                      ! should match nlist after allocation
   character(10), dimension(:), allocatable :: names,xnames
   integer, dimension(:), allocatable :: xsizes
   integer XItype, XBtype, LIMtype
!
!  array for mg member names
!
   integer infuns                    ! no. of mg members
   character(10), dimension(naxmgf) :: mgnames ! names
   integer, dimension(naxmgf) :: mgsigns ! sign factors +/-1

   character(10) aname               ! an item name
   character(64) label               ! item label
   character(32) units               ! item physical units

   integer imulti                    ! multigraph flag
   integer istype                    ! data item subtype code
   integer idims(8),irank

   character*1 cmgsign(-1:1)
!
!  array for species list...
!
   integer, parameter :: maxn=150
   character*20 zlbla(maxn)
   character*10 abray(4,maxn)
   integer n, n_species
   integer itype(maxn),ifast(maxn)
   real zz(maxn),aa(maxn)
   integer izc(maxn)

   integer :: n_thi      ! no. of "non-impurity" therm. ion species
   integer :: n_thx      ! no. of "impurity" thermal ion species
   integer :: n_bi       ! no. of "beam" ion species
   integer :: n_rfi      ! no. of "RF tail" ion species
   integer :: n_fusi     ! no. of "fusion product" ion species.

   integer :: nstate     ! number of each category of impurity species
   integer :: nion       ! number of ids ions
   integer :: nbeam, nbgr
   character*20 tmpstrng, int2strng
   character*20 tmps1, tmps2, tmps3, tmps4
   integer :: lentmps2, ktype, kerr
!
!  namelist test...
!
   integer ilun,ngmax
   real aplasm(5),backz(5)
!
! limiter data
!
   integer  nlimiter
   integer, dimension(:), allocatable ::  ilimiter
   real, dimension(:), allocatable ::  rlimiter,ylimiter
!
!  ids
!
   integer :: idsidx, shot, run, refshot, refrun
   character(len=32)::treename
!
!  work stuff
!
   real :: twopi = 2.*4.*atan(1.)
   logical :: ozero = .false.
   integer :: ideriv = 0
   integer i,j,ij,k
   integer it, ir, offset
   integer :: iion=0, iion_start
   integer :: iprofile, nprofile=1
   real, dimension(:), allocatable ::  bzxr
   real, dimension(:, :), allocatable ::  XI, XB
   real, dimension(:, :), allocatable ::  PLFLX, dXBdPLFLX
   real*8, dimension(:), allocatable ::  xibuf1, xibuf2, xbbuf1, xbbuf2, xbbuf3

   type(ezspline1_r8) :: spln1

   real :: rdum, rvdum(20), rvdum2(20)
   integer :: ivdum(20)

!----------------------------------------------------------------
!
   cmgsign(-1)='-'
   cmgsign(0)='?'
   cmgsign(1)='+'
!
!  file for messages, file for normal program output
!
   imsg=77
   call plc_msgs(imsg,'transp2imas.msgs')

   iout=66
   rpfile = '11114P04'
   call get_arg_count(k)
   do i = 1, k
      call get_arg(i,arg)
      argu = arg
      call uupper(argu)
      if(argu.eq.'STDOUT') then
         iout=6
      else
         rpfile = arg
      end if
   end do
   if (iout .eq. 66) then
      open(unit=iout,file='transp2imas.output',status='unknown')
   endif
!
!--------------------------
!  open imas data tree to write

   write(*,*) 'Open shot and write in IMAS !'

   shot =386
   run = 1
   !shot =120
   !run = 28
   refshot = 100
   refrun = 000
   treename = 'ids'
   ! get ids shot & run number
   call getids_shotid(rpfile,shot,run)

   call imas_create(treename,shot,run,refshot,refrun,idsidx)
!
!--------------------------
!  connect to run
!
   call kconnect(rpfile,rlabel,nsctime,nprtime,nxmax,nmax,ier)
   if(ier.ne.0) call transp2imas_error('kconnect',ier)

   write(iout,*) ' kconnect:  run label (rlabel) = "',rlabel,'"'
   write(iout,*) ' '
   write(iout,*) '   #scalar time pts (nsctime)  =  ',nsctime
   write(iout,*) '   #profile time pts (nprtime)  =  ',nprtime
   write(iout,*) '   max x axis size (nxmax)     =  ',nxmax
   write(iout,*) '   max data item size (nmax)   =  ',nmax
   write(iout,*) ' '
!
!--------------------------
!  scalar timebase, profile timebase
!  scalar data, profile data
!
   allocate(sctime(nsctime))
   allocate(prtime(nprtime))
   allocate(scdata(nsctime))
   allocate(prdata(nmax))

   allocate(mgdata(nmax,max(66,naxmgf)))
   allocate(mgslice(nxmax,max(66,naxmgf)))

   ! fill ids
   allocate(cp%time(nsctime))
   allocate(cp%profiles_1d(nprtime))
   allocate(eq%time(nsctime))
   allocate(eq%time_slice(nprtime))
   allocate(ct%time(nsctime))
   allocate(ct%model(1))
   allocate(ct%model(1)%profiles_1d(nprtime))
   allocate(es%time(nprtime))
   allocate(es%source(2)) ! Gas-flow (1) and recycling (2) sources, respectively
   allocate(es%source(1)%ggd(1))
   allocate(es%source(1)%ggd(1)%neutral(5)) ! H, D, T, He3, He4 are the possible options
   allocate(es%time(nsctime))
   allocate(es%source(2)%ggd(1))
   allocate(es%source(2)%ggd(1)%neutral(5))

   write(iout,*) ' '
   ilun=99
   call tr_getnl_text(ilun,ier)
   if(ier.ne.0) then
      call transp2imas_error('tr_getnl_text',ier)
   endif

   nbeam = 0
#if 1
   ! Get number of beams from namelist
   ivdum(1) = 0
   call tr_getnl_intvec('NBEAM', ivdum, 1, istat)
   write(*,*) 'istat = ', istat, ivdum(1)
   if (istat .eq. 1) nbeam = ivdum(1)
#else
! Count the number of beams
   do k = 1, 9
      write(int2strng, *) k
      tmps1 = adjustl(int2strng) ! delete the leading blanks
      !tmpstrng = 'BDC0'//tmps1(1:1)
      tmpstrng = 'PINJ0'//tmps1(1:1)
      !write(*, *) '.', trim(tmpstrng), '.'
      ! Now check if BDC01, BDC02, ..., BDC09 exist
      kerr = -99 ! if MDS+ tree is opened, leave open
      call rplabel(trim(tmpstrng), zlabel, zunits, ktype, kerr)
      !write(*, *) 'x', ktype, 'x'
      !write(*, *) 'y', kerr, 'y'
      !if (kerr.gt.0) nbeam = nbeam + 1
      if (kerr.eq.-1) nbeam = nbeam + 1
   end do
#endif
   write(*, *) 'nbeam =', nbeam

   if (nbeam .gt. 0) then
      allocate(nbi%time(nsctime))
      allocate(nbi%unit(nbeam))
      do k = 1, nbeam
         allocate(nbi%unit(k)%beamlets_group(1))
         allocate(nbi%unit(k)%beamlets_group(1)%beamlets%tangency_radii(1))
         allocate(nbi%unit(k)%beamlets_group(1)%beamlets%positions%r(1))
         allocate(nbi%unit(k)%beamlets_group(1)%beamlets%positions%z(1))
         allocate(nbi%unit(k)%beamlets_group(1)%beamlets%positions%phi(1))
         !allocate(nbi%unit(k)%beamlets_group(1)%divergence_component)
         allocate(nbi%unit(k)%power%time(nsctime))
         allocate(nbi%unit(k)%power%data(nsctime))
         allocate(nbi%unit(k)%energy%time(nsctime))
         allocate(nbi%unit(k)%energy%data(nsctime))
         allocate(nbi%unit(k)%beam_power_fraction%time(nsctime))
         allocate(nbi%unit(k)%beam_power_fraction%data(3, nsctime))
      end do
   end if

!---------------------------------
!  get the x axis info for profiles f(x,t)
!
   write(iout,*) ' transp2imas:  rpnumx, rpxname ... '
   call rpnumx(naxes)

   write(iout,*) ' '
   write(iout,*) ' ',naxes,' x axes:'
   write(iout,*) ' '
   allocate(xnames(naxes))
   allocate(xsizes(naxes))
   do i=1,naxes
      call rpxname(i,xnames(i))
      write(iout,*) ' ',i,xnames(i)
   enddo

   write(iout,*) ' '
   write(iout,*) ' transp2imas:  rpdims... '
   write(iout,*) ' '
   do i=0,naxes
      ij=i
      if(ij.eq.0) ij=-1
      idims=0
      call rpdims(ij,irank,idims,aname,ier)
      if(ier.ne.0) call transp2imas_error('rpdims',ier)
      if(aname.eq.' ') aname='scalar'
      write(iout,1001) ij,aname,irank,(idims(j),j=1,4)
1001  format(' istype=',i4,1x,a,1x,' rank=',i1,'   dims=',4(1x,i4))
      if (ij.gt.0) xsizes(ij) = idims(1)

      if(trim(aname).eq.'X') XItype=ij
      if(trim(aname).eq.'XB') XBtype=ij
      if(trim(aname).eq.'ILIM') LIMtype=ij
   enddo
   !write(iout,*) 'cj2: ', XItype, XBtype, LIMtype

   ! In theory, X (sometimes called XI) and XB could have different sizes.
   ! Let's try to keep that possible even if it's currently not used in practice.
   offset=xsizes(1); allocate(xibuf1(offset), xibuf2(offset))
   offset=xsizes(2); allocate(xbbuf1(offset), xbbuf2(offset), xbbuf3(offset))

   write(iout,*) ' '

!
!--------------------------
!  species list...
!
   call rd_nspecies(n_species,n_thi,n_thx,n_bi,n_rfi,n_fusi)
   ! return
   !   n_species=0 if no run is currently open.
   !   n_species=-1 if some other error occurred (rare).
   ! n_species  ! total no. of species including electrons
   ! n_thi      ! no. of "non-impurity" therm. ion species
   ! n_thx      ! no. of "impurity" thermal ion species
   ! n_bi       ! no. of "beam" ion species
   ! n_rfi      ! no. of "RF tail" ion species
   ! n_fusi     ! no. of "fusion product" ion species.

   if(maxn.lt.n_species) then
      write(iout,*) ' err cj1: rd_species not allocated enough space'
      write(iout,*) ' increase maxn from ', maxn, ' to ', n_species
      call transp2imas_exit(' ?? rd_nspecies not enough space')
   endif

   call rd_species(maxn,n_species,zlbla,abray,itype,ifast,zz,aa,izc)
   ! return
   !   n_species=0 if no run is currently open.
   !   n_species=-1 if some other error occurred (rare).
   ! n_species  ! total no. of species including electrons
   ! zlbla(maxn)      !  descr. label for each specie
   ! abray(4,maxn)    !  specie profile names
   !   for specie j:
   !     abray(1,j) = name of density profile-- for all species (n/cm**3).
   !     abray(2,j) = name of temperature or "average energy" profile.
   !     abray(3,j) = name of perpendicular energy density profile or
   !                  average perp energy per particle or " ".
   !     abray(4,j) = name of parallel energy density profile or
   !                  average pll energy per particle or " ".
   !        (generally abray(3:4,j) will be set for fast species only)
   ! itype       ! species type code
   !   itype(j)=-1  -- electrons
   !   itype(j)=+1  -- non-impurity thermal specie, usually H or He isotope
   !                   can be Li
   !   itype(j)=+2  -- impurity thermal specie:  Z and A are known, constant
   !   itype(j)=+3  -- impurity thermal specie:  Z and A are functions of time
   !   itype(j)=+4  -- beam ion specie
   !   itype(j)=+5  -- rf tail ion specie
   !   itype(j)=+6  -- fusion product ion specie
   ! ifast      ! thermal/fast flag
   ! zz     ! charge of specie
   ! aa     ! mass of specie (amu)
   ! izc     ! atomic number of species

   write(iout,*) ' '
   write(iout,*) ' number of plasma species in run:  ',n_species
   write(iout,*) ' associated profiles: saved in fort.131'
   do i=1,n_species
      write(131,'(1x,a,'' :: '',4(a,1x))') &
         zlbla(i),(abray(j,i),j=1,4)
      write(131, &
         '(10x,'' Z='',f10.5,'' A='',f10.5,'' IZ='',1x,i3,'// &
         ''' itype,ifast='',2(1x,i3))') zz(i),aa(i),izc(i), &
         itype(i),ifast(i)
   enddo
   write(iout,*) '       # of ion: ',n_thi
   write(iout,*) '       # of impurity: ',n_thx
   write(iout,*) '       # of beam ion: ',n_bi
   write(iout,*) '       # of rf ion: ',n_rfi
   write(iout,*) '       # of fusion ion: ',n_fusi
   write(iout,*) ' '

   ! how many categories of impurity
   ! each category of impurity specie is represented by NIMPS_,
   ! such as NIMPS_Ar NIMPS_B NIMPS_O NIMPS_MO
   ! here nlist=4
   ! example:
   ! profile:NIMPS_: category of impurity species=4
   CALL RPNLIST('profile:NIMPS_',0,nlist)
   ! get number of profile functions whose names contain
   ! "NIMP_" as a substring.  The substring test is based on
   ! a case blind comparison.
   write(iout,*) ' profile:NIMPS_: category # of impurity species=',nlist

   ! count the total number of ions
   nion=nlist+n_thi
   do it=1,nprtime
      allocate(cp%profiles_1d(it)%ion(nion))
      allocate(ct%model(1)%profiles_1d(it)%ion(nion))
      do i = 1, nion
         allocate(cp%profiles_1d(it)%ion(i)%element(1))
         allocate(cp%profiles_1d(it)%ion(i)%label(1))
         allocate(ct%model(1)%profiles_1d(it)%ion(i)%element(1))
         allocate(ct%model(1)%profiles_1d(it)%ion(i)%label(1))
      enddo
   enddo
   write(iout,*) ' There are total ', nion, ' ions', &
      nlist, ' of them are impurity'

   !
   ! work on thermal ions first ...
   ! the default iion is 0
   !

   iion=0
   do i=1,n_species
      if(itype(i) .eq. 1) then
         iion=iion+1
         iion_start=i  ! position in n_species
         write(iout,*) ' specie', i, 'th ion is ', zlbla(i)

         do it=1,nprtime
            cp%profiles_1d(it)%ion(iion)%element(1)%a=aa(i)
            cp%profiles_1d(it)%ion(iion)%z_ion=zz(i)
            cp%profiles_1d(it)%ion(iion)%element(1)%z_n=izc(i)
            cp%profiles_1d(it)%ion(iion)%label(1)=zlbla(i)
!           cp%profiles_1d(it)%ion(iion)%multiple_charge_states_flag=0
            ct%model(1)%profiles_1d(it)%ion(iion)%element(1)%a=aa(i)
            ct%model(1)%profiles_1d(it)%ion(iion)%element(1)%z_n=izc(i)
            ct%model(1)%profiles_1d(it)%ion(iion)%label(1)=zlbla(i)
         enddo

         call rplabel(abray(1,i),label,units,imulti,istype)
         if(imulti.ne.0) call transp2imas_exit(' ?? imulti.ne.0')
         !write(iout,*) ' ',i,abray(1,i),abray(2,i),&
         !'"',label,'" "',units,'"  istype=', istype

         if(istype.eq.1) then
            !get ion density
            write(iout,*) ' '
            call rprofile(trim(abray(1,i)),prdata,nprtime*xsizes(istype),iret,ier)
            if(ier.ne.0) call transp2imas_error('rprofile(abray(1,i))',ier)
            if(iret.ne.nprtime*xsizes(istype)) &
               call transp2imas_exit(' ?? '//abray(1,i)//' read error')
            call transp2imas_echo(abray(1,i),prdata,xsizes(istype),nprtime)

            offset=xsizes(istype)
            do it=1,nprtime
               allocate(cp%profiles_1d(it)%ion(iion)%density(offset))
               cp%profiles_1d(it)%ion(iion)%density(1:offset)=&
                  prdata(1+(it-1)*offset:it*offset)*1.e6
            enddo

            !get ion temperature
            write(iout,*) ' '
            call rprofile(trim(abray(2,i)),prdata,nprtime*xsizes(istype),iret,ier)
            if(ier.ne.0) call transp2imas_error('rprofile(abray(2,i))',ier)
            if(iret.ne.nprtime*xsizes(istype)) &
               call transp2imas_exit(' ?? '//abray(2,i)//' read error')
            call transp2imas_echo(abray(2,i),prdata,xsizes(istype),nprtime)

            offset=xsizes(istype)
            do it=1,nprtime
               allocate(cp%profiles_1d(it)%ion(iion)%temperature(offset))
               cp%profiles_1d(it)%ion(iion)%temperature(1:offset)=&
                  prdata(1+(it-1)*offset:it*offset)
            enddo
         endif
      endif
   enddo

   ! now handle the fast ions

   do n = 1, n_species
      if (itype(n) .ge. 4) then
         call rplabel(abray(1,n),label,units,imulti,istype)
         if (imulti.ne.0) call transp2imas_exit(' ?? imulti.ne.0')
         if (istype.eq.1) then
            !get ion density
            write(iout,*) ' '
            call rprofile(trim(abray(1,n)),prdata,nprtime*xsizes(istype),iret,ier)
            if(ier.ne.0) call transp2imas_error('rprofile(abray(1,n))',ier)
            if(iret.ne.nprtime*xsizes(istype)) &
               call transp2imas_exit(' ?? '//abray(1,n)//' read error')
            call transp2imas_echo(abray(1,n),prdata,xsizes(istype),nprtime)
            offset=xsizes(istype)
            do it = 1, nprtime
               do i = 1, iion ! Add density contribution from fast ions of same species [ion(i)]
                  if ((aa(n)  .eq. cp%profiles_1d(it)%ion(i)%element(1)%a) .and. &
                      (zz(n)  .eq. cp%profiles_1d(it)%ion(i)%z_ion)        .and. &
                      (izc(n) .eq. cp%profiles_1d(it)%ion(i)%element(1)%z_n)) then
                         if (associated(cp%profiles_1d(it)%ion(i)%density_fast)) then
                            cp%profiles_1d(it)%ion(i)%density_fast(1:offset) = &
                               cp%profiles_1d(it)%ion(i)%density_fast(1:offset) + &
                               prdata(1+(it-1)*offset:it*offset) * 1.0e6
                         else
                            allocate(cp%profiles_1d(it)%ion(i)%density_fast(offset))
                            cp%profiles_1d(it)%ion(i)%density_fast(1:offset) = &
                               prdata(1+(it-1)*offset:it*offset) * 1.0e6
                         endif
                         cp%profiles_1d(it)%ion(i)%density(1:offset) =      &
                            cp%profiles_1d(it)%ion(i)%density(1:offset) +   &
                            prdata(1+(it-1)*offset:it*offset) * 1.0e6
                  endif
               enddo
            enddo
            !get ion perpendicular energy density
            write(iout,*) ' '
            call rprofile(trim(abray(3,n)),prdata,nprtime*xsizes(istype),iret,ier)
            if(ier.ne.0) call transp2imas_error('rprofile(abray(3,n))',ier)
            if(iret.ne.nprtime*xsizes(istype)) &
               call transp2imas_exit(' ?? '//abray(3,n)//' read error')
            call transp2imas_echo(abray(3,n),prdata,xsizes(istype),nprtime)
            offset=xsizes(istype)
            do it=1,nprtime
               do i = 1, iion
                  if ((aa(n)  .eq. cp%profiles_1d(it)%ion(i)%element(1)%a) .and. &
                      (zz(n)  .eq. cp%profiles_1d(it)%ion(i)%z_ion)        .and. &
                      (izc(n) .eq. cp%profiles_1d(it)%ion(i)%element(1)%z_n)) then
                         if (associated(cp%profiles_1d(it)%ion(i)%pressure_fast_perpendicular)) then
                            cp%profiles_1d(it)%ion(i)%pressure_fast_perpendicular(1:offset) = &
                               cp%profiles_1d(it)%ion(i)%pressure_fast_perpendicular(1:offset) + &
                               prdata(1+(it-1)*offset:it*offset) * 1.0e6 ! J/cm^3 to Pa
                         else
                            allocate(cp%profiles_1d(it)%ion(i)%pressure_fast_perpendicular(offset))
                            cp%profiles_1d(it)%ion(i)%pressure_fast_perpendicular(1:offset) = &
                               prdata(1+(it-1)*offset:it*offset) * 1.0e6 ! J/cm^3 to Pa
                         endif
                     !if (it .eq. nprtime) then
                     !   write(*,*) n, abray(3,n)
                     !   write(*,*) i, cp%profiles_1d(it)%ion(i)%pressure_fast_perpendicular(1:offset)
                     !endif
                  endif
               enddo
            enddo
            !get ion parallel energy density
            write(iout,*) ' '
            call rprofile(trim(abray(4,n)),prdata,nprtime*xsizes(istype),iret,ier)
            if(ier.ne.0) call transp2imas_error('rprofile(abray(4,n))',ier)
            if(iret.ne.nprtime*xsizes(istype)) &
               call transp2imas_exit(' ?? '//abray(4,n)//' read error')
            call transp2imas_echo(abray(4,n),prdata,xsizes(istype),nprtime)
            offset=xsizes(istype)
            do it=1,nprtime
               do i = 1, iion
                  if ((aa(n)  .eq. cp%profiles_1d(it)%ion(i)%element(1)%a) .and. &
                      (zz(n)  .eq. cp%profiles_1d(it)%ion(i)%z_ion)        .and. &
                      (izc(n) .eq. cp%profiles_1d(it)%ion(i)%element(1)%z_n)) then
                         if (associated(cp%profiles_1d(it)%ion(i)%pressure_fast_parallel)) then
                            cp%profiles_1d(it)%ion(i)%pressure_fast_parallel(1:offset) = &
                               cp%profiles_1d(it)%ion(i)%pressure_fast_parallel(1:offset) + &
                               prdata(1+(it-1)*offset:it*offset) * 1.0e6 ! J/cm^3 to Pa
                         else
                            allocate(cp%profiles_1d(it)%ion(i)%pressure_fast_parallel(offset))
                            cp%profiles_1d(it)%ion(i)%pressure_fast_parallel(1:offset) = &
                               prdata(1+(it-1)*offset:it*offset) * 1.0e6 ! J/cm^3 to Pa
                         endif
                  endif
               enddo
            enddo
         endif
      endif
   enddo

   ! set back to default 0
   iion=0

   !
   ! work on impurity then ...
   !

   ! get the actual list name of "nlist" categories of impurity
   ! example:
   ! profile:NIMPS_: category of impurity species named NIMPS_AR
   ! profile:NIMPS_: category of impurity species named NIMPS_B
   ! profile:NIMPS_: category of impurity species named NIMPS_MO
   ! profile:NIMPS_: category of impurity species named NIMPS_O
   allocate(names(nlist))
   call rplist('profile:NIMPS_',0,names,nlist,ngot,ier)
   if(ier.ne.0) call transp2imas_error('rplist',ier)
   if(ngot.ne.nlist) then
      !  nlist should match due to prior rpnlist call...
      call transp2imas_exit(' ?? ngot.ne.nlist')
   endif

   ! for each list i of "nlist" categories of impurity, find
   ! 1) its istype
   ! 2) how many charge states "nstate",
   !    such as NIMP_B1, NIMP_B2, NIMP_B3, ...
   ! 3) profile data "prdata" for each charge state
   do i=1,nlist
      write(iout,*) ' '
      call rplabel(names(i),label,units,imulti,istype)
      if(imulti.ne.0) call transp2imas_exit(' ?? imulti.ne.0')
      write(iout,*) ' impurity',i,names(i),&
         '"',label,'" "',units,'"  xid=', istype, &
         'position in n_species=',iion_start+1

      iion=n_thi+i
      !do it=1,nprtime
      !cp%profiles_1d(it)%ion(iion)%multiple_charge_states_flag=1
      !enddo
      if(istype.eq.1) then
         !get impurity density
         write(iout,*) ' '
         call rprofile(trim(names(i)),prdata,nprtime*xsizes(istype),iret,ier)
         if(ier.ne.0) call transp2imas_error('rprofile(names(i))',ier)
         if(iret.ne.nprtime*xsizes(istype)) &
            call transp2imas_exit(' ?? '//names(i)//' read error')
         call transp2imas_echo(names(i),prdata,xsizes(istype),nprtime)

         offset=xsizes(istype)
         do it=1,nprtime
            allocate(cp%profiles_1d(it)%ion(iion)%density(offset))
            cp%profiles_1d(it)%ion(iion)%density(1:offset)=&
               prdata(1+(it-1)*offset:it*offset)*1.e6
         enddo

         !get impurity temperature
         write(iout,*) ' '
         call rprofile('TX',prdata,nprtime*xsizes(istype),iret,ier)
         if(ier.ne.0) call transp2imas_error('rprofile(TX)',ier)
         if(iret.ne.nprtime*xsizes(istype)) &
            call transp2imas_exit(' ?? '//'TX'//' read error')
         call transp2imas_echo('TX',prdata,xsizes(istype),nprtime)

         offset=xsizes(istype)
         do it=1,nprtime
            allocate(cp%profiles_1d(it)%ion(iion)%temperature(offset))
            cp%profiles_1d(it)%ion(iion)%temperature(1:offset)=&
               prdata(1+(it-1)*offset:it*offset)
         enddo
      endif

      write(iout,*) ' '

      ! count the number of charge states in each category of impurity
      ! by getting rid of "S" at position 5
      ! example:
      ! cj3 profile impurity profile:NIMP_AR has           18  charge states
      ! cj3 profile impurity profile:NIMP_B has            5 charge states
      ! cj3 profile impurity profile:NIMP_MO has           42 charge states
      ! cj3 profile impurity profile:NIMP_O has            8 charge states
      tmpstrng= 'profile:'//names(i)(1:4)//names(i)(6:)
      CALL RPNLIST(trim(tmpstrng),istype,nstate)
      write(iout,*) ' impurity ', trim(tmpstrng), &
         ' has ', nstate, ' charge states', &
         ' position in n_species ', iion_start+1

      do it=1,nprtime
         allocate(cp%profiles_1d(it)%ion(n_thi+i)%state(nstate))
         do j = 1, nstate
            allocate(cp%profiles_1d(it)%ion(n_thi+i)%state(j)%label(1))
         enddo
      enddo

      ! loop through each charge state for all impurity species
      do j = 1, nstate

         write(iout,*) ' '

         ! Split each TRANSP impurity species into a separate IMAS ion species for each charge state
         ! Instead of guessing charge numbers from names, we'll loop from 1 to 100 and create a species only when we get a hit
         do k = 1, 100
            write(int2strng,*) k
            tmps1=adjustl(int2strng) ! delete the leading blank
            tmps2=names(i)(6:)
            lentmps2=len_trim((tmps2)) ! delete the trailing blank
            tmpstrng= names(i)(1:4)//names(i)(6:6+lentmps2-1)//'_'//tmps1
            ! rprofile() is called with first argument 'NIMP_AR_1', which doesn't exist, but 'NIMP_AR_18' does
            call rprofile(trim(tmpstrng),prdata,nprtime*xsizes(istype),iret,ier)
            if (ier.eq.0) then
	            write(iout,*) ' cj3 ', i, 'th impurity', &
	               ' charge_state: ', trim(tmps1), &
	               ' name: ', trim(tmpstrng), &
	               ' tmps2: ', tmps2, &
	               ' lentmps2: ', lentmps2
	            if (iret.ne.nprtime*xsizes(istype)) &
             	  call transp2imas_exit(' ?? '//tmpstrng//' read error')
            	call transp2imas_echo(tmpstrng,prdata,xsizes(istype),nprtime)
            	!cp%profiles_1d(it)%ion(iion)%state(j)%n_z=prdata
    	        write(iout,*) ' ion',iion,'position in n_species ',&
        	       iion_start+j, ' named ', abray(1,iion_start+j), abray(2,iion_start+j)

        	    do it = 1, nprtime
        	       cp%profiles_1d(it)%ion(iion)%element(1)%a=aa(iion_start+j)
        	       cp%profiles_1d(it)%ion(iion)%z_ion=zz(iion_start+j)
        	       cp%profiles_1d(it)%ion(iion)%element(1)%z_n=izc(iion_start+j)
        	       cp%profiles_1d(it)%ion(iion)%label(1)=label
        	       cp%profiles_1d(it)%ion(iion)%state(j)%label(1)=zlbla(iion_start+j)
        	       ct%model(1)%profiles_1d(it)%ion(iion)%element(1)%a=aa(iion_start+j)
        	       ct%model(1)%profiles_1d(it)%ion(iion)%element(1)%z_n=izc(iion_start+j)
        	       ct%model(1)%profiles_1d(it)%ion(iion)%label(1)=label
        	    end do

	            call rplabel(abray(1,iion_start+j),label,units,imulti,istype)
	            if (imulti.ne.0) call transp2imas_exit(' ?? imulti.ne.0')
	            write(iout,*) ' ',i,abray(1,i),abray(2,i),&
    	           '"',label,'" "',units,'"  istype=', istype

    	        if (istype.eq.1) then
	               !get charge state density
	               write(iout,*) ' '
	               call rprofile(trim(abray(1,iion_start+j)),prdata,nprtime*xsizes(istype),iret,ier)
	               if (ier.ne.0) call transp2imas_error('rprofile(abray(1,iion_start+j))',ier)
	               if (iret.ne.nprtime*xsizes(istype)) &
	                  call transp2imas_exit(' ?? '//abray(1,iion_start+j)//' read error')
    	           call transp2imas_echo(abray(1,iion_start+j),prdata,xsizes(istype),nprtime)

    	           offset=xsizes(istype)
    	           do it = 1, nprtime
	                  allocate(cp%profiles_1d(it)%ion(iion)%state(j)%density(offset))
    	              cp%profiles_1d(it)%ion(iion)%state(j)%density(1:offset)=&
    	                 prdata(1+(it-1)*offset:it*offset)*1.e6
    	           end do

    	           !get charge state temperature
    	           write(iout,*) ' '
    	           call rprofile(trim(abray(2,iion_start+j)),prdata,nprtime*xsizes(istype),iret,ier)
    	           if (ier.ne.0) call transp2imas_error('rprofile(abray(2,iion_start+j))',ier)
    	           if (iret.ne.nprtime*xsizes(istype)) &
    	              call transp2imas_exit(' ?? '//abray(2,iion_start+j)//' read error')
    	           call transp2imas_echo(abray(2,iion_start+j),prdata,xsizes(istype),nprtime)

    	           offset=xsizes(istype)
    	           do it = 1, nprtime
    	              allocate(cp%profiles_1d(it)%ion(iion)%state(j)%temperature(offset))
    	              cp%profiles_1d(it)%ion(iion)%state(j)%temperature(1:offset)=&
    	                 prdata(1+(it-1)*offset:it*offset)
    	           end do
		        end if
            end if
         end do
      end do
      iion_start = iion_start+nstate ! position in n_species
   end do
   if (allocated(names)) deallocate(names)

!
!---------------------------------
! get the scalar info f(t)
!
   write(iout,*) ' '
   write(iout,*) ' transp2imas:  rpnlist, rplist, rplabel (scalars)...'
   call rpnlist('scalar',0,nlist)

   write(iout,*) ' '
   write(iout,*) ' ',nlist,' scalar functions: saved in fort.132'
   write(iout,*) ' '
   allocate(names(nlist))
   call rplist('scalar',0,names,nlist,ngot,ier)
   if(ier.ne.0) call transp2imas_error('rplist',ier)
   if(ngot.ne.nlist) then
      !  nlist should match due to prior rpnlist call...
      call transp2imas_exit(' ?? ngot.ne.nlist')
   endif

   do i=1,nlist
      call rplabel(names(i),label,units,imulti,istype)
      if(imulti.ne.0) call transp2imas_exit(' ?? imulti.ne.0')
      if(istype.ne.-1) call transp2imas_exit(' ?? istype.ne.-1')
      write(132,*) ' ',i,names(i),'"',label,'" "',units,'"'
   enddo
   if(allocated(names)) deallocate(names)

!---------------------------------
! get the profile info f(x,t)
!---------------------------------
!
! Profiles in the core_profiles IDS are evaluated at zone centers, i.e. they're functions of X.
! Profiles in the equilibrium IDS are evaluated at zone boundaries, i.e. they're functions of XB.
! When needed, we use interpolation (linear or spline) to transform between the two.

   write(iout,*) ' '
   write(iout,*) &
      ' transp2imas:  rpnlist, rplist, rplabel (profiles)...'
   call rpnlist('profile',0,nlist)

   write(iout,*) ' '
   write(iout,*) ' ',nlist,' profile functions: saved in fort.133'
   write(iout,*) ' '
   allocate(names(nlist))
   call rplist('profile',0,names,nlist,ngot,ier)
   if(ier.ne.0) call transp2imas_error('rplist',ier)
   if(ngot.ne.nlist) then
      !  nlist should match due to prior rpnlist call...
      call transp2imas_exit(' ?? ngot.ne.nlist')
   endif

   do i=1,nlist
      call rplabel(names(i),label,units,imulti,istype)
      if(imulti.ne.0) call transp2imas_exit(' ?? imulti.ne.0')
      write(133,*) ' ',i,names(i),'"',label,'" "',units,'"  xid=', &
         istype
   enddo
   if(allocated(names)) deallocate(names)

!---------------------------------
! get profile info vs. specific x axis
!
   do j=1,naxes
      call rpnlist('profile',j,nlist)
      !  if ztype='scalar', istype is ignored; for 'profile'
      !  set istype = 0 to not restrict the list by subtype.
      !  set istype = -1 for scalars only (e.g. if zytpe.eq.'multi')
      !  set istype = +N for type N profiles only

      write(iout,*) ' '
      write(iout,*) ' ',nlist,' profile functions vs. ',xnames(j), &
         'saved in fort.', 133+j
      write(133+j,*) ' '
      allocate(names(nlist))
      call rplist('profile',j,names,nlist,ngot,ier)
      if(ier.ne.0) call transp2imas_error('rplist',ier)
      if(ngot.ne.nlist) then
         !  nlist should match due to prior rpnlist call...
         call transp2imas_exit(' ?? ngot.ne.nlist')
      endif

      do i=1,nlist
         call rplabel(names(i),label,units,imulti,istype)
         if(imulti.ne.0) call transp2imas_exit(' ?? imulti.ne.0')
         if(istype.ne.j) call transp2imas_exit(' ?? istype.ne.j')
         write(133+j,*) ' ',i,names(i),'"',label,'" "',units,'"'
      enddo
      if(allocated(names)) deallocate(names)
   enddo

!---------------------------------
! get mg info
!
   write(iout,*) ' '
   write(iout,*) ' transp2imas multigraphs:'
   write(iout,*) '     rpnlist, rplist, rplabel, rpmulti...'
   do j=0,naxes
      ij=j
      if(ij.eq.0) ij=-1
      call rpnlist('multi',ij,nlist)

      write(iout,*) ' '
      if(ij.gt.0) then
         write(iout,*) ' ',nlist,' profile multigraphs vs. ', &
            xnames(j), 'saved in fort.',133+naxes+j+1
      else
         write(iout,*) ' ',nlist,' scalar multigraphs',&
            'saved in fort.',133+naxes+j+1
      endif
      if(nlist.gt.0) then
         write(iout,*) ' '
         allocate(names(nlist))
         call rplist('multi',ij,names,nlist,ngot,ier)
         if(ier.ne.0) call transp2imas_error('rplist',ier)
         if(ngot.ne.nlist) then
            !  nlist should match due to prior rpnlist call...
            call transp2imas_exit(' ?? ngot.ne.nlist')
         endif

         do i=1,nlist
            call rplabel(names(i),label,units,imulti,istype)
            !write(*,*) j, i, 'multigraph name:', names(i)
            if(imulti.ne.1) call transp2imas_exit(' ?? imulti.ne.1')
            if(istype.ne.ij) &
               call transp2imas_exit(' ?? istype.ne.ij')
            write(133+naxes+j+1,*) ' ',i,names(i),'"',label,'" "',units,'"'
            call rpmulti(names(i),istype,label,units,infuns,mgsigns, &
               mgnames,ier)
            if(ier.ne.0) call transp2imas_error('rpmulti',ier)

            do k=1,infuns
               if(mgsigns(k).eq.-1) then
                  write(133+naxes+j+1,*) '     -',mgnames(k)
               else if(mgsigns(k).eq.1) then
                  write(133+naxes+j+1,*) '     +',mgnames(k)
               else
                  call transp2imas_exit(' ?? sign factor not +/-1')
               endif
            enddo

         enddo
         if(allocated(names)) deallocate(names)
      endif
   enddo

!---------------------------------
! get actual data
!
   write(iout,*) ' '
   write(iout,*) &
      ' get data using:  rptime_s, rptime_p, rpscalar, rprofile'

!
! rptime_s get timebase for f(t)
!
   write(iout,*) ' '
   call rptime_s(sctime,nsctime,iret)
   if(iret.ne.nsctime) &
      call transp2imas_exit(' ?? sctime read error')
   call transp2imas_echo('scalar time',sctime,1,nsctime)

   ! fill ids scalar time
   write(*,*) 'fill ids sctime time'
   cp%time = sctime
   eq%time = sctime
   ct%time = sctime
   if (nbeam .gt. 0) nbi%time = sctime
   do k = 1, nbeam
      nbi%unit(k)%power%time(:) = sctime(:)
      nbi%unit(k)%energy%time = sctime
      nbi%unit(k)%beam_power_fraction%time = sctime
   end do

!
! rptime_p get timebase for f(X,t)
!
   write(iout,*) ' '
   call rptime_p(prtime,nprtime,iret)
   if(iret.ne.nprtime) &
      call transp2imas_exit(' ?? prtime read error')
   call transp2imas_echo('profile time',prtime,1,nprtime)

   ! fill ids profile time
   write(*,*) 'fill ids prtime'
   cp%profiles_1d(:)%time = prtime
   eq%time_slice(:)%time = prtime
   es%time = prtime
   do it = 1, nprtime
      ct%model(1)%profiles_1d(it)%time = prtime(it)
   enddo

!
! rpscalar get scalar data
!

   ! fill core_profiles IDS

   write(iout,*) ' '
   ! Inductance Definition for LI_3: source/doc/beta.doc
   ! same as ids li_3 definition : i.e.
   ! li_3 = 2/R0/mu0^2/Ip^2 * int(Bp^2 dV).
   call rpscalar('LI_3',scdata,nsctime,iret,ier)
   if(ier.ne.0) call transp2imas_error('rpscalar',ier)
   if(iret.ne.nsctime) &
      call transp2imas_exit(' ?? LI_3 read error')
   call transp2imas_echo('LI_3',scdata,1,nsctime)

   allocate(cp%global_quantities%li(nsctime))
   cp%global_quantities%li(:)= scdata(:)
   eq%time_slice(:)%global_quantities%li_3 = scdata(:)

   write(iout,*) ' '
   !VSUR0: measured value from VSF ufile
   !call rpscalar('VSUR0',scdata,nsctime,iret,ier)
   call rpscalar('VSURC',scdata,nsctime,iret,ier)
   if(ier.ne.0) call transp2imas_error('rpscalar',ier)
   if(iret.ne.nsctime) &
      call transp2imas_exit(' ?? VSURC read error')
   call transp2imas_echo('VSURC',scdata,1,nsctime)
   allocate(cp%global_quantities%v_loop(nsctime))
   cp%global_quantities%v_loop(:)= scdata(:)

   write(iout,*) ' '
   !ids: Poloidal beta. Defined as betap = 4 int(p dV) / [R_0 * mu_0 * Ip^2]
   call rpscalar('BETAT',scdata,nsctime,iret,ier)
   if(ier.ne.0) call transp2imas_error('rpscalar',ier)
   if(iret.ne.nsctime) &
      call transp2imas_exit(' ?? BETAT read error')
   call transp2imas_echo('BETAT',scdata,1,nsctime)

   allocate(cp%global_quantities%beta_pol(nsctime))
   cp%global_quantities%beta_pol(:)= scdata(:)

   ! fill equilibrium IDS

   write(iout,*) ' '
   call rpscalar('BZ',scdata,nsctime,iret,ier)
   if(ier.ne.0) call transp2imas_error('rpscalar',ier)
   if(iret.ne.nsctime) &
      call transp2imas_exit(' ?? BZ read error')
   call transp2imas_echo('BZ',scdata,1,nsctime)

   allocate(eq%vacuum_toroidal_field%b0(nsctime))
   eq%vacuum_toroidal_field%b0(:) = scdata(:)
   !cp%vacuum_toroidal_field%b0(:) = eq%vacuum_toroidal_field%b0(:)

   write(iout,*) ' '
   call rpscalar('BZXR',scdata,nsctime,iret,ier)
   if(ier.ne.0) call transp2imas_error('rpscalar',ier)
   if(iret.ne.nsctime) &
      call transp2imas_exit(' ?? BZXR read error')
   call transp2imas_echo('BZXR',scdata,1,nsctime)

   ! save to calculate F (also called G) function
   allocate(bzxr(nsctime))
   bzxr(:) = scdata(:)

   ! transp doesn't have R0, it needs to be calculated from BZXR/BZ
   write(iout,*) ' '
   call rpcalc('BZXR/BZ',scdata,nsctime,iret, istype, &
      iwarn,ier)
   if(ier.ne.0) call transp2imas_error('rpcalc(BZXR/BZ) ier',ier)
   if(iwarn.ne.0) call transp2imas_error('rpcalc(BZXR/BZ) iwarn',iwarn)
   write(iout,*) ' rpcalc call OK, returned istype = ',istype
   if(iret.ne.nsctime) &
      call transp2imas_exit(' ?? BZXR/BZ calc size result error')
   call transp2imas_echo('BZXR/BZ (rpcalc)',scdata,1,nsctime)
   ! take the average to fill R0 where vacuum toroidal field B0
   ! is measured
   eq%vacuum_toroidal_field%r0= &
      .5* ( scdata(1)+scdata(nsctime) )*1.e-2
   !cp%vacuum_toroidal_field%r0(:) = eq%vacuum_toroidal_field%r0(:)

   write(iout,*) ' '
   call rpscalar('PSI0_TR',scdata,nsctime,iret,ier)
   if(ier.ne.0) call transp2imas_error('rpscalar',ier)
   if(iret.ne.nsctime) &
      call transp2imas_exit(' ?? PSI0_TR read error')
   call transp2imas_echo('PSI0_TR',scdata,1,nsctime)

   eq%time_slice(:)%global_quantities%psi_axis= scdata(:) * twopi

   write(iout,*) ' '
   call rpscalar('PLFLXA',scdata,nsctime,iret,ier)
   if(ier.ne.0) call transp2imas_error('rpscalar',ier)
   if(iret.ne.nsctime) &
      call transp2imas_exit(' ?? PLFLXA read error')
   call transp2imas_echo('PLFLXA',scdata,1,nsctime)
   eq%time_slice(:)%global_quantities%psi_boundary = &
      scdata(:) * twopi

   write(iout,*) ' '
   ! PVOL, PVOLB, PVOLF, PVOL is the total volume of plasma
   call rpscalar('PVOL',scdata,nsctime,iret,ier)
   if(ier.ne.0) call transp2imas_error('rpscalar',ier)
   if(iret.ne.nsctime) &
      call transp2imas_exit(' ?? PVOL read error')
   call transp2imas_echo('PVOL',scdata,1,nsctime)

   eq%time_slice(:)%global_quantities%volume = scdata(:) * 1.e-6

   write(iout,*) ' '
   call rpscalar('PAREA',scdata,nsctime,iret,ier)
   if(ier.ne.0) call transp2imas_error('rpscalar',ier)
   if(iret.ne.nsctime) &
      call transp2imas_exit(' ?? PAREA read error')
   call transp2imas_echo('PAREA',scdata,1,nsctime)

   eq%time_slice(:)%global_quantities%area = scdata(:) * 1.e-4

   write(iout,*) ' '
   call rpscalar('RAXIS',scdata,nsctime,iret,ier)
   if(ier.ne.0) call transp2imas_error('rpscalar',ier)
   if(iret.ne.nsctime) &
      call transp2imas_exit(' ?? RAXIS read error')
   call transp2imas_echo('RAXIS',scdata,1,nsctime)

   eq%time_slice(:)%global_quantities%magnetic_axis%r = &
      scdata(:) * 1.e-2

   write(iout,*) ' '
   call rpscalar('YAXIS',scdata,nsctime,iret,ier)
   if(ier.ne.0) call transp2imas_error('rpscalar',ier)
   if(iret.ne.nsctime) &
      call transp2imas_exit(' ?? YAXIS read error')
   call transp2imas_echo('YAXIS',scdata,1,nsctime)

   eq%time_slice(:)%global_quantities%magnetic_axis%z = &
      scdata(:) * 1.e-2

   write(iout,*) ' '
   ! PCUR, PCUREQ, PCURC, which one is correct
   ! PCUR is the measured (from input data) plasma current.
   ! PCURC is computed from the plasma parameters
   call rpscalar('PCURC',scdata,nsctime,iret,ier)
   if(ier.ne.0) call transp2imas_error('rpscalar',ier)
   if(iret.ne.nsctime) &
      call transp2imas_exit(' ?? PCURC read error')
   call transp2imas_echo('PCURC',scdata,1,nsctime)
   !write(*, *) ' PCURC', nsctime, scdata(nsctime)

   eq%time_slice(:)%global_quantities%ip = scdata(1:)

   write(iout,*) ' '
   call rpscalar('Q0',scdata,nsctime,iret,ier)
   if(ier.ne.0) call transp2imas_error('rpscalar',ier)
   if(iret.ne.nsctime) &
      call transp2imas_exit(' ?? Q0 read error')
   call transp2imas_echo('Q0',scdata,1,nsctime)

   eq%time_slice(:)%global_quantities%q_axis = scdata(:)

   !? which one to choose
   !BPDM                 MAGNETICS EST. BETA(DIA)
   !BPDC                 KINETIC BETA(DIA)
   !BETAE                ELECTRON BETA (POLOIDAL)
   !BETAR                ROTATION BETA (POLOIDAL)
   !BETAI                THERMAL ION BETA POLOIDAL
   !BETAT                TOTAL BETA(POLOIDAL)
   !LI2PB                LI/2 + BETA(POLOIDAL)
   !BPEQ                 EQUILIBRIUM BETA(POLOIDAL)
   !BPDIA                DIAMAGNETIC BETA(POLOIDAL)
   !BTEQ                 EQUILIBRIUM BETA(TOROIDAL)
   !BTDIA                DIAMAGNETIC BETA(TOROIDAL)
   !BPEQ1                1D EQUILIBRIUM BETA(POLOIDAL)
   !BPDA1                1D DIAMAGNETIC BETA(POLOIDAL)
   !L2PB1                1D DEFINITION LI/2+BETA
   !BPFASTPP             TOTAL FAST ION BETA(POL) PERP
   !BPFASTPA             TOTAL FAST ION BETA(POL) PLL
   !XBFAC                MHD BETA ADJUSTMENT FACTOR

   write(iout,*) ' '
   !? Poloidal beta. Defined as betap = 4 int(p dV) / [R_0 * mu_0 * Ip^2]
   call rpscalar('BPEQ',scdata,nsctime,iret,ier)
   if(ier.ne.0) call transp2imas_error('rpscalar',ier)
   if(iret.ne.nsctime) &
      call transp2imas_exit(' ?? BPEQ read error')
   call transp2imas_echo('BPEQ',scdata,1,nsctime)

   eq%time_slice(:)%global_quantities%beta_pol = scdata(:)

   write(iout,*) ' '
   !? Toroidal beta, defined as the volume-averaged total perpendicular pressure
   !? divided by (B0^2/(2*mu0)), i.e. beta_toroidal = 2 mu0 int(p dV) / V / B0^2
   call rpscalar('BTEQ',scdata,nsctime,iret,ier)
   if(ier.ne.0) call transp2imas_error('rpscalar',ier)
   if(iret.ne.nsctime) &
      call transp2imas_exit(' ?? BTEQ read error')
   call transp2imas_echo('BTEQ',scdata,1,nsctime)

   eq%time_slice(:)%global_quantities%beta_tor = scdata(:)

!      write(iout,*) ' '
!      call rpscalar('X',scdata,nsctime,iret,ier)
!      if(ier.ne.0) call transp2imas_error('rpscalar',ier)
!      if(iret.ne.nsctime) &
!         call transp2imas_exit(' ?? X read error')
!      call transp2imas_echo('X',scdata,1,nsctime)

   ! fill es IDS

   call rpscalar('GASD',scdata,nsctime,iret,ier)
   if ((ier .eq. 0) .and. (iret .eq. nsctime)) then
      do n = 1, n_thi
         !write(*,*) 'Writing GASD: n, a(n), z_n(n) =', n, cp%profiles_1d(nsctime)%ion(n)%element(1)%a, &
         !   cp%profiles_1d(nsctime)%ion(n)%element(1)%z_n
         if ((2.0 .eq. cp%profiles_1d(nsctime)%ion(n)%element(1)%a) .and. &
             (1.0 .eq. cp%profiles_1d(nsctime)%ion(n)%element(1)%z_n)) then
            write(*,*) 'Using index', n, 'for GASD data'
            !allocate(es%source(1)%ggd(1)%neutral(n)%particles(nsctime))
            !do k = 1, nsctime
            !   allocate(es%source(1)%ggd(1)%neutral(n)%particles(k)%values(1))
            !   es%source(1)%ggd(1)%neutral(n)%particles(k)%values = scdata(k)
            !end do
            allocate(es%source(1)%ggd(1)%neutral(n)%particles(1))
            allocate(es%source(1)%ggd(1)%neutral(n)%particles(1)%values(nsctime))
            es%source(1)%ggd(1)%neutral(n)%particles(1)%values = scdata
            allocate(es%source(1)%ggd(1)%neutral(n)%label(1))
            es%source(1)%ggd(1)%neutral(n)%label = 'D'
         endif
      enddo
   endif

   call rpscalar('GAST',scdata,nsctime,iret,ier)
   if ((ier .eq. 0) .and. (iret .eq. nsctime)) then
      do n = 1, n_thi
         !write(*,*) 'Writing GAST: n, a(n), z_n(n) =', n, cp%profiles_1d(nsctime)%ion(n)%element(1)%a, &
         !   cp%profiles_1d(nsctime)%ion(n)%element(1)%z_n
         if ((3.0 .eq. cp%profiles_1d(nsctime)%ion(n)%element(1)%a) .and. &
             (1.0 .eq. cp%profiles_1d(nsctime)%ion(n)%element(1)%z_n)) then
            write(*,*) 'Using index', n, 'for GAST data'
            allocate(es%source(1)%ggd(1)%neutral(n)%particles(1))
            allocate(es%source(1)%ggd(1)%neutral(n)%particles(1)%values(nsctime))
            es%source(1)%ggd(1)%neutral(n)%particles(1)%values = scdata
            allocate(es%source(1)%ggd(1)%neutral(n)%label(1))
            es%source(1)%ggd(1)%neutral(n)%label = 'T'
         endif
      enddo
   endif

   call rpscalar('GAS4',scdata,nsctime,iret,ier)
   if ((ier .eq. 0) .and. (iret .eq. nsctime)) then
      do n = 1, n_thi
         !write(*,*) 'Writing GAS4: n, a(n), z_n(n) =', n, cp%profiles_1d(nsctime)%ion(n)%element(1)%a, &
         !   cp%profiles_1d(nsctime)%ion(n)%element(1)%z_n
         if ((4.0 .eq. cp%profiles_1d(nsctime)%ion(n)%element(1)%a) .and. &
             (2.0 .eq. cp%profiles_1d(nsctime)%ion(n)%element(1)%z_n)) then
            write(*,*) 'Using index', n, 'for GAS4 data'
            allocate(es%source(1)%ggd(1)%neutral(n)%particles(1))
            allocate(es%source(1)%ggd(1)%neutral(n)%particles(1)%values(nsctime))
            es%source(1)%ggd(1)%neutral(n)%particles(1)%values = scdata
            allocate(es%source(1)%ggd(1)%neutral(n)%label(1))
            es%source(1)%ggd(1)%neutral(n)%label = '4'
         endif
      enddo
   endif

   call rpscalar('RCYD',scdata,nsctime,iret,ier)
   if ((ier .eq. 0) .and. (iret .eq. nsctime)) then
      do n = 1, n_thi
         if ((2.0 .eq. cp%profiles_1d(nsctime)%ion(n)%element(1)%a) .and. &
             (1.0 .eq. cp%profiles_1d(nsctime)%ion(n)%element(1)%z_n)) then
            write(*,*) 'Using index', n, 'for RCYD data'
            allocate(es%source(2)%ggd(1)%neutral(n)%particles(1))
            allocate(es%source(2)%ggd(1)%neutral(n)%particles(1)%values(nsctime))
            es%source(2)%ggd(1)%neutral(n)%particles(1)%values = scdata
            allocate(es%source(2)%ggd(1)%neutral(n)%label(1))
            es%source(2)%ggd(1)%neutral(n)%label = 'D'
         endif
      enddo
   endif

   call rpscalar('RCYT',scdata,nsctime,iret,ier)
   if ((ier .eq. 0) .and. (iret .eq. nsctime)) then
      do n = 1, n_thi
         if ((3.0 .eq. cp%profiles_1d(nsctime)%ion(n)%element(1)%a) .and. &
             (1.0 .eq. cp%profiles_1d(nsctime)%ion(n)%element(1)%z_n)) then
            write(*,*) 'Using index', n, 'for RCYT data'
            allocate(es%source(2)%ggd(1)%neutral(n)%particles(1))
            allocate(es%source(2)%ggd(1)%neutral(n)%particles(1)%values(nsctime))
            es%source(2)%ggd(1)%neutral(n)%particles(1)%values = scdata
            allocate(es%source(2)%ggd(1)%neutral(n)%label(1))
            es%source(2)%ggd(1)%neutral(n)%label = 'T'
         endif
      enddo
   endif

   call rpscalar('RCY4',scdata,nsctime,iret,ier)
   if ((ier .eq. 0) .and. (iret .eq. nsctime)) then
      do n = 1, n_thi
         if ((4.0 .eq. cp%profiles_1d(nsctime)%ion(n)%element(1)%a) .and. &
             (2.0 .eq. cp%profiles_1d(nsctime)%ion(n)%element(1)%z_n)) then
            write(*,*) 'Using index', n, 'for RCY4 data'
            allocate(es%source(2)%ggd(1)%neutral(n)%particles(1))
            allocate(es%source(2)%ggd(1)%neutral(n)%particles(1)%values(nsctime))
            es%source(2)%ggd(1)%neutral(n)%particles(1)%values = scdata
            allocate(es%source(2)%ggd(1)%neutral(n)%label(1))
            es%source(2)%ggd(1)%neutral(n)%label = '4'
         endif
      enddo
   endif

   open(unit=333, file='Non_IMAS_data_for_NUBEAM.txt', status='replace')

   call rpmulti('DENS0ARC',istype,label,units,infuns,mgsigns,mgnames,ier)
   if(ier.ne.0) call transp2imas_error('rpmulti(DENS0ARC)',ier)
   call rpmg0cal(mgnames,mgsigns,infuns,mgdata,nmax,iret,istype,iwarn,ier)
   if(ier.ne.0) call transp2imas_error('rpmg0cal ier',ier)
   if(iwarn.ne.0) call transp2imas_error('rpmg0cal iwarn',iwarn)
   write(333,*) infuns, nprtime, xsizes(istype)
   write(333,*) 'DENS0ARC: ', label, ' in ', units
   do i = 1, infuns
      !write(*,*) i, mgnames(i), mgnames(i)(6:6)
      write(333,*) mgnames(i)
      do j = 1, nprtime
         write(333,*) j
         do k = 1, xsizes(istype)
            rdum = 0.0
            do n = 1, 5
               if (associated(es%source(2)%ggd(1)%neutral(n)%label)) then
                  !write(*,*) k, n, trim(es%source(2)%ggd(1)%neutral(n)%label(1)), trim(mgnames(i)(6:6))
                  if (trim(es%source(2)%ggd(1)%neutral(n)%label(1)) == trim(mgnames(i)(6:6))) then
                     if (j == 1 .and. k == 1) write(*,*) n, es%source(2)%ggd(1)%neutral(n)%label, mgnames(i)
                     rdum = es%source(2)%ggd(1)%neutral(n)%particles(1)%values(k)
                  endif
               endif
            enddo
            write(333,*) mgdata((j-1)*xsizes(istype)+k, i) / rdum * 1.0e6 ! /cm3 to /m3
            !write(333,*) mgdata((j-1)*xsizes(istype)+k, i) / rdum * 1.0e6, mgdata((j-1)*xsizes(istype)+k, i), rdum
         enddo
      enddo
   enddo

   call rpmulti('DENS0AGF',istype,label,units,infuns,mgsigns,mgnames,ier)
   if(ier.ne.0) call transp2imas_error('rpmulti(DENS0AGF)',ier)
   write(333,*) 'DENS0AGF: ', label, ' in ', units
   call rpmg0cal(mgnames,mgsigns,infuns,mgdata,nmax,iret,istype,iwarn,ier)
   if(ier.ne.0) call transp2imas_error('rpmg0cal ier',ier)
   if(iwarn.ne.0) call transp2imas_error('rpmg0cal iwarn',iwarn)
   do i = 1, infuns
      write(333,*) mgnames(i)
      do j = 1, nprtime
         write(333,*) j
         do k = 1, xsizes(istype)
            rdum = 0.0
            do n = 1, 5
               if (associated(es%source(1)%ggd(1)%neutral(n)%label)) then
                  if (trim(es%source(1)%ggd(1)%neutral(n)%label(1)) == trim(mgnames(i)(6:6))) then
                     if (j == 1 .and. k == 1) write(*,*) n, es%source(1)%ggd(1)%neutral(n)%label, mgnames(i)
                     rdum = es%source(1)%ggd(1)%neutral(n)%particles(1)%values(k)
                  endif
               endif
            enddo
            write(333,*) mgdata((j-1)*xsizes(istype)+k, i) / rdum * 1.0e6 ! /cm3 to /m3
         enddo
      enddo
   enddo

   call rpmulti('T0ARC',istype,label,units,infuns,mgsigns,mgnames,ier)
   if(ier.ne.0) call transp2imas_error('rpmulti(T0ARC)',ier)
   write(333,*) 'T0ARC: ', label, ' in ', units
   call rpmg0cal(mgnames,mgsigns,infuns,mgdata,nmax,iret,istype,iwarn,ier)
   if(ier.ne.0) call transp2imas_error('rpmg0cal ier',ier)
   if(iwarn.ne.0) call transp2imas_error('rpmg0cal iwarn',iwarn)
   do i = 1, infuns
      write(333,*) mgnames(i)
      do j = 1, nprtime
         write(333,*) j
         do k = 1, xsizes(istype)
            write(333,*) mgdata((j-1)*xsizes(istype)+k, i) * 1.0e-3 ! eV to keV
         enddo
      enddo
   enddo

   call rpmulti('T0AGF',istype,label,units,infuns,mgsigns,mgnames,ier)
   if(ier.ne.0) call transp2imas_error('rpmulti(T0AGF)',ier)
   write(333,*) 'T0AGF: ', label, ' in ', units
   call rpmg0cal(mgnames,mgsigns,infuns,mgdata,nmax,iret,istype,iwarn,ier)
   if(ier.ne.0) call transp2imas_error('rpmg0cal ier',ier)
   if(iwarn.ne.0) call transp2imas_error('rpmg0cal iwarn',iwarn)
   do i = 1, infuns
      write(333,*) mgnames(i)
      do j = 1, nprtime
         write(333,*) j
         do k = 1, xsizes(istype)
            write(333,*) mgdata((j-1)*xsizes(istype)+k, i) * 1.0e-3 ! eV to keV
         enddo
      enddo
   enddo

   call rpmulti('OMEG0ARC',istype,label,units,infuns,mgsigns,mgnames,ier)
   if(ier.ne.0) call transp2imas_error('rpmulti(OMEG0ARC)',ier)
   write(333,*) 'OMEG0ARC: ', label, ' in ', units
   call rpmg0cal(mgnames,mgsigns,infuns,mgdata,nmax,iret,istype,iwarn,ier)
   if(ier.ne.0) call transp2imas_error('rpmg0cal ier',ier)
   if(iwarn.ne.0) call transp2imas_error('rpmg0cal iwarn',iwarn)
   do i = 1, infuns
      write(333,*) mgnames(i)
      do j = 1, nprtime
         write(333,*) j
         do k = 1, xsizes(istype)
            write(333,*) mgdata((j-1)*xsizes(istype)+k, i)
         enddo
      enddo
   enddo

   call rpmulti('OMEG0AGF',istype,label,units,infuns,mgsigns,mgnames,ier)
   if(ier.ne.0) call transp2imas_error('rpmulti(OMEG0AGF)',ier)
   write(333,*) 'OMEG0AGF: ', label, ' in ', units
   call rpmg0cal(mgnames,mgsigns,infuns,mgdata,nmax,iret,istype,iwarn,ier)
   if(ier.ne.0) call transp2imas_error('rpmg0cal ier',ier)
   if(iwarn.ne.0) call transp2imas_error('rpmg0cal iwarn',iwarn)
   do i = 1, infuns
      write(333,*) mgnames(i)
      do j = 1, nprtime
         write(333,*) j
         do k = 1, xsizes(istype)
            write(333,*) mgdata((j-1)*xsizes(istype)+k, i)
         enddo
      enddo
   enddo

   close(333)

   ! fill nbi IDS

   if (nbeam.gt.0) then

      do k = 1, nbeam
         write(iout,*) ' '
         !? Beam power
         write(int2strng, *) k
         tmps1 = adjustl(int2strng) ! delete the leading blanks
         tmpstrng = 'PINJ0'//tmps1(1:1)
         call rpscalar(trim(tmpstrng),scdata,nsctime,iret,ier)
         if (ier.ne.0) call transp2imas_error('rpscalar',ier)
         if (iret.ne.nsctime) &
            call transp2imas_exit(' ?? '//trim(tmpstrng)//' read error')
         call transp2imas_echo(trim(tmpstrng),scdata,1,nsctime)
         nbi%unit(k)%power%data = scdata(:)
      end do

      do k = 1, nbeam
         write(iout,*) ' '
         !? Beam energy
         write(int2strng, *) k
         tmps1 = adjustl(int2strng) ! delete the leading blanks
         tmpstrng = 'EINJ0'//tmps1(1:1)//'_E1'
         !write(*, *) 'x', trim(tmpstrng), 'x'
         !stop
         call rpscalar(trim(tmpstrng),scdata,nsctime,iret,ier)
         if (ier.ne.0) call transp2imas_error('rpscalar',ier)
         if (iret.ne.nsctime) &
            call transp2imas_exit(' ?? '//trim(tmpstrng)//' read error')
         call transp2imas_echo(trim(tmpstrng),scdata,1,nsctime)
         nbi%unit(k)%energy%data = scdata(:)
      end do

   end if

!
! rprofile get profile data f(X,t)
!

   ! fill core_profiles ids

   ! first define the profile X coordinate
   ! X (formerly called XI) zone center
   ! XB zone boundary (excluding magnetic axis)

   write(iout,*) ' '
   call rprofile('X',prdata,nprtime*xsizes(1),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(XI)',ier)
   if(iret.ne.nprtime*xsizes(1)) &
      call transp2imas_exit(' ?? XI read error')
   call transp2imas_echo('XI',prdata,xsizes(1),nprtime)

   offset=xsizes(1)
   allocate(XI(offset,nprtime))
   do it=1,nprtime
      XI(1:offset,it) = prdata(1+(it-1)*offset:it*offset)
      allocate(cp%profiles_1d(it)%grid%rho_tor_norm(offset))
      cp%profiles_1d(it)%grid%rho_tor_norm(1:offset) = XI(1:offset,it) * XI(1:offset,it)
   enddo

   write(iout,*) ' '
   call rprofile('XB',prdata,nprtime*xsizes(2),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(XB)',ier)
   if(iret.ne.nprtime*xsizes(2)) &
      call transp2imas_exit(' ?? XB read error')
   call transp2imas_echo('XB',prdata,xsizes(2),nprtime)

   offset=xsizes(2)
   allocate(XB(offset,nprtime))
   do it=1,nprtime
      XB(1:offset,it) = prdata(1+(it-1)*offset:it*offset)
   enddo

   write(iout,*) ' '
   call rprofile('NE',prdata,nprtime*xsizes(1),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(NE)',ier)
   if(iret.ne.nprtime*xsizes(1)) &
      call transp2imas_exit(' ?? NE read error')
   call transp2imas_echo('NE',prdata,xsizes(1),nprtime)

   offset=xsizes(1)
   do it=1,nprtime
      allocate(cp%profiles_1d(it)%electrons%density(offset))
      cp%profiles_1d(it)%electrons%density(1:offset)=&
         prdata(1+(it-1)*offset:it*offset)*1.e6
   enddo

   write(iout,*) ' '
   call rprofile('TE',prdata,nprtime*xsizes(1),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(TE)',ier)
   if(iret.ne.nprtime*xsizes(1)) &
      call transp2imas_exit(' ?? TE read error')
   call transp2imas_echo('TE',prdata,xsizes(1),nprtime)

   offset=xsizes(1)
   do it=1,nprtime
      allocate(cp%profiles_1d(it)%electrons%temperature(offset))
      cp%profiles_1d(it)%electrons%temperature(1:offset)=&
         prdata(1+(it-1)*offset:it*offset)
   enddo

   write(iout,*) ' '
   call rprofile('OMEGA', prdata, nprtime * xsizes(1), iret, ier)
   !if (ier.ne.0) call transp2imas_error('rprofile(OMEGA)',ier)
   !if (iret .ne. nprtime * xsizes(1)) &
   !   call transp2imas_exit(' ?? OMEGA read error')
   if (ier .eq. 0) then
      if (iret .eq. nprtime * xsizes(1)) then
         call transp2imas_echo('OMEGA', prdata, xsizes(1), nprtime)
         offset = xsizes(1)
         do it = 1, nprtime
            allocate(cp%profiles_1d(it)%electrons%velocity_tor(offset))
            cp%profiles_1d(it)%electrons%velocity_tor(1:offset) = &
               eq%vacuum_toroidal_field%r0 * prdata(1+(it-1)*offset:it*offset)
            ! Assume ions rotate with electrons
            do iion = 1, n_thi + nlist
               allocate(cp%profiles_1d(it)%ion(iion)%velocity_tor(offset))
               cp%profiles_1d(it)%ion(iion)%velocity_tor(1:offset) = &
               cp%profiles_1d(it)%electrons%velocity_tor(1:offset)
            enddo
         enddo
      endif
   endif

   write(iout,*) ' '
   call rprofile('NI',prdata,nprtime*xsizes(1),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(NI)',ier)
   if(iret.ne.nprtime*xsizes(1)) &
      call transp2imas_exit(' ?? NI read error')
   call transp2imas_echo('NI',prdata,xsizes(1),nprtime)
   ! ids doesn't have profile for NI, only for ni_over_ne
   write(iout,*) ' '
   call rpcalc('NI/NE',prdata,nprtime*xsizes(1),iret, istype, &
      iwarn,ier)
   if(ier.ne.0) call transp2imas_error('rpcalc(NI/NE) ier',ier)
   if(iwarn.ne.0) call transp2imas_error('rpcalc(NI/NE) iwarn',iwarn)
   write(iout,*) ' rpcalc call OK, returned istype = ',istype
   if(iret.ne.nprtime*xsizes(1)) &
      call transp2imas_exit(' ?? NI/NE calc size result error')
   call transp2imas_echo('NI/NE (rpcalc)',prdata,xsizes(1),nprtime)

   offset=xsizes(1)
   do it=1,nprtime
      allocate(cp%profiles_1d(it)%n_i_total_over_n_e(offset))
      cp%profiles_1d(it)%n_i_total_over_n_e(1:offset)=&
         prdata(1+(it-1)*offset:it*offset)
   enddo

   write(iout,*) ' '
   call rprofile('TI',prdata,nprtime*xsizes(1),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(TI)',ier)
   if(iret.ne.nprtime*xsizes(1)) &
      call transp2imas_exit(' ?? TI read error')
   call transp2imas_echo('TI',prdata,xsizes(1),nprtime)

   offset=xsizes(1)
   do it=1,nprtime
      allocate(cp%profiles_1d(it)%t_i_average(offset))
      cp%profiles_1d(it)%t_i_average(1:offset)=&
         prdata(1+(it-1)*offset:it*offset)
   enddo

   write(iout,*) ' '
   ! ZEFMD                MAGDIF ZEFF PROFILE
   ! ZEFFP                PLASMA COMPOSITION ZEFF PROFILE
   ! ZEFFI                ZEFF DATA (UNCONSTRAINED)
   call rprofile('ZEFFP',prdata,nprtime*xsizes(1),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(ZEFFP)',ier)
   if(iret.ne.nprtime*xsizes(1)) &
      call transp2imas_exit(' ?? ZEFFP read error')
   call transp2imas_echo('ZEFFP',prdata,xsizes(1),nprtime)

   offset=xsizes(1)
   do it=1,nprtime
      allocate(cp%profiles_1d(it)%zeff(offset))
      cp%profiles_1d(it)%zeff(1:offset)=&
         prdata(1+(it-1)*offset:it*offset)
   enddo

   write(iout,*) ' '
   call rprofile('Q',prdata,nprtime*xsizes(2),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(Q)',ier)
   if(iret.ne.nprtime*xsizes(1)) &
      call transp2imas_exit(' ?? Q read error')
   call transp2imas_echo('Q',prdata,xsizes(1),nprtime)

   ! Use linear interpolation to transform Q(XB) -> Q(X),
   ! from zone boundary (XB) to zone center (X)
   offset = xsizes(1)
   do it = 1, nprtime
      allocate(cp%profiles_1d(it)%q(offset))
      cp%profiles_1d(it)%q(1:offset) = &
         0.5 * prdata(1+(it-1)*offset:it*offset)
      do ir = offset, 2, -1
         cp%profiles_1d(it)%q(ir) = &
         cp%profiles_1d(it)%q(ir) + &
         cp%profiles_1d(it)%q(ir-1)
      enddo
      cp%profiles_1d(it)%q(1) = &
      cp%profiles_1d(it)%q(1) + &
      0.5 * eq%time_slice(it)%global_quantities%q_axis
   enddo

   write(iout,*) ' '
   call rprofile('CURBS',prdata,nprtime*xsizes(1),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(CURBS)',ier)
   if(iret.ne.nprtime*xsizes(1)) &
      call transp2imas_exit(' ?? CURBS read error')
   call transp2imas_echo('CURBS',prdata,xsizes(1),nprtime)

   offset=xsizes(1)
   do it=1,nprtime
      allocate(cp%profiles_1d(it)%j_bootstrap(offset))
      cp%profiles_1d(it)%j_bootstrap(1:offset)=&
         prdata(1+(it-1)*offset:it*offset) * 1.e4
   enddo

   write(iout,*) ' '
   call rprofile('CUROH',prdata,nprtime*xsizes(1),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(CUROH)',ier)
   if(iret.ne.nprtime*xsizes(1)) &
      call transp2imas_exit(' ?? CUROH read error')
   call transp2imas_echo('CUROH',prdata,xsizes(1),nprtime)

   offset=xsizes(1)
   do it=1,nprtime
      allocate(cp%profiles_1d(it)%j_ohmic(offset))
      cp%profiles_1d(it)%j_ohmic(1:offset)=&
         prdata(1+(it-1)*offset:it*offset) * 1.e4
   enddo
#if 1
   write(iout, *) ' '
   call rprofile('VRPOT', prdata, nprtime * xsizes(2), iret, ier)
   !if (ier .ne. 0) call transp2imas_error('rprofile(VRPOT)', ier)
   !if (iret .ne. nprtime * xsizes(2)) &
   !   call transp2imas_exit(' ?? VRPOT read error')
   if (ier .eq. 0) then
      if (iret .eq. nprtime * xsizes(2)) then
         call transp2imas_echo('VRPOT', prdata, xsizes(2), nprtime)
         ! Use linear differentiation to transform VRPOT(XB) -> "ERAD(X)"
         offset = xsizes(1)
         do it = 1, nprtime
            allocate(cp%profiles_1d(it)%e_field%radial(offset))
            do ir = 2, offset
               cp%profiles_1d(it)%e_field%radial(ir) = &
               prdata(ir+(it-1)*offset) - prdata(ir-1+(it-1)*offset)
            enddo
         enddo
         write(iout,*) ' '
         call rprofile('RBOUN', prdata, nprtime * xsizes(2), iret, ier)
         if (ier .ne. 0) call transp2imas_error('rprofile(RBOUN)', ier)
         if (iret .ne. nprtime * xsizes(2)) &
            call transp2imas_exit(' ?? RBOUN read error')
         call transp2imas_echo('RBOUN', prdata, xsizes(2), nprtime)
         offset = xsizes(1)
         do it = 1, nprtime
            do ir = 2, offset
               cp%profiles_1d(it)%e_field%radial(ir) = &
               cp%profiles_1d(it)%e_field%radial(ir) / &
               (prdata(ir+(it-1)*offset) - prdata(ir-1+(it-1)*offset)) * &
               1.0e2 ! V/cm -> V/m
            enddo
            ! Assume on-axis radial electric field is zero
            cp%profiles_1d(it)%e_field%radial(1) = &
            0.5 * cp%profiles_1d(it)%e_field%radial(2)
         enddo
      endif
   endif
#endif
   write(iout,*) ' '
   !? correct?
   ! CUR(X)            TOTAL PLASMA CURRENT         AMPS/CM2
   ! CURGP(X)          GRAD(P) TOROIDAL CUR         AMPS/CM2
   ! PLCURPLL(XB)      POLOIDAL CUR (J PLL)         AMPS
   ! PLCURPRP(XB)      POLOIDAL CUR (J PERP)        AMPS
   ! PLCURTOT(XB)      TOTAL POLOIDAL CUR TO WALL   AMPS
   !call rprofile('CURGP',prdata,nprtime*xsizes(1),iret,ier)
   call rprofile('CUR',prdata,nprtime*xsizes(1),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(CUR)',ier)
   if(iret.ne.nprtime*xsizes(1)) &
      call transp2imas_exit(' ?? CUR read error')
   call transp2imas_echo('CUR',prdata,xsizes(1),nprtime)

   offset=xsizes(1)
   do it=1,nprtime
      allocate(cp%profiles_1d(it)%j_tor(offset))
      cp%profiles_1d(it)%j_tor(1:offset)=&
         prdata(1+(it-1)*offset:it*offset) * 1.e4
   enddo

   write(iout,*) ' '
   call rprofile('SHAT',prdata,nprtime*xsizes(1),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(SHAT)',ier)
   if(iret.ne.nprtime*xsizes(1)) &
      call transp2imas_exit(' ?? SHAT read error')
   call transp2imas_echo('SHAT',prdata,xsizes(1),nprtime)

   offset=xsizes(1)
   do it=1,nprtime
      allocate(cp%profiles_1d(it)%magnetic_shear(offset))
      cp%profiles_1d(it)%magnetic_shear(1:offset)=&
         prdata(1+(it-1)*offset:it*offset)
   enddo

   write(iout,*) ' '
   call rprofile('PLFLX2PI',prdata,nprtime*xsizes(2),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(PLFLX2PI)',ier)
   if(iret.ne.nprtime*xsizes(2)) &
      call transp2imas_exit(' ?? PLFLX2PI read error')
   call transp2imas_echo('PLFLX2PI',prdata,xsizes(2),nprtime)

   ! Use linear interpolation to transform PLFLX2PI(XB) -> PLFLX2PI(X),
   ! from zone boundary (XB) to zone center (X)
   offset = xsizes(1)
   do it = 1, nprtime
      allocate(cp%profiles_1d(it)%grid%psi(offset))
      cp%profiles_1d(it)%grid%psi(1:offset) = &
         0.5 * prdata(1+(it-1)*offset:it*offset)
      do ir = offset, 2, -1
         cp%profiles_1d(it)%grid%psi(ir) = &
         cp%profiles_1d(it)%grid%psi(ir) + &
         cp%profiles_1d(it)%grid%psi(ir-1)
      enddo
      cp%profiles_1d(it)%grid%psi(1) = &
      cp%profiles_1d(it)%grid%psi(1) + &
      0.5 * eq%time_slice(it)%global_quantities%psi_axis
   enddo

   write(iout,*) ' '
   call rprofile('TRFLX',prdata,nprtime*xsizes(2),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(TRFLX)',ier)
   if(iret.ne.nprtime*xsizes(2)) &
      call transp2imas_exit(' ?? TRFLX read error')
   call transp2imas_echo('TRFLX',prdata,xsizes(2),nprtime)

   ! Use linear interpolation to transform TRFLX(XB) -> TRFLX(X),
   ! from zone boundary (XB) to zone center (X)
   offset = xsizes(1)
   do it = 1, nprtime
      allocate(cp%profiles_1d(it)%grid%rho_tor(offset))
      cp%profiles_1d(it)%grid%rho_tor(1:offset) = &
         0.5 * prdata(1+(it-1)*offset:it*offset)
      do ir = offset, 2, -1
         cp%profiles_1d(it)%grid%rho_tor(ir) = &
         cp%profiles_1d(it)%grid%rho_tor(ir) + &
         cp%profiles_1d(it)%grid%rho_tor(ir-1)
      enddo
      cp%profiles_1d(it)%grid%rho_tor(1) = &
      cp%profiles_1d(it)%grid%rho_tor(1) + &
      0.5 * 0.0
   enddo

   write(iout,*) ' '
   call rprofile('DVOL',prdata,nprtime*xsizes(1),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(DVOL)',ier)
   if(iret.ne.nprtime*xsizes(1)) &
      call transp2imas_exit(' ?? DVOL read error')
   call transp2imas_echo('DVOL',prdata,xsizes(1),nprtime)

   offset=xsizes(1)
   do it=1,nprtime
      allocate(cp%profiles_1d(it)%grid%volume(offset))
      cp%profiles_1d(it)%grid%volume(1:offset)= &
         prdata(1+(it-1)*offset:it*offset) * 1.e-6
   enddo

   write(iout,*) ' '
   call rprofile('DAREA',prdata,nprtime*xsizes(1),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(DAREA)',ier)
   if(iret.ne.nprtime*xsizes(1)) &
      call transp2imas_exit(' ?? DAREA read error')
   call transp2imas_echo('DAREA',prdata,xsizes(1),nprtime)

   offset=xsizes(1)
   do it=1,nprtime
      allocate(cp%profiles_1d(it)%grid%area(offset))
      cp%profiles_1d(it)%grid%area(1:offset)= &
         prdata(1+(it-1)*offset:it*offset) * 1.e-4
   enddo

   !write(iout,*) ' '
   !call rprofile('NT',prdata,nprtime*xsizes(1),iret,ier)
   !if(ier.ne.0) call transp2imas_error('rprofile(NT)',ier)
   !if(iret.ne.nprtime*xsizes(1)) &
   !   call transp2imas_exit(' ?? NT read error')
   !call transp2imas_echo('NT',prdata,xsizes(1),nprtime)

   !write(iout,*) ' '
   !call rprofile('ND',prdata,nprtime*xsizes(1),iret,ier)
   !if(ier.ne.0) call transp2imas_error('rprofile(ND)',ier)
   !if(iret.ne.nprtime*xsizes(1)) &
   !   call transp2imas_exit(' ?? ND read error')
   !call transp2imas_echo('ND',prdata,xsizes(1),nprtime)

   !write(iout,*) ' '
   !call rprofile('NH',prdata,nprtime*xsizes(1),iret,ier)
   !if(ier.ne.0) call transp2imas_error('rprofile(NH)',ier)
   !if(iret.ne.nprtime*xsizes(1)) &
   !   call transp2imas_exit(' ?? NH read error')
   !call transp2imas_echo('NH',prdata,xsizes(1),nprtime)

   !write(iout,*) ' '
   !call rprofile('NHE3',prdata,nprtime*xsizes(1),iret,ier)
   !if(ier.ne.0) call transp2imas_error('rprofile(NHE3)',ier)
   !if(iret.ne.nprtime*xsizes(1)) &
   !   call transp2imas_exit(' ?? NHE3 read error')
   !call transp2imas_echo('NHE3',prdata,xsizes(1),nprtime)

   !write(iout,*) ' '
   !call rprofile('NHE4',prdata,nprtime*xsizes(1),iret,ier)
   !if(ier.ne.0) call transp2imas_error('rprofile(NHE4)',ier)
   !if(iret.ne.nprtime*xsizes(1)) &
   !   call transp2imas_exit(' ?? NHE4 read error')
   !call transp2imas_echo('NHE4',prdata,xsizes(1),nprtime)

   !write(iout,*) ' '
   !call rprofile('NLITH',prdata,nprtime*xsizes(1),iret,ier)
   !if(ier.ne.0) call transp2imas_error('rprofile(NLITH)',ier)
   !if(iret.ne.nprtime*xsizes(1)) &
   !   call transp2imas_exit(' ?? NLITH read error')
   !call transp2imas_echo('NLITH',prdata,xsizes(1),nprtime)

   ! First set all ion diffusivities to zero
   do iion = 1, nion
      do it = 1, nprtime
         allocate(ct%model(1)%profiles_1d(it)%ion(iion)%particles%d(offset))
         ct%model(1)%profiles_1d(it)%ion(iion)%particles%d(1:offset) = 0.0
      enddo
   enddo

   ! Next look for actual ion-diffusivity data
   ! First for beam deuterons...
   write(iout,*) ' '
   call rprofile('BDIFBX_D', prdata, nprtime * xsizes(2), iret, ier)
   if (ier .eq. 0) then
      if (iret .ne. nprtime * xsizes(2)) &
         call transp2imas_exit(' ?? BDIFBX_D read error')
      call transp2imas_echo('BDIFBX_D', prdata, xsizes(2), nprtime)
      offset = xsizes(2)
      do it = 1, nprtime
            do iion = 1, nion
               !write(*,*) nion, iion, cp%profiles_1d(it)%ion(iion)%element(1)%a, &
               !cp%profiles_1d(it)%ion(iion)%z_ion, cp%profiles_1d(it)%ion(iion)%element(1)%z_n
               if ((2.0 .eq. cp%profiles_1d(it)%ion(iion)%element(1)%a) .and. &
                   (1.0 .eq. cp%profiles_1d(it)%ion(iion)%z_ion)        .and. &
                   (1.0 .eq. cp%profiles_1d(it)%ion(iion)%element(1)%z_n)) &
                      ct%model(1)%profiles_1d(it)%ion(iion)%particles%d(1:offset) = &
                      prdata(1+(it-1)*offset:it*offset) * 1.0e-4 ! cm^2/s -> m^2/s
            enddo
      enddo
   endif
   ! Then for fusion alphas...
   write(iout,*) ' '
   call rprofile('FDIFBX_4', prdata, nprtime * xsizes(2), iret, ier)
   if (ier .eq. 0) then
      if (iret .ne. nprtime * xsizes(2)) &
         call transp2imas_exit(' ?? FDIFBX_4 read error')
      call transp2imas_echo('FDIFBX_4', prdata, xsizes(2), nprtime)
      offset = xsizes(2)
      do it = 1, nprtime
            do iion = 1, nion
               if ((4.0 .eq. cp%profiles_1d(it)%ion(iion)%element(1)%a) .and. &
                   (1.0 .eq. cp%profiles_1d(it)%ion(iion)%z_ion)        .and. &
                   (2.0 .eq. cp%profiles_1d(it)%ion(iion)%element(1)%z_n)) &
                      ct%model(1)%profiles_1d(it)%ion(iion)%particles%d(1:offset) = &
                      prdata(1+(it-1)*offset:it*offset) * 1.0e-4 ! cm^2/s -> m^2/s
            enddo
      enddo
   endif

   ! fill ids_equilibrium

   write(iout,*) ' '
   call rprofile('PLFLX2PI',prdata,nprtime*xsizes(2),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(PLFLX2PI)',ier)
   if(iret.ne.nprtime*xsizes(2)) &
      call transp2imas_exit(' ?? PLFLX2PI read error')
   call transp2imas_echo('PLFLX2PI',prdata,xsizes(2),nprtime)

   offset=xsizes(2)
   ! Calculate dXB/dPLFLX, it will be used later to calculate p' and FF'
   allocate(PLFLX(offset, nprtime), dXBdPLFLX(offset, nprtime))
   do it = 1, nprtime
      PLFLX(:, it) = prdata(1+(it-1)*offset:it*offset)
      allocate(eq%time_slice(it)%profiles_1d%psi(offset))
      eq%time_slice(it)%profiles_1d%psi(1:offset) = PLFLX(:, it)
      ! Make psi a STRICTLY monotonic profile
      do i = 1, offset
         eq%time_slice(it)%profiles_1d%psi(i) = &
         eq%time_slice(it)%profiles_1d%psi(i) * (1.0 + 1.0e-6 * (i - offset))
      enddo

      xbbuf1(:) = PLFLX(:, it)
      xbbuf2(:) = XB(:, it)

      call ezspline_init(spln1, offset, (/0, 0/), ier)
      call ezspline_error(ier)

      !spln1%bcval1min = 0.
      !spln1%isHermite = 1
      spln1%x1 = xbbuf1(1:offset)

      call ezspline_setup(spln1, xbbuf2(1:offset), ier)
      call ezspline_error(ier)

      ! This call might be needed to init derivative calculation below?
      !call ezspline_interp(spln1, offset, xbbuf1, xibuf2, ier)
      !call ezspline_error(ier)

      call ezspline_derivative(spln1, 1, offset, xbbuf1, xbbuf3, ier)
      call ezspline_error(ier)

      dXBdPLFLX(:, it) = xbbuf3(:)

      call ezspline_free(spln1, ier)
      call ezspline_error(ier)

      !write(312,*) 'poloidal flux at', it
      !write(312,*) 'size', size(xbbuf2), size(xibuf2), size(xibuf1), size(xbbuf1)
      !do ir=1,offset
      !   write(312,*) xbbuf1(ir), xibuf1(ir), xbbuf2(ir), xibuf2(ir), xibuf22(ir)
         !write(312,*) PLFLX(ir,it), dXBdPLFLX(ir,it)
      !enddo
   enddo

   write(iout,*) ' '
   call rprofile('TRFLX',prdata,nprtime*xsizes(2),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(TRFLX)',ier)
   if(iret.ne.nprtime*xsizes(2)) &
      call transp2imas_exit(' ?? TRFLX read error')
   call transp2imas_echo('TRFLX',prdata,xsizes(2),nprtime)

   offset = xsizes(2)
   do it = 1, nprtime
      allocate(eq%time_slice(it)%profiles_1d%phi(offset))
      eq%time_slice(it)%profiles_1d%phi(1:offset) = prdata(1+(it-1)*offset:it*offset)
   enddo

   write(iout,*) ' '
   call rprofile('PMHD_IN',prdata,nprtime*xsizes(1),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(PMHD_IN)',ier)
   if(iret.ne.nprtime*xsizes(1)) &
      call transp2imas_exit(' ?? PMHD_IN read error')
   call transp2imas_echo('PMHD_IN',prdata,xsizes(1),nprtime)

   offset = xsizes(1)

   ! Calculate p and p' = dp/dPLFLX on zone boundaries
   do it = 1, nprtime
      xibuf1(:) = XI(:, it)
      xibuf2(:) = prdata(1+(it-1)*offset:it*offset)
      xbbuf1(:) = XB(:, it)

      call ezspline_init(spln1, offset, (/0, 0/), ier)
      call ezspline_error(ier)

      !spln1%bcval1min = 0.
      !spln1%isHermite = 1
      spln1%x1 = xibuf1(1:offset)

      call ezspline_setup(spln1, xibuf2(1:offset), ier)
      call ezspline_error(ier)

      call ezspline_interp(spln1, offset, xbbuf1, xbbuf2, ier)
      ! xbbuf2 now has p on zone boundaries
      call ezspline_error(ier)
      allocate(eq%time_slice(it)%profiles_1d%pressure(offset))
      eq%time_slice(it)%profiles_1d%pressure(1:offset) = xbbuf2(1:offset)

      call ezspline_free(spln1, ier)
      call ezspline_error(ier)

      ! Calculate dp/dXB
      call ezspline_init(spln1, offset, (/0, 0/), ier)
      call ezspline_error(ier)

      spln1%x1 = xbbuf1(1:offset)

      call ezspline_setup(spln1, xbbuf2(1:offset), ier)
      call ezspline_error(ier)

      ! this call only made to init derivate calculation below (might be unnecessary)
      !call ezspline_interp(spln1, offset, xbbuf1, xibuf2, ier)
      !call ezspline_error(ier)

      call ezspline_derivative(spln1, 1, offset, xbbuf1, xbbuf3, ier)
      ! xbbuf3 now contains dp/dXB
      call ezspline_error(ier)

      ! dp/dPLFLX = dp/dXB * dXB/dPLFLX
      allocate(eq%time_slice(it)%profiles_1d%dpressure_dpsi(offset))
      eq%time_slice(it)%profiles_1d%dpressure_dpsi(1:offset) = xbbuf3(:) * dXBdPLFLX(:, it)
#if 0
      do ir = 2, offset-1
         eq%time_slice(it)%profiles_1d%dpressure_dpsi(ir) = &
         (eq%time_slice(it)%profiles_1d%pressure(ir+1) - eq%time_slice(it)%profiles_1d%pressure(ir-1)) / &
         (eq%time_slice(it)%profiles_1d%psi(ir+1) - eq%time_slice(it)%profiles_1d%psi(ir-1))
      enddo
      eq%time_slice(it)%profiles_1d%dpressure_dpsi(1) = eq%time_slice(it)%profiles_1d%dpressure_dpsi(2)
      eq%time_slice(it)%profiles_1d%dpressure_dpsi(offset) = eq%time_slice(it)%profiles_1d%dpressure_dpsi(offset-1)
#endif

      call ezspline_free(spln1, ier)
      call ezspline_error(ier)

      write(313, *) 'foobar, p, pprime at ', it
      do ir = 1, offset
         write(313, *) eq%time_slice(it)%profiles_1d%pressure(ir), &
                       eq%time_slice(it)%profiles_1d%dpressure_dpsi(ir)
      enddo
   enddo

   write(iout,*) ' '
   call rprofile('ELONG',prdata,nprtime*xsizes(2),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(ELONG)',ier)
   if(iret.ne.nprtime*xsizes(1)) &
      call transp2imas_exit(' ?? ELONG read error')
   call transp2imas_echo('ELONG',prdata,xsizes(1),nprtime)

   offset = xsizes(2)
   do it = 1, nprtime
      allocate(eq%time_slice(it)%profiles_1d%elongation(offset))
      eq%time_slice(it)%profiles_1d%elongation(1:offset) = prdata(1+(it-1)*offset:it*offset)
   enddo

   write(iout,*) ' '
   call rprofile('TRIANGU',prdata,nprtime*xsizes(2),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(TRIANGU)',ier)
   if(iret.ne.nprtime*xsizes(2)) &
      call transp2imas_exit(' ?? TRIANGU read error')
   call transp2imas_echo('TRIANGU',prdata,xsizes(2),nprtime)

   offset = xsizes(2)
   do it = 1, nprtime
      allocate(eq%time_slice(it)%profiles_1d%triangularity_upper(offset))
      eq%time_slice(it)%profiles_1d%triangularity_upper(1:offset) = prdata(1+(it-1)*offset:it*offset)
   enddo

   write(iout,*) ' '
   call rprofile('TRIANGL',prdata,nprtime*xsizes(2),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(TRIANGL)',ier)
   if(iret.ne.nprtime*xsizes(2)) &
      call transp2imas_exit(' ?? TRIANGL read error')
   call transp2imas_echo('TRIANGL',prdata,xsizes(2),nprtime)

   offset = xsizes(2)
   do it = 1, nprtime
      allocate(eq%time_slice(it)%profiles_1d%triangularity_lower(offset))
      eq%time_slice(it)%profiles_1d%triangularity_lower(1:offset) = prdata(1+(it-1)*offset:it*offset)
   enddo

   write(iout,*) ' '
   call rprofile('DVOL',prdata,nprtime*xsizes(1),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(DVOL)',ier)
   if(iret.ne.nprtime*xsizes(1)) &
      call transp2imas_exit(' ?? DVOL read error')
   call transp2imas_echo('DVOL',prdata,xsizes(1),nprtime)

   offset=xsizes(1)
   do it=1,nprtime
      allocate(eq%time_slice(it)%profiles_1d%volume(offset))
      eq%time_slice(it)%profiles_1d%volume(1:offset)= &
         prdata(1+(it-1)*offset:it*offset) * 1.e-6
   enddo

   write(iout,*) ' '
   call rprofile('DAREA',prdata,nprtime*xsizes(1),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(DAREA)',ier)
   if(iret.ne.nprtime*xsizes(1)) &
      call transp2imas_exit(' ?? DAREA read error')
   call transp2imas_echo('DAREA',prdata,xsizes(1),nprtime)

   offset=xsizes(1)
   do it=1,nprtime
      allocate(eq%time_slice(it)%profiles_1d%area(offset))
      eq%time_slice(it)%profiles_1d%area(1:offset)= &
         prdata(1+(it-1)*offset:it*offset) * 1.e-4
   enddo

   write(iout,*) ' '
   call rprofile('SURF',prdata,nprtime*xsizes(2),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(SURF)',ier)
   if(iret.ne.nprtime*xsizes(1)) &
      call transp2imas_exit(' ?? SURF read error')
   call transp2imas_echo('SURF',prdata,xsizes(1),nprtime)

   offset = xsizes(2)
   do it = 1, nprtime
      allocate(eq%time_slice(it)%profiles_1d%surface(offset))
      eq%time_slice(it)%profiles_1d%surface(1:offset) = prdata(1+(it-1)*offset:it*offset) * 1.0e-4
   enddo

   write(iout,*) ' '
   call rprofile('Q',prdata,nprtime*xsizes(2),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(Q)',ier)
   if(iret.ne.nprtime*xsizes(2)) &
      call transp2imas_exit(' ?? Q read error')
   call transp2imas_echo('Q',prdata,xsizes(2),nprtime)

   offset = xsizes(2)
   do it = 1, nprtime
      allocate(eq%time_slice(it)%profiles_1d%q(offset))
      eq%time_slice(it)%profiles_1d%q(1:offset) = prdata(1+(it-1)*offset:it*offset)
   enddo

   write(iout,*) ' '
   call rprofile('SHAT',prdata,nprtime*xsizes(1),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(SHAT)',ier)
   if(iret.ne.nprtime*xsizes(1)) &
      call transp2imas_exit(' ?? SHAT read error')
   call transp2imas_echo('SHAT',prdata,xsizes(1),nprtime)

   ! Use linear interpolation to transform SHAT(X) -> SHAT(XB),
   ! from zone center (X) to zone boundary (XB)
   offset = xsizes(2)
   do it = 1, nprtime
      allocate(eq%time_slice(it)%profiles_1d%magnetic_shear(offset))
      eq%time_slice(it)%profiles_1d%magnetic_shear(1:offset) = &
         0.5 * prdata(1+(it-1)*offset:it*offset)
      ! Extrapolate SHAT(XB(offset)) and save for later
      rdum = 3 * eq%time_slice(it)%profiles_1d%magnetic_shear(offset) - &
                 eq%time_slice(it)%profiles_1d%magnetic_shear(offset-1)
      do ir = 1, offset-1
         eq%time_slice(it)%profiles_1d%magnetic_shear(ir) = &
         eq%time_slice(it)%profiles_1d%magnetic_shear(ir) + &
         eq%time_slice(it)%profiles_1d%magnetic_shear(ir+1)
      enddo
      ! Now write extrapolated boundary value
      eq%time_slice(it)%profiles_1d%magnetic_shear(offset) = rdum
   enddo

   write(iout,*) ' '
   call rprofile('GRI',prdata,nprtime*xsizes(1),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(GRI)',ier)
   if(iret.ne.nprtime*xsizes(1)) &
      call transp2imas_exit(' ?? GRI read error')
   call transp2imas_echo('GRI',prdata,xsizes(1),nprtime)

   offset=xsizes(1)
   do it=1,nprtime
      allocate(eq%time_slice(it)%profiles_1d%gm9(offset))
      eq%time_slice(it)%profiles_1d%gm9(1:offset)= &
         prdata(1+(it-1)*offset:it*offset) * 1.e2
   enddo

   write(iout,*) ' '
   call rprofile('GR2I',prdata,nprtime*xsizes(1),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(GR2I)',ier)
   if(iret.ne.nprtime*xsizes(1)) &
      call transp2imas_exit(' ?? GR2I read error')
   call transp2imas_echo('GR2I',prdata,xsizes(1),nprtime)

   offset=xsizes(1)
   do it=1,nprtime
      allocate(eq%time_slice(it)%profiles_1d%gm1(offset))
      eq%time_slice(it)%profiles_1d%gm1(1:offset)= &
         prdata(1+(it-1)*offset:it*offset) * 1.e4
   enddo

   write(iout,*) ' '
   !? correct. gm2 is m^-2, but GX2R2I is cm^-4
   ! <GRAD(XI)**2/R**2> FLX.SURF.AVG.
   ! The grad has terms like dxi/dR,dxi/dZ
   ! ==>> diffetent from gm2, redo it!!!!
   call rprofile('GX2R2I',prdata,nprtime*xsizes(1),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(GX2R2I)',ier)
   if(iret.ne.nprtime*xsizes(1)) &
      call transp2imas_exit(' ?? GX2R2I read error')
   call transp2imas_echo('GX2R2I',prdata,xsizes(1),nprtime)

   offset=xsizes(1)
   do it=1,nprtime
      allocate(eq%time_slice(it)%profiles_1d%gm2(offset))
      eq%time_slice(it)%profiles_1d%gm2(1:offset)= &
         prdata(1+(it-1)*offset:it*offset) * 1.e8
   enddo

   write(iout,*) ' '
   !? correct. gm3 is dimensonless, but GXI2 is cm^-2
   ! The grad has terms like dxi/dR,dxi/dZ
   ! ==>> diffetent from gm2, redo it!!!!
   call rprofile('GXI2',prdata,nprtime*xsizes(1),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(GXI2)',ier)
   if(iret.ne.nprtime*xsizes(1)) &
      call transp2imas_exit(' ?? GXI2 read error')
   call transp2imas_echo('GXI2',prdata,xsizes(1),nprtime)

   offset=xsizes(1)
   do it=1,nprtime
      allocate(eq%time_slice(it)%profiles_1d%gm3(offset))
      eq%time_slice(it)%profiles_1d%gm3(1:offset)= &
         prdata(1+(it-1)*offset:it*offset) * 1.e4
   enddo

   write(iout,*) ' '
   !? correct. gm7 is dimensonless, but GXI is cm^-1
   ! The grad has terms like dxi/dR,dxi/dZ
   ! ==>> diffetent from gm2, redo it!!!!
   call rprofile('GXI',prdata,nprtime*xsizes(1),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(GXI)',ier)
   if(iret.ne.nprtime*xsizes(1)) &
      call transp2imas_exit(' ?? GXI read error')
   call transp2imas_echo('GXI',prdata,xsizes(1),nprtime)

   offset=xsizes(1)
   do it=1,nprtime
      allocate(eq%time_slice(it)%profiles_1d%gm7(offset))
      eq%time_slice(it)%profiles_1d%gm7(1:offset)= &
         prdata(1+(it-1)*offset:it*offset) * 1.e2
   enddo

   write(iout,*) ' '
   call rprofile('GB2I',prdata,nprtime*xsizes(1),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(GB2I)',ier)
   if(iret.ne.nprtime*xsizes(1)) &
      call transp2imas_exit(' ?? GB2I read error')
   call transp2imas_echo('GB2I',prdata,xsizes(1),nprtime)

   offset=xsizes(1)
   do it=1,nprtime
      allocate(eq%time_slice(it)%profiles_1d%gm4(offset))
      eq%time_slice(it)%profiles_1d%gm4(1:offset)= &
         prdata(1+(it-1)*offset:it*offset)
   enddo

   write(iout,*) ' '
   call rprofile('GB2',prdata,nprtime*xsizes(1),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(GB2)',ier)
   if(iret.ne.nprtime*xsizes(1)) &
      call transp2imas_exit(' ?? GB2 read error')
   call transp2imas_echo('GB2',prdata,xsizes(1),nprtime)

   offset=xsizes(1)
   do it=1,nprtime
      allocate(eq%time_slice(it)%profiles_1d%gm5(offset))
      eq%time_slice(it)%profiles_1d%gm5(1:offset)= &
         prdata(1+(it-1)*offset:it*offset)
   enddo

   write(iout,*) ' '
   call rprofile('GB1',prdata,nprtime*xsizes(1),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(GB1)',ier)
   if(iret.ne.nprtime*xsizes(1)) &
      call transp2imas_exit(' ?? GB1 read error')
   call transp2imas_echo('GB1',prdata,xsizes(1),nprtime)

   offset=xsizes(1)
   do it=1,nprtime
      allocate(eq%time_slice(it)%profiles_1d%b_average(offset))
      eq%time_slice(it)%profiles_1d%b_average(1:offset)= &
         prdata(1+(it-1)*offset:it*offset)
   enddo

   write(iout,*) ' '
   ! calculate Fprime on zone boundaries
   ! GFUN(XB) is para/dia-magnetic factor, already on zone boundaries
   call rprofile('GFUN',prdata,nprtime*xsizes(2),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(GFUN)',ier)
   if(iret.ne.nprtime*xsizes(2)) &
      call transp2imas_exit(' ?? GFUN read error')
   call transp2imas_echo('GFUN',prdata,xsizes(2),nprtime)

   offset = xsizes(2)
   ! Calculate dF/dPLFLX
   do it = 1, nprtime
      ! Multiply to get F...
      xbbuf2(:) = prdata(1+(it-1)*offset:it*offset) * bzxr(it) * 1.0e-2 ! T * cm -> T * m
      allocate(eq%time_slice(it)%profiles_1d%f(offset))
      eq%time_slice(it)%profiles_1d%f(1:offset) = xbbuf2(:)
      
      ! Calculate dF/dXB and put result in xbbuf3.
      xbbuf1(:) = XB(:, it)

      call ezspline_init(spln1, offset, (/0, 0/), ier)
      call ezspline_error(ier)

      !spln1%bcval1min = 0.
      !spln1%isHermite = 1
      spln1%x1 = xbbuf1(1:offset)

      call ezspline_setup(spln1, xbbuf2(1:offset), ier)
      call ezspline_error(ier)

      ! this call only made to init derivate calculation below (xibuf2 should be identical to xbbuf2)
      !call ezspline_interp(spln1, offset, xbbuf1, xibuf2, ier)
      !call ezspline_error(ier)

      call ezspline_derivative(spln1, 1, offset, xbbuf1, xbbuf3, ier)
      call ezspline_error(ier)

      ! F' = dF/dPLFLX = dF/dXB * dXB/dPLFLX => FF' = F * dF/dXB * dXB/dPLFLX
      allocate(eq%time_slice(it)%profiles_1d%f_df_dpsi(offset))
      eq%time_slice(it)%profiles_1d%f_df_dpsi(1:offset) = xbbuf2(:) * xbbuf3(:) * dXBdPLFLX(:, it)
#if 0
      do ir = 2, offset-1
         eq%time_slice(it)%profiles_1d%f_df_dpsi(ir) = eq%time_slice(it)%profiles_1d%f(ir) * &
         (eq%time_slice(it)%profiles_1d%f(ir+1) - eq%time_slice(it)%profiles_1d%f(ir-1)) / &
         (eq%time_slice(it)%profiles_1d%psi(ir+1) - eq%time_slice(it)%profiles_1d%psi(ir-1))
      enddo
      eq%time_slice(it)%profiles_1d%f_df_dpsi(1) = eq%time_slice(it)%profiles_1d%f_df_dpsi(2)
      eq%time_slice(it)%profiles_1d%f_df_dpsi(offset) = eq%time_slice(it)%profiles_1d%f_df_dpsi(offset-1)
#endif

      call ezspline_free(spln1, ier)
      call ezspline_error(ier)

      !write(314,*) 'F func at', it
      !write(314,*) 'size', size(xbbuf2), size(xibuf2), size(xbbuf1)
      !do ir=1,offset
      !   write(314,*) xbbuf1(ir), xbbuf2(ir), xibuf2(ir), xibuf22(ir)
      !enddo
      write(314, *) 'barbaz, F, FFprime at ', it
      do ir = 1, offset
         write(314, *) eq%time_slice(it)%profiles_1d%f(ir), &
                       eq%time_slice(it)%profiles_1d%f_df_dpsi(ir)
      enddo
   enddo
! TRANSP does not have data on LCFS. Need to think about how to provide it. Johan 12/21/16
#if 0
   offset=xsizes(2)
   do it = 1, nprtime
      allocate(eq%time_slice(it)%boundary%lcfs%r(offset))
#endif
   ! fetch limiter data and write it into file
   ! or can get it from dumped g-eqdsk file 'eqdskin'
   write(iout,*) ' '
   call rprofile('ILIM',prdata,nprtime*xsizes(LIMtype),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(ILIM)',ier)
   if(iret.ne.nprtime*xsizes(LIMtype)) &
      call transp2imas_exit(' ?? ILIM read error')
   call transp2imas_echo('ILIM',prdata,xsizes(LIMtype),nprtime)

   offset=xsizes(LIMtype)
   nlimiter=offset
   allocate(ilimiter(offset))
   ilimiter(1:offset)= prdata(1:offset)

   write(iout,*) ' '
   call rprofile('RLIM',prdata,nprtime*xsizes(LIMtype),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(RLIM)',ier)
   if(iret.ne.nprtime*xsizes(LIMtype)) &
      call transp2imas_exit(' ?? RLIM read error')
   call transp2imas_echo('RLIM',prdata,xsizes(LIMtype),nprtime)

   offset=xsizes(LIMtype)
   allocate(rlimiter(offset))
   rlimiter(1:offset)= prdata(1:offset)*1.e-2

   write(iout,*) ' '
   call rprofile('YLIM',prdata,nprtime*xsizes(LIMtype),iret,ier)
   if(ier.ne.0) call transp2imas_error('rprofile(YLIM)',ier)
   if(iret.ne.nprtime*xsizes(LIMtype)) &
      call transp2imas_exit(' ?? YLIM read error')
   call transp2imas_echo('YLIM',prdata,xsizes(LIMtype),nprtime)

   offset=xsizes(LIMtype)
   allocate(ylimiter(offset))
   ylimiter(1:offset)= prdata(1:offset)*1.e-2

   if(allocated(ilimiter) .and. &
      allocated(rlimiter) .and. &
      allocated(ylimiter)) then
      write(130,'(1x,i4,x,a)') nlimiter, "limiter points"
      write(130,'(5x,a,5x,a)') "R, m", "Z, m"
      do ir=1,nlimiter
         write(130,'(4x,f7.4,2x,f7.4)') rlimiter(ir),ylimiter(ir)
      enddo
      deallocate(ilimiter, rlimiter, ylimiter)
   endif
   write(iout,*) ' LIMITER DATA are saved in fort.130'

#ifdef CJDEBUG
!
!  calculator -- simple test
!
   write(iout,*) ' '
   call rpcalc('NE*TE',prdata,nprtime*xsizes(1),iret, istype, &
      iwarn,ier)
   if(ier.ne.0) call transp2imas_error('rpcalc(NE*TE) ier',ier)
   if(iwarn.ne.0) call transp2imas_error('rpcalc(NE*TE) iwarn',iwarn)
   write(iout,*) ' rpcalc call OK, returned istype = ',istype
   if(iret.ne.nprtime*xsizes(1)) &
      call transp2imas_exit(' ?? NE*TE calc size result error')
   call transp2imas_echo('NE*TE (rpcalc)',prdata,xsizes(1),nprtime)

   write(iout,*) ' '
   write(iout,*) ' call RPCAL0("TMP=NE*TE", iwarn,ier)'
   call RPCAL0('TMP=NE*TE', iwarn,ier)
   if(ier.ne.0) call transp2imas_error('rpcal0(NE*TE) ier',ier)
   if(iwarn.ne.0) call transp2imas_error('rpcal0(NE*TE) iwarn',iwarn)

   write(iout,*) ' call RPCALC("TMP+NI*TI", ..., iwarn,ier)'
   call rpcalc('TMP + NI*TI',prdata,nprtime*xsizes(1),iret, istype, &
      iwarn,ier)
   if(ier.ne.0) call transp2imas_error('rpcalc(...) ier',ier)
   if(iwarn.ne.0) call transp2imas_error('rpcalc(...) iwarn',iwarn)
   write(iout,*) ' rpcalc call OK, returned istype = ',istype
   if(iret.ne.nprtime*xsizes(1)) &
      call transp2imas_exit(' ?? TMP+NI*TI calc size result error')
   call transp2imas_echo('NE*TE+NI*TI (rpcal0,rpcalc)', &
      prdata,xsizes(1),nprtime)

!
!  multigraph access
!
   write(iout,*) ' '
   call rpmulti('EEBAL',istype,label,units,infuns,mgsigns, &
      mgnames,ier)
   if(ier.ne.0) call transp2imas_error('rpmulti(EEBAL)',ier)
   write(iout,*) 'EEBAL ',label,' ',units,' istype=',istype

   write(iout,*) 'read the data...'
   call rpmg0cal(mgnames,mgsigns,infuns,mgdata,nmax,iret,istype, &
      iwarn,ier)
   if(ier.ne.0) call transp2imas_error('rpmg0cal ier',ier)
   if(iwarn.ne.0) call transp2imas_error('rpmg0cal iwarn',iwarn)

   write(iout,*) 'EEBAL istype=',istype,':  read results...'
   do i=1,infuns
      write(iout,*) ' '
      call transp2imas_echo(cmgsign(mgsigns(i))//mgnames(i), &
         mgdata(1,i),xsizes(istype),nprtime)
   enddo

   write(iout,*) 'read volint(EEBAL members)'
   call rpmg1cal(mgnames,mgsigns,infuns,'VOLINT(@)',mgdata,nmax, &
      iret,istype,iwarn,ier)
   if(ier.ne.0) call transp2imas_error('rpmg1cal ier',ier)
   if(iwarn.ne.0) call transp2imas_error('rpmg1cal iwarn',iwarn)

   write(iout,*) 'VOLINT(EEBAL) istype=',istype,':  read results...'
   do i=1,infuns
      write(iout,*) ' '
      call transp2imas_echo( &
         'VOLINT('//cmgsign(mgsigns(i))//mgnames(i)//')', &
         mgdata(1,i),xsizes(istype),nprtime)
   enddo
!
!  time slice... the next two calls to rptimav should produce the
!  same numbers, since the linear interpolatin half way btw two
!  time points is the same as the time average btw the two time points
!  of the linearly interpolated data...
!
   ztime=0.5*(prtime(1)+prtime(2))
   call rptimav(mgdata,nmax,infuns,istype,ztime,0.0, &
      mgslice,nxmax, ier)
   if(ier.ne.0) call transp2imas_error('rptimav(0)',ier)
   write(iout,*) ' '
   write(iout,*) ' rptimav test (no time averaging)'
   do i=1,infuns
      write(iout,*) ' '
      call transp2imas_echo( &
         'VOLINT('//cmgsign(mgsigns(i))//mgnames(i)//')', &
         mgslice(1,i),xsizes(istype),1)
   enddo

   dt_avg=0.5*(prtime(2)-prtime(1))
   call rptimav(mgdata,nmax,infuns,istype,ztime,dt_avg, &
      mgslice,nxmax, ier)
   if(ier.ne.0) call transp2imas_error('rptimav(0)',ier)
   write(iout,*) ' '
   write(iout,*) ' rptimav test (time averaging 1st time zone)'
   do i=1,infuns
      write(iout,*) ' '
      call transp2imas_echo( &
         'VOLINT('//cmgsign(mgsigns(i))//mgnames(i)//')', &
         mgslice(1,i),xsizes(istype),1)
   enddo

   call t1scalar('RAXIS',label,units,ztime,dt_avg,zanswer,ier)
   if(ier.ne.0) call transp2imas_error('t1scalar',ier)
   write(iout,'(1x,a,1x,1pe12.5)') &
      ' t1scalar(RAXIS,...) zanswer = ',zanswer

   call t1profil('TE',label,units,ztime,dt_avg,istype,mgslice, &
      nxmax,iret,ier)
   if(ier.ne.0) call transp2imas_error('t1scalar',ier)
   write(iout,*) ' '
   write(iout,*) ' t1profil test...'
   call transp2imas_echo('TE',mgslice(1,1),xsizes(istype),1)
#endif

!---------------------------------
   do it=1,nprtime
      allocate(eq%time_slice(it)%profiles_2d(nprofile))
      r8ztime=prtime(it)

      write(iout,*) ' '
      write(iout,*) 'call t1mhdgeq at time ...', it, r8ztime

      do iprofile=1,nprofile
         whichtimeslice=it
         whichprofile=iprofile
         timeofinterest=r8ztime

         call t1mhdeq_geq(trim(rpfile),r8ztime)
      enddo
   enddo
   !stop
   !   call t1mhdeq_test(ztime,dt_avg,xsizes(1)+1,50)
   !   t1mhdeq_geq
   !      t1mhdeq_compgeq
   !         xplasma_wr_geqdsk  ==>> t1mhdeq_xplasma2geq, somehow
   !         didn't work, module name change also won't work. ????
!---------------------------------
#if 1
!  read namelist file

   write(iout,*) ' '
   ilun=99
   call tr_getnl_text(ilun,ier)
   if(ier.ne.0) then
      call transp2imas_error('tr_getnl_text',ier)
   endif
#endif
!  get number of plasma species

   call tr_getnl_intvec('NGMAX',ngmax,1,istat)
   if(istat.ne.1) then
      call transp2imas_exit('NGMAX not in namelist.')
   else if(ngmax.gt.5) then
      call transp2imas_exit('NGMAX.gt.5')
   endif
   write(iout,*) ' NAMELIST read OK:  NGMAX=',ngmax

!  get A of plasma species

   call tr_getnl_r4vec('APLASM',aplasm,ngmax,istat)
   if(istat.ne.ngmax) then
      call transp2imas_exit('APLASM missing or wrong size')
   endif
   write(iout,3001) (aplasm(i),i=1,ngmax)

!  get Z of plasma species

   call tr_getnl_r4vec('BACKZ',backz,ngmax,istat)
   if(istat.ne.ngmax) then
      call transp2imas_exit('BACKZ missing or wrong size')
   endif
   write(iout,3002) (backz(i),i=1,ngmax)

3001 format(' OK APLASM:  ',5(1x,1pe12.5))
3002 format(' OK BACKZ:   ',5(1x,1pe12.5))

!  get beam data from namelist

   if (nbeam .gt. 0) then
      call tr_getnl_r4vec('FFULLA', rvdum, nbeam, istat)
      call tr_getnl_r4vec('FHALFA', rvdum2, nbeam, istat)
      do i = 1, nbeam
         do k = 1, nsctime
            if (nbi%unit(i)%power%data(k) .gt. 0.0) then ! Is beam on?
               rdum = & ! Normalization constant
               nbi%unit(i)%energy%data(k) * rvdum(i) + &
               nbi%unit(i)%energy%data(k) / 2.0 * rvdum2(i) + &
               nbi%unit(i)%energy%data(k) / 3.0 * (1.0 - rvdum(i) - rvdum2(i))
               nbi%unit(i)%beam_power_fraction%data(1, k) = &
               (nbi%unit(i)%energy%data(k) * rvdum(i)) / rdum ! Fraction of power at full energy
               nbi%unit(i)%beam_power_fraction%data(2, k) = &
               (nbi%unit(i)%energy%data(k) / 2.0 * rvdum2(i)) / rdum ! Fraction of power at 1/2 energy
               nbi%unit(i)%beam_power_fraction%data(3, k) = &
               (nbi%unit(i)%energy%data(k) / 3.0 * (1.0 - rvdum(i) - rvdum2(i))) / rdum ! Fraction of power at 1/3 energy
            else ! Beam is off
               nbi%unit(i)%beam_power_fraction%data(1, k) = 0.0
               nbi%unit(i)%beam_power_fraction%data(2, k) = 0.0
               nbi%unit(i)%beam_power_fraction%data(3, k) = 0.0
            endif
         enddo
      enddo
      ! Distances are converted from cm to m below
      call tr_getnl_r4vec('RTCENA', rvdum, nbeam, istat)
      do i = 1, nbeam
         nbi%unit(i)%beamlets_group(1)%beamlets%tangency_radii(1) = 1.0e-2 * rvdum(i)
         nbi%unit(i)%beamlets_group(1)%tangency_radius = nbi%unit(i)%beamlets_group(1)%beamlets%tangency_radii(1)
      enddo
      call tr_getnl_r4vec('XYBSCA', rvdum, nbeam, istat)
      do i = 1, nbeam
         nbi%unit(i)%beamlets_group(1)%beamlets%positions%z(1) = 1.0e-2 * rvdum(i)
         nbi%unit(i)%beamlets_group(1)%position%z = nbi%unit(i)%beamlets_group(1)%beamlets%positions%z(1)
      enddo
      call tr_getnl_r4vec('XLBTNA', rvdum, nbeam, istat)
      do i = 1, nbeam
         nbi%unit(i)%beamlets_group(1)%beamlets%positions%r(1) = sqrt(1.0e-4 * rvdum(i)**2 - &
         nbi%unit(i)%beamlets_group(1)%beamlets%positions%z(1)**2)
            nbi%unit(i)%beamlets_group(1)%position%r = nbi%unit(i)%beamlets_group(1)%beamlets%positions%r(1)
      enddo
      ! Direction must be calculated from sign of plasma current?
      do i = 1, nbeam
         nbi%unit(i)%beamlets_group(1)%direction = 1
      enddo
      call tr_getnl_r4vec('ABEAMA', rvdum, nbeam, istat)
      do i = 1, nbeam
         nbi%unit(i)%species%a = rvdum(i)
         nbi%unit(i)%species%z_n = 1.0
      enddo
      ! Toroidal angle in degrees (the default is zero)
      call tr_getnl_r4vec('XBZETA', rvdum, nbeam, istat)
      if (istat.eq.nbeam) then
         do i = 1, nbeam ! Switch sign and convert to radians (Johan 01/25/2017)
            nbi%unit(i)%beamlets_group(1)%beamlets%positions%phi(1) = -twopi * rvdum(i) / 360.
            nbi%unit(i)%beamlets_group(1)%position%phi = nbi%unit(i)%beamlets_group(1)%beamlets%positions%phi(1)
         enddo
      else
         do i = 1, nbeam ! Set to zero if XBZETA was not defined
            nbi%unit(i)%beamlets_group(1)%beamlets%positions%phi(1) = 0.
            nbi%unit(i)%beamlets_group(1)%position%phi = 0.
         enddo
      endif
      call tr_getnl_r4vec('FOCLR', rvdum, nbeam, istat)
      do i = 1, nbeam
         nbi%unit(i)%beamlets_group(1)%focus%focal_length_horizontal = 1.0e-2 * rvdum(1) ! Not an index bug
      enddo
      call tr_getnl_r4vec('FOCLZ', rvdum, nbeam, istat)
      do i = 1, nbeam
         nbi%unit(i)%beamlets_group(1)%focus%focal_length_vertical = 1.0e-2 * rvdum(1) ! Not an index bug
      enddo
      call tr_getnl_r4vec('DIVR', rvdum, nbeam, istat)
      do i = 1, nbeam
         !nbi%unit(i)%beamlets_group(1)%focus%width_min_horizontal = 1.0e-2 * rvdum(i)
         nbi%unit(i)%beamlets_group(1)%focus%width_min_horizontal = 0.
      enddo
      call tr_getnl_r4vec('DIVZ', rvdum, nbeam, istat)
      do i = 1, nbeam
         !nbi%unit(i)%beamlets_group(1)%focus%width_min_vertical = 1.0e-2 * rvdum(i)
         nbi%unit(i)%beamlets_group(1)%focus%width_min_vertical = 0.
      enddo
      call tr_getnl_r4vec('BMWIDRA', rvdum, nbeam, istat)
      do i = 1, nbeam
         nbi%unit(i)%beamlets_group(1)%width_horizontal = 1.0e-2 * rvdum(i)
      enddo
      call tr_getnl_r4vec('BMWIDZA', rvdum, nbeam, istat)
      do i = 1, nbeam
         nbi%unit(i)%beamlets_group(1)%width_vertical = 1.0e-2 * rvdum(i)
      enddo
      call tr_getnl_r4vec('XYBSCA', rvdum, nbeam, istat)
      call tr_getnl_r4vec('XLBTNA', rvdum2, nbeam, istat)
      do i = 1, nbeam
         nbi%unit(i)%beamlets_group(1)%angle = asin(rvdum(i) / rvdum2(i))
      enddo
      ! I don't know how to get stuff in nbi%unit(i)%beamlets_group(1)%divergence_component allocated
#if 0
      write(*, *) 'Here -->', nbi%unit(1)%beamlets_group(1)%divergence_component%horizontal
      stop
      do i = 1, nbeam
         ! Set both components to 7.e-3 / (math.sqrt(2.) * 180. / math.pi)
         nbi%unit(i)%beamlets_group(1)%divergence_component%horizontal = 8.638939046419044e-05
         nbi%unit(i)%beamlets_group(1)%divergence_component%vertical = &
            nbi%unit(i)%beamlets_group(1)%divergence_component%horizontal
         nbi%unit(i)%beamlets_group(1)%divergence_component%particles_fraction = 1.0
      enddo
#endif
      call tr_getnl_r4vec('DN0OUT', rvdum, 10, istat)
      offset = size(cp%profiles_1d(1)%electrons%density)
      do it = 1, nprtime
         allocate(cp%profiles_1d(it)%neutral(1))
         allocate(cp%profiles_1d(it)%neutral(1)%density(offset))
         cp%profiles_1d(it)%neutral(1)%density(:) = 0.0
         if (istat .eq. 1) &
            cp%profiles_1d(it)%neutral(1)%density(:) = rvdum(1) * 1.0e6 ! /cm3 to /m3
      enddo
   endif
!
!  save ids data
!
   write(*,*) ' transp2imas: save eq ids'
   call ids_put(idsidx, "equilibrium", eq)
   write(*,*) ' transp2imas: save cp ids'
   call ids_put(idsidx, "core_profiles", cp)
   write(*,*) ' transp2imas: save ct ids'
   call ids_put(idsidx, "core_transport", ct)
   write(*,*) ' transp2imas: save es ids'
   es%source(1)%ggd(1)%time = 0.
   es%source(2)%ggd(1)%time = 0.
   call ids_put(idsidx, "edge_sources", es)
   write(*,*) ' transp2imas: save nbi ids'
   call ids_put(idsidx, "nbi", nbi)

   write(*,*) 'Close shot in IMAS!'

   call ids_deallocate(eq)
   call ids_deallocate(cp)
   call ids_deallocate(ct)
   call ids_deallocate(es)
   call ids_deallocate(nbi)
   call imas_close(idsidx)

!
!---------------------------------
!  exit
   deallocate(xnames)
   deallocate(xsizes)
!
   deallocate(sctime)
   deallocate(prtime)
   deallocate(scdata)
   deallocate(prdata)
   deallocate(mgdata)
   deallocate(mgslice)
!
   if (allocated(bzxr)) deallocate(bzxr)
   if (allocated(PLFLX)) deallocate(PLFLX)
   if (allocated(dXBdPLFLX)) deallocate(dXBdPLFLX)
   if (allocated(XI)) deallocate(XI)
   if (allocated(XB)) deallocate(XB)
   if (allocated(xibuf1)) deallocate(xibuf1)
   if (allocated(xibuf2)) deallocate(xibuf2)
   if (allocated(xbbuf1)) deallocate(xbbuf1)
   if (allocated(xbbuf2)) deallocate(xbbuf2)
   if (allocated(xbbuf3)) deallocate(xbbuf3)

   stop

end

!----------------------------------------------------
!  test t1mhdeq -- mhd timeslice reader
!
subroutine t1mhdeq_test(ztime,zdelta,nsmax,ntheta)

   implicit NONE

! passed:
   real ztime                        ! time of interest
   real zdelta                       ! +/- delta(t)
   integer nsmax                     ! max no. surfaces
   integer ntheta                    ! no. of theta pts

! local:
   integer nsgot

   real theta(ntheta)
   real Rarr(ntheta,nsmax)
   real Zarr(ntheta,nsmax)
   real rho(nsmax)
   real psi(nsmax)
   real pmhd(nsmax)
   real qmhd(nsmax)
   real gmhd(nsmax)
   real tflux
   real pcur

   integer ier

!----------------------
   integer i
   real zpi,z2pi

   integer imsg,iout                 ! i/o unit numbers

   common/msgs/ imsg,iout
!
!----------------------------
!
   zpi=acos(-1.0)
   z2pi=2.0*zpi
!
   do i=1,ntheta
      theta(i)= -zpi + (i-1)*z2pi/(ntheta-1)
   enddo
!
!***
!
   call t1mhdeq(ztime,zdelta,nsmax,ntheta,nsgot,'MKS', &
      theta,Rarr,Zarr, &
      rho,psi,pmhd,qmhd,gmhd, &
      tflux,pcur, ier)
!
!***
!
   if(ier.ne.0) call transp2imas_error('t1scalar',ier)

   write(iout,*) ' '
   write(iout,*) ' t1mhdeq output in MKS units.'

   write(iout,*) ' '
   call transp2imas_echo('Rarr on axis',Rarr,ntheta,1)
   write(iout,*) ' '
   call transp2imas_echo('Zarr on axis',Zarr,ntheta,1)
   write(iout,*) ' '
   write(iout,*) ' Rarr @edge:'
   write(iout,'(6(1x,1pe12.5))') (Rarr(i,nsgot),i=1,ntheta)

   write(iout,*) ' '
   write(iout,*) ' Zarr @edge:'
   write(iout,'(6(1x,1pe12.5))') (Zarr(i,nsgot),i=1,ntheta)

   write(iout,*) ' '
   call transp2imas_echo('rho(t1mhdeq)',rho,nsgot,1)
   write(iout,*) ' '
   call transp2imas_echo('psi(t1mhdeq)',psi,nsgot,1)
   write(iout,*) ' '
   call transp2imas_echo('pmhd(t1mhdeq)',pmhd,nsgot,1)
   write(iout,*) ' '
   call transp2imas_echo('gmhd(t1mhdeq)',gmhd,nsgot,1)
   write(iout,*) ' '
   call transp2imas_echo('qmhd(t1mhdeq)',qmhd,nsgot,1)

   write(iout,*) ' '
   call transp2imas_echo('tflux(t1mhdeq)',tflux,1,1)
   write(iout,*) ' '
   call transp2imas_echo('pcur(t1mhdeq)',pcur,1,1)

   return
end

!----------------------------------------------------
subroutine transp2imas_error(srname,ier)
!
!  write message w/status code and exit
!
   character*(*) srname
   integer ier

   integer imsg,iout                 ! i/o unit numbers

   common/msgs/ imsg,iout
!
!----------
!
   write(iout,*) ' ?? subroutine ',srname,' returned ier=',ier
   write(imsg,*) ' ?? subroutine ',srname,' returned ier=',ier

   call bad_exit
   stop
end

!----------------------------------------------------
subroutine transp2imas_exit(msg)
!
!  write message and exit
!
   character*(*) msg

   integer imsg,iout                 ! i/o unit numbers

   common/msgs/ imsg,iout
!
!----------
!
   write(iout,*) ' transp2imas error exit:  ',msg
   write(imsg,*) ' transp2imas error exit:  ',msg

   call bad_exit
   stop
end

!----------------------------------------------------
subroutine transp2imas_echo(dname,ditem,idimx,idimt)
!
!  write a few words of the data...
!
   character*(*) dname               ! data name or label
   integer idimx,idimt               ! dimensioning info
   real ditem(idimx,idimt)           ! the data
!
!----------
!
   integer imsg,iout                 ! i/o unit numbers

   common/msgs/ imsg,iout
!
!----------
!
   write(iout,*) ' data name:  ',dname
   if(idimx.lt.4) then
      if(idimx.gt.1) write(6,*) ' nx = ',idimx,'; @x(1):'
      if(idimt.gt.1) then
         do it=1,4
            jt=it
            if(jt.eq.4) then
               write(iout,*) '         ... '
               jt=idimt
            endif
            write(iout,1000) jt,ditem(1,jt)
1000        format(' @time(',i4,'):  ',1pe12.5)
         enddo
      else
         write(iout,2000) ditem(1,1)
2000     format(' @0.5*(time(1)+time(2)): ',1pe12.5)
      endif
   else
      write(iout,*) ' nx = ',idimx
      if(idimt.gt.1) then
         do it=1,4
            jt=it
            if(jt.eq.4) then
               write(iout,*) '         ... '
               jt=idimt
            endif
            write(iout,1001) jt,ditem(1,jt),ditem(2,jt),ditem(3,jt), &
               ditem(idimx,jt)
1001        format(' @time(',i4,'): ',3(1x,1pe12.5),' ... ',1pe12.5)
         enddo
      else
         write(iout,2001) ditem(1,1),ditem(2,1),ditem(3,1), &
            ditem(idimx,1)
2001     format(' @0.5*(time(1)+time(2)):',3(1x,1pe12.5), &
            ' ... ',1pe12.5)
      endif
   endif

   return
end

!----------------------------------------------------

SUBROUTINE getids_shotid(runid,shot,run)
   implicit none
   CHARACTER*(*),INTENT(IN) :: runid  !transp runid
   INTEGER :: shot,run  !ids shot run numnber
   !INTEGER,INTENT(OUT) :: shot,run  !ids shot run numnber

   CHARACTER(26) :: delim='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
   CHARACTER(10) :: instring,string1,string2
   INTEGER :: index=0,itmp,stat

   write(*,*) ' '
   instring = TRIM(runid)

   !string1 will be parsed to ids shot and run numbers
   !string2 is transp user specified runid
   index = SCAN(instring,delim)
   string1 = instring(1:index-1)
   string2 = instring(index+1:)

   write(*,*) '       transp runid: ', trim(runid)
   !write(*,*) 'string1=',string1
   !write(*,*) 'string2=',string2
   !write(*,*) 'index=',index

   if(index.eq.6) then
      ! there are 5 digits in the transp runid
      read(string1,'(i5)') itmp
      !read(string1,*,iostat=stat) itmp
      shot=itmp/100
      run=itmp-shot*100
   else if(index.eq.7) then
      ! there are 6 digits in the transp runid
      read(string1,'(i6)') itmp
      shot=itmp/1000
      run=itmp-shot*1000
   else
      write(*,*) 'transp runid must have 5 or 6 digits'
      call transp2imas_exit(' ?? transp runid is wrong')
   endif
   write(*,*) '       ids shot: ',shot
   write(*,*) '       ids run: ',run
   write(*,*) ' '
END SUBROUTINE getids_shotid

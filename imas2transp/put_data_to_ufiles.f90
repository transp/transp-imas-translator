
!                                                                            !
! prepare data for TRANSP from ITER IDS                                      !
! TRANSP UFILES LIB is called                                                !
! The correct order of calls: UFSETR, UFOPRD/UFOPWR, append comment, UFCLOS, !
! 2015-march-27 jinc chen                                                    !
!                                                                            !
!............................................................................!

subroutine put_data_to_ufiles(ilun, &
                              prefix,suffix,disk,directory, &
                              ishot, &
                              tdev,ndim,shdate, &
                              comment, &
                        t_uf0d,t_uf1d,t_uf2d,t_uf3d, &
                              ierr)
use transp_ufiles_0d
use transp_ufiles_1d
use transp_ufiles_2d
use transp_ufiles_3d
implicit none

    integer,intent(out) :: ierr
    ! IERR: completion code.
    ! IERR=0 success success.
    ! IERR=1 err from file open
    ! IERR=2 err from write scalar
    ! IERR=3 err from write 1d array
    ! IERR=4 err from write 2d array
    ! IERR=5 err from write 3d array
    ! IERR=6 err from write comment
    ! IERR=7 err from file close

!1. associate a FORTRAN logical unit number with a filename
    integer,intent(in) :: ilun
    ! ilun - integer, logic unit number

    CHARACTER(len=*),intent(in) :: prefix, suffix, disk, directory

    ! prefix - character, max length 16
    ! suffix - character, max length 16
    ! disk - character, max length 16, here blank
    ! directory - character, max length 64

!2. open file with a shot number
    integer,intent(in) :: ishot  !(6-digit shot number)

!.
    CHARACTER(len=*) TDEV !(tokamak id)
    !(0,1,2,3. 0 for scalar, 1 for 1d array, 2 for 2d array, 3 for 3darray)
    integer,intent(in) :: ndim
    CHARACTER(len=*),intent(in) :: SHDATE !(shot date)

!3.1 write scalar
    !PARAMETER (NSCMAX= maximum number of scalars; array dimension ndim=0)
    !REAL SCVAL (NSCMAX) ! (scalar data array)

    ! SCVAL- VALUE(S) OF ASSOCIATED SCALARS
    ! SCLAB- ASSOCIATED SCALARS KEYWORDS, LABELS, PHYS. UNITS.

!3.2 write 1d array
    !PARAMETER (NXMAX= size of arrays to hold x,f(x))
    !REAL X(NXMAX), F(NXMAX) !( 1d arrays to hold x, f(x))
    !CHARACTER*10 XLAB(3), FLAB(3) !( label arrays for X and F data) ! 
    !integer iproc  ! process code-- historical, use a value of zero.

    ! X: array containing independent variable values to write 
    ! F: array containing dependent variable values to write.
    ! NX: number of values in X, F to write
    ! IPROC: process code-- historical, use a value of zero.
    ! XLAB: 30 character string containing: characters 1-10: the name of "x" (e.g. "POSITION"); 
    !       characters 11-20: blank; characters 21-30: physical units of "x"  (e.g. "METERS   ").
    ! FLAB: 30 character string containing: characters 1-20: the name of "f" (e.g. "ELECTRON DENSITY   "); 
    !       characters 21-30: physical units of "f" (e.g. "N/CM**3  ").
    ! IPROC: process code-- see description, UFILES 2-d file RECORD 2*NSC+7.

!3.3 write 2d array
    !PARAMETER (NXMAX= size of array to hold x)
    !PARAMETER (NYMAX= size of array to hold y)
    !PARAMETER (NFMAX= NXMAX*NYMAX)  ! size of array to hold f (optional)
    !PARAMETER (NF1= NXMAX) 
    !REAL X(NXMAX), Y(NYMAX) !(arrays to hold x,y)
    !REAL F(NXMAX,NYMAX) !or REAL F(NFMAX) (array to hold f. The dual 
    !CHARACTER*10 XLAB(3), YLAB(3), FLAB(3) ! (label arrays for X, Y, and  F)

    ! F: array containing dependent variable data.
    ! NF1: size of the first dimension of F, as it will be declared inside the 
    ! subroutine UF2DWR. If the data in F is stored contiguously then NF1=NX should be set. 
    ! REAL F(NF1,NY) is to be written with a statement of the form
    ! WRITE(ILUN,format) ((F(I,J),I=1,NX),J=1,NY).
    ! X: array containing 1st independent variable data
    ! NX: number of X values to write
    ! Y: array containing 2nd independent variable data.
    ! NY: number of Y values to write
    ! YLAB: 30 character string containing the name of "y" characters 1-10) and 
    ! the physical units of "y" (characters 21-30) For example: "TIME       SECONDS   ".
    ! FLAB: 30 character string containing: characters 1--20: the name of "f" 
    ! (e.g. "ELECTRON DENSITY "); characters 21--30: physical units of "f" (e.g.  "N/CM**3 ").

!3.4 write 3d array
    !PARAMETER (NXMAX= size of array to hold x)
    !PARAMETER (NYMAX= size of array to hold y)
    !PARAMETER (NZMAX= size of array to hold z)
    !PARAMETER (NFMAX= NXMAX*NYMAX*NZMAX)  ! size of array to hold f (optional)
    !PARAMETER (NF1= NXMAX) 
    !PARAMETER (NF2= NYMAX) 
    !REAL X(NXMAX), Y(NYMAX), Z(NZMAX) !(arrays to hold x,y,z)
    !REAL F(NXMAX,NYMAX) !or REAL F(NFMAX) (array to hold f)
    !CHARACTER*10 XLAB(3), YLAB(3), ZLAB(3), FLAB(3) ! (label arrays for X, Y, X, and F)

    ! X(NX),F(NF1,NF2,NZ),Y(NY),Z(NZ)-- REAL NUMERIC DATA TO BE WRITTEN
    ! ARRAYS AND DIMENSIONS
    ! IPROC- PROCESS CODE 0= "RAW", 1= AVERAGED 2= SMOOTHED, 3=AVERAGED AND SMOOTHED
    ! XLAB-- 30 CHARACTER LABEL FOR X FIRST 20 CHARACTERS-- NAME OF X
    ! LAST 10 CHARACTERS--  PHYSICAL UNITS OF X E.G. XLAB = 'MAJOR RADIUS        CM        '
    ! YLAB-- 30 CHARACTER LABEL FOR Y - SAME FORMAT AS XLAB
    ! ZLAB-- 30 CHARACTER LABEL FOR Z - SAME FORMAT AS XLAB
    ! FLAB-- 30 CHARACTER LABEL FOR F - SAME FORMAT AS XLAB
 
!4. append comments to the end of the file after the standard UFILES data write is complete.
    CHARACTER(len=*),intent(in) :: comment
 

    ! ... local ... !
    ! label X (coor) , Y (coor), Z (coor), and F (function of X,Y,Z)
    CHARACTER*30 XLAB, YLAB, ZLAB, FLAB
    type(transp_ufiles_0d_data) :: t_uf0d
    type(transp_ufiles_1d_data) :: t_uf1d
    type(transp_ufiles_2d_data) :: t_uf2d
    type(transp_ufiles_3d_data) :: t_uf3d
    integer :: is



               ! ............ code start here ............ !

!1. The association of a FORTRAN logical unit number (ILUN)
!   with a filename form, device , and a directory path.
!   call UFSETR(lun,prefix,suffix,disk,directory)
!
    !force uppercase filename prefix and suffix (default)
    call UFNCAP(1)

    !set filename prefix, suffix, and directory
    call UFSETR(ilun,prefix,suffix,'',directory)

    !enhance the default ascii format writing
    CALL UFCMPR(ILUN,0)


!2. A shot (experiment) number is supplied (this completes the filename whose 
!   form was specified in step 1), and the UFILES file is opened for write access. 
!   call UFOPWR(lun,ishot,ierr)
!
!   ISHOT - the 6 digit shot number between 1000 and 999999

    CALL UFOPWR(ilun,ishot,ierr) 


!3. write


   if( ndim==0 ) then

!3.1 write scalar
!   CALL UF0DWR(ILUN,TDEV,SHDATE,SCVAL,SCLAB,NSC,IER)
!
!   RECORD 1: DATA ID [1X,I6,1X,A4,1X,I1]: 6 digit shot number, 4 character 
!source tokamak id, and maximum dimensionality of data in the file (0, 1, 2, or 3. 
!0 for a scalar data file, 1 for a one-dimensional data file, 
!2 for a two-dimensional data file, and 3 for a three-dimensional data file.
!
!   RECORD 2: SHOT DATE [1X,A10]: 10 character shot date, i.e. date when the
!experiment was run which produced the data in the file. UFILES does not 
!require that the date be in any particular format.
!
!   RECORD 3: NUMBER OF SCALARS IN FILE [1X,I3]: number of scalars, NSC, in 
!the data file. In 1d and 2d UFILES files NSC=0 is allowed; in 0d (scalar) 
!files NSC must be greater than or equal to 1.
!
!   RECORDS  4,6,  ... 3+2*NSC-1:   SCALAR QUANTITIES: [1X,1PE13.6]:   scalar 
!data in floating point format with 5 digit precision.
!
!   RECORDS 5,7,  ... , 3+2*NSC: SCALAR LABELS: [1X,3A10]:   30 character 
!labels for each scalar quantity. the first 10 characters are reserved for the
!keyword which must end in a colon leading blanks, and the character 
!preceding the colon must not be a blank. The second 10 characters are
!reserved for a descriptive label (arbitrary format), and the last 10 
!characters are reserved for a physical units label (arbitrary format).

   !nsc=NSCMAX
   write(*,*) "number of scalars=",t_uf0d%nsc
   allocate( t_uf0d%SCLAB( t_uf0d%nsc ) )
   do is=1,t_uf0d%nsc
      write(t_uf0d%SCLAB(is),101) t_uf0d%labels(is),t_uf0d%unitss(is)
      !write(*,*) "scalar val=", is, t_uf0d%SCLAB(is)
   enddo
   CALL R8_UF0DWR(ILUN,TDEV,SHDATE, &
                  t_uf0d%SCVAL,t_uf0d%SCLAB,t_uf0d%nsc, &
                  IERR)
   deallocate( t_uf0d%SCLAB )


   else if( ndim==1 ) then

!3.2 write 1d array
!   CALL UF1DWR(ILUN,TDEV,SHDATE,X,F,NX,IPROC,XLAB,FLAB,NSC,SCVAL,SCLAB,IER)
!
!     RECORDS 1, ... , 2*NSC+3: shot number, tokamak id, dimensionality
! (=1), shot date, number of scalars (NSC), and scalar data. Format 
!descriptions are given in subsection 2.2 above - these records are identical
!to corresponding records in UFILES scalar data files. However, in 
!one-dimensional data file NSC=0 is allowed.
!
!     RECORD 2*NSC+4:   INDEPENDENT COORDINATE (X) LABEL:[1X, 3A10]: 30
!character label for the one-dimensional data x coordinate. The first 20
!characters are reserved for the coordinate name; the last 10 characters for 
!physical units. If the data is to "concatenated" with other data of a similar 
!form to create a function of 2 independent coordinates (cf. subsection 2.4 
!below), it may be desirable to leave the 11th thru 20th character positions 
!in the name blank, because these will be inaccessible to graphics routines 
!utilized by UFILES-based utilities for two-dimensional data.
!
!     RECORD 2*NSC+5:  DEPENDENT COORDINATE (F) LABEL: [1X, 3A10]: 30 
!character label for the one dimensional data function (f) coordinate. The 
!first 20 characters are reserved for the name; the last 10 characters for 
!physical units.
!
!     RECORD 2*NSC+6: PROCESS CODE [I1]: Integer code used by UFILES-based 
!utility programs. =0 to indicate unprocessed data, =1 to indicate averaged 
!data (formed from more than one shot), =2 to indicate smoothed data, =3 to 
!indicate averaged and smoothed data. Thus if the smoothing utility GSMOO1 is 
!used on data with the processing code set to 2, a warning message to the 
!effect that the data has already been smoothed, will be generated. But the 
!user will not be prevented from applying further smoothing to the data.
!
!     RECORD 2*NSC+7: NUMBER OF POINTS (NX): [1X,I10]: The number of X and 
!F(X) data points in the file.
!
!     RECORD 2*NSC+8: X ARRAY: [6(1X,1PE13.6) repeating]: The X data array 
!which is always written out explicitly. Thus an unevenly spaced or even 
!non-monotonic sequence of X values may be supplied. As many lines are written 
!as needed to specify NX points 6 per line.
!
!     RECORD 2*NSC+9: F(X) ARRAY: [6(1X,1PE13.6) repeating]: The F(X) data 
!array. As many lines are written as needed to specify NX points 6 per line.

   !nsc=0
   !nx=NXMAX
   !iproc=0
   if(t_uf1d%nsc > 0) then
      write(*,*) "number of scalars=",t_uf1d%nsc
      allocate( t_uf1d%SCLAB( t_uf1d%nsc ) )
      do is=1,t_uf1d%nsc
         write(t_uf1d%SCLAB(is),101) t_uf1d%labels(is),t_uf1d%unitss(is)
         !write(*,*) "scalar val=", is, t_uf1d%SCLAB(is)
      enddo
    endif
   write(XLAB,101) t_uf1d%labelx, t_uf1d%unitsx
   write(FLAB,101) t_uf1d%labelf, t_uf1d%unitsf
   CALL R8_UF1DWR(ILUN,TDEV,SHDATE, &
                  t_uf1d%X,t_uf1d%F, &
                  t_uf1d%nx, &
                  t_uf1d%iproc, &
                  XLAB,FLAB, &
                  t_uf1d%nsc,t_uf1d%SCVAL,t_uf1d%SCLAB, &
                  IERR)
   if(t_uf1d%nsc > 0) then
      deallocate( t_uf1d%SCLAB )
   endif

   else if( ndim==2 ) then

!3.3 write 2d array
!   CALL UF2DWR(ILUN, TDEV,SHDATE,F,NF1,X,NX,Y,NY,IPROC,FLAB,XLAB,YLAB,NSC,SCVAL,SCLAB,IER)
!
!   RECORDS 1, ... , 2*NSC+3: shot number, tokamak id, dimensionality (=2),
!   shot date, number of scalars (NSC), and scalar data. 
!   In two- dimensional data files NSC=0 is allowed.
!
!   RECORD  2*NSC+4:   INDEPENDENT   COORDINATE (X) LABEL: [1X,3A10]: 30 
!   character label for the two-dimensional data x coordinate. The first 20 
!   characters are reserved for the coordinate name; the last 10 characters for 
!   physical units. However, graphics subroutines utilized by UFILE-based 
!   utilities often only allow 10 characters for the coordinate name, therefore 
!   it is preferable to leave character positions 11-20 blank if possible.
!
!   RECORD  2*NSC+5:   INDEPENDENT COORDINATE  (Y) LABEL: [1X,3A10]:  30 
!characters for the second independent coordinate of the two dimensional data 
!function.  Format the same as the X label.
!
!   RECORD 2*NSC+6:  DEPENDENT COORDINATE (F) LABEL:  [1X,3A10]: 30 character 
!label for the two dimensional data function (f) coordinate. The first 20 
!characters are reserved for the name; the last 10 characters for physical units.
!
!   RECORD 2*NSC+7: PROCESS CODE [I1]: Integer code used by UFILES-based 
!utility programs . =0 to indicate unprocessed data, =1 to indicate averaged 
!data (formed from more than one shot),  =2 to indicate smoothed data, =3 to 
!indicate averaged and smoothed data. Thus if the smoothing utility GSMOO2 is 
!used on data with processing code set to 2, a warning message to the effect 
!that the data has already been smoothed, will be generated. But the user will 
!not be prevented from applying further smoothing to the data.
!
!   RECORD 2*NSC+8: NUMBER OF X POINTS (NX):  [1X,I10]: The number of X points 
!i.e. the length of the first dimension of the data function.
!
!   RECORD 2*NSC+9: NUMBER OF Y POINTS (NY):  [1X,I10]: The number of Y points 
!i.e. the length of the second dimension of the data function.
!
!   RECORD 2*NSC+10: X ARRAY: [6(1X,1PE13.6) repeating]: The data array which 
!is always written out explicitly. Thus an unevenly spaced or even 
!non-monotonic sequence of Y values may be supplied. As many lines are written 
!as needed to specify NY points 6 per line.
!
!   RECORD 2*NSC+11: YARRAY: [6((1X,1PE13.6) repeating]: The data array which 
!is always written out explicitly. Thus an unevenly space or non-monotonic 
!sequence of Y values may be supplied. As many lines are written as needed to 
!specify NY points per line.
!
!   RECORD 2*NSC+12: F(X,Y) ARRAY: [6(1X,1PE13.6) repeating]:
!The F(X,Y) data array. As many lines are written as needed to specify NX*NY 
!points 6 per line. NOTE that the X-variation of the data function is stored 
!contiguously, i.e. the order of the data is all of f vs. x at the first y, 
!then all of f vs. x at the second y, etc. This is consistent with the FORTRAN 
!convention for the storage of data in multiply-subscripted arrays.


   !ndim=2
   !nsc=0
   !nx=NXMAX
   !ny=NYMAX
   !nf1=NXMAX
   !iproc=0
   if(t_uf2d%nsc > 0) then
      write(*,*) "number of scalars=",t_uf2d%nsc
      allocate( t_uf2d%SCLAB( t_uf2d%nsc ) )
      do is=1,t_uf2d%nsc
         write(t_uf2d%SCLAB(is),101) t_uf2d%labels(is),t_uf2d%unitss(is)
      enddo
    endif
   write(XLAB,101) t_uf2d%labelx, t_uf2d%unitsx
   write(YLAB,101) t_uf2d%labely, t_uf2d%unitsy
   write(FLAB,101) t_uf2d%labelf, t_uf2d%unitsf
   CALL R8_UF2DWR(ILUN, TDEV,SHDATE, &
                  t_uf2d%F,t_uf2d%NF1, &
                  t_uf2d%X,t_uf2d%nx, &
                  t_uf2d%Y,t_uf2d%ny, &
                  t_uf2d%IPROC, &
                  FLAB,XLAB,YLAB, &
                  t_uf2d%NSC,t_uf2d%SCVAL,t_uf2d%SCLAB, &
                  IERR)
   if(t_uf2d%nsc > 0) then
      deallocate( t_uf2d%SCLAB )
   endif

   else if( ndim==3 ) then

!3.4 write 3d array
!   SUBROUTINE UF3DWR(ILUN,IDEV,ZSDATE,F,NF1,NF2,X,NX,Y,NY,Z,NZ,IPROC,
!   FLAB,XLAB,YLAB,ZLAB,NSC,SCVAL,SCLAB,IER)
 
   !ndim=3
   !nsc=0
   !nx=NXMAX
   !ny=NYMAX
   !nz=NZMAX
   !nf1=NXMAX
   !nf2=NYMAX
   !iproc=0
   if(t_uf3d%nsc > 0) then
      write(*,*) "number of scalars=",t_uf3d%nsc
      allocate( t_uf3d%SCLAB( t_uf3d%nsc ) )
      do is=1,t_uf3d%nsc
         write(t_uf3d%SCLAB(is),101) t_uf3d%labels(is),t_uf3d%unitss(is)
      enddo
    endif
   write(XLAB,101) t_uf3d%labelx, t_uf3d%unitsx
   write(YLAB,101) t_uf3d%labely, t_uf3d%unitsy
   write(ZLAB,101) t_uf3d%labelz, t_uf3d%unitsz
   write(FLAB,101) t_uf3d%labelf, t_uf3d%unitsf
   call R8_UF3DWR(ILUN,TDEV,shdate, &
                  t_uf3d%F,t_uf3d%nf1,t_uf3d%nf2, &
                  t_uf3d%X,t_uf3d%nx, &
                  t_uf3d%Y,t_uf3d%ny, &
                  t_uf3d%Z,t_uf3d%nz, &
                  t_uf3d%iproc, &
                  FLAB,XLAB,YLAB,ZLAB, &
                  t_uf3d%nsc,t_uf3d%SCVAL,t_uf3d%SCLAB, &
                  IERR)
   if(t_uf3d%nsc > 0) then
      deallocate( t_uf3d%SCLAB )
   endif
 
   endif

!4. This step is optional, but may be important: 
!   From the point of view of the UFILES system, any 
!   non-standard data is a "comment". 

!all UFILES files are divided into two parts:
!a data section at the top and a comments section at the bottom,
!separated by the line:
!  ;---END-OF-DATA---------------------COMMENTS:---------
!This line is the last line written by any UFILES write subroutine (UF0DWR, UF1DWR, or UF2DWR).
!so that on exit from these routines the file being written is left open and the next line is the first line of comments.

!Comment lines are never accessed by UFILES subroutines, so the format of these lines is free
!(but still formatted ascii data with a limit of 160 characters per line). 

    !comment=";---END-OF-DATA---------------------COMMENTS:---------"
    WRITE(ILUN,8001) comment
    8001 FORMAT (a160)

!C USE A TEMPORARY *ASCII* FILE OPENED BY UFOPCF SUBROUTINE 
!   CALL UFOPCF (LUNT, IER) 
!   WRITE(LUNT,8001) comment
!   8001 FORMAT 'a160'
!C SPOOL CONTENTS OF TEMPORARY FILE
!   CALL UCSEND(LUNT, ILUN)


!5. Close the file.
   CALL UFCLOS(ILUN) 

      101 format ( a20,a10 )

end subroutine put_data_to_ufiles

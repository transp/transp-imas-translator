!
! prepare data for TRANSP from ITER ids
! TRANSP UFILES lib is called
! The correct order of calls: UFSETR, UFOPRD/UFOPWR, append comment, UFCLOS,
! check ierr for errors, ierr=0 indicates success
! 2015-march-27 jinc chen
!

module transp_ufiles_0d
  implicit none

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

!   CHARACTER*16 prefix, suffix, disk, directory

!   integer ishot  !(6-digit shot number)
!   CHARACTER*4 TDEV !(tokamak id)
!   integer ndim !(0,1,2,3. 0 for scalar, 1 for 1d array, 2 for 2d array, 3 for 3darray)

!   CHARACTER*10 SHDATE !(shot date)

  type transp_ufiles_0d_data
    integer nsc   !number of scalars 

    REAL,allocatable,dimension(:) :: SCVAL ! (scalar data array)

    CHARACTER(len=30),allocatable,dimension(:) :: SCLAB ! (scalar labels array)
    character(len=20),allocatable,dimension(:) :: labels
    character(len=10),allocatable,dimension(:) :: unitss

  end type transp_ufiles_0d_data

!   CHARACTER*160 comment
end module transp_ufiles_0d


module transp_ufiles_1d
  implicit none

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
!   CHARACTER*16 prefix, suffix, disk, directory

!   integer ishot  !(6-digit shot number)
!   CHARACTER*4 TDEV !(tokamak id)
!   integer ndim !(0,1,2,3. 0 for scalar, 1 for 1d array, 2 for 2d array, 3 for 3darray)

!   CHARACTER*10 SHDATE !(shot date)

  type transp_ufiles_1d_data
    integer nsc   !number of scalars 

    REAL,allocatable,dimension(:) :: SCVAL ! (scalar data array)

    CHARACTER(len=30),allocatable,dimension(:) :: SCLAB ! (scalar labels array)
    character(len=20),allocatable,dimension(:) :: labels
    character(len=10),allocatable,dimension(:) :: unitss

    CHARACTER(len=30) :: XLAB, FLAB !( label arrays for X and F data) ! character*30 XLAB, FLAB
    character(len=20) :: labelx,labelf
    character(len=10) :: unitsx,unitsf

    integer iproc  ! process code-- historical, use a value of zero.

    integer nx   ! actual number of values in X, F to write

    REAL,allocatable,dimension(:) :: X, F !( 1d arrays to hold X, F(X))
  end type transp_ufiles_1d_data

!   CHARACTER*160 comment
end module transp_ufiles_1d

module transp_ufiles_2d
  implicit none

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

!   CHARACTER*16 prefix, suffix, disk, directory

!   integer ishot  !(6-digit shot number)
!   CHARACTER*4 TDEV !(tokamak id)
!   integer ndim !(0,1,2,3. 0 for scalar, 1 for 1d array, 2 for 2d array, 3 for 3darray)

!   CHARACTER*10 SHDATE !(shot date)

  type transp_ufiles_2d_data
    integer nsc   !number of scalars 

    REAL,allocatable,dimension(:) :: SCVAL ! (scalar data array)

    CHARACTER(len=30),allocatable,dimension(:) :: SCLAB ! (scalar labels array)
    character(len=20),allocatable,dimension(:) :: labels
    character(len=10),allocatable,dimension(:) :: unitss

    CHARACTER(len=30) :: XLAB, YLAB, FLAB ! (label arrays for X, Y, and  F) ( CHARACTER*30 XLAB,FLAB )
    character(len=20) :: labelx,labely,labelf
    character(len=10) :: unitsx,unitsy,unitsf

    integer iproc  ! process code-- historical, use a value of zero.

    integer nx, ny, nf1   ! actual number of values in X and Y , first dimension of F

    REAL,allocatable,dimension(:) :: X, Y !(arrays to hold x,y)

    REAL,allocatable,dimension(:,:) :: F !(array to hold F)
  end type transp_ufiles_2d_data

!   CHARACTER*160 comment
end module transp_ufiles_2d

module transp_ufiles_3d
  implicit none

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
!   RECORD  2*NSC+6:   INDEPENDENT COORDINATE  (Z) LABEL: [1X,3A10]:  30 
!characters for the second independent coordinate of the two dimensional data 
!function.  Format the same as the X label.
!
!   RECORD 2*NSC+7:  DEPENDENT COORDINATE (F) LABEL:  [1X,3A10]: 30 character 
!label for the two dimensional data function (f) coordinate. The first 20 
!characters are reserved for the name; the last 10 characters for physical units.
!
!   RECORD 2*NSC+8: PROCESS CODE [I1]: Integer code used by UFILES-based 
!utility programs . =0 to indicate unprocessed data, =1 to indicate averaged 
!data (formed from more than one shot),  =2 to indicate smoothed data, =3 to 
!indicate averaged and smoothed data. Thus if the smoothing utility GSMOO2 is 
!used on data with processing code set to 2, a warning message to the effect 
!that the data has already been smoothed, will be generated. But the user will 
!not be prevented from applying further smoothing to the data.
!
!   RECORD 2*NSC+9: NUMBER OF X POINTS (NX):  [1X,I10]: The number of X points 
!i.e. the length of the first dimension of the data function.
!
!   RECORD 2*NSC+10: NUMBER OF Y POINTS (NY):  [1X,I10]: The number of Y points 
!i.e. the length of the second dimension of the data function.
!
!   RECORD 2*NSC+11: NUMBER OF Z POINTS (NZ):  [1X,I10]: The number of Y points 
!i.e. the length of the second dimension of the data function.
!
!   RECORD 2*NSC+12: X ARRAY: [6(1X,1PE13.6) repeating]: The data array which 
!is always written out explicitly. Thus an unevenly spaced or even 
!non-monotonic sequence of Y values may be supplied. As many lines are written 
!as needed to specify NY points 6 per line.
!
!   RECORD 2*NSC+13: YARRAY: [6((1X,1PE13.6) repeating]: The data array which 
!is always written out explicitly. Thus an unevenly space or non-monotonic 
!sequence of Y values may be supplied. As many lines are written as needed to 
!specify NY points per line.
!
!   RECORD 2*NSC+14: ZARRAY: [6((1X,1PE13.6) repeating]: The data array which 
!is always written out explicitly. Thus an unevenly space or non-monotonic 
!sequence of Y values may be supplied. As many lines are written as needed to 
!specify NY points per line.
!
!   RECORD 2*NSC+15: F(X,Y,Z) ARRAY: [6(1X,1PE13.6) repeating]:
!The F(X,Y) data array. As many lines are written as needed to specify NX*NY 
!points 6 per line. NOTE that the X-variation of the data function is stored 
!contiguously, i.e. the order of the data is all of f vs. x at the first y, 
!then all of f vs. x at the second y, etc. This is consistent with the FORTRAN 
!convention for the storage of data in multiply-subscripted arrays.

!   CHARACTER*16 prefix, suffix, disk, directory

!   integer ishot  !(6-digit shot number)
!   CHARACTER*4 TDEV !(tokamak id)
!   integer ndim !(0,1,2,3. 0 for scalar, 1 for 1d array, 2 for 2d array, 3 for 3darray)

!   CHARACTER*10 SHDATE !(shot date)

  type transp_ufiles_3d_data
    integer nsc   !number of scalars 

    REAL,allocatable,dimension(:) :: SCVAL ! (scalar data array)

    CHARACTER(len=30),allocatable,dimension(:) :: SCLAB ! (scalar labels array)
    character(len=20),allocatable,dimension(:) :: labels
    character(len=10),allocatable,dimension(:) :: unitss


    CHARACTER(len=30) :: XLAB, YLAB, ZLAB, FLAB ! (label arrays for X, Y, and  F) ( CHARACTER*30 XLAB,FLAB )
    character(len=20) :: labelx,labely,labelz,labelf
    character(len=10) :: unitsx,unitsy,unitsz,unitsf

    integer iproc  ! process code-- historical, use a value of zero.

    integer nx, ny, nz, nf1, nf2   ! actual number of values in X and Y and Z, first and second dimension of F

    REAL,allocatable,dimension(:) :: X, Y, Z !(arrays to hold x,y,z)

    REAL,allocatable,dimension(:,:,:) :: F !(array to hold F)
  end type transp_ufiles_3d_data

!   CHARACTER*160 comment
end module transp_ufiles_3d

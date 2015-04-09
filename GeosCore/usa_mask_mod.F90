!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: usa_mask
!
! !DESCRIPTION: Module USA_MASK contains the mask for NEI08 and NE11
!\\
!\\
! !INTERFACE: 
!
      MODULE USA_MASK_MOD
! 
! !USES:
!
      IMPLICIT NONE
#     include "netcdf.inc"
      PRIVATE
!
! !PUBLIC DATA MEMBERS:
!
      REAL*8, PUBLIC, ALLOCATABLE :: USA_MASK(:,:)
      REAL*8, PUBLIC, ALLOCATABLE, TARGET :: USA_MASK_TMP(:,:)

      REAL*8,  ALLOCATABLE :: XEDGE_NEI11(:)
      REAL*8,  ALLOCATABLE :: YEDGE_NEI11(:)
      REAL*8,  ALLOCATABLE :: XEDGE_MODELG(:)
      REAL*8,  ALLOCATABLE :: YEDGE_MODELG(:)
!
! !PUBLIC MEMBER FUNCTIONS:
!

      PUBLIC  :: GET_MASK_FORFIRE
      PRIVATE :: READ_USA_MASK
!
! !REMARKS:
!     
! !REVISION HISTORY:
!  12 Feb 2013 - K. Travis   - initial version, adapted from Aaron von 
!                              Donkelaar's NEI05
!  28 Jun 2013 - R. Yantosca - Now read data from global data paths
!EOP

      CONTAINS
!BOP
!------------------------------------------------------------------------------
! !IROUTINE: GET_MASK_FORFIRE
! !DESCRIPTION: Subroutine GET_MASK_FORFIRE initializes the mask for the fire injection routine

      SUBROUTINE GET_MASK_FORFIRE( Input_Opt, RC )
! !INPUT PARAMETERS:
!
        USE GIGC_Input_Opt_Mod, ONLY : OptInput
        USE GIGC_ErrCode_Mod
        USE CMN_SIZE_MOD    ! Size parameters
        USE ERROR_MOD,   ONLY : ALLOC_ERR
        USE GLOBAL_GRID_MOD, ONLY : GET_IIIPAR, GET_JJJPAR
       
        TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
        INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?!
        LOGICAL,        SAVE          :: FIRST = .TRUE.
        INTEGER                       :: IIIPAR0, JJJPAR0

      ! Get global longitude/latitude extent [# of boxes]
#if   defined( GRID05x0666 ) || defined( GRID025x03125 )

      ! Nested grids utilize global longitude and latidude extent
      ! parameters IIIPAR and JJJPAR (from global_grid_mod)
        IIIPAR0 = GET_IIIPAR()
        JJJPAR0 = GET_JJJPAR()

#else

      ! Global grids utilize window longitude and latidude extent
      ! parameters IIPAR and JJPAR (from CMN_SIZE_mod)
        IIIPAR0 = IIPAR
        JJJPAR0 = JJPAR

#endif

       ! First-time initialization
        IF ( FIRST ) THEN
         ! Assume success
           RC        =  GIGC_SUCCESS
           ! allocate and read USA Mask
           ALLOCATE( USA_MASK( IIPAR, JJPAR ), STAT=RC )
           IF ( RC /= 0 ) CALL ALLOC_ERR( 'USA_MASK' )
           USA_MASK = 0d0
           ALLOCATE( USA_MASK_TMP( IIIPAR0, JJJPAR0 ), STAT=RC )
           IF ( RC /= 0 ) CALL ALLOC_ERR( 'USA_MASK' )
           USA_MASK_TMP = 0d0
           
           ! Allocate array to hold NEI11 grid box lon edges
           ALLOCATE( XEDGE_NEI11( I01x01+1 ), STAT=RC )
           IF ( RC /= 0 ) CALL ALLOC_ERR( 'XEDGE_NEI11' )
           XEDGE_NEI11 = 0.d0
           
           ! Allocate array to hold NEI11 grid box lat edges
           ALLOCATE( YEDGE_NEI11( J01x01+1 ), STAT=RC )
           IF ( RC /= 0 ) CALL ALLOC_ERR( 'YEDGE_NEI11' )
           YEDGE_NEI11 = 0.d0
           
           ! Allocate array to hold GEOS-Chem grid box lon edges
           ALLOCATE( XEDGE_MODELG( IIIPAR0+1 ), STAT=RC )
           IF ( RC /= 0 ) CALL ALLOC_ERR( 'XEDGE_MODELG' )
           XEDGE_MODELG = 0.d0
           
           ! Allocate array to hold GEOS-Chem grid box lat edges
           ALLOCATE( YEDGE_MODELG( JJJPAR0+1 ), STAT=RC )
           IF ( RC /= 0 ) CALL ALLOC_ERR( 'YEDGE_MODELG' )
           YEDGE_MODELG = 0.d0
           
           CALL READ_USA_MASK( Input_Opt )
           FIRST = .FALSE.
        ENDIF
        ! Return to calling program
      END SUBROUTINE GET_MASK_FORFIRE
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_usa_mask
!
! !DESCRIPTION: Subroutine READ\_USA\_MASK reads the mask for NEI data  
!\\
!\\
! !INTERFACE:
      
      SUBROUTINE READ_USA_MASK( Input_Opt )
!
! !USES:
!     
      ! Reference to F90 modules
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE BPCH2_MOD,       ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE LOGICAL_MOD,     ONLY : LCAC,            LBRAVO
      USE DIRECTORY_MOD,   ONLY : DATA_DIR_1x1
      USE REGRID_A2A_MOD,  ONLY : DO_REGRID_A2A, MAP_A2A
      USE TRANSFER_MOD,    ONLY : TRANSFER_2D
      USE GRID_MOD,        ONLY : GET_XOFFSET
      USE GRID_MOD,        ONLY : GET_YOFFSET
      USE GLOBAL_GRID_MOD, ONLY : GET_XEDGE_G, GET_YEDGE_G
      USE GRID_MOD,        ONLY : GET_XEDGE, GET_YEDGE
      USE GLOBAL_GRID_MOD, ONLY : GET_IIIPAR, GET_JJJPAR

      USE CMN_SIZE_MOD         ! Size parameters

      USE m_netcdf_io_open     
      USE m_netcdf_io_read
      USE m_netcdf_io_readattr
      USE m_netcdf_io_close
      USE m_netcdf_io_get_dimlen
      
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !REMARKS:
!     
! !REVISION HISTORY: 
!  20 Oct 2009 - P. Le Sager - init
!  26 Oct 2009 - P. Le Sager - new masks
!  13 Mar 2012 - M. Cooper   - Changed regrid algorithm to map_a2a
!  24 May 2012 - R. Yantosca - Fixed minor bugs in map_a2a implementation
!  15 Aug 2012 - M. Payer    - Fixed minor bugs in regridding of mask; Also
!                              set mask to 1 if greater than 0 (L. Murray)
!  24 Aug 2012 - R. Yantosca - DO_REGRID_A2A now reads netCDF input file
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*4             :: ARRAY2(900,400)
      REAL*8, TARGET     :: GEOS_MASK(I01x01,J01x01)
      CHARACTER(LEN=255) :: FILENAME1, FILENAME2
      CHARACTER(LEN=255) :: LLFILENAME
      REAL*8, POINTER    :: INGRID(:,:) => NULL()
      INTEGER            :: st2d(2), ct2d(2)
      INTEGER            :: fId1, I0, J0, I, J
      LOGICAL,  SAVE     :: FIRST = .TRUE.
      INTEGER            :: IIIPAR0, JJJPAR0
      REAL*4             :: DEG2RAD
      !=================================================================
      ! Mask specific to NEI2011 and NEI2008 data
      !=================================================================
      IF ( FIRST ) THEN
         CALL INIT_USA_MASK( Input_Opt )
      ! Get global longitude/latitude extent [# of boxes]
#if   defined( GRID05x0666 ) || defined( GRID025x03125 )
      ! Nested grids utilize global longitude and latidude extent
      ! parameters IIIPAR and JJJPAR (from global_grid_mod)
        IIIPAR0 = GET_IIIPAR()
        JJJPAR0 = GET_JJJPAR()

#else
      ! Global grids utilize window longitude and latidude extent
      ! parameters IIPAR and JJPAR (from CMN_SIZE_mod)
        IIIPAR0 = IIPAR
        JJJPAR0 = JJPAR
#endif

         DEG2RAD = (4. * ATAN(1.) ) /180.

         ! Define NEI11 grid box lat and lon edges
         XEDGE_NEI11( 1 ) = -180.d0
         DO I = 2,I01x01 +1
            XEDGE_NEI11( I ) = XEDGE_NEI11( I-1 ) + 1.d-1
         END DO

         YEDGE_NEI11( 1 ) = -89.975d0
         DO J = 2, J01x01+1
            YEDGE_NEI11( J ) = YEDGE_NEI11( J-1 ) + 1.d-1
         END DO

         DO J = 1,J01x01+1
            YEDGE_NEI11( J ) = SIN( YEDGE_NEI11( J ) * DEG2RAD)
         END DO

         ! Define global grid box lat and lon edges at model resolution
#if   defined( GRID05x0666 ) || defined( GRID025x03125 )

         DO I = 1,IIIPAR0+1
            XEDGE_MODELG( I ) = GET_XEDGE_G ( I )
         END DO

         DO J = 1,JJJPAR0+1
            YEDGE_MODELG( J ) = GET_YEDGE_G ( J )
         END DO

         DO J = 1,JJJPAR0+1
            YEDGE_MODELG( J ) = SIN( YEDGE_MODELG( J ) * DEG2RAD)
         END DO

#else

         DO I = 1,IIIPAR0+1
            XEDGE_MODELG( I ) = GET_XEDGE( I, 1, 1 )
         END DO

         DO J = 1,JJJPAR0+1
            YEDGE_MODELG( J ) = GET_YEDGE( 1, J, 1 )
         END DO

         DO J = 1,JJJPAR0+1
            YEDGE_MODELG( J ) = SIN( YEDGE_MODELG( J ) * DEG2RAD)
         END DO

#endif
     
         FIRST = .FALSE.
      ENDIF            

      FILENAME2 = '/as/scratch/krt/NEI11/VERYNESTED/USA_LANDMASK_NEI2011_0.1x0.1.nc'

      ! Echo info
      WRITE( 6, 200 ) TRIM( FILENAME2 )
200   FORMAT( '     - READ_USA_MASK: Reading ', a )
      WRITE(*,*) 'WARNING - NEI11/NEI08 CONTAINS EMISSIONS IN CANADA, MEXICO, and OVER WATER'
      WRITE(*,*) 'TO GET JUST U.S. TOTALS, USE A MASK JUST FOR THE U.S.'

      I0    = GET_XOFFSET( GLOBAL=.TRUE. )
      J0    = GET_YOFFSET( GLOBAL=.TRUE. )

      ! Allocate start and count arrays
      st2d = (/1, 1/)
      ct2d = (/900, 400/)
      ! Open and read model_ready data from netCDF file - wkday
      CALL Ncop_Rd(fId1, TRIM(FILENAME2))
      Call NcRd(ARRAY2, fId1, 'LANDMASK', st2d,  ct2d )        !Start andCount lat/lon
      ! Close netCDF file
      CALL NcCl( fId1 )

      GEOS_MASK = 0.0
      ! Cast to REAL*8 before regridding
      GEOS_MASK(402:1301,1101:1500) = ARRAY2(:,:)
      ! File with lat/lon edges for regridding
      LLFILENAME = TRIM( DATA_DIR_1x1) // &
           'MAP_A2A_Regrid_201203/MAP_A2A_latlon_generic01x01.nc'
      ! Regrid from GEOS 0.1x0.1 --> current model resolution [unitless]
      INGRID => GEOS_MASK(:,:)

      CALL MAP_A2A( I01x01, J01x01, XEDGE_NEI11, &
           YEDGE_NEI11, INGRID, IIIPAR0,  &
           JJJPAR0,     XEDGE_MODELG, YEDGE_MODELG, &
           USA_MASK_TMP,         0,            0 )
      ! Free pointer
      NULLIFY( INGRID )
       
      ! Add offset
!$OMP PARALLEL DO        &
!$OMP DEFAULT( SHARED )  &
!$OMP PRIVATE( I, J )
      DO J=1,JJPAR
         DO I=1,IIPAR
            USA_MASK(I,J)=USA_MASK_TMP(I+I0,J+J0)
         END DO
      END DO
!$OMP END PARALLEL DO

      WRITE(*,*) 'READ USA MASK!'
      WHERE ( USA_MASK > 0D0 ) USA_MASK = 1D0
      ! Return to calling program
      END SUBROUTINE READ_USA_MASK
!------------------------------------------------------------------------------
!EOC


!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_nei2011_anthro
!
! !DESCRIPTION: Subroutine INIT\_NEI2011\_ANTHRO allocates and zeroes all 
!  module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_USA_MASK(  Input_Opt )
!
! !USES:
! 
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE TIME_MOD,    ONLY : ITS_A_NEW_MONTH
      USE CMN_SIZE_MOD    ! Size parameters
!
! !INPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !REVISION HISTORY:
!  02 Mar 2012 - R. Yantosca - Remove A_CM2 array
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt,  RC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      !=================================================================
      ! INIT_USA_MASK begins here!
      !=================================================================            
      !--------------------------------------------------
      ! Allocate and zero arrays for emissions
      !--------------------------------------------------
      IF ( .not. Input_Opt%LBIOMASS ) THEN
      ! allocate and read USA Mask
         ALLOCATE( USA_MASK( IIPAR, JJPAR ) )
         !IF ( RC /= 0 ) CALL ALLOC_ERR( 'USA_MASK' )
         USA_MASK = 0d0
         
         CALL READ_USA_MASK( Input_Opt )
      ENDIF
    END SUBROUTINE INIT_USA_MASK

 
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_nei2011anthro
!
! !DESCRIPTION: Subroutine CLEANUP\_NEI2011\_ANTHRO deallocates all module 
!  arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_USA_MASK
!
! !REVISION HISTORY: 
!  01 Mar 2012 - R. Yantosca - Remove reference to A_CM2 array
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_NEIO2011_ANTHRO begins here!
      !=================================================================
      ! USA mask
      IF ( ALLOCATED( USA_MASK_TMP) ) DEALLOCATE( USA_MASK_TMP )
      IF ( ALLOCATED( USA_MASK) ) DEALLOCATE( USA_MASK )
      
      END SUBROUTINE CLEANUP_USA_MASK
!EOC
      END MODULE USA_MASK_MOD

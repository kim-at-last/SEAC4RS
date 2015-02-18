!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: nei2011_anthro_mod
!
! !DESCRIPTION: Module NEI2011\_ANTHRO\_MOD contains variables and routines to 
!  read the NEI2011 anthropogenic emissions.   
!\\
!\\
! !INTERFACE: 
!
      MODULE NEI2011_ANTHRO_MOD
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
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: CLEANUP_NEI2011_ANTHRO
      PUBLIC  :: EMISS_NEI2011_ANTHRO
      PUBLIC  :: GET_NEI2011_ANTHRO
      PUBLIC  :: GET_MASK_FORFIRE2
      !--------------------------------------
      ! Leave for future use (bmy, 12/3/09)
      !PUBLIC  :: GET_NEI2005_MASK
      !--------------------------------------
!
! !PRIVATE MEMBER FUNCTIONS:
! No longer need scaling except for future emissions
      PRIVATE :: NEI2011_SCALE_FUTURE
      PRIVATE :: INIT_NEI2011_ANTHRO
      PRIVATE :: TOTAL_ANTHRO_TG
      PRIVATE :: READ_NEI2011_MASK
!
! !REMARKS:
!     
! !REVISION HISTORY:
!  12 Feb 2013 - K. Travis   - initial version, adapted from Aaron von 
!                              Donkelaar's NEI05
!  28 Jun 2013 - R. Yantosca - Now read data from global data paths
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      ! Arrays for emissions (lat/lon/lev/hrs)
      REAL*8,  ALLOCATABLE, TARGET :: TMP(:,:,:,:)
      REAL*8,  ALLOCATABLE, TARGET :: TMPARR(:,:)
      REAL*8,  ALLOCATABLE         :: T_CO_MON(:)
      REAL*8,  ALLOCATABLE         :: T_NO_MON(:)
      REAL*8,  ALLOCATABLE         :: T_NO2_MON(:)
      REAL*8,  ALLOCATABLE         :: T_HNO2_MON(:)
      REAL*8,  ALLOCATABLE         :: T_SO2_MON(:)
      REAL*8,  ALLOCATABLE         :: T_SO4_MON(:)
      REAL*8,  ALLOCATABLE         :: T_OC_MON(:)
      REAL*8,  ALLOCATABLE         :: T_BC_MON(:)
      REAL*8,  ALLOCATABLE         :: T_NH3_MON(:)
      REAL*8,  ALLOCATABLE         :: T_C2H6_MON(:)
!      REAL*8,  ALLOCATABLE         :: T_C2H4_MON(:)
      REAL*8,  ALLOCATABLE         :: T_C3H8_MON(:)
      REAL*8,  ALLOCATABLE         :: T_ALD2_MON(:)
      REAL*8,  ALLOCATABLE         :: T_RCHO_MON(:)
      REAL*8,  ALLOCATABLE         :: T_CH2O_MON(:)
      REAL*8,  ALLOCATABLE         :: T_ALK4_MON(:)
      REAL*8,  ALLOCATABLE         :: T_ACET_MON(:)
      REAL*8,  ALLOCATABLE         :: T_MACR_MON(:)
      REAL*8,  ALLOCATABLE         :: T_MEK_MON(:)
      REAL*8,  ALLOCATABLE         :: T_PRPE_MON(:)
      REAL*8,  ALLOCATABLE         :: T_TOLU_MON(:)
      REAL*8,  ALLOCATABLE         :: T_BENZ_MON(:)
      REAL*8,  ALLOCATABLE         :: T_XYLE_MON(:)
      REAL*8,  ALLOCATABLE         :: T_EOH_MON(:)
      REAL*8,  ALLOCATABLE         :: T_MOH_MON(:)

      REAL*8,  ALLOCATABLE        :: CO(:,:,:,:) 
      REAL*8,  ALLOCATABLE        :: NO(:,:,:,:) 
      REAL*8,  ALLOCATABLE        :: NO2(:,:,:,:)
      REAL*8,  ALLOCATABLE        :: HNO2(:,:,:,:)
      REAL*8,  ALLOCATABLE        :: SO2a(:,:,:,:) 
      REAL*8,  ALLOCATABLE        :: SO2b(:,:,:,:) 
      REAL*8,  ALLOCATABLE        :: NH3(:,:,:,:) 
      REAL*8,  ALLOCATABLE        :: ALD2(:,:,:,:)
      REAL*8,  ALLOCATABLE        :: RCHO(:,:,:,:)
      REAL*8,  ALLOCATABLE        :: BENZ(:,:,:,:)
      !REAL*8,  ALLOCATABLE        :: CH4(:,:,:,:)
      REAL*8,  ALLOCATABLE        :: C2H6(:,:,:,:) 
      REAL*8,  ALLOCATABLE        :: PRPEa(:,:,:,:)
      REAL*8,  ALLOCATABLE        :: PRPEb(:,:,:,:)
      REAL*8,  ALLOCATABLE        :: ALK4(:,:,:,:)
      REAL*8,  ALLOCATABLE        :: TOLU(:,:,:,:)
      REAL*8,  ALLOCATABLE        :: XYLE(:,:,:,:)
      !REAL*8,  ALLOCATABLE        :: C2H4(:,:,:,:)
      REAL*8,  ALLOCATABLE        :: MOH(:,:,:,:) 
      REAL*8,  ALLOCATABLE        :: EOH(:,:,:,:) 
      REAL*8,  ALLOCATABLE        :: CH2O(:,:,:,:) 
      REAL*8,  ALLOCATABLE        :: OCPO(:,:,:,:) 
      REAL*8,  ALLOCATABLE        :: BCPO(:,:,:,:)
      REAL*8,  ALLOCATABLE        :: SO4(:,:,:,:) 
      REAL*8,  ALLOCATABLE        :: MACR(:,:,:,:) 
      REAL*8,  ALLOCATABLE        :: MEK(:,:,:,:) 
      REAL*8,  ALLOCATABLE        :: C3H8(:,:,:,:) 
      REAL*8,  ALLOCATABLE        :: ACET(:,:,:,:) 

      ! Shadow logical variables from Input_Opt
      LOGICAL                      :: LBRAVO
      LOGICAL                      :: LCAC
      LOGICAL                      :: LFUTURE
      LOGICAL                      :: LNEI11
!
!
! !DEFINED PARAMETERS:
!
      REAL*8,  PARAMETER   :: SEC_IN_YEAR  = 86400d0 * 365.25d0

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_nei2011_anthro
!
! !DESCRIPTION: Function GET\_NEI2011\_ANTHRO returns the NEI2011
!  emission for GEOS-Chem grid box (I,J) and tracer N and hour IH.  
!  Emissions can be returned in units of [kg/s] or [molec/cm2/s].
!\\ (krt, 2/10/13), now need IH
!\\
! !INTERFACE:
!
      FUNCTION GET_NEI2011_ANTHRO( I, J, L, IH, N ) RESULT( VALUE )
!
! !USES:
!
      USE GRID_MOD,     ONLY : GET_AREA_M2
      USE TRACERID_MOD, ONLY : IDTCO, IDTNO,IDTNO2, IDTHNO2 
      USE TRACERID_MOD, ONLY : IDTSO2, IDTNH3, IDTMACR
      USE TRACERID_MOD, ONLY : IDTALD2, IDTRCHO, IDTC2H6
      USE TRACERID_MOD, ONLY : IDTPRPE, IDTALK4!, IDTC2H4
      USE TRACERID_MOD, ONLY : IDTBENZ, IDTTOLU, IDTXYLE
      USE TRACERID_MOD, ONLY : IDTSO4, IDTCH2O
      USE TRACERID_MOD, ONLY : IDTOCPO,  IDTBCPO, IDTMEK
      USE TRACERID_MOD, ONLY : IDTMOH, IDTEOH!, IDTCH4
      USE TRACERID_MOD, ONLY : IDTMEK, IDTC3H8, IDTACET !added 9/24/14, krt 
!
! !INPUT PARAMETERS: 
!
      ! Longitude, latitude, hour, and tracer indices
      INTEGER, INTENT(IN)           :: I, J, L, IH, N
 
      ! ALWAYS returns emissions in [molec/cm2/s]

! !RETURN VALUE:
!     
      ! Emissions output
      REAL*8                        :: VALUE  

! !REMARKS:
!  SET LEVEL to 0   
!
! !REVISION HISTORY: 
!  12 Feb 2013 - K. Travis - initial version  
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL :: DO_KGS, DO_MCS

      !=================================================================
      ! GET_NEI2011_ANTHRO begins here!
      !=================================================================

      ! Initialize
      IF ( N == IDTCO ) THEN ! Second IF
         ! CO [kg/m2/s] to [molec/cm2/s]
         VALUE = CO(I,J,L,IH)
      ELSE IF ( N == IDTNO ) THEN
         ! NO[molec/cm2/s]
         VALUE = NO(I,J,L,IH)   
      ELSE IF ( N == IDTNO2 ) THEN
         ! NO2[molec/cm2/s]
         VALUE = NO2(I,J,L,IH)
      ELSE IF ( N == IDTHNO2 ) THEN
         ! HNO2[molec/cm2/s]
         VALUE = HNO2(I,J,L,IH)
      ELSE IF ( N == IDTSO2 ) THEN
         ! SO2 [molec/cm2/s]
         VALUE = SO2a(I,J,L,IH) + SO2b(I,J,L,IH)
      ELSE IF ( N == IDTNH3 ) THEN
         ! NH3 [molec/cm2/s]
         VALUE = NH3(I,J,L,IH)
      ELSE IF ( N == IDTALD2 ) THEN
         ! [molec/cm2/s]
         VALUE = ALD2(I,J,L,IH)
      ELSE IF ( N == IDTRCHO ) THEN
         ! [molec/cm2/s]
         VALUE = RCHO(I,J,L,IH)
      ELSE IF ( N == IDTBENZ ) THEN
         ! [molec/cm2/s]
         VALUE = BENZ(I,J,L,IH)
      ELSE IF ( N == IDTC2H6 ) THEN
         ! [molec/cm2/s]
         VALUE = C2H6(I,J,L,IH)
      ELSE IF ( N == IDTPRPE ) THEN
         ! [molec/cm2/s]
         VALUE = PRPEa(I,J,L,IH)+PRPEb(I,J,L,IH)
      ELSE IF ( N == IDTALK4 ) THEN
            ! [molec/cm2/s]
         VALUE = ALK4(I,J,L,IH)
      ELSE IF ( N == IDTC3H8 ) THEN
            ! [molec/cm2/s]
         VALUE = C3H8(I,J,L,IH)
      ELSE IF ( N == IDTACET ) THEN
            ! [molec/cm2/s]
         VALUE = ACET(I,J,L,IH)
      ELSE IF ( N == IDTMEK ) THEN
            ! [molec/cm2/s]
         VALUE = MEK(I,J,L,IH)
      ! NOT EMITTING MACR AT THIS TIME
      ELSE IF ( N == IDTTOLU ) THEN
         ! [molec/cm2/s]
         VALUE = TOLU(I,J,L,IH)
      ELSE IF ( N == IDTXYLE ) THEN
         ! [molec/cm2/s]
         VALUE = XYLE(I,J,L,IH)
      ELSE IF ( N == IDTCH2O ) THEN
         ! [molec/cm2/s]
         VALUE = CH2O(I,J,L,IH)
      ELSE IF ( N == IDTOCPO ) THEN
            ! [g/cm2/s]
         VALUE = OCPO(I,J,L,IH)
      ELSE IF ( N == IDTBCPO ) THEN
         ! [g/cm2/s]
         VALUE = BCPO(I,J,L,IH)
      ELSE IF ( N == IDTSO4 ) THEN
         ! [g/cm2/s]
         VALUE = SO4(I,J,L,IH)
      ELSE IF ( N == IDTEOH ) THEN
         ! [molec/cm2/s]
         VALUE = EOH(I,J,L,IH)
      ELSE IF ( N == IDTMOH ) THEN
            ! [molec/cm2/s]
         VALUE = MOH(I,J,L,IH)
      ELSE
         ! Otherwise return a negative value to indicate
         ! that there are no NEI2011 emissions for tracer N
         VALUE = -1d0
         RETURN
         
      ENDIF ! END Second IF

      ! Return to calling program
      END FUNCTION GET_NEI2011_ANTHRO
!EOC
!BOP
!------------------------------------------------------------------------------
! !IROUTINE: GET_MASK_FORFIRE2
! !DESCRIPTION: Subroutine GET_MASK_FORFIRE initializes the mask for the fire injection routine

      SUBROUTINE GET_MASK_FORFIRE2(  RC )
! !INPUT PARAMETERS:
!
        USE GIGC_ErrCode_Mod
        USE CMN_SIZE_MOD    ! Size parameters
        USE ERROR_MOD,   ONLY : ALLOC_ERR

        INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?!
        LOGICAL,        SAVE          :: FIRST = .TRUE.

       ! First-time initialization
        IF ( FIRST ) THEN
         ! Assume success
           RC        =  GIGC_SUCCESS
           ! allocate and read USA Mask
           ALLOCATE( USA_MASK( IIPAR, JJPAR ), STAT=RC )
           IF ( RC /= 0 ) CALL ALLOC_ERR( 'USA_MASK' )
           USA_MASK = 0d0

           CALL READ_NEI2011_MASK
           FIRST = .FALSE.
        ENDIF
    ! Return to calling program
      END SUBROUTINE GET_MASK_FORFIRE2
!EOC
!------------------------------------------------------------------------------
!       Adopted from NEI05 from
!       Dalhousie University Atmospheric Compositional Analysis Group         !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emiss_nei2011_anthro
!
! !DESCRIPTION: Subroutine EMISS\_NEI2011\_ANTHRO reads the NEI2011
!  emission fields at .25 x 0.3125 resolution
!  Designed to work with IIPAR and JJPAR as long as emissions are on the
!  same nested grid.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMISS_NEI2011_ANTHRO( am_I_Root, Input_Opt, &
                                             State_Chm, RC         )
!
! !USES:
!
      USE DIRECTORY_MOD,     ONLY : DATA_DIR
      USE DIRECTORY_MOD,     ONLY : DATA_DIR_1x1
      USE GRID_MOD,          ONLY : GET_XOFFSET
      USE GRID_MOD,          ONLY : GET_YOFFSET
      USE REGRID_A2A_MOD,    ONLY : DO_REGRID_A2A
      USE LOGICAL_MOD,       ONLY : LFUTURE
      USE TIME_MOD,          ONLY : GET_YEAR, GET_MONTH, GET_DAY
      USE TIME_MOD,          ONLY : GET_HOUR
      USE TRACER_MOD,        ONLY : XNUMOL

      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE NCDF_MOD,          ONLY : NC_READ
      USE TRACERID_MOD,      ONLY : IDTCO, IDTNO, IDTHNO2, IDTNO2 
      USE TRACERID_MOD,      ONLY : IDTSO2, IDTNH3, IDTMACR, IDTACET
      USE TRACERID_MOD,      ONLY : IDTALD2, IDTRCHO, IDTC2H6, IDTMEK
      USE TRACERID_MOD,      ONLY : IDTPRPE, IDTALK4, IDTC3H8!, IDTC2H4
      USE TRACERID_MOD,      ONLY : IDTBENZ, IDTTOLU, IDTXYLE
      USE TRACERID_MOD,      ONLY : IDTSO4, IDTCH2O, IDTOCPO, IDTBCPO
      USE TRACERID_MOD,      ONLY : IDTEOH, IDTMOH

      USE CMN_SIZE_MOD          ! Size parameters
      USE CMN_O3_MOD            ! FSCALYR
      
      USE m_netcdf_io_open     
      USE m_netcdf_io_read
      USE m_netcdf_io_readattr
      USE m_netcdf_io_close
      USE m_netcdf_io_get_dimlen
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?

!
! !REVISION HISTORY:
!  16 Feb 2013 - K. Travis   - initial version
!  28 Jun 2013 - R. Yantosca - Now reads data from global data path
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE              :: FIRST = .TRUE.
      INTEGER                    :: I, J, IH, I0, J0, THISYEAR, THISMONTH, SNo
      INTEGER                    :: THISDAY
      INTEGER                    :: L, HH, KLM, SPECIES_ID(23)
      INTEGER                    :: st3d(3), ct3d(3)
      INTEGER                    :: st4d(4), ct4d(4)
      INTEGER                    :: fId1, fId1b, fId1c, fId1d, fId1e
      INTEGER                    :: fId1f, fId1g, fId1h
      REAL*8                     :: ARRAY(225,202,24)
      REAL*8                     :: ARRAYOTH(225,202,6,24)
      REAL*8                     :: ARRAYPTN(225,202,6,24)
      REAL*8                     :: ARRAYOIL(225,202,6,24)
      REAL*8                     :: ARRAYEGU(225,202,6,24)
      REAL*8                     :: ARRAYEGUPK(225,202,6,24)
      REAL*8                     :: ARRAYC3(225,202,6,24)
      REAL*4                     :: ARRAY_NH3ag(225,202,24)
      REAL*8, TARGET             :: GEOS_NATIVE_NH3ag(225,202,24)
      REAL*8, TARGET             :: GEOS_NATIVE(I025x03125,J025x03125,6,24)
      REAL*4                     :: ScCO, ScNOx, ScSO2, ScNH3, ScPM10
      REAL*4                     :: ScNH3_Ag, ScNH3_NonAg
      REAL*4                     :: ScPM25, ScVOC
      CHARACTER(LEN=255)         :: DATA_DIR_NEI
      CHARACTER(LEN=255)         :: FILENAME, FILENAMEOTH, FILENAMEOIL
      CHARACTER(LEN=255)         :: FILENAMEPTN, FILENAMEC3, FILENAMEEGU
      CHARACTER(LEN=255)         :: FILENAMEEGUPK, LLFILENAME
      CHARACTER(LEN=24)          :: SPCLIST(23)
      CHARACTER(LEN=8)           :: SId
      CHARACTER(LEN=5)           :: SNAME
      CHARACTER(LEN=3)           :: TTMON
      CHARACTER(LEN=2)           :: FDAY, FMON
      REAL*8, POINTER            :: INGRID(:,:) => NULL()
      REAL*8, POINTER            :: OUTGRID(:,:) => NULL()

      ! For scaling NH3 agricultural emissions (jaf, 12/12/13)
      REAL*4, POINTER            :: NCARR(:,:,:) => NULL()
      REAL*8                     :: ScAgNH3_MASAGE(225,202)
      LOGICAL                    :: LSCALE2MASAGE
      CHARACTER(LEN=255)         :: DATA_DIR_NH3_ag
      CHARACTER(LEN=255)         :: FILENAME_NH3ag
      CHARACTER(LEN=255)         :: FILENAME_ScAg

      !=================================================================
      ! EMISS_NEI2011_ANTHRO begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_NEI2011_ANTHRO( am_I_Root, Input_Opt, RC )
         FIRST = .FALSE.
      ENDIF

      ! Copy values from Input_Opt
      LFUTURE   = Input_Opt%LFUTURE
      LSCALE2MASAGE = Input_Opt%LSCALE2MASAGE

      ! Get emissions year
      THISYEAR = GET_YEAR()

      ! Initialize scaling factors
      ScCO        = 1.0
      ScNOx       = 1.0
      ScPM10      = 1.0
      ScPM25      = 1.0
      ScSO2       = 1.0
      ScVOC       = 1.0
      ScNH3       = 1.0
     
      ! Apply annual scalar factor.
      ! Using EPA's National Tier1 CAPS (http://www.epa.gov/ttnchie1/trends/)
      IF ( THISYEAR .eq. 2012 ) THEN ! scale based on 2010
         ScCO   = 0.981
         ScNOx  = 0.940
         ScPM10 = 0.999
         ScPM25 = 0.995
         ScSO2  = 0.802
         ScVOC  = 0.986
         ScNH3_NonAg = 0.999
         ScNH3_Ag    = 1.0
      ELSE IF ( THISYEAR .ge. 2013 ) THEN 
         ScCO   = 0.962
         ScNOx  = 0.899
         ScPM10 = 0.998
         ScPM25 = 0.991
         ScSO2  = 0.785
         ScVOC  = 0.971
         ScNH3_NonAg = 0.998
         ScNH3_Ag    = 1.0
      ENDIF
        
      I0    = GET_XOFFSET( GLOBAL=.TRUE. )
      J0    = GET_YOFFSET( GLOBAL=.TRUE. )

      SPECIES_ID = (/ IDTCO,   IDTNO,   IDTNO2,  IDTHNO2, IDTSO2, IDTSO2, IDTNH3, &
           IDTMACR, IDTALD2, IDTRCHO, IDTC2H6, IDTCH2O, IDTPRPE,  IDTPRPE, &
           IDTALK4, IDTTOLU, IDTXYLE, IDTOCPO, IDTBCPO, IDTSO4, IDTBENZ, IDTEOH, IDTMOH /)

      SPCLIST =    (/ 'CO','NO','NO2','HONO','SO2','SULF', 'NH3', &
           'ACROLEIN','ALD2', 'ALDX', 'ETHA', 'FORM','IOLE', 'OLE', &
           'PAR', 'TOL', 'XYL', 'POC', 'PEC','PSO4','BENZENE','ETOH','MEOH'/)

      ! File with lat/lon edges for regridding
      LLFILENAME = TRIM( DATA_DIR_1x1) // &
                  'MAP_A2A_Regrid_201203/MAP_A2A_latlon_geos025x03125.nc'

      ! Base data directory
      DATA_DIR_NEI =  '/as/cache/2015-02/krt/'// &
           'NEI11_025x03125_2011'
      ! For NH3 -- files with agricultural emissions only (jaf, 12/12/13)
      ! Eventually these files will move to the data directory
      DATA_DIR_NH3_ag = '/as/cache/2015-02/krt/ag/'

     
      ! Get month
      THISMONTH = GET_MONTH()
      ! Get day
      THISDAY   = GET_DAY()
      WRITE (FDAY, '(I2.2)') THISDAY
      WRITE (FMON, '(I2.2)')  THISMONTH
      
      ! surface
      FILENAME    = TRIM( DATA_DIR_NEI ) //                        &
           TRIM(FMON) // TRIM(FDAY) // '.nc'
      ! othpt
      FILENAMEOTH  = TRIM( DATA_DIR_NEI ) //  &
           TRIM(FMON) // TRIM(FDAY) // '_othpt.nc'
      ! ptnonipm
      FILENAMEPTN = TRIM( DATA_DIR_NEI ) // &
           TRIM(FMON) // TRIM(FDAY)//  '_ptnonipm.nc'
      ! oilgas
      FILENAMEOIL  = TRIM( DATA_DIR_NEI ) //  &
           TRIM(FMON) // TRIM(FDAY) //  '_oilgas.nc'
      ! ptegu
      FILENAMEEGU  = TRIM( DATA_DIR_NEI ) //  &
           TRIM(FMON) // TRIM(FDAY) //  '_egu.nc'
      ! pteguok
      FILENAMEEGUPK  = TRIM( DATA_DIR_NEI ) //  &
           TRIM(FMON) // TRIM(FDAY) //  '_egupk.nc'
      ! c3marine
      FILENAMEC3  = TRIM( DATA_DIR_NEI ) //  &
           TRIM(FMON) // TRIM(FDAY)//  '_c3marine.nc'

      ! Allocate start and count arrays
      st3d = (/1, 1, 1/)            !Start lat/lon/time
      ct3d = (/225, 202, 24/)     !Count lat/lon/time
      st4d = (/1, 1, 1, 1/)         !Start lat/lon/time/lev
      ct4d= (/225, 202, 6, 24/)  !Count lat/lon/time/lev 
 

      ! Open netCDF files for reading
      CALL Ncop_Rd(fId1,  TRIM(FILENAME))
      CALL Ncop_Rd(fId1b, TRIM(FILENAMEOTH))     ! ptipm
      CALL Ncop_Rd(fId1c, TRIM(FILENAMEPTN))    ! ptnonipm
      CALL Ncop_Rd(fId1d, TRIM(FILENAMEC3))     ! c3marine
      CALL Ncop_Rd(fId1e, TRIM(FILENAMEOIL))     ! oilgas
      CALL Ncop_Rd(fId1f, TRIM(FILENAMEEGU))     ! egu
      CALL Ncop_Rd(fId1g, TRIM(FILENAMEEGUPK))     ! egupk

      ! Open NH3 ag files (only avail at 025x03125)
      IF ( LSCALE2MASAGE ) THEN
         FILENAME_NH3ag = TRIM(DATA_DIR_NH3_ag) // &
            'NEI11_025x03125_2011'              // &
            TRIM(FMON) // TRIM(FDAY) //  '_ag.nc'
         FILENAME_ScAg  = TRIM(DATA_DIR_NH3_ag) // &
            'MASAGE_NEI11_Ratio.geos.025x03125.nc'
         CALL Ncop_Rd(fId1h, TRIM(FILENAME_NH3ag)) ! NH3ag
      ENDIF

 100  FORMAT( '     - EMISS_NEI2011_ANTHRO:  Reading ', &
                      a, ' -> ', a )

      ! Loop over species
      DO KLM = 1, SIZE( SPECIES_ID )

         SId = SPCLIST( KLM )
         SNo = SPECIES_ID( KLM ) 

         ! Skip undefined tracers
         !IF ( SNo == 0 ) CYCLE

         ! Read variable from netCDF files
         ! Units are in kg/m2/s
         WRITE( 6, 100 )  TRIM( FILENAME ), SID
         Call NcRd(ARRAY,       fId1,  TRIM(SId), st3d, ct3d )
         Call NcRd(ARRAYEGU,    fId1f, TRIM(SId), st4d, ct4d )
         Call NcRd(ARRAYEGUPK,  fId1g, TRIM(SId), st4d, ct4d )
         ! no ship emissions of acrolein or ammonia
         IF ( TRIM(SId) .ne. 'ACROLEIN' .and. TRIM(SId) .ne. 'NH3') THEN
            Call NcRd(ARRAYC3, fId1d, TRIM(SId), st4d, ct4d )
         ELSE
            ARRAYC3  = ARRAYC3  * 0d0
         ENDIF
         IF ( TRIM(SId) .ne. 'ACROLEIN' ) THEN
            Call NcRd(ARRAYOIL,    fId1e, TRIM(SId), st4d, ct4d )
            Call NcRd(ARRAYOTH,    fId1b, TRIM(SId), st4d, ct4d )
            Call NcRd(ARRAYPTN,    fId1c, TRIM(SId), st4d, ct4d )
         ELSE
            ARRAYOIL = ARRAYOIL * 0d0
            ARRAYOTH = ARRAYOTH * 0d0
            ARRAYPTN = ARRAYPTN * 0d0
         ENDIF

         GEOS_NATIVE = 0.0d0

         ! Cast to REAL*8 before regridding
         GEOS_NATIVE(160:384,399:600,1,:) = ARRAY(:,:,:) + ARRAYOTH(:,:,1,:)  &
              + ARRAYPTN(:,:,1,:) + ARRAYC3(:,:,1,:) + ARRAYOIL(:,:,1,:) &
              + ARRAYEGU(:,:,1,:) + ARRAYEGUPK(:,:,1,:)
         GEOS_NATIVE(160:384,399:600,2,:) = ARRAYOTH(:,:,2,:)  &
              + ARRAYPTN(:,:,2,:) + ARRAYC3(:,:,2,:)+ ARRAYOIL(:,:,2,:) &
              + ARRAYEGU(:,:,2,:) + ARRAYEGUPK(:,:,2,:)
         GEOS_NATIVE(160:384,399:600,3,:) =  ARRAYOTH(:,:,3,:)  &
              + ARRAYPTN(:,:,3,:) + ARRAYC3(:,:,3,:)+ ARRAYOIL(:,:,3,:) &
              + ARRAYEGU(:,:,3,:) + ARRAYEGUPK(:,:,3,:)
         GEOS_NATIVE(160:384,399:600,4,:) =  ARRAYOTH(:,:,4,:)  &
              + ARRAYPTN(:,:,4,:) + ARRAYC3(:,:,4,:)+ ARRAYOIL(:,:,4,:) &
              + ARRAYEGU(:,:,4,:) + ARRAYEGUPK(:,:,4,:)
         GEOS_NATIVE(160:384,399:600,5,:) =  ARRAYOTH(:,:,5,:)  &
              + ARRAYPTN(:,:,5,:) + ARRAYC3(:,:,5,:)+ ARRAYOIL(:,:,5,:) &
              + ARRAYEGU(:,:,5,:) + ARRAYEGUPK(:,:,5,:)
         GEOS_NATIVE(160:384,399:600,6,:) =  ARRAYOTH(:,:,6,:)  &
              + ARRAYPTN(:,:,6,:) + ARRAYC3(:,:,6,:)+ ARRAYOIL(:,:,6,:) &
              + ARRAYEGU(:,:,6,:) + ARRAYEGUPK(:,:,6,:)
                 

         ! Special case for NH3 emissions -- scale agricultural
         ! component based on MASAGE monthly gridded values from Paulot
         ! et al., 2013 (jaf, 12/10/13)
         IF ( LSCALE2MASAGE .and. TRIM(SId) == 'NH3') THEN

            ! Read/close ag files
            CALL NcRd(ARRAY_NH3ag, fId1h, TRIM(SId), st3d, ct3d )
            CALL NcCl( fId1h )

            ! Cast to REAL*8
            GEOS_NATIVE_NH3ag = ARRAY_NH3ag

            ! Subtract agricultural component from total
            ! This is a global array so only put in the US segment,
            ! as above
            GEOS_NATIVE(160:384,399:600,1,:) = &
               GEOS_NATIVE(160:384,399:600,1,:) - &
               GEOS_NATIVE_NH3ag(:,:,:)

            ! Read scaling factor (ratio of MASAGE to NEI08
            CALL NC_READ( NC_PATH = TRIM(FILENAME_ScAg),     &
                          PARA = 'ratio', ARRAY = NCARR,     &
                          YEAR = 2011,    MONTH = THISMONTH, &
                          DAY = 01,       VERBOSE = .FALSE.    )

            ! Cast to REAL*8
            ScAgNH3_MASAGE = NCARR(:,:,1)

            ! Deallocate ncdf-array
            IF ( ASSOCIATED ( NCARR ) ) DEALLOCATE ( NCARR )

            ! Scale agricultural component to MASAGE monthly totals
            DO HH = 1, 24
               GEOS_NATIVE_NH3ag(:,:,HH) =              &
                  GEOS_NATIVE_NH3ag(:,:,HH) * ScAgNH3_MASAGE
            ENDDO

            ! Add scaled agricultural component back to total and
            ! apply interannual scaling factors
            ! This is a global array so only put in the US segment,
            ! as above
            GEOS_NATIVE(160:384,399:600,1,:) = &
               GEOS_NATIVE(160:384,399:600,1,:) * ScNH3_NonAg + &
               GEOS_NATIVE_NH3ag(:,:,:) * ScNH3_Ag

            ELSE
               ! If we can't separate out the agricultural component
               ! then just apply the single
               ! annual scaling factor to NH3 emissions.
               GEOS_NATIVE = GEOS_NATIVE * ScNH3_NonAg
            ENDIF

         DO L=1,6
            DO HH=1,24
               ! Point to array slices
               INGRID  => GEOS_NATIVE(:,:,L,HH)
               OUTGRID => TMPARR(:,:)
               ! Regrid
               CALL DO_REGRID_A2A( LLFILENAME, I025x03125,  J025x03125, &
                    INGRID,     OUTGRID, IS_MASS=0, &
                    netCDF=.TRUE.                   )
               ! Free pointers
               NULLIFY( INGRID, OUTGRID )
                  ! Add offset
!$OMP PARALLEL DO        &
!$OMP DEFAULT( SHARED )  &
!$OMP PRIVATE( I, J )
               DO J=1,JJPAR
                  DO I=1,IIPAR
                     TMP(I,J,L,HH)=TMPARR(I+I0,J+J0)
                  END DO
               END DO
!$OMP END PARALLEL DO
            ENDDO
         ENDDO

         ! Begin loopthrough tracers
         ![kg/m2/s] to [molec/cm2/s]
         IF ( TRIM(SId) == 'CO') THEN !CO
            CO = ( TMP * XNUMOL(IDTCO) / 1E4 )  * ScCO
            ! USE VALUES FROM BORBON ET AL, 2013 (JGR) FOR UNKNOWN VOCS
            !WRITE(*,*) 'USING BORBON ET AL, 2013 for C3H8, ACET'
            !C3H8 = CO * 1.12E-02 * 3.0 ! FOR 3 CARBONS
            !ACET = CO * 1.18E-02 * 3.0 ! FOR 3 CARBONS
            !MVK = CO * 2.40E-04
         ELSEIF ( TRIM(SId) == 'NO') THEN !NO  !convert from kg/m2/s of NO
            NO  = ( TMP * XNUMOL(IDTNO) / 1E4 )  * ScNOx
         ELSEIF ( TRIM(SId) == 'NO2') THEN !NO2
            NO2 = ( TMP * XNUMOL(IDTNO2) / 1E4 ) * ScNOx
         ELSEIF( TRIM(SId) == 'HONO') THEN !HNO2
            HNO2 = ( TMP * XNUMOL(IDTHNO2)/ 1E4 ) * ScNOx
         ELSEIF ( TRIM(SId) == 'SO2') THEN !SO2
            SO2a = ( TMP * XNUMOL(IDTSO2) / 1E4 ) * ScSO2
         ELSEIF ( TRIM(SId) == 'SULF') THEN !SO2
            SO2b = ( TMP * XNUMOL(IDTSO2) / 1E4 ) * ScSO2
         ELSEIF ( TRIM(SId) == 'NH3') THEN !NH3
            NH3 = ( TMP * XNUMOL(IDTNH3) / 1E4 ) !Scaling above
         ELSEIF ( TRIM(SId) == 'ALD2') THEN !ALD2
            ALD2 = ( TMP * XNUMOL(IDTALD2)/ 1E4 ) * ScVOC
         ELSEIF ( TRIM(SId) == 'ALDX') THEN !RCHO
            RCHO = ( TMP * XNUMOL(IDTRCHO)/ 1E4 ) * ScVOC
         ELSEIF ( TRIM(SId) == 'ACROLEIN') THEN !MACR
            MACR = ( TMP * XNUMOL(IDTMACR)/ 1E4 ) * ScVOC
         ELSEIF ( TRIM(SId) == 'ETHA') THEN !C2H6
            C2H6 = ( TMP * XNUMOL(IDTC2H6)/1E4 ) * ScVOC
         ELSEIF ( TRIM(SId) == 'IOLE' ) THEN !PRPE
            PRPEa = ( TMP * XNUMOL(IDTPRPE)/1E4 ) * ScVOC
         ELSEIF ( TRIM(SId) == 'OLE' ) THEN !PRPE
            PRPEb = ( TMP * XNUMOL(IDTPRPE)/1E4 ) * ScVOC
         ELSEIF ( TRIM(SId) == 'PAR' ) THEN 
            ! I used the graph of the top 50 VOCs from the following paper:
            ! http://www.epa.gov/ttnchie1/software/speciate/atmospheric.pdf
            ! And determined the PAR weighting from the following report:
            ! http://www.camx.com/publ/pdfs/cb05_final_report_120805.pdf
            ! I then determined the fraction by carbon of PAR for each of the
            !  following species.
            ! Note that it will not add to 100% since some species are double
            ! counted (benzene) (krt, 2/3/15)
            ACET = ( TMP * 0.06 * XNUMOL(IDTACET)/1E4 ) * ScVOC
            C3H8 = ( TMP * 0.03 * XNUMOL(IDTC3H8)/1E4 ) * ScVOC
            MEK  = ( TMP * 0.02 * XNUMOL(IDTMEK)/1E4  ) * ScVOC
            ALK4 = ( TMP * 0.87 * XNUMOL(IDTALK4)/1E4 ) * ScVOC
         ELSEIF ( TRIM(SId) == 'BENZENE' ) THEN !BENZ
            BENZ =  ( TMP * XNUMOL(IDTALD2)*3 / 1E4) *ScVOC
         ELSEIF ( TRIM(SId) == 'TOL' ) THEN !TOLU
            TOLU =  ( TMP * XNUMOL(IDTALD2)*3.5 / 1E4) *ScVOC
         ELSEIF ( TRIM(SId) == 'XYL') THEN !XYLE
            XYLE  =  ( TMP * XNUMOL(IDTALD2)*4 / 1E4) *ScVOC
         ELSEIF ( TRIM(SId) == 'FORM') THEN !CH2O
            CH2O =  ( TMP * XNUMOL(IDTCH2O) / 1E4 ) * ScVOC
         ELSEIF ( TRIM(SId) == 'PSO4') THEN !SO4 - scale to SO2
            SO4 =  ( TMP * 1E3 / 1E4 ) *ScSO2
         ELSEIF ( TRIM(SId) == 'POC') THEN !OCPO
            OCPO =  ( TMP * 1E3/ 1E4 ) * ScPM25
         ELSEIF ( TRIM(SId) == 'PEC') THEN !BCPO
            BCPO  =  ( TMP * 1E3 / 1E4 ) * ScPM25
            ! Species not currently used
!         ELSEIF ( TRIM(SId) == 'C2H4') THEN !C2H4 - no scaling
!            C2H4       = TMP
         ELSEIF ( TRIM(SId) == 'MEOH') THEN !MOH - no scaling
            MOH   =  ( TMP * XNUMOL(IDTALD2)*0.5/ 1E4 ) * ScVOC
         ELSEIF ( TRIM(SId) == 'ETOH') THEN !EOH - no scaling
            EOH   = ( TMP * XNUMOL(IDTALD2)/ 1E4 ) * ScVOC
         ENDIF
      ENDDO

      ! Close netCDF file
      CALL NcCl( fId1  )
      CALL NcCl( fId1b )
      CALL NcCl( fId1c )
      CALL NcCl( fId1d )
      CALL NcCl( fId1e  )
      CALL NcCl( fId1f ) 
      CALL NcCl( fId1g ) 

      !--------------------------
      ! Compute future emissions
      !--------------------------
      IF ( LFUTURE ) THEN
         CALL NEI2011_SCALE_FUTURE
      ENDIF

      !--------------------------
      ! Print emission totals
      !--------------------------

      CALL TOTAL_ANTHRO_Tg( THISMONTH, THISDAY )

      ! Return to calling program
      END SUBROUTINE EMISS_NEI2011_ANTHRO

!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_nei2011_mask
!
! !DESCRIPTION: Subroutine READ\_NEI2011\_MASK reads the mask for NEI data  
!\\
!\\
! !INTERFACE:
      
      SUBROUTINE READ_NEI2011_MASK
!
! !USES:
!     
      ! Reference to F90 modules
      USE BPCH2_MOD,      ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE LOGICAL_MOD,    ONLY : LCAC,            LBRAVO
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1
      USE REGRID_A2A_MOD, ONLY : DO_REGRID_A2A
      USE TRANSFER_MOD,   ONLY : TRANSFER_2D

      USE CMN_SIZE_MOD         ! Size parameters

      USE m_netcdf_io_open     
      USE m_netcdf_io_read
      USE m_netcdf_io_readattr
      USE m_netcdf_io_close
      USE m_netcdf_io_get_dimlen
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
      REAL*4             :: ARRAY2(I1x1,J1x1)
      REAL*8, TARGET     :: GEOS_1x1(I1x1,J1x1)
      CHARACTER(LEN=255) :: FILENAME
      CHARACTER(LEN=255) :: LLFILENAME
      REAL*8, POINTER    :: INGRID(:,:) => NULL()
      INTEGER            :: st2d(2), ct2d(2)
      INTEGER            :: fId1
      !=================================================================
      ! Mask specific to NEI2011 data
      !=================================================================
      
      !SNAME = 'usa.'

      ! NEI2008 covers CANADA if we do not use CAC     
      !IF ( .NOT. LCAC ) SNAME = TRIM( SNAME ) // 'can.'

      ! NEI2008 covers Mexico if we do not use BRAVO      
      !IF ( .NOT. LBRAVO ) SNAME = TRIM( SNAME ) // 'mex.'

      
!%%%      FILENAME  = '/as/home/ktravis/' // &     
      !FILENAME  = TRIM( DATA_DIR_1x1) // &
           !'NEI2008_201307/usa.mask.nei2005.geos.1x1.nc' !kyu, 7Feb2014
      FILENAME = '/as/home/kyu/' // &
            'regrid/nei2005masktest.nc'

      ! Echo info
      WRITE( 6, 200 ) TRIM( FILENAME )
200   FORMAT( '     - READ_NEI2011_MASK: Reading ', a )
     
      ! Allocate start and count arrays
      st2d = (/1, 1/)
      ct2d = (/I1x1, J1x1/)
      ! Open and read model_ready data from netCDF file - wkday
      CALL Ncop_Rd(fId1, TRIM(FILENAME))
      Call NcRd(ARRAY2, fId1, 'MASK',   &
           st2d,  ct2d )        !Start andCount lat/lon
      ! Close netCDF file
      CALL NcCl( fId1 )
      ! Cast to REAL*8 before regridding
      GEOS_1x1(:,:) = ARRAY2(:,:)
     
      ! File with lat/lon edges for regridding
      LLFILENAME = TRIM( DATA_DIR_1x1) // &
           'MAP_A2A_Regrid_201203/MAP_A2A_latlon_geos1x1.nc'
      ! Regrid from GEOS 1x1 --> current model resolution [unitless]
      INGRID => GEOS_1x1(:,:)
      CALL DO_REGRID_A2A( LLFILENAME, I1x1,     J1x1, &
           INGRID,     USA_MASK, IS_MASS=0, &
           netCDF=.TRUE.                   )

      ! Free pointer
      NULLIFY( INGRID )

      WHERE ( USA_MASK > 0D0 ) USA_MASK = 1D0
      WRITE(*,*) 'READ NEI2011 MASK!'
      ! Return to calling program
      END SUBROUTINE READ_NEI2011_MASK
!------------------------------------------------------------------------------
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: nei2011_scale_future
!
! !DESCRIPTION: Subroutine NEI2011\_SCALE\_FUTURE applies the IPCC future 
!  scale factors to the NEI2011 anthropogenic emissions.
!\\
!\\
! !INTERFACE:

      SUBROUTINE NEI2011_SCALE_FUTURE
!
! !USES:
! 
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_COff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NH3an 
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NOxff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_SO2ff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_OCff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_BCff

      USE CMN_SIZE_MOD             ! Size parameters
!
! !REMARKS:
!    VOC are not scaled, however scale factors are available (see
!     epa_nei_mod.f for procedure)
!     
! !REVISION HISTORY: 
!    7 Oct 2009 - A. van Donkelaar - initial version
!   20 Oct 2009 - P. Le Sager - set L OpenMP private, put L loop first  
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                       :: I, J, L, HH

      !=================================================================
      ! NEI2011_SCALE_FUTURE begins here!
      !=================================================================

!$OMP PARALLEL DO &
!$OMP DEFAULT( SHARED ) &
!$OMP PRIVATE( I, J, L, HH )

      DO HH=1,24
      DO L = 1,3
      DO J = 1, JJPAR
      DO I = 1, IIPAR
          ! Future NO2 [molec/cm2/s]
          NO2(I,J,L,HH) = NO2(I,J,L,HH) * GET_FUTURE_SCALE_NOxff( I, J )

          ! Future CO  [molec/cm2/s]
          CO(I,J,L,HH) = CO(I,J,L,HH)  * GET_FUTURE_SCALE_COff(  I, J )

          ! Future SO2 [molec/cm2/s] 
          SO2a(I,J,L,HH) = SO2a(I,J,L,HH) * GET_FUTURE_SCALE_SO2ff( I, J )

          ! Future SO2 [molec/cm2/s] 
          SO2b(I,J,L,HH) = SO2b(I,J,L,HH) * GET_FUTURE_SCALE_SO2ff( I, J )

          ! Future SO4 [molec/cm2/s]
          SO4(I,J,L,HH)  = SO4(I,J,L,HH) * GET_FUTURE_SCALE_SO2ff( I, J )

          ! Future NH3 [molec/cm2/s] 
          NH3(I,J,L,HH)  = NH3(I,J,L,HH) * GET_FUTURE_SCALE_NH3an( I, J )

          ! Future OC [molec/cm2/s]
          OCPO(I,J,L,HH)  = OCPO(I,J,L,HH) * GET_FUTURE_SCALE_OCff( I, J )

          ! Future BC [molec/cm2/s]
          BCPO(I,J,L,HH) = BCPO(I,J,L,HH) * GET_FUTURE_SCALE_BCff( I, J )

      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE NEI2011_SCALE_FUTURE
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: total_anthro_Tg
!
! !DESCRIPTION: Subroutine TOTAL\_ANTHRO\_TG prints the totals for the 
!  anthropogenic emissions of NOx, CO, SO2 and NH3.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TOTAL_ANTHRO_TG( THISMONTH, THISDAY )
!
! !USES:
! 
      USE CMN_SIZE_MOD            ! Size parameters
      USE TIME_MOD,     ONLY : ITS_A_NEW_MONTH
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TRACER_MOD,   ONLY : XNUMOL
      USE TRACERID_MOD, ONLY : IDTCO, IDTNO,IDTNO2, IDTHNO2 
      USE TRACERID_MOD, ONLY : IDTSO2, IDTNH3, IDTSO4, IDTMACR
      USE TRACERID_MOD, ONLY : IDTALD2, IDTRCHO, IDTC2H6, IDTC3H8
      USE TRACERID_MOD, ONLY : IDTPRPE, IDTALK4, IDTACET
      USE TRACERID_MOD, ONLY : IDTBENZ, IDTTOLU, IDTXYLE, IDTMEK
      USE TRACERID_MOD, ONLY : IDTSO4, IDTCH2O, IDTOCPO,  IDTBCPO 
      USE TRACERID_MOD, ONLY : IDTMOH, IDTEOH

! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: THISMONTH   ! Month of data to compute totals
      INTEGER, INTENT(IN) :: THISDAY   ! Month of data to compute totals
!
! !REVISION HISTORY: 
!   7 Oct 2009 - A. van Donkelaar - initial version
!   13 Dec 2014 - K. Travis - revised for NEI2011
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE       :: FIRST = .TRUE.
      INTEGER             :: II, JJ, IH, LL
      REAL*8              :: T_CO, T_NO, T_NO2, T_HNO2, T_SO2, T_NH3
      REAL*8              :: T_ALD2,  T_RCHO, T_C2H6, T_C3H8, T_MEK
      REAL*8              :: T_PRPE, T_ALK4, T_TOLU, T_XYLE, T_ACET
      REAL*8              :: T_CH2O,T_BC, T_OC, T_SO4
      REAL*8              :: T_BENZ, T_EOH, T_MOH!, T_C2H4
      REAL*8              :: tmpArea(IIPAR, JJPAR,6)
  
      CHARACTER(LEN=3)    :: UNIT
      REAL*8,  PARAMETER  :: SEC_IN_HOUR  = 3600d0! * 365.25d0
                                 
      !=================================================================
      ! TOTAL_ANTHRO_TG begins here!
      !=================================================================

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, 100  )
 100  FORMAT( 'N. E. I. 2011  U. S. A.   E M I S S I O N S', / )
      
      IF ( ITS_A_NEW_MONTH() ) THEN
         T_CO_MON = 0d0 
         T_NO_MON = 0d0
         T_NO2_MON = 0d0
         T_HNO2_MON = 0d0
         T_SO2_MON = 0d0
         T_NH3_MON = 0d0
         T_ALD2_MON = 0d0
         T_RCHO_MON = 0d0
         T_C2H6_MON = 0d0
         T_C3H8_MON = 0d0
         T_MEK_MON = 0d0
         T_PRPE_MON = 0d0
         T_ALK4_MON = 0d0
         T_BENZ_MON = 0d0
         T_TOLU_MON = 0d0
         T_XYLE_MON = 0d0
         T_ACET_MON = 0d0
         T_CH2O_MON = 0d0
         T_BC_MON = 0d0
         T_OC_MON = 0d0
         T_SO4_MON = 0d0
         T_EOH_MON = 0d0
         T_MOH_MON = 0d0
      ENDIF

      tmpArea = 0d0
      DO II = 1, IIPAR
         DO JJ = 1, JJPAR
            DO LL=1, 6
               tmpArea(II,JJ,LL) = GET_AREA_CM2(II,JJ,LL)
            ENDDO
         ENDDO
      ENDDO
         
      ! Units are in [molec/cm2/s]
      ! Total CO  [Tg CO]
      IF ( IDTCO .NE. 0 ) &
           ! Total CO [Tg]
           T_CO  = SUM(SUM(CO,4)*tmpArea) * &
           (SEC_IN_HOUR *1.0d-9/XNUMOL(IDTCO))
           T_CO_MON = T_CO_MON + T_CO
   
      ! Total NO  [Tg NO]
      IF ( IDTNO .NE. 0 ) &
           ! Total NO [Tg N]
           T_NO =   SUM(SUM(NO, 4)*tmpArea ) *14d0/30d0 * &
           SEC_IN_HOUR *1.0d-9/XNUMOL(IDTNO)
           T_NO_MON = T_NO_MON + T_NO

      ! Total NO2  [Tg N]
      IF ( IDTNO2 .NE. 0 ) &
           ! Total NO [Tg N]
           T_NO2 =   SUM(SUM(NO2, 4)*tmpArea ) *14d0/46d0 * &
           SEC_IN_HOUR *1.0d-9/XNUMOL(IDTNO2)
           T_NO2_MON = T_NO2_MON + T_NO2

      ! Total HNO2  [Tg HNO2]
      IF ( IDTHNO2 .NE. 0 ) &
           ! Total HNO2 [Tg]
           T_HNO2 =   SUM(SUM(HNO2, 4)*tmpArea ) *14d0/47d0 * &
           SEC_IN_HOUR *1.0d-9/XNUMOL(IDTHNO2)
           T_HNO2_MON = T_HNO2_MON + T_HNO2

      ! Total SO2  [Tg SO2]
      IF ( IDTSO2 .NE. 0 ) &
           ! Total SO2 [Tg]
           T_SO2 =   SUM(SUM(SO2a, 4)*tmpArea ) *32d0/64d0 * &
           SEC_IN_HOUR *1.0d-9/XNUMOL(IDTSO2) + &
           SUM(SUM(SO2b, 4)*tmpArea ) *32d0/64d0*  &
           SEC_IN_HOUR *1.0d-9/XNUMOL(IDTSO2)
           T_SO2_MON = T_SO2_MON + T_SO2

      ! Total SO4  [Tg SO4]
      IF ( IDTSO4 .NE. 0 ) &
           ! Total SO4 [Tg]
           T_SO4 =   SUM(SUM(SO4, 4)*tmpArea ) * 32d0/96d0 * &
           SEC_IN_HOUR *1.0d-12
           T_SO4_MON = T_SO4_MON + T_SO4

      ! Total NH3  [Tg NH3]
      IF ( IDTNH3 .NE. 0 ) &
           ! Total NH3 [Tg]
           T_NH3 =   SUM(SUM(NH3, 4)*tmpArea ) * &
           SEC_IN_HOUR *1.0d-9/XNUMOL(IDTNH3) 
           T_NH3_MON = T_NH3_MON + T_NH3

      ! Total OC  [Tg OC]
      IF ( IDTOCPO .NE. 0 ) &
           ! Total OC [Tg] 
           T_OC =   SUM(SUM(OCPO, 4)*tmpArea ) * &
           SEC_IN_HOUR *1.0d-12
           T_OC_MON = T_OC_MON + T_OC

      ! Total BC  [Tg BC]
      IF ( IDTBCPO .NE. 0 ) &
           ! Total BC [Tg]
           T_BC =   SUM(SUM(BCPO, 4)*tmpArea ) * &
           SEC_IN_HOUR *1.0d-12
           T_BC_MON = T_BC_MON + T_BC

      ! Total MEK  [Tg MEK]
      IF ( IDTMEK .NE. 0 ) &
           ! Total MEK [Tg]
           T_MEK =   SUM(SUM(MEK, 4)*tmpArea ) * &
           SEC_IN_HOUR *1.0d-9/XNUMOL(IDTMEK)
           T_MEK_MON = T_MEK_MON + T_MEK

      ! Total RCHO  [Tg C]
      IF ( IDTRCHO .NE. 0 ) &
           ! Total RCHO [Tg]
           T_RCHO =   SUM(SUM(RCHO, 4)*tmpArea ) * &
           SEC_IN_HOUR *1.0d-9/XNUMOL(IDTRCHO)
           T_RCHO_MON = T_RCHO_MON + T_RCHO

      ! Total CH2O  [Tg]
      IF ( IDTCH2O .NE. 0 ) &
           ! Total CH2O [Tg]
           T_CH2O =   SUM(SUM(CH2O, 4)*tmpArea ) * &
           SEC_IN_HOUR *1.0d-9/XNUMOL(IDTCH2O)
           T_CH2O_MON = T_CH2O_MON + T_CH2O

      ! Total PRPE  [Tg C]
      IF ( IDTPRPE .NE. 0 ) &
           ! Total PRPE [Tg C ]
           T_PRPE =   SUM(SUM(PRPEa, 4)*tmpArea ) * &
           SEC_IN_HOUR *1.0d-9/XNUMOL(IDTPRPE) + &
           SUM(SUM(PRPEb, 4)*tmpArea ) * &
           SEC_IN_HOUR *1.0d-9/XNUMOL(IDTPRPE)
           T_PRPE_MON = T_PRPE_MON + T_PRPE

      ! Total ACET  [Tg C]
      IF ( IDTACET .NE. 0 ) &
           ! Total ACET [Tg C ]
           T_ACET =   SUM(SUM(ACET, 4)*tmpArea ) * &
           SEC_IN_HOUR *1.0d-9/XNUMOL(IDTACET)
           T_ACET_MON = T_ACET_MON + T_ACET

      ! Total ALK4 [Tg C]
      IF ( IDTALK4 .NE. 0 ) &
           ! Total ALK4 [Tg C ]
           T_ALK4 =   SUM(SUM(ALK4, 4)*tmpArea ) * &
           SEC_IN_HOUR *1.0d-9/XNUMOL(IDTALK4)
           T_ALK4_MON = T_ALK4_MON + T_ALK4

     ! Total ALD2 [Tg C]
      IF ( IDTALD2 .NE. 0 ) &
           ! Total ALD2 [Tg C ]
           T_ALD2 =   SUM(SUM(ALD2, 4)*tmpArea ) * &
           SEC_IN_HOUR *1.0d-9/XNUMOL(IDTALD2)
           T_ALD2_MON = T_ALD2_MON + T_ALD2

      ! Total C3H8  [Tg]
      IF ( IDTC3H8 .NE. 0 ) &
           ! Total C3H8 [Tg]
           T_C3H8 =   SUM(SUM(C3H8, 4)*tmpArea ) * &
           SEC_IN_HOUR *1.0d-9/XNUMOL(IDTC3H8)
           T_C3H8_MON = T_C3H8_MON + T_C3H8

      ! Total C2H6  [Tg]
      IF ( IDTC2H6 .NE. 0 ) &
           ! Total C2H6 [Tg]
           T_C2H6 =   SUM(SUM(C2H6, 4)*tmpArea ) * &
           SEC_IN_HOUR *1.0d-9/XNUMOL(IDTC2H6)
           T_C2H6_MON = T_C2H6_MON + T_C2H6
      
      ! Total EOH  [Tg C]
      IF ( IDTEOH .NE. 0 ) &
           ! Total EOH [Tg]
           T_EOH =   SUM(SUM(EOH, 4)*tmpArea ) * &
           SEC_IN_HOUR *1.0d-9/XNUMOL(IDTALD2)
           T_EOH_MON = T_EOH_MON + T_EOH

      ! Total MOH  [Tg]
      IF ( IDTMOH .NE. 0 ) &
           ! Total MOH [Tg]
           T_MOH =   SUM(SUM(MOH, 4)*tmpArea ) * &
           SEC_IN_HOUR *1.0d-9/(XNUMOL(IDTALD2)*0.5)
           T_MOH_MON = T_MOH_MON + T_MOH

      ! Total BENZ  [Tg]
!      IF ( IDTBENZ .NE. 0 ) &
           ! Total BENZ [ 
           T_BENZ =   SUM(SUM(BENZ, 4)*tmpArea ) * &
           SEC_IN_HOUR *1.0d-9/(XNUMOL(IDTALD2)*3)
           T_BENZ_MON = T_BENZ_MON + T_BENZ

           T_TOLU =   SUM(SUM(TOLU, 4)*tmpArea ) * &
           SEC_IN_HOUR *1.0d-9/(XNUMOL(IDTALD2)*3.5)
           T_TOLU_MON = T_TOLU_MON + T_TOLU
           
           T_XYLE =   SUM(SUM(XYLE, 4)*tmpArea ) * &
           SEC_IN_HOUR *1.0d-9/(XNUMOL(IDTALD2)*4)
           T_XYLE_MON = T_XYLE_MON + T_XYLE

      ! Format statement
      WRITE(*,*) 'NEI2011 anthro for month', &
             THISMONTH, ' and day ', THISDAY

     ! Print totals in [Tg]
      WRITE( 6, 110 ) 'CO   ', THISMONTH, T_CO_MON, THISDAY, T_CO, '[ Tg CO ]'
      WRITE( 6, 110 ) 'NO   ', THISMONTH, T_NO_MON, THISDAY, T_NO, '[ Tg N ]'
      WRITE( 6, 110 ) 'NO2  ', THISMONTH,  T_NO2_MON,THISDAY, T_NO2, '[ Tg N ]'
      WRITE( 6, 110 ) 'HNO2  ', THISMONTH, T_HNO2_MON, THISDAY, T_HNO2, '[ Tg N ]'
      WRITE( 6, 110 ) 'CH2O  ', THISMONTH,T_CH2O_MON, THISDAY, T_CH2O, '[ Tg ]'
      WRITE( 6, 110 ) 'MEK  ', THISMONTH,T_MEK_MON, THISDAY, T_MEK, '[ Tg C ]'
      WRITE( 6, 110 ) 'RCHO  ', THISMONTH,T_RCHO_MON, THISDAY, T_RCHO, '[ Tg C ]'
      WRITE( 6, 110 ) 'PRPE  ', THISMONTH, T_PRPE_MON, THISDAY, T_PRPE, '[ Tg C ]'
      WRITE( 6, 110 ) 'ACET  ', THISMONTH,T_ACET_MON, THISDAY, T_ACET, '[ Tg C ]'
      WRITE( 6, 110 ) 'ALK4  ', THISMONTH,T_ALK4_MON, THISDAY, T_ALK4, '[ Tg C ]'
      WRITE( 6, 110 ) 'ALD2  ', THISMONTH,T_ALD2_MON, THISDAY, T_ALD2, '[ Tg C ]'
      WRITE( 6, 110 ) 'C3H8  ', THISMONTH, T_C3H8_MON, THISDAY, T_C3H8, '[ Tg ]'
      WRITE( 6, 110 ) 'C2H6  ', THISMONTH,T_C2H6_MON, THISDAY, T_C2H6, '[ Tg ]'
      WRITE( 6, 110 ) 'SO2  ', THISMONTH,T_SO2_MON, THISDAY, T_SO2, '[ Tg S ]'
      WRITE( 6, 110 ) 'SO4  ', THISMONTH, T_SO4_MON, THISDAY, T_SO4, '[ Tg S ]'
      WRITE( 6, 110 ) 'NH3  ', THISMONTH, T_NH3_MON,THISDAY, T_NH3, '[ Tg ]'
      WRITE( 6, 110 ) 'OC  ', THISMONTH,T_OC_MON, THISDAY, T_OC, '[ Tg ]'
      WRITE( 6, 110 ) 'BC  ', THISMONTH,T_BC_MON, THISDAY, T_BC, '[ Tg ]'
      WRITE( 6, 110 ) 'EOH  ', THISMONTH,T_EOH_MON, THISDAY, T_EOH, '[ Tg C ]'
      WRITE( 6, 110 ) 'MOH  ', THISMONTH,T_MOH_MON, THISDAY, T_MOH, '[ Tg C ]'
      WRITE( 6, 110 ) 'BENZ  ', THISMONTH,T_BENZ_MON, THISDAY, T_BENZ, '[ Tg C ]'
      WRITE( 6, 110 ) 'TOLU  ', THISMONTH,T_TOLU_MON, THISDAY, T_TOLU, '[ Tg C ]'
      WRITE( 6, 110 ) 'XYLE  ', THISMONTH,T_XYLE_MON, THISDAY, T_XYLE, '[ Tg C ]'

      ! Format statement
 110  FORMAT( 'NEI2011 anthro ', a5, &
             'for cumulative month ', i2,': ' f11.5,', and just day ', i2, ': ', f11.5, 1x, a9 )
      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      
      ! Return to calling program
      END SUBROUTINE TOTAL_ANTHRO_Tg
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
      SUBROUTINE INIT_NEI2011_ANTHRO( am_I_Root, Input_Opt, RC )
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
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
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
      ! INIT_NEI2011_ANTHRO begins here!
      !=================================================================
      ! Assume success
      RC        =  GIGC_SUCCESS
      
      ! Return if LNEI11 is false
      IF ( .not. Input_Opt%LNEI11 ) RETURN
      
      !--------------------------------------------------
      ! Allocate and zero arrays for emissions
      !--------------------------------------------------
      IF ( .not. Input_Opt%LBIOMASS ) THEN
      ! allocate and read USA Mask
         ALLOCATE( USA_MASK( IIPAR, JJPAR ), STAT=RC )
         IF ( RC /= 0 ) CALL ALLOC_ERR( 'USA_MASK' )
         USA_MASK = 0d0
         
         CALL READ_NEI2011_MASK
      ENDIF

      ALLOCATE( TMP( IIPAR, JJPAR, 6, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'TMP' )
      TMP = 0d0

      ALLOCATE( TMPARR( IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'TMPARR' )
      TMPARR = 0d0

      IF ( ITS_A_NEW_MONTH() ) THEN
         ALLOCATE( T_CO_MON( 1 ), STAT=RC )
         IF ( RC /= 0 ) CALL ALLOC_ERR( 'T_CO_MON' )
         T_CO_MON = 0d0

         ALLOCATE( T_NO_MON( 1 ), STAT=RC )
         IF ( RC /= 0 ) CALL ALLOC_ERR( 'T_NO_MON' )
         T_NO_MON = 0d0

         ALLOCATE( T_NO2_MON( 1 ), STAT=RC )
         IF ( RC /= 0 ) CALL ALLOC_ERR( 'T_NO2_MON' )
         T_NO2_MON = 0d0

         ALLOCATE( T_HNO2_MON( 1 ), STAT=RC )
         IF ( RC /= 0 ) CALL ALLOC_ERR( 'T_HNO2_MON' )
         T_HNO2_MON = 0d0

         ALLOCATE( T_NH3_MON( 1 ), STAT=RC )
         IF ( RC /= 0 ) CALL ALLOC_ERR( 'T_NH3_MON' )
         T_NH3_MON = 0d0
         
         ALLOCATE( T_SO2_MON( 1 ), STAT=RC )
         IF ( RC /= 0 ) CALL ALLOC_ERR( 'T_SO2_MON' )
         T_SO2_MON = 0d0
         
         ALLOCATE( T_SO4_MON( 1 ), STAT=RC )
         IF ( RC /= 0 ) CALL ALLOC_ERR( 'T_SO4_MON' )
         T_SO4_MON = 0d0
         
         ALLOCATE( T_OC_MON( 1 ), STAT=RC )
         IF ( RC /= 0 ) CALL ALLOC_ERR( 'T_OC_MON' )
         T_OC_MON = 0d0
         
         ALLOCATE( T_BC_MON( 1 ), STAT=RC )
         IF ( RC /= 0 ) CALL ALLOC_ERR( 'T_BC_MON' )
         T_BC_MON = 0d0
         
         ALLOCATE( T_EOH_MON( 1 ), STAT=RC )
         IF ( RC /= 0 ) CALL ALLOC_ERR( 'T_EOH_MON' )
         T_EOH_MON = 0d0
         
         ALLOCATE( T_MOH_MON( 1 ), STAT=RC )
         IF ( RC /= 0 ) CALL ALLOC_ERR( 'T_MOH_MON' )
         T_MOH_MON = 0d0
         
         ALLOCATE( T_C3H8_MON( 1 ), STAT=RC )
         IF ( RC /= 0 ) CALL ALLOC_ERR( 'T_C3H8_MON' )
         T_C3H8_MON = 0d0
         
         ALLOCATE( T_C2H6_MON( 1 ), STAT=RC )
         IF ( RC /= 0 ) CALL ALLOC_ERR( 'T_C2H6_MON' )
         T_C2H6_MON = 0d0
         
         ALLOCATE( T_RCHO_MON( 1 ), STAT=RC )
         IF ( RC /= 0 ) CALL ALLOC_ERR( 'T_RCHO_MON' )
         T_RCHO_MON = 0d0
         
         ALLOCATE( T_ALD2_MON( 1 ), STAT=RC )
         IF ( RC /= 0 ) CALL ALLOC_ERR( 'T_ALD2_MON' )
         T_ALD2_MON = 0d0
         
         ALLOCATE( T_CH2O_MON( 1 ), STAT=RC )
         IF ( RC /= 0 ) CALL ALLOC_ERR( 'T_CH2O_MON' )
         T_CH2O_MON = 0d0
         
         ALLOCATE( T_ALK4_MON( 1 ), STAT=RC )
         IF ( RC /= 0 ) CALL ALLOC_ERR( 'T_ALK4_MON' )
         T_ALK4_MON = 0d0
         
         ALLOCATE( T_ACET_MON( 1 ), STAT=RC )
         IF ( RC /= 0 ) CALL ALLOC_ERR( 'T_ACET_MON' )
         T_ACET_MON = 0d0
         
         ALLOCATE( T_MACR_MON( 1 ), STAT=RC )
         IF ( RC /= 0 ) CALL ALLOC_ERR( 'T_MACR_MON' )
         T_MACR_MON = 0d0

         ALLOCATE( T_MEK_MON( 1 ), STAT=RC )
         IF ( RC /= 0 ) CALL ALLOC_ERR( 'T_MEK_MON' )
         T_MEK_MON = 0d0
         
         ALLOCATE( T_PRPE_MON( 1 ), STAT=RC )
         IF ( RC /= 0 ) CALL ALLOC_ERR( 'T_PRPE_MON' )
         T_PRPE_MON = 0d0
         
         ALLOCATE( T_TOLU_MON( 1 ), STAT=RC )
         IF ( RC /= 0 ) CALL ALLOC_ERR( 'T_TOLU_MON' )
         T_TOLU_MON = 0d0
         
         ALLOCATE( T_XYLE_MON( 1 ), STAT=RC )
         IF ( RC /= 0 ) CALL ALLOC_ERR( 'T_XYLE_MON' )
         T_XYLE_MON = 0d0

         ALLOCATE( T_BENZ_MON( 1 ), STAT=RC )
         IF ( RC /= 0 ) CALL ALLOC_ERR( 'T_BENZ_MON' )
         T_BENZ_MON = 0d0
      
      ENDIF

      ALLOCATE( CO( IIPAR, JJPAR, 6, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'CO' )
      CO = 0d0

      ALLOCATE( NO( IIPAR, JJPAR, 6, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'NO' )
      NO = 0d0

      ALLOCATE( NO2( IIPAR, JJPAR, 6, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'NO2' )
      NO2 = 0d0
      
      ALLOCATE( HNO2( IIPAR, JJPAR, 6, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'HNO2' )
      HNO2 = 0d0

      ALLOCATE( SO2a( IIPAR, JJPAR, 6, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'SO2a' )
      SO2a = 0d0

      ALLOCATE( SO2b( IIPAR, JJPAR, 6, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'SO2b' )
      SO2b = 0d0

      ALLOCATE( NH3( IIPAR, JJPAR, 6, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'NH3' )
      NH3 = 0d0 

      ALLOCATE( ALD2( IIPAR, JJPAR, 6, 24), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'ALD2' )
      ALD2 = 0d0 

      ALLOCATE( RCHO( IIPAR, JJPAR, 6, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'RCHO' )
      RCHO = 0d0

      ALLOCATE( BENZ( IIPAR, JJPAR, 6, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'BENZ' )
      BENZ = 0d0

      ALLOCATE( C2H6( IIPAR, JJPAR, 6, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'C2H6' )
      C2H6 = 0d0

      ALLOCATE( PRPEa( IIPAR, JJPAR, 6, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'PRPEa' )
      PRPEa = 0d0

      ALLOCATE( PRPEb( IIPAR, JJPAR, 6, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'PRPEb' )
      PRPEb = 0d0

      ALLOCATE( ALK4( IIPAR, JJPAR, 6, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'ALK4' )
      ALK4 = 0d0 

      ALLOCATE( TOLU( IIPAR, JJPAR, 6, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'TOLU' )
      TOLU = 0d0 

      ALLOCATE( XYLE( IIPAR, JJPAR, 6, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'XYLE' )
      XYLE = 0d0 
     
      !ALLOCATE( C2H4( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      !IF ( RC /= 0 ) CALL ALLOC_ERR( 'C2H4' )
      !C2H4 = 0d0 

      ALLOCATE( CH2O( IIPAR, JJPAR, 6, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'CH2O' )
      CH2O = 0d0

      ALLOCATE( C3H8( IIPAR, JJPAR, 6, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'C3H8' )
      C3H8 = 0d0
      
      ALLOCATE( BCPO( IIPAR, JJPAR, 6, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'BCPO' )
      BCPO = 0d0
      
      ALLOCATE( OCPO( IIPAR, JJPAR, 6, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'OCPO' )
      OCPO = 0d0

      ALLOCATE( SO4( IIPAR, JJPAR, 6, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'SO4' )
      SO4 = 0d0

      ALLOCATE( MACR( IIPAR, JJPAR, 6, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'MACR' )
      MACR = 0d0
      
      ALLOCATE( MEK( IIPAR, JJPAR, 6, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'MEK' )
      MEK = 0d0

      ALLOCATE( ACET( IIPAR, JJPAR, 6, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'ACET' )
      ACET = 0d0
      
      ALLOCATE( EOH( IIPAR, JJPAR, 6, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'EOH' )
      EOH = 0d0

      ALLOCATE( MOH( IIPAR, JJPAR, 6, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'MOH' )
      MOH = 0d0

      !ALLOCATE( CH4( IIPAR, JJPAR, 6, 24 ), STAT=RC )
      !IF ( RC /= 0 ) CALL ALLOC_ERR( 'CH4' )
      !CH4 = 0d0
      
      END SUBROUTINE INIT_NEI2011_ANTHRO
!EOC
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
      SUBROUTINE CLEANUP_NEI2011_ANTHRO
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
      IF ( ALLOCATED( USA_MASK) ) DEALLOCATE( USA_MASK )
      IF ( ALLOCATED( CO     ) ) DEALLOCATE( CO    )
      IF ( ALLOCATED( T_CO_MON ) ) DEALLOCATE( T_CO_MON    )
      IF ( ALLOCATED( T_NO_MON ) ) DEALLOCATE( T_NO_MON    )
      IF ( ALLOCATED( T_NO2_MON ) ) DEALLOCATE( T_NO2_MON   )
      IF ( ALLOCATED( T_HNO2_MON ) ) DEALLOCATE( T_HNO2_MON )
      IF ( ALLOCATED( T_SO2_MON ) ) DEALLOCATE( T_SO2_MON    )
      IF ( ALLOCATED( T_SO4_MON ) ) DEALLOCATE( T_SO4_MON    )
      IF ( ALLOCATED( T_OC_MON ) ) DEALLOCATE( T_OC_MON    )
      IF ( ALLOCATED( T_BC_MON ) ) DEALLOCATE( T_BC_MON    )
      IF ( ALLOCATED( T_NH3_MON ) ) DEALLOCATE( T_NH3_MON  )
      IF ( ALLOCATED( T_C3H8_MON ) ) DEALLOCATE( T_C3H8_MON  )
      IF ( ALLOCATED( T_C2H6_MON ) ) DEALLOCATE( T_C2H6_MON  )
      IF ( ALLOCATED( T_ALD2_MON ) ) DEALLOCATE( T_ALD2_MON  )
      IF ( ALLOCATED( T_RCHO_MON ) ) DEALLOCATE( T_RCHO_MON  )
      IF ( ALLOCATED( T_ALK4_MON ) ) DEALLOCATE( T_ALK4_MON  )
      IF ( ALLOCATED( T_ACET_MON ) ) DEALLOCATE( T_ACET_MON  )
      IF ( ALLOCATED( T_MACR_MON ) ) DEALLOCATE( T_MACR_MON  )
      IF ( ALLOCATED( T_MEK_MON ) ) DEALLOCATE( T_MEK_MON  )
      IF ( ALLOCATED( T_PRPE_MON ) ) DEALLOCATE( T_PRPE_MON  )
      IF ( ALLOCATED( T_TOLU_MON ) ) DEALLOCATE( T_TOLU_MON  )
      IF ( ALLOCATED( T_XYLE_MON ) ) DEALLOCATE( T_XYLE_MON  )
      IF ( ALLOCATED( T_BENZ_MON ) ) DEALLOCATE( T_BENZ_MON  )
      IF ( ALLOCATED( T_EOH_MON ) ) DEALLOCATE( T_EOH_MON  )
      IF ( ALLOCATED( T_MOH_MON ) ) DEALLOCATE( T_MOH_MON  )

      IF ( ALLOCATED( NO     ) ) DEALLOCATE( NO    )
      IF ( ALLOCATED( NO2    ) ) DEALLOCATE( NO2   )
      IF ( ALLOCATED( HNO2   ) ) DEALLOCATE( HNO2  )
      IF ( ALLOCATED( SO2a    ) ) DEALLOCATE( SO2a   )
      IF ( ALLOCATED( SO2b    ) ) DEALLOCATE( SO2b   )
      IF ( ALLOCATED( NH3    ) ) DEALLOCATE( NH3   )
      IF ( ALLOCATED( ALD2   ) ) DEALLOCATE( ALD2  )
      IF ( ALLOCATED( RCHO   ) ) DEALLOCATE( RCHO  )
      IF ( ALLOCATED( BENZ   ) ) DEALLOCATE( BENZ  )
      IF ( ALLOCATED( C2H6   ) ) DEALLOCATE( C2H6  )
      IF ( ALLOCATED( PRPEa   ) ) DEALLOCATE( PRPEa  )
      IF ( ALLOCATED( PRPEb   ) ) DEALLOCATE( PRPEb  )
      IF ( ALLOCATED( ALK4   ) ) DEALLOCATE( ALK4  )
      IF ( ALLOCATED( TOLU   ) ) DEALLOCATE( TOLU  )
      IF ( ALLOCATED( XYLE   ) ) DEALLOCATE( XYLE  )
      !IF ( ALLOCATED( C2H4   ) ) DEALLOCATE( C2H4  )
      IF ( ALLOCATED( CH2O   ) ) DEALLOCATE( CH2O  )
      IF ( ALLOCATED( BCPO   ) ) DEALLOCATE( BCPO  )
      IF ( ALLOCATED( OCPO   ) ) DEALLOCATE( OCPO  )
      IF ( ALLOCATED( SO4    ) ) DEALLOCATE( SO4   )
      IF ( ALLOCATED( EOH    ) ) DEALLOCATE( EOH   )
      IF ( ALLOCATED( MOH    ) ) DEALLOCATE( MOH   )
      !IF ( ALLOCATED( CH4    ) ) DEALLOCATE( CH4   )
      
      END SUBROUTINE CLEANUP_NEI2011_ANTHRO
!EOC
      END MODULE NEI2011_ANTHRO_MOD

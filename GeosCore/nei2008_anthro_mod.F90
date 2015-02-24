!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: nei2008_anthro_mod
!
! !DESCRIPTION: Module NEI2008\_ANTHRO\_MOD contains variables and routines to 
!  read the NEI2008 anthropogenic emissions.   
!\\
!\\
! !INTERFACE: 
!
      MODULE NEI2008_ANTHRO_MOD
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
      PUBLIC  :: CLEANUP_NEI2008_ANTHRO
      PUBLIC  :: EMISS_NEI2008_ANTHRO
      PUBLIC  :: GET_NEI2008_ANTHRO
      PUBLIC  :: GET_MASK_FORFIRE
      !--------------------------------------
      ! Leave for future use (bmy, 12/3/09)
      !PUBLIC  :: GET_NEI2005_MASK
      !--------------------------------------
!
! !PRIVATE MEMBER FUNCTIONS:
! No longer need scaling except for future emissions
      PRIVATE :: NEI2008_SCALE_FUTURE
      PRIVATE :: INIT_NEI2008_ANTHRO
      PRIVATE :: TOTAL_ANTHRO_TG
      PRIVATE :: READ_NEI2008_MASK
!
! !REMARKS:
!  Note that NEI2008 does not have MEK, ACET, or C3H8
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
      REAL*8,  ALLOCATABLE, TARGET :: TMP_WD(:,:,:,:)
      REAL*8,  ALLOCATABLE, TARGET :: TMP_WE(:,:,:,:)
      REAL*8,  ALLOCATABLE, TARGET :: TMPARR_WD(:,:)
      REAL*8,  ALLOCATABLE, TARGET :: TMPARR_WE(:,:)

      REAL*8,  ALLOCATABLE         :: CO(:,:,:,:) 
      REAL*8,  ALLOCATABLE         :: NO(:,:,:,:) 
      REAL*8,  ALLOCATABLE         :: NO2(:,:,:,:)
      REAL*8,  ALLOCATABLE         :: HNO2(:,:,:,:)
      REAL*8,  ALLOCATABLE         :: SO2(:,:,:,:) 
      REAL*8,  ALLOCATABLE         :: NH3(:,:,:,:) 
      REAL*8,  ALLOCATABLE         :: ALD2(:,:,:,:)
      REAL*8,  ALLOCATABLE         :: RCHO(:,:,:,:)
      REAL*8,  ALLOCATABLE         :: BENZ(:,:,:,:)
      !REAL*8,  ALLOCATABLE         :: CH4(:,:,:,:)
      REAL*8,  ALLOCATABLE         :: C2H6(:,:,:,:) 
      REAL*8,  ALLOCATABLE         :: PRPE(:,:,:,:)
      REAL*8,  ALLOCATABLE         :: ALK4(:,:,:,:)
      REAL*8,  ALLOCATABLE         :: C3H8(:,:,:,:)
      REAL*8,  ALLOCATABLE         :: MEK(:,:,:,:)
      REAL*8,  ALLOCATABLE         :: ACET(:,:,:,:)
      REAL*8,  ALLOCATABLE         :: TOLU(:,:,:,:)
      REAL*8,  ALLOCATABLE         :: XYLE(:,:,:,:)
      !REAL*8,  ALLOCATABLE         :: C2H4(:,:,:,:)
      !REAL*8,  ALLOCATABLE         :: MOH(:,:,:,:) 
      !REAL*8,  ALLOCATABLE         :: EOH(:,:,:,:) 
      REAL*8,  ALLOCATABLE         :: CH2O(:,:,:,:) 
      REAL*8,  ALLOCATABLE         :: OCPO(:,:,:,:) 
      REAL*8,  ALLOCATABLE         :: BCPO(:,:,:,:)
      REAL*8,  ALLOCATABLE         :: SO4(:,:,:,:) 

      ! No longer need NOx family variables (skim, 6/26/13)
      !REAL*8,  ALLOCATABLE         :: NOX_WKEND(:,:,:,:) 
      REAL*8,  ALLOCATABLE         :: CO_WKEND(:,:,:,:) 
      REAL*8,  ALLOCATABLE         :: NO_WKEND(:,:,:,:) 
      REAL*8,  ALLOCATABLE         :: NO2_WKEND(:,:,:,:) 
      REAL*8,  ALLOCATABLE         :: HNO2_WKEND(:,:,:,:)
      REAL*8,  ALLOCATABLE         :: SO2_WKEND(:,:,:,:) 
      REAL*8,  ALLOCATABLE         :: NH3_WKEND(:,:,:,:) 
      REAL*8,  ALLOCATABLE         :: ALD2_WKEND(:,:,:,:) 
      REAL*8,  ALLOCATABLE         :: RCHO_WKEND(:,:,:,:)
      REAL*8,  ALLOCATABLE         :: BENZ_WKEND(:,:,:,:) 
      !REAL*8,  ALLOCATABLE         :: CH4_WKEND(:,:,:,:)
      REAL*8,  ALLOCATABLE         :: C2H6_WKEND(:,:,:,:)  
      REAL*8,  ALLOCATABLE         :: PRPE_WKEND(:,:,:,:) 
      REAL*8,  ALLOCATABLE         :: ALK4_WKEND(:,:,:,:) 
      REAL*8,  ALLOCATABLE         :: C3H8_WKEND(:,:,:,:) 
      REAL*8,  ALLOCATABLE         :: MEK_WKEND(:,:,:,:) 
      REAL*8,  ALLOCATABLE         :: ACET_WKEND(:,:,:,:) 
      REAL*8,  ALLOCATABLE         :: TOLU_WKEND(:,:,:,:) 
      REAL*8,  ALLOCATABLE         :: XYLE_WKEND(:,:,:,:) 
      !REAL*8,  ALLOCATABLE         :: C2H4_WKEND(:,:,:,:)
      !REAL*8,  ALLOCATABLE         :: MOH_WKEND(:,:,:,:) 
      !REAL*8,  ALLOCATABLE         :: EOH_WKEND(:,:,:,:) 
      REAL*8,  ALLOCATABLE         :: CH2O_WKEND(:,:,:,:)
      REAL*8,  ALLOCATABLE         :: OCPO_WKEND(:,:,:,:)
      REAL*8,  ALLOCATABLE         :: BCPO_WKEND(:,:,:,:)
      REAL*8,  ALLOCATABLE         :: SO4_WKEND(:,:,:,:)

      ! Shadow logical variables from Input_Opt
      LOGICAL                      :: LBRAVO
      LOGICAL                      :: LCAC
      LOGICAL                      :: LFUTURE
      LOGICAL                      :: LNEI08
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
! !IROUTINE: get_nei2008_anthro
!
! !DESCRIPTION: Function GET\_NEI2008\_ANTHRO returns the NEI2008
!  emission for GEOS-Chem grid box (I,J) and tracer N and hour IH.  
!  Emissions can be returned in units of [kg/s] or [molec/cm2/s].
!\\ (krt, 2/10/13), now need IH
!\\
! !INTERFACE:
!
      FUNCTION GET_NEI2008_ANTHRO( I, J, L, IH, N, WEEKDAY ) RESULT( VALUE )
!
! !USES:
!
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TRACER_MOD,   ONLY : XNUMOL
      USE TRACERID_MOD, ONLY : IDTCO, IDTNO,IDTNO2, IDTHNO2 
      USE TRACERID_MOD, ONLY : IDTSO2, IDTNH3
      USE TRACERID_MOD, ONLY : IDTALD2, IDTRCHO, IDTC2H6
      USE TRACERID_MOD, ONLY : IDTPRPE, IDTALK4!, IDTC2H4
      USE TRACERID_MOD, ONLY : IDTBENZ, IDTTOLU, IDTXYLE
      USE TRACERID_MOD, ONLY : IDTSO4, IDTCH2O
      USE TRACERID_MOD, ONLY : IDTOCPO,  IDTBCPO
      !USE TRACERID_MOD, ONLY : IDTMOH, IDTEOH!, IDTCH4
      USE TRACERID_MOD, ONLY : IDTMEK, IDTC3H8, IDTACET !added 9/24/14, krt 
!
! !INPUT PARAMETERS: 
!
      ! Longitude, latitude, hour, and tracer indices
      INTEGER, INTENT(IN)           :: I, J, L, IH, N
      ! OPTIONAL -- WEEKEND/WEEKDAY
      LOGICAL, INTENT(IN), OPTIONAL :: WEEKDAY
 
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
      ! GET_NEI2008_ANTHRO begins here!
      !=================================================================

      ! Initialize

      IF ( WEEKDAY ) THEN ! First IF
         IF ( N == IDTCO ) THEN ! Second IF

            ! CO [molec/cm2/s]
            VALUE = CO(I,J,L,IH)
         ELSE IF ( N == IDTNO ) THEN

            ! NOX[molec/cm2/s]
          !  NOX(I,J,L,IH) = NO(I,J,L,IH)!+NO2(I,J,L,IH)+HNO2(I,J,L,IH)
            VALUE = NO(I,J,L,IH)
            
         ELSE IF ( N == IDTNO2 ) THEN

            ! NO2[molec/cm2/s]
            VALUE = NO2(I,J,L,IH)

         ELSE IF ( N == IDTHNO2 ) THEN

            ! HNO2[molec/cm2/s]
            VALUE = HNO2(I,J,L,IH)

         ELSE IF ( N == IDTSO2 ) THEN

            ! SO2 [molec/cm2/s]
            VALUE = SO2(I,J,L,IH)

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
            VALUE = PRPE(I,J,L,IH)

            ! IN NEI08, ALK4, C3H8, ACET and MEK 
            ! are lumped into ALK4 - split here into
            ! components based on NEI05
         ELSE IF ( N == IDTALK4 ) THEN

            ! [molec/cm2/s]
            VALUE = ALK4(I,J,L,IH)!*0.67
         ELSE IF ( N == IDTC3H8 ) THEN

            ! [molec/cm2/s]
            VALUE = C3H8(I,J,L,IH)!*0.05*5

         ELSE IF ( N == IDTACET ) THEN

            ! [molec/cm2/s]
            VALUE = ACET(I,J,L,IH)

         ELSE IF ( N == IDTMEK ) THEN

            ! [molec/cm2/s]
            VALUE = MEK(I,J,L,IH)

        ELSE IF ( N == IDTTOLU ) THEN

            ! [molec/cm2/s]
            VALUE = TOLU(I,J,L,IH)
         
         ELSE IF ( N == IDTXYLE ) THEN

            ! [molec/cm2/s]
            VALUE = XYLE(I,J,L,IH)
  
         !ELSE IF ( N == IDTC2H4 ) THEN

            ! [molec/cm2/s]
          !  VALUE = C2H4(I,J,L,IH)

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
         
        ! ELSE IF ( N == IDTEOH ) THEN

            ! [molec/cm2/s]
         !   VALUE = EOH(I,J,L,IH)
         
         !ELSE IF ( N == IDTMOH ) THEN

            ! [molec/cm2/s]
          !  VALUE = MOH(I,J,L,IH)
         
         !ELSE IF ( N == IDTCH4 ) THEN

            ! [molec/cm2/s]
          !  VALUE = CH4(I,J,L,IH)
         ELSE
            ! Otherwise return a negative value to indicate
            ! that there are no NEI2008 emissions for tracer N
            VALUE = -1d0
            RETURN

         ENDIF ! END Second IF

      ELSE

         IF ( N == IDTCO ) THEN ! NEW SECOND IF

            ! CO [molec/cm2/s]
            VALUE = CO_WKEND(I,J,L,IH)

         ELSE IF ( N == IDTNO ) THEN

            ! NO [molec/cm2/s]
            !NO_WKEND(I,J,L,IH) = NO_WKEND(I,J,L,IH)+NO2_WKEND(I,J,L,IH) & 
            !     + HNO2_WKEND(I,J,L,IH)
            VALUE = NO_WKEND(I,J,L,IH)
            
         ELSE IF ( N == IDTNO2 ) THEN

            ! NO2 [molec/cm2/s]
            VALUE = NO2_WKEND(I,J,L,IH)
        
         ELSE IF ( N == IDTHNO2 ) THEN

            ! HNO2 [molec/cm2/s]
            VALUE = HNO2_WKEND(I,J,L,IH)

         ELSE IF ( N == IDTSO2 ) THEN

            ! SO2 [molec/cm2/s]
            VALUE = SO2_WKEND(I,J,L,IH)

         ELSE IF ( N == IDTNH3 ) THEN

            ! NH3 [molec/cm2/s]
            VALUE = NH3_WKEND(I,J,L,IH)

         ELSE IF ( N == IDTALD2 ) THEN

            ! [molec/cm2/s]
            VALUE = ALD2_WKEND(I,J,L,IH)

         ELSE IF ( N == IDTRCHO ) THEN

            ! [molec/cm2/s]
            VALUE = RCHO_WKEND(I,J,L,IH)

         ELSE IF ( N == IDTBENZ ) THEN

            ! [molec/cm2/s]
            VALUE = BENZ_WKEND(I,J,L,IH)

         ELSE IF ( N == IDTC2H6 ) THEN

            ! [molec/cm2/s]
            VALUE = C2H6_WKEND(I,J,L,IH)

         ELSE IF ( N == IDTPRPE ) THEN

            ! [molec/cm2/s]
            VALUE = PRPE_WKEND(I,J,L,IH)

         ELSE IF ( N == IDTALK4 ) THEN

            ! [molec/cm2/s]
            VALUE = ALK4_WKEND(I,J,L,IH)*0.67
         
         ELSE IF ( N == IDTC3H8 ) THEN

            ! [molec/cm2/s]
            VALUE = C3H8_WKEND(I,J,L,IH)!*0.05*5

         ELSE IF ( N == IDTACET ) THEN

            ! [molec/cm2/s]
            VALUE = ACET_WKEND(I,J,L,IH)!*0.06

         ELSE IF ( N == IDTMEK ) THEN

            ! [molec/cm2/s]
            VALUE = MEK_WKEND(I,J,L,IH)!*0.02

        ELSE IF ( N == IDTTOLU ) THEN

            ! [molec/cm2/s]
            VALUE = TOLU_WKEND(I,J,L,IH)
         
         ELSE IF ( N == IDTXYLE ) THEN

            ! [molec/cm2/s]
            VALUE = XYLE_WKEND(I,J,L,IH)
  
        ! ELSE IF ( N == IDTC2H4 ) THEN

            ! [molec/cm2/s]
         !   VALUE = C2H4_WKEND(I,J,L,IH)

         ELSE IF ( N == IDTCH2O ) THEN

            ! [molec/cm2/s]
            VALUE = CH2O_WKEND(I,J,L,IH)
      
         ELSE IF ( N == IDTOCPO ) THEN

            ! [g/cm2/s]
            VALUE = OCPO_WKEND(I,J,L,IH)
        
         ELSE IF ( N == IDTBCPO ) THEN

            ! [g/cm2/s]
            VALUE = BCPO_WKEND(I,J,L,IH)

         ELSE IF ( N == IDTSO4 ) THEN

            ! [g/cm2/s]
            VALUE = SO4_WKEND(I,J,L,IH)

         !ELSE IF ( N == IDTEOH ) THEN

            ! [molec/cm2/s]
          !  VALUE = EOH_WKEND(I,J,L,IH)

         !ELSE IF ( N == IDTMOH ) THEN

            ! [molec/cm2/s]
          !  VALUE = MOH_WKEND(I,J,L,IH)

         !ELSE IF ( N == IDTCH4 ) THEN

            ! [molec/cm2/s]
          !  VALUE = CH4_WKEND(I,J,L,IH)

         ELSE
            ! Otherwise return a negative value to indicate
            ! that there are no NEI2008 emissions for tracer N
            VALUE = -1d0
            RETURN
         ENDIF !END SECOND IF
      ENDIF  !END FIRST IF

      ! Return to calling program
      END FUNCTION GET_NEI2008_ANTHRO
!EOC
!BOP
!------------------------------------------------------------------------------
! !IROUTINE: GET_MASK_FORFIRE
! !DESCRIPTION: Subroutine GET_MASK_FORFIRE initializes the mask for the fire injection routine

      SUBROUTINE GET_MASK_FORFIRE(  RC )
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

           CALL READ_NEI2008_MASK
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
! !IROUTINE: emiss_nei2008_anthro
!
! !DESCRIPTION: Subroutine EMISS\_NEI2008\_ANTHRO reads the NEI2008
!  emission fields at 1x1 resolution and regrids them to the 
!  current model resolution.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMISS_NEI2008_ANTHRO(  am_I_Root, Input_Opt, &
                                      State_Chm, RC )
!
! !USES:
! 
      USE DIRECTORY_MOD,     ONLY : DATA_DIR_NATIVE
      USE LOGICAL_MOD,       ONLY : LFUTURE
      USE CMN_O3_MOD
      USE CMN_SIZE_MOD
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE NCDF_MOD,          ONLY : NC_READ
      USE GRID_MOD,          ONLY : GET_XOFFSET
      USE GRID_MOD,          ONLY : GET_YOFFSET
      USE REGRID_A2A_MOD,    ONLY : DO_REGRID_A2A
      USE TIME_MOD,          ONLY : GET_YEAR, GET_MONTH, GET_DAY
      USE TIME_MOD,          ONLY : GET_DAY_OF_WEEK, GET_HOUR
      USE TRACERID_MOD,      ONLY : IDTCO, IDTNO, IDTNO2, IDTHNO2
      USE TRACERID_MOD,      ONLY : IDTSO2, IDTNH3
      USE TRACERID_MOD,      ONLY : IDTALD2, IDTRCHO, IDTC2H6
      USE TRACERID_MOD,      ONLY : IDTPRPE, IDTALK4!, IDTC2H4
      USE TRACERID_MOD,      ONLY : IDTBENZ, IDTTOLU, IDTXYLE
      USE TRACERID_MOD,      ONLY : IDTSO4, IDTCH2O
      USE TRACERID_MOD,      ONLY : IDTOCPO, IDTBCPO
      !USE TRACERID_MOD,      ONLY : IDTCH4, IDTEOH, IDTMOH

      USE CMN_SIZE_MOD            ! Size parameters
      USE CMN_O3_MOD              ! FSCALYR

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
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?!
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
      INTEGER                    :: I, J,  IH,  THISYEAR, THISMONTH, SNo
      INTEGER                    :: L, HH, KLM, SPECIES_ID(18)
      INTEGER                    :: OFFLINE_ID(15)
      INTEGER                    :: I0, J0
      INTEGER                    :: st3d(3), ct3d(3)
      INTEGER                    :: st4d(4), ct4da(4), ct4db(4)
      INTEGER                    :: fId1, fId1b, fId1c, fId1d, fId1e
      INTEGER                    :: fId2, fId2b, fId2c, fId2d, fId2e
      REAL*4                     :: ARRAYWD(225,202,24)
      REAL*4                     :: ARRAYWE(225,202,24) 
      REAL*4                     :: ARRAYWDPT(2,225,202,24)
      REAL*4                     :: ARRAYWEPT(2,225,202,24)
      REAL*4                     :: ARRAYWDPTN(3,225,202,24)
      REAL*4                     :: ARRAYWEPTN(3,225,202,24)
      REAL*4                     :: ARRAYWDC3(225,202,24)
      REAL*4                     :: ARRAYWEC3(225,202,24)
      REAL*4                     :: ARRAYWD_NH3ag(225,202,24)
      REAL*4                     :: ARRAYWE_NH3ag(225,202,24)
      REAL*8, TARGET             :: GEOS_NATIVEWD_NH3ag(I025x03125,J025x03125,24)
      REAL*8, TARGET             :: GEOS_NATIVEWE_NH3ag(I025x03125,J025x03125,24)
      REAL*8, TARGET             :: GEOS_NATIVEWD(I025x03125,J025x03125,3,24)
      REAL*8, TARGET             :: GEOS_NATIVEWE(I025x03125,J025x03125,3,24)
      REAL*4                     :: ScCO, ScNOx, ScPM10, ScPM25
      REAL*4                     :: ScVOC, ScNH3, ScSO2
      REAL*4                     :: ScNH3_Ag, ScNH3_NonAg
      CHARACTER(LEN=255)         :: DATA_DIR_NEI
      CHARACTER(LEN=255)         :: FILENAMEWD, FILENAMEWE
      CHARACTER(LEN=255)         :: FILENAMEWDPT, FILENAMEWEPT
      CHARACTER(LEN=255)         :: FILENAMEWDPTN, FILENAMEWEPTN
      CHARACTER(LEN=255)         :: FILENAMEWDC3, FILENAMEWEC3
      CHARACTER(LEN=4)           :: SId
      CHARACTER(LEN=255)         :: LLFILENAME
      CHARACTER(LEN=3)           :: TTMON
      CHARACTER(LEN=24)          :: SPCLIST(18)
      REAL*8, POINTER            :: OUTGRID(:,:) => NULL()
      REAL*8, POINTER            :: INGRID(:,:) => NULL()

      ! For scaling NH3 agricultural emissions (jaf, 12/10/13)
      REAL*4, POINTER            :: NCARR(:,:,:) => NULL()
      REAL*8                     :: ScAgNH3_MASAGE(225,202)
      REAL*8                     :: ScAgNH3_MASAGE2(I025x03125,J025x03125)
      LOGICAL                    :: LSCALE2MASAGE
      CHARACTER(LEN=255)         :: DATA_DIR_NH3_ag
      CHARACTER(LEN=255)         :: FILENAMEWD_NH3ag, FILENAMEWE_NH3ag
      CHARACTER(LEN=255)         :: FILENAME_ScAg

      !=================================================================
      ! EMISS_NEI2008_ANTHRO begins here!
      !=================================================================
      ! Assume success
      RC        =  GIGC_SUCCESS

      ! Copy values from Input_Opt
      LFUTURE   = Input_Opt%LFUTURE
      LSCALE2MASAGE = Input_Opt%LSCALE2MASAGE
      
      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_NEI2008_ANTHRO( am_I_Root, Input_Opt, RC )
         FIRST = .FALSE.
      ENDIF
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
      ScNH3_Ag    = 1.0
      ScNH3_NonAg = 1.0

      ! Apply annual scalar factor.
      ! Using EPA's National Tier1 CAPS (http://www.epa.gov/ttnchie1/trends/)
      ! Based on NEI08 v2
      IF ( THISYEAR == 2011 ) THEN  ! scale based on 2010
         ScCO        = 0.916
         ScNOx       = 0.897
         ScPM10      = 0.998
         ScPM25      = 0.988
         ScSO2       = 0.905
         ScVOC       = 0.955
         ScNH3       = 0.996
         ScNH3_Ag    = 1.000
         ScNH3_NonAg = 0.952
      ELSEIF ( THISYEAR .ge. 2012 ) THEN ! scale based on 2010
         ! To better model SEAC4RS reduce mobile NOx (Anderson et al etc)
         !  Mobile NOx is 38% of total according to NEI for 2013
         !  Scale down 50% of mobile = 0.38 & 0.5 = 0.19 (krt, 9/19/14)
         !  To better reconcile inventory, scale CO and NOx both by 50%
         ScNOx       = 0.794 * 0.5
         ScCO        = 0.820 * 0.5
         ScPM10      = 0.995
         ScPM25      = 0.975
         ! Scale down on the basis of surface station & wet deposition data
         ! skim, 9/18/14
         ScSO2       = 0.700
         !ScSO2       = 0.899
         ScVOC       = 0.905
         ScNH3       = 0.991
         ScNH3_Ag    = 1.000
         ScNH3_NonAg = 0.901
      ENDIF

      SPECIES_ID = (/ IDTCO,   IDTNO,  IDTNO2, IDTHNO2,           &
                      IDTSO2,  IDTNH3, IDTALD2, IDTRCHO, IDTC2H6, &
                      IDTPRPE, IDTALK4, IDTSO4, IDTCH2O, IDTOCPO, &
                      IDTBCPO, IDTTOLU, IDTXYLE,  IDTBENZ /)!, IDTC2H4/)
                      !IDTMOH, IDTEOH, IDTCH4/)

      SPCLIST =    (/ 'CO',   'NO',   'NO2',  'HNO2', 'SO2',  'NH3',  &
                      'ALD2', 'RCHO', 'C2H6', 'PRPE', 'ALK4', 'SO4',  &
                      'CH2O', 'OC',   'BC',   'TOLU', 'XYLE', 'BENZ'/)!, &
                      !'C2H4'/)!,'MOH', 'EOH','CH4' /)

     ! ID #'s for that are not tied to IDTxxxx flags
      OFFLINE_ID = (/ 2, 1, 64, 66, 26, 30, 11, 12, &
                     21, 18, 5, 27, 20, 36, 37    /)

      ! File with lat/lon edges for regridding
      LLFILENAME = TRIM( DATA_DIR_1x1) // &
                  'MAP_A2A_Regrid_201203/MAP_A2A_latlon_geos025x03125.nc'

      ! DataDir for year
      ! model ready
      DATA_DIR_NEI = TRIM( DATA_DIR_NATIVE ) // &
           'NEI08_2010_25x3125_'
   
      ! For NH3 -- files with agricultural emissions only (jaf, 12/10/13)
      ! Eventually these files will move to the data directory
!!!      DATA_DIR_NH3_ag = TRIM( Input_Opt%DATA_DIR_1x1 ) // &
!!!                        'NEI2008_201307/'
      DATA_DIR_NH3_ag = '/home/jaf/emissions/NH3/NH3_NEI08/ag_only/'

      ! Get month
      THISMONTH = GET_MONTH()

      ! GET NEI2008 FILES! 1 for wday, 1 for wkend
      IF (THISMONTH == 1)      THEN
         TTMON = 'Jan'
      ELSEIF (THISMONTH == 2)  THEN
         TTMON = 'Feb'
      ELSEIF (THISMONTH == 3)  THEN
         TTMON = 'Mar'
      ELSEIF (THISMONTH == 4)  THEN
         TTMON = 'Apr'
      ELSEIF (THISMONTH == 5)  THEN
         TTMON = 'May'
      ELSEIF (THISMONTH == 6)  THEN
         TTMON = 'Jun'
      ELSEIF (THISMONTH == 7)  THEN
         TTMON = 'Jul'
      ELSEIF (THISMONTH == 8)  THEN
         TTMON = 'Aug'
      ELSEIF (THISMONTH == 9)  THEN
         TTMON = 'Sep'
      ELSEIF (THISMONTH == 10) THEN
         TTMON = 'Oct'
      ELSEIF (THISMONTH == 11) THEN
         TTMON = 'Nov'
      ELSEIF (THISMONTH == 12) THEN
         TTMON = 'Dec'
      ENDIF

      ! model ready
      FILENAMEWD    = TRIM( DATA_DIR_NEI ) //                         &
                      TRIM( TTMON        ) // '_wkday_regrid.nc'
      FILENAMEWE    = TRIM( DATA_DIR_NEI ) //                        &
                      TRIM( TTMON        ) // '_wkend_regrid.nc'

      ! ptipm
      FILENAMEWDPT  = TRIM( DATA_DIR_NEI ) //  'ptipm_'           // &
                      TRIM( TTMON        ) // '_wkday_regrid.nc'
      FILENAMEWEPT  = TRIM( DATA_DIR_NEI ) // 'ptipm_'            // &
                      TRIM( TTMON        ) // '_wkend_regrid.nc'

      ! ptnonipm
      FILENAMEWDPTN = TRIM( DATA_DIR_NEI ) // 'ptnonipm_'         // &
                      TRIM( TTMON        ) //  '_wkday_regrid.nc'
      FILENAMEWEPTN = TRIM( DATA_DIR_NEI ) // 'ptnonipm_'         // &
                      TRIM( TTMON        ) //  '_wkend_regrid.nc'

      ! c3marine
      FILENAMEWDC3  = TRIM( DATA_DIR_NEI ) // 'c3marine_'         // &
                      TRIM( TTMON        ) //  '_wkday_regrid.nc'
      FILENAMEWEC3  = TRIM( DATA_DIR_NEI ) //  'c3marine_'        // &
                      TRIM( TTMON        ) // '_wkend_regrid.nc'

      ! Allocate start and count arrays
      st3d = (/1, 1, 1/)            !Start lat/lon/time
      st4d = (/1, 1, 1, 1/)         !Start lat/lon/time/lev
      ct3d = (/225, 202, 24/)     !Count lat/lon/time
      ct4da= (/2, 225, 202, 24/)  !Count lat/lon/time/lev - pt
      ct4db= (/3, 225, 202, 24/)  !Count lat/lon/time/lev - ptn

      ! Open weekday netCDF files for reading
      CALL Ncop_Rd(fId1,  TRIM(FILENAMEWD))
      CALL Ncop_Rd(fId1b, TRIM(FILENAMEWDPT))     ! ptipm
      CALL Ncop_Rd(fId1c, TRIM(FILENAMEWDPTN))    ! ptnonipm
      CALL Ncop_Rd(fId1d, TRIM(FILENAMEWDC3))     ! c3marine

      ! Open weekend netCDF files for reading
      CALL Ncop_Rd(fId2,  TRIM(FILENAMEWE))
      CALL Ncop_Rd(fId2b, TRIM(FILENAMEWEPT))     ! ptipm
      CALL Ncop_Rd(fId2c, TRIM(FILENAMEWEPTN))    ! ptnonipm
      CALL Ncop_Rd(fId2d, TRIM(FILENAMEWEC3))     ! c3marine

      ! Open NH3 ag files (only avail at 025x03125)
      IF ( LSCALE2MASAGE ) THEN
         FILENAMEWD_NH3ag = TRIM(DATA_DIR_NH3_ag) // &
                            'NEI08_2010_25x3125_' // &
                            TRIM( TTMON        )  // '_wkday_regrid.nc'
         FILENAMEWE_NH3ag = TRIM(DATA_DIR_NH3_ag) // &
                            'NEI08_2010_25x3125_' // &
                            TRIM( TTMON        )  // '_wkend_regrid.nc'
!        Kludge for SEAC4RS (skim, 3/6/14)
        FILENAME_ScAg    = TRIM(DATA_DIR_NH3_ag) // &
!       FILENAME_ScAg    = '/home/jaf/emissions/NH3/NH3_MASAGE/' // &
                           'MASAGE_NEI08_Ratio.geos.025x03125.nc'
         CALL Ncop_Rd(fId1e, TRIM(FILENAMEWD_NH3ag)) ! NH3ag weekday
         CALL Ncop_Rd(fId2e, TRIM(FILENAMEWE_NH3ag)) ! NH3ag weekend
      ENDIF

 100  FORMAT( '     - EMISS_NEI2008_ANTHRO:  Reading ', &
                      a, ' -> ', a )

      ! Loop over species
      DO KLM = 1, SIZE( SPCLIST )
 
         IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
            SId = SPCLIST( KLM )
            SNo = SPECIES_ID( KLM ) 
         ELSE
            SNo = OFFLINE_ID( KLM )
         ENDIF

         ! Skip undefined tracers
         IF ( SNo == 0 ) CYCLE

         GEOS_NATIVEWD = 0d0
         GEOS_NATIVEWE = 0d0
         GEOS_NATIVEWD_NH3ag = 0d0
         GEOS_NATIVEWE_NH3ag = 0d0

         ! Read variable from weekday netCDF files
         WRITE( 6, 100 )  TRIM( FILENAMEWD ), SID
         Call NcRd(ARRAYWD,       fId1,  TRIM(SId), st3d, ct3d )
         Call NcRd(ARRAYWDPT,     fId1b, TRIM(SId), st4d, ct4da)
         Call NcRd(ARRAYWDPTN,    fId1c, TRIM(SId), st4d, ct4db)
         Call NcRd(ARRAYWDC3,     fId1d, TRIM(SId), st3d, ct3d )

         GEOS_NATIVEWD(160:384,399:600,1,:) = ARRAYWD(:,:,:) + &
              ARRAYWDPT(1,:,:,:) + ARRAYWDPTN(1,:,:,:) + ARRAYWDC3(:,:,:)
         GEOS_NATIVEWD(160:384,399:600,2,:) = ARRAYWDPT(2,:,:,:) + &
              ARRAYWDPTN(2,:,:,:)
         GEOS_NATIVEWD(160:384,399:600,3,:) = ARRAYWDPTN(3,:,:,:)

         ! Read variable from weekend netCDF files
         WRITE( 6, 100 ) TRIM( FILENAMEWE ), SId
         Call NcRd(ARRAYWE,       fId2,  TRIM(SId), st3d, ct3d )
         Call NcRd(ARRAYWEPT,     fId2b, TRIM(SId), st4d, ct4da)
         Call NcRd(ARRAYWEPTN,    fId2c, TRIM(SId), st4d, ct4db)
         Call NcRd(ARRAYWEC3,     fId2d, TRIM(SId), st3d, ct3d )

         GEOS_NATIVEWE(160:384,399:600,1,:) = ARRAYWE(:,:,:) + &
              ARRAYWEPT(1,:,:,:) + ARRAYWEPTN(1,:,:,:) + ARRAYWEC3(:,:,:)
         GEOS_NATIVEWE(160:384,399:600,2,:) = ARRAYWEPT(2,:,:,:) + &
              ARRAYWEPTN(2,:,:,:)
         GEOS_NATIVEWE(160:384,399:600,3,:) = ARRAYWEPTN(3,:,:,:)

         IF ( LSCALE2MASAGE .and. SId == 'NH3' ) THEN
            CALL NcRd(ARRAYWD_NH3ag, fId1e, TRIM(SId), st3d, ct3d )
            CALL NcRd(ARRAYWE_NH3ag, fId2e, TRIM(SId), st3d, ct3d )
            ! Special case for NH3 emissions -- scale agricultural
            ! component based on MASAGE monthly gridded values from Paulot
            ! et al., 2013 (jaf, 12/10/13)
            ! Read ag files
            
            ! Cast to REAL*8
            GEOS_NATIVEWD_NH3ag(160:384,399:600,:) = ARRAYWD_NH3ag(:,:,:)
            GEOS_NATIVEWE_NH3ag(160:384,399:600,:) = ARRAYWE_NH3ag(:,:,:)

            ! Subtract agricultural component from total
            GEOS_NATIVEWD(:,:,1,:) = GEOS_NATIVEWD(:,:,1,:) - &
                 GEOS_NATIVEWD_NH3ag(:,:,:)
            GEOS_NATIVEWE(:,:,1,:) = GEOS_NATIVEWE(:,:,1,:) - &
                 GEOS_NATIVEWE_NH3ag(:,:,:)

            ! Read scaling factor (ratio of MASAGE to NEI08
            CALL NC_READ( NC_PATH = TRIM(FILENAME_ScAg),     &
                 PARA = 'ratio', ARRAY = NCARR,     & 
                 YEAR = 2010,    MONTH = THISMONTH, &
                 DAY = 01,       VERBOSE = .FALSE.    )
                     
            ! Cast to REAL*8
            ScAgNH3_MASAGE = NCARR(:,:,1)
            ScAgNH3_MASAGE2(160:384,399:600) = ScAgNH3_MASAGE(:,:)
            ! Deallocate ncdf-array
            IF ( ASSOCIATED ( NCARR ) ) DEALLOCATE ( NCARR )
            ! Scale agricultural component to MASAGE monthly totals
            DO HH = 1, 24
               GEOS_NATIVEWD_NH3ag(:,:,HH) = &
                    GEOS_NATIVEWD_NH3ag(:,:,HH) * ScAgNH3_MASAGE2
               GEOS_NATIVEWE_NH3ag(:,:,HH) = &
                    GEOS_NATIVEWE_NH3ag(:,:,HH) * ScAgNH3_MASAGE2
            ENDDO
            ! Add scaled agricultural component back to total and apply
            ! interannual scaling factors
            GEOS_NATIVEWD(:,:,1,:) = GEOS_NATIVEWD(:,:,1,:) * ScNH3_NonAg + &
                          GEOS_NATIVEWD_NH3ag(:,:,:) * ScNH3_Ag
            GEOS_NATIVEWE(:,:,1,:) = GEOS_NATIVEWE(:,:,1,:) * ScNH3_NonAg + &
                 GEOS_NATIVEWE_NH3ag(:,:,:) * ScNH3_Ag
         ENDIF
         
         ! Regrid from GEOS 0.25x0.3125 --> model resolution [molec/cm2/2]
         DO L=1,3
            DO HH=1,24
               !-------WEEKDAY------------
               ! Point to array slices
               INGRID  => GEOS_NATIVEWD(:,:,L,HH)
               OUTGRID => TMPARR_WD(:,:)
               ! Regrid
               CALL DO_REGRID_A2A( LLFILENAME, I025x03125,  J025x03125, &
                    INGRID,     OUTGRID, IS_MASS=0, &
                    netCDF=.TRUE.                   )
               ! Free pointers
               NULLIFY( INGRID, OUTGRID )
               !-------WEEKEND------------
               ! Point to array slices
               INGRID  => GEOS_NATIVEWE(:,:,L,HH)
               OUTGRID => TMPARR_WE(:,:)
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
                     TMP_WD(I,J,L,HH)=TMPARR_WD(I+I0,J+J0)
                     TMP_WE(I,J,L,HH)=TMPARR_WE(I+I0,J+J0)
                  END DO
               END DO
!$OMP END PARALLEL DO

               IF ( SId == 'CO' ) THEN
                  !-----------------
                  ! Apply scaling factor and masks
                  CO(:,:,L,HH)       = TMP_WD(:,:,L,HH)*ScCO * USA_MASK
                  CO_WKEND(:,:,L,HH) = TMP_WE(:,:,L,HH)*ScCO * USA_MASK
               ELSEIF ( SId == 'NO' ) THEN
                  ! Apply scaling factor and masks
                  NO(:,:,L,HH)       = TMP_WD(:,:,L,HH)*ScNOx * USA_MASK
                  NO_WKEND(:,:,L,HH) = TMP_WE(:,:,L,HH)*ScNOx * USA_MASK
               ELSEIF ( TRIM(SId) == 'NO2' ) THEN
                  ! Apply scaling factor and masks
                  NO2(:,:,L,HH)       = TMP_WD(:,:,L,HH)*ScNOx * USA_MASK
                  NO2_WKEND(:,:,L,HH) = TMP_WE(:,:,L,HH)*ScNOx * USA_MASK
               ELSEIF ( TRIM(SId) == 'HNO2' ) THEN
                  ! Apply scaling factor and masks
                  HNO2(:,:,L,HH)       = TMP_WD(:,:,L,HH)*ScNOx * USA_MASK
                  HNO2_WKEND(:,:,L,HH) = TMP_WE(:,:,L,HH)*ScNOx * USA_MASK
               ELSEIF ( TRIM(SId) == 'ALD2' ) THEN
                  ! Apply scaling factor and masks
                  ALD2(:,:,L,HH)       = TMP_WD(:,:,L,HH)*ScVOC * USA_MASK
                  ALD2_WKEND(:,:,L,HH) = TMP_WE(:,:,L,HH)*ScVOC * USA_MASK
               ELSEIF ( TRIM(SId) == 'RCHO' ) THEN
                  ! Apply scaling factor and masks
                  RCHO(:,:,L,HH)       = TMP_WD(:,:,L,HH)*ScVOC * USA_MASK
                  RCHO_WKEND(:,:,L,HH) = TMP_WE(:,:,L,HH)*ScVOC * USA_MASK
               ELSE IF ( TRIM(SId) == 'BENZ' ) THEN
                  ! Apply scaling factor and masks
                  BENZ(:,:,L,HH)       = TMP_WD(:,:,L,HH)*ScVOC * USA_MASK
                  BENZ_WKEND(:,:,L,HH) = TMP_WE(:,:,L,HH)*ScVOC * USA_MASK
               ELSE IF ( TRIM(SId) == 'TOLU' ) THEN
                  ! Apply scaling factor and masks
                  TOLU(:,:,L,HH)       = TMP_WD(:,:,L,HH)*ScVOC * USA_MASK
                  TOLU_WKEND(:,:,L,HH) = TMP_WE(:,:,L,HH)*ScVOC * USA_MASK
               ELSE IF ( TRIM(SId) == 'XYLE' ) THEN
                  ! Apply scaling factor and masks
                  XYLE(:,:,L,HH)       = TMP_WD(:,:,L,HH)*ScVOC * USA_MASK
                  XYLE_WKEND(:,:,L,HH) = TMP_WE(:,:,L,HH)*ScVOC * USA_MASK
               ELSE IF ( TRIM(SId) == 'C2H6' ) THEN
                  ! Apply scaling factor and masks
                  C2H6(:,:,L,HH)       = TMP_WD(:,:,L,HH)*ScVOC * USA_MASK
                  C2H6_WKEND(:,:,L,HH) = TMP_WE(:,:,L,HH)*ScVOC * USA_MASK
               ELSE IF ( TRIM(SId) == 'PRPE' ) THEN
                  ! Apply scaling factor and masks
                  PRPE(:,:,L,HH)       = TMP_WD(:,:,L,HH)*ScVOC * USA_MASK
                  PRPE_WKEND(:,:,L,HH) = TMP_WE(:,:,L,HH)*ScVOC * USA_MASK
               ELSE IF ( TRIM(SId) == 'ALK4' ) THEN
                  ! Apply scaling factor and masks
                  ALK4(:,:,L,HH)       = TMP_WD(:,:,L,HH)* 0.87 * ScVOC * USA_MASK
                  ALK4_WKEND(:,:,L,HH) = TMP_WE(:,:,L,HH)* 0.87 * ScVOC * USA_MASK
                  C3H8(:,:,L,HH)       = TMP_WD(:,:,L,HH)* 0.03 * ScVOC * USA_MASK
                  C3H8_WKEND(:,:,L,HH) = TMP_WE(:,:,L,HH)* 0.03 * ScVOC * USA_MASK
                  ACET(:,:,L,HH)       = TMP_WD(:,:,L,HH)* 0.06 * ScVOC * USA_MASK
                  ACET_WKEND(:,:,L,HH) = TMP_WE(:,:,L,HH)* 0.06 * ScVOC * USA_MASK
                  MEK(:,:,L,HH)        = TMP_WD(:,:,L,HH)* 0.02 * ScVOC * USA_MASK
                  MEK_WKEND(:,:,L,HH)  = TMP_WE(:,:,L,HH)* 0.02 * ScVOC * USA_MASK
                  ! I used the graph of the top 50 VOCs from the following paper:
                  ! http://www.epa.gov/ttnchie1/software/speciate/atmospheric.pdf
                  ! And determined the PAR weighting from the following report:
                  ! http://www.camx.com/publ/pdfs/cb05_final_report_120805.pdf
                  ! I then determined the fraction by carbon of PAR for each of the
                  !  following species.
                  ! Note that it will not add to 100% since some species are double
                  ! counted (benzene) (krt, 2/3/15)
               ELSE IF ( TRIM(SId) == 'CH2O' ) THEN
                  ! Apply scaling factor and masks
                  CH2O(:,:,L,HH)       = TMP_WD(:,:,L,HH)*ScVOC * USA_MASK
                  CH2O_WKEND(:,:,L,HH) = TMP_WE(:,:,L,HH)*ScVOC * USA_MASK
               ELSE IF ( TRIM(SId) == 'BC' ) THEN
                  ! Scale down by 30% on the basis of comparisons to
                  ! SEAC4RS (skim, 5/23/14)
                  ! Apply scaling factor and masks
                  BCPO(:,:,L,HH)       = TMP_WD(:,:,L,HH)*ScPM25 * USA_MASK *0.7
                  BCPO_WKEND(:,:,L,HH) = TMP_WE(:,:,L,HH)*ScPM25 * USA_MASK *0.7
               ELSE IF ( TRIM(SId) == 'OC' ) THEN
                  ! Apply scaling factor and masks
                  OCPO(:,:,L,HH)       = TMP_WD(:,:,L,HH)*ScPM25 * USA_MASK
                  OCPO_WKEND(:,:,L,HH) = TMP_WE(:,:,L,HH)*ScPM25 * USA_MASK
               ELSE IF ( TRIM(SId) == 'SO2' ) THEN
                  ! Apply scaling factor and masks
                  SO2(:,:,L,HH)       = TMP_WD(:,:,L,HH)*ScSO2 * USA_MASK
                  SO2_WKEND(:,:,L,HH) = TMP_WE(:,:,L,HH)*ScSO2 * USA_MASK
               ELSE IF ( TRIM(SId) == 'SO4' ) THEN
                  ! Apply scaling factor and masks
                  SO4(:,:,L,HH)       = TMP_WD(:,:,L,HH)*ScSO2 * USA_MASK
                  SO4_WKEND(:,:,L,HH) = TMP_WE(:,:,L,HH)*ScSO2 * USA_MASK
                  ! ELSE IF ( TRIM(SId) == 'C2H4' ) THEN
                  ! ELSE IF ( TRIM(SId) == 'CH4' ) THEN
                  ! ELSE IF ( TRIM(SId) == 'EOH' ) THEN
                  ! ELSE IF ( TRIM(SId) == 'MOH' ) THEN
               ELSEIF ( TRIM(SId) == 'NH3' ) THEN
                  IF ( .not. LSCALE2MASAGE ) THEN
                     ! If we can't separate out the agricultural component
                     ! (e.g. for the 05x0667 grid), then just apply the single
                     ! annual scaling factor to NH3 emissions.
                     NH3(:,:,L,HH)  = TMP_WD(:,:,L,HH)* ScNH3 * USA_MASK
                     NH3_WKEND(:,:,L,HH)  = TMP_WE(:,:,L,HH)* ScNH3 * USA_MASK
                  ELSE
                     ! Apply masks
                     NH3(:,:,L,HH)       = TMP_WD(:,:,L,HH) * USA_MASK
                     NH3(:,:,L,HH)       = TMP_WD(:,:,L,HH) * USA_MASK
                  ENDIF  ! MASAGE scaling                 
               ENDIF
            ENDDO
         ENDDO

      ENDDO ! END SPECIES

      ! Close netCDF files
      CALL NcCl( fId1  )
      CALL NcCl( fId1b )
      CALL NcCl( fId1c )
      CALL NcCl( fId1d )
      CALL NcCl( fId2  )
      CALL NcCl( fId2b )
      CALL NcCl( fId2c )
      CALL NcCl( fId2d )

      ! Close ag files
      IF ( LSCALE2MASAGE ) THEN
         CALL NcCl( fId1e )
         CALL NcCl( fId2e )
      ENDIF

      !--------------------------
      ! Compute future emissions
      !--------------------------
      IF ( LFUTURE ) THEN 
         CALL NEI2008_SCALE_FUTURE
      ENDIF

      !--------------------------
      ! Print emission totals
      !--------------------------

      CALL TOTAL_ANTHRO_Tg( THISMONTH )

      ! Return to calling program
      END SUBROUTINE EMISS_NEI2008_ANTHRO
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_nei2008_mask
!
! !DESCRIPTION: Subroutine READ\_NEI2008\_MASK reads the mask for NEI data  
!\\
!\\
! !INTERFACE:
      
      SUBROUTINE READ_NEI2008_MASK
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
      ! Mask specific to NEI2008 data
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
200   FORMAT( '     - READ_NEI2008_MASK: Reading ', a )
     
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
      WRITE(*,*) 'READ NEI2008 MASK!'
      ! Return to calling program
      END SUBROUTINE READ_NEI2008_MASK
!------------------------------------------------------------------------------
! Prior to 12/3/09:
! Leave for future use (bmy, 12/3/09)
!!EOC
!!------------------------------------------------------------------------------
!!          Harvard University Atmospheric Chemistry Modeling Group            !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: get_nei2005_mask
!!
!! !DESCRIPTION: Subroutine GET\_NEI2005\_MASK returns the value of the 
!!  NEI 2005 mask to the calling program.  Values of 1 denote grid boxes 
!!  within the EPA/NEI2005 emission region.!  
!!\\
!!\\
!! !INTERFACE:
!      
!      FUNCTION GET_NEI2005_MASK( I, J ) RESULT ( USA )
!!
!! !INPUT PARAMETERS:
!!     
!      INTEGER, INTENT(IN) :: I, J   ! GEOS-Chem lon & lat indices
!!
!! !RETURN VALUE:
!!
!      REAL*8              :: USA    ! Value of the mask  
!!
!! !REMARKS:
!!  This is entended to encapsulate the USA_MASK variable.
!!     
!! !REVISION HISTORY: 
!!  02 Dec 2009 - R. Yantosca - Initial version
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!      USA = USA_MASK(I,J)
!
!      END FUNCTION GET_NEI2005_MASK
!------------------------------------------------------------------------------
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: nei2008_scale_future
!
! !DESCRIPTION: Subroutine NEI2008\_SCALE\_FUTURE applies the IPCC future 
!  scale factors to the NEI2008 anthropogenic emissions.
!\\
!\\
! !INTERFACE:

      SUBROUTINE NEI2008_SCALE_FUTURE
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
      ! NEI2008_SCALE_FUTURE begins here!
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
          SO2(I,J,L,HH) = SO2(I,J,L,HH) * GET_FUTURE_SCALE_SO2ff( I, J )

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
      END SUBROUTINE NEI2008_SCALE_FUTURE
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
      SUBROUTINE TOTAL_ANTHRO_TG( MONTH )
!
! !USES:
! 
      USE CMN_SIZE_MOD            ! Size parameters
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TRACER_MOD,   ONLY : XNUMOL
      USE TRACERID_MOD, ONLY : IDTCO, IDTNO,IDTNO2, IDTHNO2 
      USE TRACERID_MOD, ONLY : IDTSO2, IDTNH3, IDTMEK
      USE TRACERID_MOD, ONLY : IDTALD2, IDTRCHO, IDTC2H6
      USE TRACERID_MOD, ONLY : IDTPRPE, IDTALK4!, IDTC2H4
      USE TRACERID_MOD, ONLY : IDTBENZ, IDTTOLU, IDTXYLE
      USE TRACERID_MOD, ONLY : IDTSO4, IDTCH2O, IDTACET, IDTC3H8
      USE TRACERID_MOD, ONLY : IDTOCPO,  IDTBCPO 
      !USE TRACERID_MOD, ONLY : IDTMOH, IDTEOH, IDTCH4

! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: MONTH   ! Month of data to compute totals
!
! !REVISION HISTORY: 
!   7 Oct 2009 - A. van Donkelaar - initial version
!   9 May 2013 - K. Travis - revised for NEI2008
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: II, JJ, IH, LL
      INTEGER             :: DAY_LIST(12), DL
      REAL*8              :: T_CO, T_NO, T_NO2, T_HNO2, T_SO2, T_NH3
      REAL*8              :: T_ALD2,  T_RCHO, T_C2H6, T_C3H8, T_ACET
      REAL*8              :: T_PRPE, T_ALK4, T_TOLU, T_XYLE, T_MEK
      REAL*8              :: T_CH2O,T_BC, T_OC, T_SO4
      REAL*8              :: T_BENZ!, T_C2H4
      REAL*8              :: tmpArea(IIPAR, JJPAR,3)
      REAL*4              :: WDFRAC, WEFRAC
      CHARACTER(LEN=3)    :: UNIT
      REAL*8,  PARAMETER  :: SEC_IN_HOUR  = 3600d0! * 365.25d0
                                 
      !=================================================================
      ! TOTAL_ANTHRO_TG begins here!
      !=================================================================

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, 100  )
 100  FORMAT( 'N. E. I. 2008 U. S. A.   E M I S S I O N S', / )
      
      DO II = 1, IIPAR
         DO JJ = 1, JJPAR
            DO LL=1, 3
               tmpArea(II,JJ,LL) = GET_AREA_CM2(II,JJ,LL)
            ENDDO
         ENDDO
      ENDDO
      
      !           J,  F,  M,  A,  Ma,  Ju, J, Au, Se, Oc, No, Dec
      !DAY_LIST = (/31,28,31,30,31,30,31,31,30,31,30,31/)
      !DL = DAY_LIST(MONTH)
      ! Annual average weekends and weekdays
      WDFRAC = 21.7d0
      WEFRAC = 8.7d0
      
      WRITE(6,101) WDFRAC
 101  FORMAT('WEEKDAY FRACTION = ', f11.4)

      ! Total CO  [Tg CO]
      IF ( IDTCO .NE. 0 ) &
      T_CO  = SUM(SUM( CO,4)  * tmpArea ) * &
              SEC_IN_HOUR *1d-9/XNUMOL(IDTCO)* WDFRAC + &
              SUM(SUM( CO_WKEND,4)  * tmpArea ) * &
              SEC_IN_HOUR *1d-9/XNUMOL(IDTCO) * WEFRAC

      IF ( IDTNO .NE. 0 ) &
      ! Total NO [Tg N]
      T_NO =   SUM(SUM(NO, 4)*tmpArea ) * &
               SEC_IN_HOUR *1d-9/XNUMOL(IDTNO)*14d0/30d0 * WDFRAC + &
               SUM(SUM(NO_WKEND, 4)*tmpArea ) * &
               SEC_IN_HOUR *1d-9/XNUMOL(IDTNO)*14d0/30d0 * WEFRAC
      
      IF ( IDTNO2 .NE. 0 ) &
      ! Total NO2 [Tg N]
      T_NO2 =  SUM(SUM( NO2, 4 ) * tmpArea) * &
               SEC_IN_HOUR *1d-9/XNUMOL(IDTNO2)*14d0/46d0 * WDFRAC + &
               SUM(SUM( NO2_WKEND, 4 ) * tmpArea) * &
               SEC_IN_HOUR *1d-9/XNUMOL(IDTNO2)*14d0/46 * WEFRAC

      IF ( IDTHNO2 .NE. 0 ) &
      ! Total HNO2 [Tg N]
      T_HNO2 = SUM(SUM( HNO2,4 )  * tmpArea) * &
               SEC_IN_HOUR *1d-9/XNUMOL(IDTHNO2)*14d0/47d0 * WDFRAC + &
               SUM(SUM( HNO2_WKEND,4 )  * tmpArea) * &
               SEC_IN_HOUR *1d-9/XNUMOL(IDTHNO2)*14d0/47d0 * WEFRAC
     
      ! Total SO2 [Tg S]
      IF ( IDTSO2 .NE. 0 ) &
      T_SO2 =  SUM( SUM( SO2,4)  * tmpArea )* &
               SEC_IN_HOUR *1d-9/XNUMOL(IDTSO2)*32d0/64d0 * WDFRAC + &
               SUM(SUM( SO2_WKEND,4 )  * tmpArea) * &
               SEC_IN_HOUR *1d-9/XNUMOL(IDTSO2)*32d0/64d0 * WEFRAC

      ! Total NH3 [Tg NH3]
      IF ( IDTNH3 .NE. 0 ) &
      T_NH3 =  SUM( SUM( NH3,4)  * tmpArea) * &
               SEC_IN_HOUR *1d-9/XNUMOL(IDTNH3) * WDFRAC + &
               SUM(SUM( NH3_WKEND,4 )  * tmpArea) * &
               SEC_IN_HOUR *1d-9/XNUMOL(IDTNH3) * WEFRAC

      ! Total ALD2 [Tg C]
      IF ( IDTALD2 .NE. 0 ) &
      T_ALD2 = SUM( SUM( ALD2,4)  * tmpArea) * &
               SEC_IN_HOUR *1d-9/XNUMOL(IDTALD2) * WDFRAC + &
               SUM(SUM( ALD2_WKEND,4 )  * tmpArea) * &
               SEC_IN_HOUR *1d-9/XNUMOL(IDTALD2) * WEFRAC

      ! Total RCHO [Tg C]
      IF ( IDTRCHO .NE. 0 ) &
      T_RCHO = SUM( SUM( RCHO,4)  * tmpArea) * &
               SEC_IN_HOUR *1d-9/XNUMOL(IDTRCHO) * WDFRAC + &
               SUM(SUM( RCHO_WKEND,4 )  * tmpArea) * &
               SEC_IN_HOUR *1d-9/XNUMOL(IDTRCHO) * WEFRAC

      ! Total BENZ [Tg C]
      IF ( IDTBENZ .NE. 0 ) &
      T_BENZ = SUM( SUM( BENZ,4)  * tmpArea )* &
               SEC_IN_HOUR *1d-9/XNUMOL(IDTBENZ) * WDFRAC + &
               SUM(SUM( BENZ_WKEND,4 )  * tmpArea) * &
               SEC_IN_HOUR *1d-9/XNUMOL(IDTBENZ) * WEFRAC

      ! Total C2H6 [Tg C]
      IF ( IDTC2H6 .NE. 0 ) &
      T_C2H6 = SUM( SUM( C2H6,4)  * tmpArea )* &
               SEC_IN_HOUR *1d-9/XNUMOL(IDTC2H6) * WDFRAC + &
               SUM(SUM( C2H6_WKEND,4 )  * tmpArea) * &
               SEC_IN_HOUR *1d-9/XNUMOL(IDTC2H6) * WEFRAC

      ! Total PRPE [Tg C]
      IF ( IDTPRPE .NE. 0 ) &
      T_PRPE = SUM( SUM( PRPE,4)  * tmpArea )* &
               SEC_IN_HOUR *1d-9/XNUMOL(IDTPRPE) * WDFRAC + &
               SUM(SUM( PRPE_WKEND,4 )  * tmpArea) * &
               SEC_IN_HOUR *1d-9/XNUMOL(IDTPRPE) * WEFRAC

      ! Total ALK4 [Tg C]
      IF ( IDTALK4 .NE. 0 ) &
      T_ALK4 = SUM( SUM( ALK4,4)  * tmpArea )* &
               SEC_IN_HOUR *1d-9/XNUMOL(IDTALK4) * WDFRAC + &
               SUM(SUM( ALK4_WKEND,4 )  * tmpArea) * &
               SEC_IN_HOUR *1d-9/XNUMOL(IDTALK4)* WEFRAC

      ! Total C3H8 [Tg C]
      IF ( IDTC3H8 .NE. 0 ) &
      T_C3H8 = SUM( SUM( C3H8,4)  * tmpArea )* &
               SEC_IN_HOUR *1d-9/XNUMOL(IDTC3H8) * WDFRAC + &
               SUM(SUM( C3H8_WKEND,4 )  * tmpArea) * &
               SEC_IN_HOUR *1d-9/XNUMOL(IDTC3H8)* WEFRAC

      ! Total MEK [Tg C]
      IF ( IDTMEK .NE. 0 ) &
      T_MEK = SUM( SUM( MEK,4)  * tmpArea )* &
               SEC_IN_HOUR *1d-9/XNUMOL(IDTMEK) * WDFRAC + &
               SUM(SUM( MEK_WKEND,4 )  * tmpArea) * &
               SEC_IN_HOUR *1d-9/XNUMOL(IDTMEK)* WEFRAC

      ! Total ACET [Tg C]
      IF ( IDTACET .NE. 0 ) &
      T_ACET = SUM( SUM( ACET,4)  * tmpArea )* &
               SEC_IN_HOUR *1d-9/XNUMOL(IDTACET) * WDFRAC + &
               SUM(SUM( ACET_WKEND,4 )  * tmpArea) * &
               SEC_IN_HOUR *1d-9/XNUMOL(IDTACET)* WEFRAC

      ! Total TOLU [Tg C]
      IF ( IDTTOLU .NE. 0 ) &
      T_TOLU = SUM( SUM( TOLU,4)  *tmpArea )* &
               SEC_IN_HOUR *1d-9/XNUMOL(IDTTOLU) * WDFRAC + &
               SUM(SUM( TOLU_WKEND,4 )  * tmpArea) * &
               SEC_IN_HOUR *1d-9/XNUMOL(IDTTOLU) * WEFRAC

      ! Total XYLE [Tg C]
      IF ( IDTXYLE .NE. 0 ) &
      T_XYLE = SUM( SUM( XYLE,4)  * tmpArea )* &
               SEC_IN_HOUR *1d-9/XNUMOL(IDTXYLE) * WDFRAC + &
               SUM(SUM( XYLE_WKEND,4 )  * tmpArea) * &
               SEC_IN_HOUR *1d-9/XNUMOL(IDTXYLE) * WEFRAC

      ! Total CH2O [Tg CH2O]
      IF ( IDTCH2O .NE. 0 ) &
      T_CH2O = SUM( SUM( CH2O,4)  * tmpArea) * &
               SEC_IN_HOUR *1d-9/XNUMOL(IDTCH2O) * WDFRAC + &
               SUM(SUM( CH2O_WKEND,4 )  * tmpArea) * &
               SEC_IN_HOUR *1d-9/XNUMOL(IDTCH2O) * WEFRAC
      
      ! Total BC [Tg]
      IF ( IDTBCPO .NE. 0 ) &
      T_BC =  SUM( SUM( BCPO,4)   * tmpArea) * &
               SEC_IN_HOUR *1d-12 * WDFRAC + &
               SUM(SUM( BCPO_WKEND,4 )  * tmpArea) * &
               SEC_IN_HOUR *1d-12 * WEFRAC

      ! Total OC [Tg]
      IF ( IDTOCPO .NE. 0 ) &
      T_OC =  SUM( SUM( OCPO,4)  *  tmpArea )* &
              SEC_IN_HOUR *1d-12 * WDFRAC + &
               SUM(SUM( OCPO_WKEND,4 )  * tmpArea) * &
               SEC_IN_HOUR *1d-12 * WEFRAC
      
      ! Total SO4 [Tg S]
      IF ( IDTSO4 .NE. 0 ) &
      T_SO4 =  SUM( SUM( SO4,4)  *  tmpArea) * &
               SEC_IN_HOUR *1d-12 * 32d0/96d0 * WDFRAC + &
               SUM(SUM( SO4_WKEND,4 )  * tmpArea) * &
               SEC_IN_HOUR *1d-12 * 32d0/96d0 * WEFRAC

      ! Total C2H4 [Tg C]
      !T_C2H4 = SUM( C2H4 * tmpArea )* &
      !         SEC_IN_HOUR *1d-9/XNUMOL(IDTC2H4) * WDFRAC + &
      !         SUM(SUM( C2H4_WKEND,4 )  * tmpArea) * &
      !         SEC_IN_HOUR *1d-12/XNUMOL(IDTC2H4) * WEFRAC

      ! Print totals in [Tg]
      WRITE( 6, 110 ) 'CO   ', MONTH, T_CO,   '[Tg CO ]'
      WRITE( 6, 110 ) 'NO   ', MONTH, T_NO,   '[Tg N ]'
      WRITE( 6, 110 ) 'NO2  ', MONTH, T_NO2,  '[Tg N ]'
      WRITE( 6, 110 ) 'HNO2 ', MONTH, T_HNO2, '[Tg N ]'
      WRITE( 6, 110 ) 'SO2  ', MONTH, T_SO2,  '[Tg S]'
      WRITE( 6, 110 ) 'SO4  ', MONTH, T_SO4,  '[Tg S]'
      WRITE( 6, 110 ) 'NH3  ', MONTH, T_NH3,  '[Tg NH3]'
      WRITE( 6, 110 ) 'ALD2 ', MONTH, T_ALD2, '[Tg C]'
      WRITE( 6, 110 ) 'RCHO ', MONTH, T_RCHO, '[Tg C]'
      WRITE( 6, 110 ) 'BENZ ', MONTH, T_BENZ, '[Tg C]'
      WRITE( 6, 110 ) 'C2H6 ', MONTH, T_C2H6, '[Tg C]'
      WRITE( 6, 110 ) 'PRPE ', MONTH, T_PRPE, '[Tg C]'
      WRITE( 6, 110 ) 'ALK4 ', MONTH, T_ALK4, '[Tg C]'
      WRITE( 6, 110 ) 'C3H8 ', MONTH, T_C3H8, '[Tg C]'
      WRITE( 6, 110 ) 'MEK ',  MONTH, T_MEK, '[Tg C]'
      WRITE( 6, 110 ) 'ACET ', MONTH, T_ACET, '[Tg C]'
      WRITE( 6, 110 ) 'TOLU ', MONTH, T_TOLU, '[Tg C]'
      WRITE( 6, 110 ) 'XYLE ', MONTH, T_XYLE, '[Tg C]'
      WRITE( 6, 110 ) 'CH2O ', MONTH, T_CH2O, '[Tg CH2O]'
      WRITE( 6, 110 ) 'BC   ', MONTH, T_BC,   '[Tg ]'
      WRITE( 6, 110 ) 'OC   ', MONTH, T_OC,   '[Tg ]'
      !WRITE( 6, 110 ) 'C2H4 ', MONTH, T_C2H4, '[Tg C]'

      ! Format statement
 110  FORMAT( 'NEI2008 anthro ', a5, &
             'for month', i4, ': ', f11.4, 1x, a8 )

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
! !IROUTINE: init_nei2008_anthro
!
! !DESCRIPTION: Subroutine INIT\_NEI2008\_ANTHRO allocates and zeroes all 
!  module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_NEI2008_ANTHRO( am_I_Root, Input_Opt, RC )
!
! !USES:
! 
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE ERROR_MOD,   ONLY : ALLOC_ERR

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
      ! INIT_NEI2008_ANTHRO begins here!
      !=================================================================
      ! Assume success
      RC        =  GIGC_SUCCESS
      
      ! Return if LNEI08 is false
      IF ( .not. Input_Opt%LNEI08 ) RETURN
      
      !--------------------------------------------------
      ! Allocate and zero arrays for emissions
      !--------------------------------------------------
      
      IF ( .not. Input_Opt%LBIOMASS ) THEN
      ! allocate and read USA Mask
         ALLOCATE( USA_MASK( IIPAR, JJPAR ), STAT=RC )
         IF ( RC /= 0 ) CALL ALLOC_ERR( 'USA_MASK' )
         USA_MASK = 0d0
         
         CALL READ_NEI2008_MASK
      ENDIF

      ALLOCATE( TMP_WD( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'TMP_WD' )
      TMP_WD = 0d0

      ALLOCATE( TMPARR_WD( IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'TMPARR_WD' )
      TMPARR_WD = 0d0

      ALLOCATE( TMP_WE( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'TMP_WE' )
      TMP_WE = 0d0

      ALLOCATE( TMPARR_WE( IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'TMPARR_WE' )
      TMPARR_WE = 0d0

      ALLOCATE( CO( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'CO' )
      CO = 0d0

      !ALLOCATE( NOX( IIPAR, JJPAR, 3, 24 ), STAT=RC)
      !IF ( RC /= 0 ) CALL ALLOC_ERR( 'NOX' )
      !NOX = 0d0

      ALLOCATE( NO( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'NO' )
      NO = 0d0

      ALLOCATE( NO2( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'NO2' )
      NO2 = 0d0
      
      ALLOCATE( HNO2( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'HNO2' )
      HNO2 = 0d0

      ALLOCATE( SO2( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'SO2' )
      SO2 = 0d0

      ALLOCATE( NH3( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'NH3' )
      NH3 = 0d0 

      ALLOCATE( ALD2( IIPAR, JJPAR, 3, 24), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'ALD2' )
      ALD2 = 0d0 

      ALLOCATE( RCHO( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'RCHO' )
      RCHO = 0d0

      ALLOCATE( BENZ( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'BENZ' )
      BENZ = 0d0

      ALLOCATE( C2H6( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'C2H6' )
      C2H6 = 0d0

      ALLOCATE( PRPE( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'PRPE' )
      PRPE = 0d0

      ALLOCATE( ALK4( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'ALK4' )
      ALK4 = 0d0 

      ALLOCATE( C3H8( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'C3H8' )
      C3H8 = 0d0 

      ALLOCATE( MEK( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'MEK' )
      MEK = 0d0 

      ALLOCATE( ACET( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'ACET' )
      ACET = 0d0 

      ALLOCATE( TOLU( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'TOLU' )
      TOLU = 0d0 

      ALLOCATE( XYLE( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'XYLE' )
      XYLE = 0d0 
     
      !ALLOCATE( C2H4( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      !IF ( RC /= 0 ) CALL ALLOC_ERR( 'C2H4' )
      !C2H4 = 0d0 

      ALLOCATE( CH2O( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'CH2O' )
      CH2O = 0d0
      
      ALLOCATE( BCPO( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'BCPO' )
      BCPO = 0d0
      
      ALLOCATE( OCPO( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'OCPO' )
      OCPO = 0d0

      ALLOCATE( SO4( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'SO4' )
      SO4 = 0d0
      
      !ALLOCATE( EOH( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      !IF ( RC /= 0 ) CALL ALLOC_ERR( 'EOH' )
      !EOH = 0d0

      !ALLOCATE( MOH( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      !IF ( RC /= 0 ) CALL ALLOC_ERR( 'MOH' )
      !MOH = 0d0

      !ALLOCATE( CH4( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      !IF ( RC /= 0 ) CALL ALLOC_ERR( 'CH4' )
      !CH4 = 0d0
      
!     Weekend

      ALLOCATE( CO_WKEND( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'CO_WKEND' )
      CO_WKEND = 0d0

      !ALLOCATE( NOX_WKEND( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      !IF ( RC /= 0 ) CALL ALLOC_ERR( 'NOX_WKEND' )
      !NOX_WKEND = 0d0
      
      ALLOCATE( NO_WKEND( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'NO_WKEND' )
      NO_WKEND = 0d0

      ALLOCATE( NO2_WKEND( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'NO2_WKEND' )
      NO2_WKEND = 0d0

      ALLOCATE( HNO2_WKEND( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'HNO2_WKEND' )
      HNO2_WKEND = 0d0

      ALLOCATE( SO2_WKEND( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'SO2_WKEND' )
      SO2_WKEND = 0d0

      ALLOCATE( NH3_WKEND( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'NH3_WKEND' )
      NH3_WKEND = 0d0 

      ALLOCATE( ALD2_WKEND( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'ALD2_WKEND' )
      ALD2_WKEND = 0d0 

      ALLOCATE( RCHO_WKEND( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'RCHO_WKEND' )
      RCHO_WKEND = 0d0

      ALLOCATE( BENZ_WKEND( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'BENZ_WKEND' )
      BENZ_WKEND = 0d0

      ALLOCATE( C2H6_WKEND( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'C2H6_WKEND' )
      C2H6_WKEND = 0d0

      ALLOCATE( PRPE_WKEND( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'PRPE_WKEND' )
      PRPE_WKEND = 0d0

      ALLOCATE( ALK4_WKEND( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'ALK4_WKEND' )
      ALK4_WKEND = 0d0 

      ALLOCATE( C3H8_WKEND( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'C3H8_WKEND' )
      C3H8_WKEND = 0d0 

      ALLOCATE( MEK_WKEND( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'MEK_WKEND' )
      MEK_WKEND = 0d0 

      ALLOCATE( ACET_WKEND( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'ACET_WKEND' )
      ACET_WKEND = 0d0 

      ALLOCATE( TOLU_WKEND( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'TOLU_WKEND' )
      TOLU_WKEND = 0d0 

      ALLOCATE( XYLE_WKEND( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'XYLE_WKEND' )
      XYLE_WKEND = 0d0 
     
      !ALLOCATE( C2H4_WKEND( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      !IF ( RC /= 0 ) CALL ALLOC_ERR( 'C2H4_WKEND' )
      !C2H4_WKEND = 0d0 

      ALLOCATE( CH2O_WKEND( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'CH2O_WKEND' )
      CH2O_WKEND = 0d0
      
      ALLOCATE( BCPO_WKEND( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'BCPO_WKEND' )
      BCPO_WKEND = 0d0
     
      ALLOCATE( OCPO_WKEND( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'OCPO_WKEND' )
      OCPO_WKEND = 0d0

      ALLOCATE( SO4_WKEND( IIPAR, JJPAR, 3, 24  ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'SO4_WKEND' )
      SO4_WKEND = 0d0

      !ALLOCATE( EOH_WKEND( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      !IF ( RC /= 0 ) CALL ALLOC_ERR( 'EOH_WKEND' )
      !EOH_WKEND = 0d0

      !ALLOCATE( MOH_WKEND( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      !IF ( RC /= 0 ) CALL ALLOC_ERR( 'MOH_WKEND' )
      !MOH_WKEND = 0d0
      
      !ALLOCATE( CH4_WKEND( IIPAR, JJPAR, 3, 24 ), STAT=RC )
      !IF ( RC /= 0 ) CALL ALLOC_ERR( 'CH4_WKEND' )
      !CH4_WKEND = 0d0

      ! Return to calling program
      END SUBROUTINE INIT_NEI2008_ANTHRO
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_nei2008anthro
!
! !DESCRIPTION: Subroutine CLEANUP\_NEI2008\_ANTHRO deallocates all module 
!  arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_NEI2008_ANTHRO
!
! !REVISION HISTORY: 
!  01 Mar 2012 - R. Yantosca - Remove reference to A_CM2 array
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_NEIO2008_ANTHRO begins here!
      !=================================================================
      ! USA mask
      IF ( ALLOCATED( USA_MASK) ) DEALLOCATE( USA_MASK )
      IF ( ALLOCATED( TMP_WD ) ) DEALLOCATE( TMP_WD )
      IF ( ALLOCATED( TMP_WE ) ) DEALLOCATE( TMP_WE )
      IF ( ALLOCATED( TMPARR_WD ) ) DEALLOCATE( TMPARR_WD )
      IF ( ALLOCATED( TMPARR_WE ) ) DEALLOCATE( TMPARR_WE )
      !IF ( ALLOCATED( NOX    ) ) DEALLOCATE( NOX   )
      IF ( ALLOCATED( CO     ) ) DEALLOCATE( CO    )
      IF ( ALLOCATED( NO     ) ) DEALLOCATE( NO    )
      IF ( ALLOCATED( NO2    ) ) DEALLOCATE( NO2   )
      IF ( ALLOCATED( HNO2   ) ) DEALLOCATE( HNO2  )
      IF ( ALLOCATED( SO2    ) ) DEALLOCATE( SO2   )
      IF ( ALLOCATED( NH3    ) ) DEALLOCATE( NH3   )
      IF ( ALLOCATED( ALD2   ) ) DEALLOCATE( ALD2  )
      IF ( ALLOCATED( RCHO   ) ) DEALLOCATE( RCHO  )
      IF ( ALLOCATED( BENZ   ) ) DEALLOCATE( BENZ  )
      IF ( ALLOCATED( C2H6   ) ) DEALLOCATE( C2H6  )
      IF ( ALLOCATED( PRPE   ) ) DEALLOCATE( PRPE  )
      IF ( ALLOCATED( ALK4   ) ) DEALLOCATE( ALK4  )
      IF ( ALLOCATED( C3H8   ) ) DEALLOCATE( C3H8  )
      IF ( ALLOCATED( MEK    ) ) DEALLOCATE( MEK  )
      IF ( ALLOCATED( ACET   ) ) DEALLOCATE( ACET  )
      IF ( ALLOCATED( TOLU   ) ) DEALLOCATE( TOLU  )
      IF ( ALLOCATED( XYLE   ) ) DEALLOCATE( XYLE  )
      !IF ( ALLOCATED( C2H4   ) ) DEALLOCATE( C2H4  )
      IF ( ALLOCATED( CH2O   ) ) DEALLOCATE( CH2O  )
      IF ( ALLOCATED( BCPO   ) ) DEALLOCATE( BCPO  )
      IF ( ALLOCATED( OCPO   ) ) DEALLOCATE( OCPO  )
      IF ( ALLOCATED( SO4    ) ) DEALLOCATE( SO4   )
      !IF ( ALLOCATED( EOH    ) ) DEALLOCATE( EOH   )
      !IF ( ALLOCATED( MOH    ) ) DEALLOCATE( MOH   )
      !IF ( ALLOCATED( CH4    ) ) DEALLOCATE( CH4   )

      IF ( ALLOCATED( CO_WKEND    ) ) DEALLOCATE( CO_WKEND   )
      !IF ( ALLOCATED( NOX_WKEND   ) ) DEALLOCATE( NOX_WKEND  )
      IF ( ALLOCATED( NO_WKEND    ) ) DEALLOCATE( NO_WKEND   )
      IF ( ALLOCATED( NO2_WKEND   ) ) DEALLOCATE( NO2_WKEND  )
      IF ( ALLOCATED( HNO2_WKEND  ) ) DEALLOCATE( HNO2_WKEND )
      IF ( ALLOCATED( SO2_WKEND   ) ) DEALLOCATE( SO2_WKEND  )
      IF ( ALLOCATED( NH3_WKEND   ) ) DEALLOCATE( NH3_WKEND  )
      IF ( ALLOCATED( ALD2_WKEND  ) ) DEALLOCATE( ALD2_WKEND )
      IF ( ALLOCATED( RCHO_WKEND  ) ) DEALLOCATE( RCHO_WKEND )
      IF ( ALLOCATED( BENZ_WKEND  ) ) DEALLOCATE( BENZ_WKEND )
      IF ( ALLOCATED( C2H6_WKEND  ) ) DEALLOCATE( C2H6_WKEND )
      IF ( ALLOCATED( PRPE_WKEND  ) ) DEALLOCATE( PRPE_WKEND )
      IF ( ALLOCATED( ALK4_WKEND  ) ) DEALLOCATE( ALK4_WKEND )
      IF ( ALLOCATED( C3H8_WKEND  ) ) DEALLOCATE( C3H8_WKEND )
      IF ( ALLOCATED( MEK_WKEND   ) ) DEALLOCATE( MEK_WKEND  )
      IF ( ALLOCATED( ACET_WKEND  ) ) DEALLOCATE( ACET_WKEND )
      IF ( ALLOCATED( TOLU_WKEND  ) ) DEALLOCATE( TOLU_WKEND )
      IF ( ALLOCATED( XYLE_WKEND  ) ) DEALLOCATE( XYLE_WKEND )
      !IF ( ALLOCATED( C2H4_WKEND  ) ) DEALLOCATE( C2H4_WKEND )
      IF ( ALLOCATED( CH2O_WKEND  ) ) DEALLOCATE( CH2O_WKEND )
      IF ( ALLOCATED( BCPO_WKEND  ) ) DEALLOCATE( BCPO_WKEND )
      IF ( ALLOCATED( OCPO_WKEND  ) ) DEALLOCATE( OCPO_WKEND )
      IF ( ALLOCATED( SO4_WKEND   ) ) DEALLOCATE( SO4_WKEND  ) 
      !IF ( ALLOCATED( EOH_WKEND   ) ) DEALLOCATE( EOH_WKEND  )
      !IF ( ALLOCATED( MOH_WKEND   ) ) DEALLOCATE( MOH_WKEND  )
      !IF ( ALLOCATED( CH4_WKEND   ) ) DEALLOCATE( CH4_WKEND  )
      
      END SUBROUTINE CLEANUP_NEI2008_ANTHRO
!EOC
      END MODULE NEI2008_ANTHRO_MOD

!------------------------------------------------------------------------------
!          J. Mao (NOAA/GFDL) & J.A. Fisher (U. Wollongong)                   !
!------------------------------------------------------------------------------
!BOP
!     
! !MODULE: qfed_biomass_mod
!     
! !DESCRIPTION: Module QFED\_BIOMASS\_MOD contains routines and variables 
!  used to incorporate QFED emissions into GEOS-Chem
!\\   
!\\   
! !INTERFACE: 
!
      MODULE QFED_BIOMASS_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE

!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: QFED_COMPUTE_BIOMASS
      PUBLIC  :: CLEANUP_QFED_BIOMASS
      PUBLIC  :: QFED_IS_NEW
!
! PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: CHECK_QFED
      PRIVATE :: QFED_AVAILABLE
      PRIVATE :: READ_QFED
      PRIVATE :: QFED_SCALE_FUTURE
      PRIVATE :: QFED_TOTAL_Tg
      PRIVATE :: INIT_QFED_BIOMASS
      PRIVATE :: REARRANGE_BIOM
!
! !REMARKS:
!  Emissions of biomass burning species are read at daily
!  resolution. Note: no emission factors are used here - emissions of
!  individual species are given in input files. Emissions on the QFED 0.1x0.1 
!  degree grid are regridded to the current GEOS-Chem or GCAP grid.
!                                                                             .
!  QFED biomass burning emissions are computed for the following gas-phase 
!  and aerosol-phase species:
!                                                                             .
!     (1 ) NO   [  molec/cm2/s]     (13) BC   [atoms C/cm2/s]
!     (2 ) CO   [  molec/cm2/s]     (14) OC   [atoms C/cm2/s]                 
!     (3 ) ALK4 [atoms C/cm2/s]     (15) CO2  [  molec/cm2/s]
!     (4 ) ACET [atoms C/cm2/s]     (16) CH4  [  molec/cm2/s]
!     (5 ) MEK  [atoms C/cm2/s]     
!     (6 ) ALD2 [atoms C/cm2/s]     
!     (7 ) PRPE [atoms C/cm2/s]     
!     (8 ) C3H8 [atoms C/cm2/s]     
!     (9 ) CH2O [  molec/cm2/s]    
!     (10) C2H6 [atoms C/cm2/s]     
!     (11) SO2  [  molec/cm2/s]     
!     (12) NH3  [  molec/cm2/s]     
!                                                                             .
!  References:
!  ============================================================================
!  (1 ) QFED info and files:
!  Darmenov, A., and da Silva, A. M.: The Quick Fire Emissions Dataset (QFED) -
!  Documentation of versions 2.1, 2.2 and 2.4, NASA TM-2013-104606, Vol. 32, 
!  (http://gmao.gsfc.nasa.gov/pubs/tm/), 183 pp, 2013.
!  (2 ) You can get 0.1 deg gridded data here:
!   ftp://iesa@ftp.nccs.nasa.gov/aerosol/emissions/QFED/v2.4r6/0.1/
!
! !REVISION HISTORY: 
!  02 Sep 2014 - K. Travis - Initial version, based on GFED3
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      !=================================================================
      ! MODULE PARAMETERS
      !
      ! N_SPEC  : Number of species
      !=================================================================
      INTEGER,          PARAMETER   :: N_SPEC  = 16

!
! PRIVATE TYPES:
!
      !=================================================================
      ! MODULE VARIABLES:
      !
      ! Scalars
      !
      ! UPDATED         : flag to indicate if QFED emissions are updated
      ! UPDATED_MON     : flag to indicate if new month
      ! UPDATED_DAY     : flag to indicate if new day
      !                   - only set to true if daily emissions are used
      ! SECONDS         : Number of seconds in the current month
      !
      ! Arrays
      !
      ! QFED_SPEC_NAME  : Array for QFED biomass species names
      ! QFED_SPEC_MOLWT : Array for QFED biomass species molecular wts
      ! QFED_SPEC_UNIT  : Array for QFED biomass species emissions units
      ! BIOMASS_MODEL   : Array for QFED species emissions on model grid
      ! BIO_SAVE        : Index array to store IDBxxx values
      !=================================================================

      ! Scalars
      LOGICAL                       :: UPDATED
      LOGICAL                       :: UPDATED_MON
      LOGICAL                       :: UPDATED_DAY
      REAL*8                        :: SECONDS

      ! Arrays
      CHARACTER(LEN=7), ALLOCATABLE :: QFED_SPEC_NAME(:)
      REAL*8,           ALLOCATABLE :: QFED_SPEC_MOLWT(:)
      CHARACTER(LEN=6), ALLOCATABLE :: QFED_SPEC_UNIT(:)
      REAL*8,           ALLOCATABLE :: BIOMASS_MODEL(:,:,:)
      INTEGER,          ALLOCATABLE :: BIO_SAVE(:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: qfed_is_new
!
! !DESCRIPTION: Function QFED\_IS\_NEW returns TRUE if QFED emissions 
!  have been updated.
!\\
!\\
! !INTERFACE:
!
      FUNCTION QFED_IS_NEW( ) RESULT( IS_UPDATED )
!
! !RETURN VALUE:
!
      LOGICAL :: IS_UPDATED    ! =T if QFED is updated; =F otherwise
!
! !REMARKS:
!  Called from carbon_mod.f and sulfate_mod.f
!
! !REVISION HISTORY: 
!  02 Sep 2014 - K. Travis - Initial version, based on GFED3
!EOP
!------------------------------------------------------------------------------
!BOC
      IS_UPDATED = UPDATED      

      END FUNCTION QFED_IS_NEW
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_qfed
!
! !DESCRIPTION: Subroutine CHECK\_QFED checks if we entered a new QFED period
!  since last emission timestep (ie, last call). The result depends
!  on the emissions time step, and the QFED time period used, as well
!  as MMDD at beginning of the GEOS-Chem run
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHECK_QFED( DOY )
!
! !USES:
!
      USE TIME_MOD,    ONLY : ITS_A_NEW_MONTH
      USE TIME_MOD,    ONLY : ITS_A_NEW_DAY
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: DOY   ! Day of year (0-365 or 0-366 leap years)
!
! !REMARKS:
!  The routine computes the DOY at start of the 1-day period we are in, 
!  if the 1-day QFED. option is on. Result is compared to previous value
!  to indicate if new data should be read.
! 
! !REVISION HISTORY: 
!  02 Sep 2014 - K. Travis - Initial version, based on GFED3
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Reset to default
      UPDATED     = .FALSE.
      UPDATED_MON = .FALSE.
      UPDATED_DAY = .FALSE.

      ! Check if it is a new month
      IF ( ITS_A_NEW_MONTH() ) THEN
            UPDATED     = .TRUE.
            UPDATED_MON = .TRUE.
      ENDIF
      
      ! Check if it is a new day
      ! Use the daily emissions, hardcoded (skim, 7/27/13)
      IF ( ITS_A_NEW_DAY() ) THEN
         UPDATED     = .TRUE.
         UPDATED_DAY = .TRUE.
      ENDIF
      

      END SUBROUTINE CHECK_QFED
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: qfed_available
!
! !DESCRIPTION: Function QFED\_AVAILABLE checks an input YYYY year and MM
!  month against the available data dates.  If the requested YYYY and MM
!  lie outside of the valid range of dates, then QFED\_AVAILABLE will return 
!  the last valid YYYY and MM.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE QFED_AVAILABLE( YYYY, YMIN, YMAX, MM, MMIN, MMAX )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)              :: YMIN, YMAX   ! Min & max years 
      INTEGER, INTENT(IN),    OPTIONAL :: MMIN, MMAX   ! Min & max months
!
! !INPUT/OUTPUT PARAMETERS: 
!
      INTEGER, INTENT(INOUT)           :: YYYY         ! Year of QFED data
      INTEGER, INTENT(INOUT), OPTIONAL :: MM           ! Month of QFED data
! 
! !REVISION HISTORY: 
!  02 Sep 2014 - K. Travis - Initial version, based on GFED3
!EOP
!------------------------------------------------------------------------------
!BOC

      !!!!! Kludge - also have 2012 emissions at the moment (jaf)
      IF ( YYYY == 2012 ) RETURN

      !!!!! Kludge - also have 2013 emissions for the SEAC4RS mission
      IF ( YYYY == 2013 ) RETURN

      ! Check year
      IF ( YYYY > YMAX .OR. YYYY < YMIN ) THEN
         
         YYYY = MAX( YMIN, MIN( YYYY, YMAX) )
         
         WRITE( 6, 100 ) YMAX, YMIN, YYYY
 100     FORMAT( 'YEAR > ', i4, ' or YEAR < ', i4, 
     $        '. Using QFED biomass for ', i4)
      ENDIF
      

      ! Check month
      IF ( PRESENT( MM ) ) THEN 
         IF ( MM > MMAX .OR. MM < MMIN ) THEN

            MM = MAX( MMIN, MIN( MM, MMAX) )
            
            WRITE( 6, 200 ) MMIN, MMAX, MM
 200        FORMAT( ' ** WARNING ** : MONTH is not within ', i2,'-',
     $              i2, '. Using QFED biomass for month #', i2)
         ENDIF
      ENDIF

      END SUBROUTINE QFED_AVAILABLE
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: qfed_compute_biomass
!
! !DESCRIPTION: Subroutine QFED\_COMPUTE\_BIOMASS computes the daily
!  QFED biomass burning emissions for a given year and month. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE QFED_COMPUTE_BIOMASS( THIS_YYYY, THIS_MM, BIOM_OUT )
!
! !USES:
!
      USE BPCH2_MOD,       ONLY : GET_TAU0
      USE DIRECTORY_MOD,   ONLY : DATA_DIR_NATIVE => DATA_DIR_1x1
      USE JULDAY_MOD,      ONLY : JULDAY
      USE JULDAY_MOD,      ONLY : CALDATE
      USE LOGICAL_MOD,     ONLY : LFUTURE
      USE LOGICAL_MOD,     ONLY : LQFED, LFINN, LQFED_FINN
      USE TIME_MOD,        ONLY : EXPAND_DATE
      USE TIME_MOD,        ONLY : GET_DAY
      USE TIME_MOD,        ONLY : GET_DAY_OF_YEAR
      USE TIME_MOD,        ONLY : YMD_EXTRACT

      USE CMN_SIZE_MOD                                       ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)  :: THIS_YYYY                      ! Current year   
      INTEGER, INTENT(IN)  :: THIS_MM                        ! Current month
!
! !OUTPUT PARAMETERS:
!
      REAL*8,  INTENT(OUT) :: BIOM_OUT(IIPAR,JJPAR,NBIOMAX)  ! BB emissions
                                                             ! [molec/cm2/s]
!
! !REVISION HISTORY: 
!  02 Sep 2014 - K. Travis - Initial verson, based on GFED3
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE       :: FIRST = .TRUE.
      INTEGER             :: N
      INTEGER             :: YYYY, MM, MM1, YYYY1
      INTEGER             :: YYYYMMDD, YYYYMMDD1, HHMMSS
      REAL*8              :: ARRAY_QFED(IIPAR, JJPAR)
      REAL*8              :: TAU0, TAU1
      REAL*8              :: JD,   JD1
      INTEGER             :: DD,   DD1,  DOY
      REAL*4              :: TMP
      CHARACTER(LEN=255)  :: FILENAME,  QFED_DIR
      CHARACTER(LEN=7  )  :: TIME_STR
      CHARACTER(LEN=4  )  :: SPEC_NAME

      !=================================================================
      ! QFED_COMPUTE_BIOMASS begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_QFED_BIOMASS
         FIRST = .FALSE.

      ENDIF

      ! Save in local variables
      YYYY = THIS_YYYY
      MM   = THIS_MM
      DD   = GET_DAY()
      DOY  = GET_DAY_OF_YEAR()
      
      ! Check if we need to update QFED 
      CALL CHECK_QFED( DOY )
      
      ! If no updating is needed, module variable BIOMASS_MODEL
      ! from last update can be used
      IF ( .not. UPDATED ) THEN
         CALL REARRANGE_BIOM(BIOMASS_MODEL,BIOM_OUT)
         RETURN
      ENDIF

      ! Echo info
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, '(a)' ) 
     &  'Q F E D   B I O M A S S   B U R N I N G   E M I S S I O N S'

      !=================================================================
      ! Set-up for reading QFED daily emissions
      !=================================================================
      
      ! DAILY data
      !-------------------------------

      ! Check availability
      CALL QFED_AVAILABLE( YYYY, 2012, 2013 )

      ! Create YYYYMMDD integer values
      YYYYMMDD = YYYY*10000 + MM*100 + DD

      ! Use Julian dates to automatically increase year/month
      JD  = JULDAY( YYYY, MM, DBLE(DD) )
      JD1 = JD + 1d0

      CALL CALDATE( JD1, YYYYMMDD1, HHMMSS )
      CALL YMD_EXTRACT( YYYYMMDD1, YYYY1, MM1, DD1 )

      ! TAU value at start of this day
      TAU0     = GET_TAU0( MM, DD, YYYY )
      ! TAU value at start of next day
      TAU1     = GET_TAU0( MM1, DD1, YYYY1 )

      ! Directory with QFED emissions
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         QFED_DIR = '/as/data/geos/ExtData/HEMCO/QFED/v2014-09/'
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         TIME_STR =  'YYYY/MM'
      !=================================================================
      ! Calculate biomass species emissions on QFED grid
      ! and regrid to model grid
      !
      ! Emission are in [molec/cm2/s] and converted here to
      ! [molec/cm2/timeperiod] or [atoms C/cm2/timeperiod]
      !
      ! Units:
      !  [  molec/cm2/month] : NO,   CO,   CH2O, SO2,  NH3,  CO2
      !  [atoms C/cm2/month] : ALK4, ACET, MEK,  ALD2, PRPE, C3H8,
      !                        C2H6, BC,   OC, C2H6
      !=================================================================

      ! Loop over biomass species 
      DO N = 1, N_SPEC

         ! Initialize temporary array
         ARRAY_QFED = 0d0

         ! Set filename for this species
         SPEC_NAME = QFED_SPEC_NAME(N)
!
         FILENAME = TRIM( QFED_DIR ) //
     &        TRIM( TIME_STR ) // '/qfed2.emis_' // TRIM( SPEC_NAME ) //
     &              '.005.YYYYMMDD.nc4'

         CALL EXPAND_DATE( FILENAME, YYYYMMDD, 000000 )

         ! Read QFED emissions for species in kg/m2/2 and convert to molecules/cm2/s
         CALL READ_QFED(FILENAME,'biomass', ARRAY_QFED, YYYY, MM, DD, N)

         ! Add to biomass array
         BIOMASS_MODEL(:,:,N) = ARRAY_QFED

      ENDDO

      ! Compute future biomass emissions (if necessary)
      IF ( LFUTURE ) THEN
         CALL QFED_SCALE_FUTURE( BIOMASS_MODEL )
      ENDIF

      ! Print totals in Tg/time period
      WRITE( 6, 411 ) YYYY, MM, DD
 411  FORMAT( 'QFED daily emissions for year, month, day: ',
     &     i4, '/', i2.2, '/', i2.2, / )
      
      ! Number of seconds in this period (day or month)
      ! (NOTE: its value will be saved until the next call)
      SECONDS  = ( TAU1 - TAU0 ) * 3600d0

      CALL QFED_TOTAL_Tg( SECONDS )
      
      ! Rearrange the species to the same order as in the IDBxxx (fp, 6/09)
      ! BIOMASS_MODEL is indexed as QFED
      ! BIOM_OUT      is indexed as IDBs
      CALL REARRANGE_BIOM( BIOMASS_MODEL, BIOM_OUT )

      ! Echo info
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      END SUBROUTINE QFED_COMPUTE_BIOMASS
!EOC
!------------------------------------------------------------------------------
!         J. Fisher - U. Wollongong                                           !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_qfed
!
! !DESCRIPTION: Subroutine READ\_QFED reads the QFED BB emissions at 0.1x0.1
!  degree resolution and regrids them to the current model resolution
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_QFED( FILENAME, PARNAME, ARR, YYYY, MM, DD, N )
!
! !USES:
!
      USE DIRECTORY_MOD,  ONLY: DATA_DIR_NATIVE => DATA_DIR_1x1
      USE NCDF_MOD,       ONLY: NC_READ
      USE REGRID_A2A_MOD, ONLY: DO_REGRID_A2A
      USE ERROR_MOD,      ONLY: IT_IS_NAN

      USE CMN_SIZE_MOD          ! Size parameters

!
! !INPUT PARAMETERS:
!
      ! File and parameter to read
      CHARACTER(LEN=255),  INTENT(IN) :: FILENAME
      CHARACTER(LEN=7),    INTENT(IN) :: PARNAME
      INTEGER,             INTENT(IN) :: YYYY
      INTEGER,             INTENT(IN) :: MM
      INTEGER,             INTENT(IN) :: DD
      INTEGER,             INTENT(IN) :: N
! 
! !OUTPUT PARAMETERS:
!
      ! Emissions output
      REAL*8,            INTENT(OUT) :: ARR(IIPAR,JJPAR)
! 
! !REVISION HISTORY:
!  02 Jan 2013 - J. Fisher - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER, PARAMETER    :: I_NATIVE = 3600 ! 0.1 x 0.1 resolution
      INTEGER, PARAMETER    :: J_NATIVE = 1800 ! 0.1 x 0.1 resolution
      REAL*8                :: GENERIC_NATIVE(I_NATIVE,J_NATIVE,1)
      CHARACTER(LEN=255)    :: LLFILENAME
      REAL*8                :: INGRID(I_NATIVE,J_NATIVE)
      ! for ncdf-reading
      REAL*4, POINTER       :: NCARR(:,:,:) => NULL()
      INTEGER               :: I, J

      !=================================================================
      ! READ_QFED begins here!
      !=================================================================

      ! Echo info
      WRITE( 6, 100 ) TRIM( QFED_SPEC_NAME(N)), TRIM( FILENAME )
 100  FORMAT( '     - READ_QFED : Reading ', a , ' from ', a ) 

      ARR   = 0D0 

      CALL NC_READ ( NC_PATH = TRIM(FILENAME) ,
     &               PARA    = TRIM(PARNAME)  ,
     &               ARRAY   = NCARR, 
     &               VERBOSE = .FALSE.       )

      ! Cast to REAL*8 before regridding
      GENERIC_NATIVE(:,:,1) = NCARR(:,:,1)
      
      ! Convert from kg/s/m2 to molecules/cm2/s
      GENERIC_NATIVE = GENERIC_NATIVE * 6.022d23 / QFED_SPEC_MOLWT(N)
     &     / 1.0d4

      ! Deallocate ncdf-array
      IF ( ASSOCIATED ( NCARR ) ) DEALLOCATE ( NCARR )

!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, J )
      DO J = 1, J_NATIVE
      DO I = 1, I_NATIVE
         IF (IT_IS_NAN(GENERIC_NATIVE(I,J,1))) THEN
            GENERIC_NATIVE(I,J,1) = 0d0
            !WRITE(*,*) 'NAN', I, J
         ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! File with lat/lon edges for regridding
      LLFILENAME = TRIM( DATA_DIR_NATIVE) //
     &             'MAP_A2A_Regrid_201203/' //
     &             'MAP_A2A_latlon_generic01x01.nc'

      ! Regrid to current model resolution
      INGRID = GENERIC_NATIVE(:,:,1)
      CALL DO_REGRID_A2A( LLFILENAME, I_NATIVE, J_NATIVE,
     &                    INGRID,     ARR,      IS_MASS=0,
     &                    netCDF=.TRUE.                    )

      ! Return to calling program
      END SUBROUTINE READ_QFED
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: qfed_scale_future
!
! !DESCRIPTION: Subroutine QFED\_SCALE\_FUTURE applies the IPCC future 
!  emissions scale factors to the QFED biomass burning emisisons in order 
!  to compute the future emissions of biomass burning for NOx, CO, and VOC's.  
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE QFED_SCALE_FUTURE( BB )
!
! !USES:
!
      USE FUTURE_EMISSIONS_MOD,   ONLY : GET_FUTURE_SCALE_BCbb
      USE FUTURE_EMISSIONS_MOD,   ONLY : GET_FUTURE_SCALE_CObb
      USE FUTURE_EMISSIONS_MOD,   ONLY : GET_FUTURE_SCALE_NH3bb
      USE FUTURE_EMISSIONS_MOD,   ONLY : GET_FUTURE_SCALE_NOxbb
      USE FUTURE_EMISSIONS_MOD,   ONLY : GET_FUTURE_SCALE_OCbb
      USE FUTURE_EMISSIONS_MOD,   ONLY : GET_FUTURE_SCALE_SO2bb
      USE FUTURE_EMISSIONS_MOD,   ONLY : GET_FUTURE_SCALE_VOCbb
      USE TRACER_MOD,             ONLY : ITS_A_CO2_SIM 
      USE TRACER_MOD,             ONLY : ITS_A_CH4_SIM
      USE TRACERID_MOD,           ONLY : IDBNO,   IDBCO,   IDBSO2 
      USE TRACERID_MOD,           ONLY : IDBNH3,  IDBBC,   IDBOC 

      USE CMN_SIZE_MOD                 ! Size parameters
!
! !OUTPUT PARAMETERS:
!
      !  Array w/ biomass burning emisisons [molec/cm2]
      REAL*8, INTENT(INOUT) :: BB(IIPAR,JJPAR,N_SPEC)
! 
! !REVISION HISTORY: 
!  02 Jan 2013 - J. Fisher - Initial version, based on GFED3
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL :: ITS_CO2
      LOGICAL :: ITS_CH4
      INTEGER :: I, J, N
      
      !=================================================================
      ! QFED_SCALE_FUTURE begins here!
      !=================================================================

      ! Test if it's a CO2 simulation outside of the loop
      ITS_CO2 = ITS_A_CO2_SIM()
      ITS_CH4 = ITS_A_CH4_SIM()

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, N )

      ! Loop over species and grid boxes
      DO N = 1, N_SPEC
      DO J = 1, JJPAR
      DO I = 1, IIPAR 

         ! Scale each species to IPCC future scenario
         IF ( BIO_SAVE(N) == IDBNO ) THEN

            ! Future biomass NOx [molec/cm2]
            BB(I,J,N) = BB(I,J,N) * GET_FUTURE_SCALE_NOxbb( I, J )

         ELSE IF ( BIO_SAVE(N) == IDBCO ) THEN

            ! Future biomass CO [molec/cm2]
            BB(I,J,N) = BB(I,J,N) * GET_FUTURE_SCALE_CObb( I, J )

         ELSE IF ( BIO_SAVE(N) == IDBSO2 ) THEN

            ! Future biomass SO2 [molec/cm2]
            BB(I,J,N) = BB(I,J,N) * GET_FUTURE_SCALE_SO2bb( I, J )

         ELSE IF ( BIO_SAVE(N) == IDBNH3 ) THEN

            ! Future biomass NH3 [molec/cm2]
            BB(I,J,N) = BB(I,J,N) * GET_FUTURE_SCALE_NH3bb( I, J )

         ELSE IF ( BIO_SAVE(N) == IDBBC ) THEN

            ! Future biomass BC [molec/cm2]
            BB(I,J,N) = BB(I,J,N) * GET_FUTURE_SCALE_BCbb( I, J )

         ELSE IF ( BIO_SAVE(N) == IDBOC ) THEN

            ! Future biomass OC [molec/cm2]
            BB(I,J,N) = BB(I,J,N) * GET_FUTURE_SCALE_OCbb( I, J )

         ! Don't scale future emissions if CO2 or CH4
         ELSE IF ( ITS_CO2 .OR. ITS_CH4 ) THEN

            ! Nothing

         ELSE

            ! Future biomass Hydrocarbons [atoms C/cm2]
            BB(I,J,N) = BB(I,J,N) * GET_FUTURE_SCALE_VOCbb( I, J )

         ENDIF
         
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE QFED_SCALE_FUTURE
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: qfed_total_Tg
!
! !DESCRIPTION: Subroutine QFED\_TOTAL\_Tg prints the amount of biomass 
!  burning emissions that are emitted each month/day in Tg or Tg C. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE QFED_TOTAL_Tg( SECONDS )
!
! !USES:
!
      USE GRID_MOD,    ONLY : GET_AREA_CM2

      USE CMN_SIZE_MOD      ! Size parameters
!
! !INPUT PARAMETERS:
!
      ! Seconds in this time period
      REAL*8,  INTENT(IN)  :: SECONDS
!
! !REVISION HISTORY: 
!  02 Sep 2014 - K. Travis - Initial version, based on GFED3
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER          :: I,    J,     N
      REAL*8           :: CONV, MOLWT, TOTAL
      CHARACTER(LEN=4) :: NAME
      CHARACTER(LEN=6) :: UNIT

      !=================================================================
      ! QFED_TOTAL_Tg begins here!
      !=================================================================

      ! Loop over biomass species
      DO N = 1, N_SPEC

         ! Initialize
         NAME  = QFED_SPEC_NAME(N)
         MOLWT = QFED_SPEC_MOLWT(N)
         UNIT  = QFED_SPEC_UNIT(N)
         TOTAL = 0d0

         ! Loop over surface boxes
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Convert to [Tg/qfed-period] (or [Tg C/qfed-period] for HC's)
            CONV  = GET_AREA_CM2( I, J, 1 ) 
     &            * ( MOLWT / 6.023d23 ) * 1d-9 * SECONDS

            TOTAL = TOTAL + ( BIOMASS_MODEL(I,J,N) * CONV )

         ENDDO
         ENDDO
     
         ! Write totals
         WRITE( 6, 110 ) NAME, TOTAL, UNIT
 110     FORMAT( 'Sum Biomass ', a4, 1x, ': ', e12.5, 1x, a6 )
      ENDDO

      END SUBROUTINE QFED_TOTAL_Tg
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_qfed_biomass
!
! !DESCRIPTION: Subroutine INIT\_QFED\_BIOMASS allocates all module arrays.  
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_QFED_BIOMASS
!
! !USES:
!
      USE DIRECTORY_MOD,   ONLY : DATA_DIR_NATIVE => DATA_DIR_1x1
      USE ERROR_MOD,       ONLY : ALLOC_ERR
      USE FILE_MOD,        ONLY : IOERROR
      USE LOGICAL_MOD,     ONLY : LDICARB
      USE TRACERID_MOD,    ONLY : IDBNO,   IDBCO,   IDBALK4
      USE TRACERID_MOD,    ONLY : IDBACET, IDBMEK,  IDBALD2
      USE TRACERID_MOD,    ONLY : IDBPRPE, IDBC3H8, IDBCH2O
      USE TRACERID_MOD,    ONLY : IDBC2H6, IDBBC,   IDBOC
      USE TRACERID_MOD,    ONLY : IDBSO2,  IDBNH3,  IDBCO2
      USE TRACERID_MOD,    ONLY : IDBCH4
      USE CMN_SIZE_MOD          ! Size parameters
! 
! !REVISION HISTORY: 
!  02 Sep 2014 - K. Travis - Initial version, based on GFED3
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: IU_FILE
      INTEGER            :: AS, IOS, M, N, NDUM
      CHARACTER(LEN=255) :: FILENAME
      
      !=================================================================
      ! INIT_QFED_BIOMASS begins here!
      !=================================================================

      ! Allocate array to hold QFED species emissions on model grid
      ALLOCATE( BIOMASS_MODEL( IIPAR, JJPAR, N_SPEC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOMASS_MODEL' )
      BIOMASS_MODEL = 0d0

      ! Allocate array for species molecular weight
      ALLOCATE( QFED_SPEC_MOLWT( N_SPEC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'QFED_SPEC_MOLWT' )
      QFED_SPEC_MOLWT = 0d0

      ! Allocate array for species name
      ALLOCATE( QFED_SPEC_NAME( N_SPEC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'QFED_SPEC_NAME' )
      QFED_SPEC_NAME = ''

      ! Allocate array for QFED biomass buning species mass units
      ALLOCATE( QFED_SPEC_UNIT( N_SPEC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'QFED_SPEC_UNIT' )
      QFED_SPEC_UNIT = ''

      !IDBs are now the same as the ones in TRACERID AND BIOMASS_MOD
      !BIOSAVE INDEX IS THE LOCATION OF THE EMISSION IN THE QFED FILE
      !(fp)
      ALLOCATE( BIO_SAVE( N_SPEC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIO_SAVE' )
      BIO_SAVE = 0

      !=================================================================
      ! Define QFED species names
      !=================================================================
      QFED_SPEC_NAME(1)  = 'no'
      QFED_SPEC_NAME(2)  = 'co'
      QFED_SPEC_NAME(3)  = 'alk4'
      QFED_SPEC_NAME(4)  = 'acet'
      QFED_SPEC_NAME(5)  = 'mek'
      QFED_SPEC_NAME(6)  = 'ald2'
      QFED_SPEC_NAME(7)  = 'c3h6'
      QFED_SPEC_NAME(8)  = 'c3h8'
      QFED_SPEC_NAME(9)  = 'ch2o'
      QFED_SPEC_NAME(10) = 'c2h6'
      QFED_SPEC_NAME(11) = 'bc'
      QFED_SPEC_NAME(12) = 'oc'
      QFED_SPEC_NAME(13) = 'so2'
      QFED_SPEC_NAME(14) = 'nh3'
      QFED_SPEC_NAME(15) = 'co2'
      QFED_SPEC_NAME(16) = 'ch4'

      !=================================================================
      ! Define local ID flags and arrays for the names, units, 
      ! and molecular weights of the QFED biomass species
      !=================================================================
      ! Save correspondance between QFED species order (N) and 
      ! species order of the simulation (IDBxxxs).(ccc, 2/4/10)
      ! and also initialize arrays for mol wts and units
      DO N = 1, N_SPEC
         SELECT CASE ( TRIM( QFED_SPEC_NAME(N) ) ) 
            CASE( 'no'  )
               BIO_SAVE(N)         = IDBNO
               QFED_SPEC_MOLWT(N) = 14d-3
               QFED_SPEC_UNIT(N)  = '[Tg N]'
            CASE( 'co'   )
               BIO_SAVE(N)         = IDBCO
               QFED_SPEC_MOLWT(N) = 28d-3
               QFED_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'alk4' )
               BIO_SAVE(N)         = IDBALK4
               QFED_SPEC_MOLWT(N) = 12d-3
               QFED_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'acet' )
               BIO_SAVE(N)         = IDBACET
               QFED_SPEC_MOLWT(N) = 12d-3
               QFED_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'mek'  )
               BIO_SAVE(N)         = IDBMEK
               QFED_SPEC_MOLWT(N) = 12d-3
               QFED_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'ald2' )
               BIO_SAVE(N)         = IDBALD2
               QFED_SPEC_MOLWT(N) = 12d-3
               QFED_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'c3h6' )
               BIO_SAVE(N)         = IDBPRPE
               QFED_SPEC_MOLWT(N) = 12d-3
               QFED_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'c3h8' )
               BIO_SAVE(N)         = IDBC3H8
               QFED_SPEC_MOLWT(N) = 12d-3
               QFED_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'ch2o' )
               BIO_SAVE(N)         = IDBCH2O
               QFED_SPEC_MOLWT(N) = 30d-3
               QFED_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'c2h6' )
               BIO_SAVE(N)         = IDBC2H6
               QFED_SPEC_MOLWT(N) = 12d-3
               QFED_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'so2'  )
               BIO_SAVE(N)         = IDBSO2
               QFED_SPEC_MOLWT(N) = 64d-3
               QFED_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'nh3'  )
               BIO_SAVE(N)         = IDBNH3
               QFED_SPEC_MOLWT(N) = 17d-3
               QFED_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'bc'   )
               BIO_SAVE(N)         = IDBBC
               QFED_SPEC_MOLWT(N) = 12d-3
               QFED_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'oc'   )
               BIO_SAVE(N)         = IDBOC
               QFED_SPEC_MOLWT(N) = 12d-3
               QFED_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'co2'  )
               BIO_SAVE(N)         = IDBCO2
               QFED_SPEC_MOLWT(N) = 44d-3
               QFED_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'ch4' )
               BIO_SAVE(N)         = IDBCH4
               QFED_SPEC_MOLWT(N) = 16d-3
               QFED_SPEC_UNIT(N)  = '[Tg  ]'
            CASE DEFAULT
               ! Nothing
               BIO_SAVE(N)         = 0

            WRITE(*,*) 'NAME',TRIM( QFED_SPEC_NAME(N) )
         END SELECT
      ENDDO

      END SUBROUTINE INIT_QFED_BIOMASS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rearrange_biom
!
! !DESCRIPTION: Subroutine REARRANGE\_BIOM takes QFED emissions (which have 
!  their own, unique ID numbers and associates them with the IDBxxxs of 
!  tracerid\_mod.F.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE REARRANGE_BIOM( BIOM_OUT, BIOM_OUTM )

!
! !USES:
!
      USE CMN_SIZE_MOD     ! Size parameters
!
! !INPUT PARAMETERS: 
!
      REAL*8, INTENT(IN)  :: BIOM_OUT (IIPAR,JJPAR,N_SPEC)
!
! !OUTPUT PARAMETERS:
!
      REAL*8, INTENT(OUT) :: BIOM_OUTM(IIPAR,JJPAR,NBIOMAX)
! 
! !REVISION HISTORY: 
!  02 Sep 2014 - K. Travis - Initial version, based on GFED3
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: N

      ! Loop over QFED species
      DO N = 1, N_SPEC

         ! Save into array w/ proper ordering for GEOS-Chem
         IF ( BIO_SAVE(N) .GT. 0 ) THEN
            BIOM_OUTM(:,:,BIO_SAVE(N)) = BIOM_OUT(:,:,N)
         ENDIF

      ENDDO

      END SUBROUTINE REARRANGE_BIOM
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_qfed_biomass
!
! !DESCRIPTION: Subroutine CLEANUP\_QFED\_BIOMASS deallocates all module 
!  arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_QFED_BIOMASS
! 
! !REVISION HISTORY: 
!  02 Sep 2014 - K. Travis - Initial version, based on GFED3
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_QFED_BIOMASS begins here!
      !=================================================================
      IF ( ALLOCATED( QFED_SPEC_MOLWT ) ) DEALLOCATE( QFED_SPEC_MOLWT )
      IF ( ALLOCATED( QFED_SPEC_NAME  ) ) DEALLOCATE( QFED_SPEC_NAME  )
      IF ( ALLOCATED( BIOMASS_MODEL   ) ) DEALLOCATE( BIOMASS_MODEL   )
      
      END SUBROUTINE CLEANUP_QFED_BIOMASS
!EOC

      END MODULE QFED_BIOMASS_MOD

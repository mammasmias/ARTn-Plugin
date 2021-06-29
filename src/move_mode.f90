SUBROUTINE move_mode(nat, dlanc, force, &
                     vel, alpha_init, dt, &
                     iperp, push, &
                     mode, prfx, tmpdir )
  !
  ! translate specified move to appropriate force and set FIRE parameters accordingly  
  !
  USE artn_params, ONLY: DP, AMU_RY 
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN)                       :: nat
  REAL(DP), INTENT(IN)                      :: dlanc  
  ! REAL(DP), DIMENSION(3,nat), INTENT(IN) :: v1
  REAL(DP), DIMENSION(3,nat), INTENT(INOUT) :: force
  REAL(DP), DIMENSION(3,nat), INTENT(INOUT) :: vel
  REAL(DP), INTENT(IN)                      :: alpha_init, dt
  INTEGER, INTENT(IN)                       :: iperp
  REAL(DP), DIMENSION(3,nat), INTENT(IN)    :: push
  CHARACTER(LEN=4), INTENT(IN)              :: mode
  CHARACTER(LEN=255), INTENT(IN)            :: tmpdir, prfx
  ! 
  REAL(DP), EXTERNAL :: ddot,dnrm2
  ! variables read from the FIRE minimization algorithm
  INTEGER :: nsteppos
  INTEGER :: ios
  REAL(DP) :: dt_curr, alpha, etot
  LOGICAL :: file_exists
  CHARACTER(len=256) :: filnam
  ! etot is set to zero to ensure that a value is written down in the initial push 
  etot = 0.D0
  !
  ! do things depending on mode of the move
  ! NOTE force units of Ry/a.u. are assumed ... 
  !
  filnam = trim(tmpdir) // '/' // trim(prfx) // '.' //'fire'
  INQUIRE( file = filnam, exist = file_exists )
  OPEN( unit = 4, file = filnam, form = 'formatted', status = 'unknown', iostat = ios)
  !
  IF (file_exists ) THEN
     ! if file exists read the data, otherwise just close it 
     READ( UNIT = 4, FMT = * ) etot, nsteppos, dt_curr, alpha
     CLOSE( UNIT = 4, STATUS = 'KEEP' )
  ELSE
     CLOSE( UNIT = 4, STATUS = 'DELETE')
  ENDIF


  !
  ! WARNING
  ! iperp is incremented after the call move_mode so should be init at 0
  ! But now no...
  !



  SELECT CASE( TRIM(mode) )

  CASE( 'perp' )
     !
     IF( iperp .eq. 0 ) THEN
        ! for the first step forget previous velocity (prevent P < 0)
        etot = 0.D0
        vel(:,:) = 0.D0
        alpha = alpha_init
        dt_curr = dt
     ELSE
        ! subtract the components that are parallel
        vel(:,:) = vel(:,:) - ddot(3*nat,vel, 1, push, 1 )*push(:,:)/ddot(3*nat,push(:,:),1, push(:,:),1)
     ENDIF
        !
  CASE( 'lanc' )
     !
     ! set the velocity and acceleration and alpha of previous step to move correctly
     !
     etot = 0.D0
     vel(:,:) = 0.D0
     dt_curr = dt
     alpha = 0.D0
     nsteppos = 0
     ! the step performed should be like this now translate it into the correct force
     force(:,:) = force(:,:)*dlanc*amu_ry/dt_curr**2
     !
  CASE( 'eign' )
     !
     etot = 0.D0
     vel(:,:) = 0.D0
     alpha = 0.0_DP
     dt_curr = dt
     nsteppos = 0
     force(:,:) = force(:,:)*amu_ry/dt_curr**2
     !
  CASE default
     write(*,*) 'Problem with move_mode!'
  END SELECT
  !
  ! write the FIRE parameters to its scratch file
  ! 
  OPEN( unit = 4, file = filnam, form = 'formatted', status = 'unknown', iostat = ios)
  WRITE( UNIT = 4, FMT = * ) etot, nsteppos, dt_curr, alpha
  !
  CLOSE( UNIT = 4, STATUS = 'KEEP' )
  !
  !
 

END SUBROUTINE move_mode




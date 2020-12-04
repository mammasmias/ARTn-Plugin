SUBROUTINE move_mode(nat, alat ,dlanc, v1, force, &
                     vel, acc, alpha_init, dt, &
                     istepperp, push, &
                     mode, prfx, tmpdir )
  !
  ! translate specified move to appropriate force and set FIRE parameters accordingly  
  !
  USE artn_params, ONLY: DP, AMU_RY 
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nat
  REAL(DP), INTENT(IN) :: dlanc, alat 
  REAL(DP), DIMENSION(3,nat), INTENT(IN) :: v1
  REAL(DP), DIMENSION(3,nat), INTENT(INOUT) :: force
  REAL(DP), DIMENSION(3,nat), INTENT(INOUT) :: vel, acc
  REAL(DP), INTENT(IN) :: alpha_init, dt
  INTEGER, INTENT(IN) :: istepperp
  REAL(DP), DIMENSION(3,nat), INTENT(IN) :: push
  CHARACTER(LEN=4), INTENT(IN) :: mode
  CHARACTER(LEN=255), INTENT(IN) :: tmpdir, prfx
  ! 
  REAL(DP), EXTERNAL :: ddot,dnrm2
  ! variables read from the FIRE minimization algorithm
  INTEGER :: nsteppos
  INTEGER :: ios
  REAL(DP) :: dt_curr, alpha, etot
  LOGICAL :: file_exists
  CHARACTER(len=256) :: filnam

  !
  ! write(*,*) 'in MOVE:'
  ! write(*,*) 'mode:',trim(mode)
  ! write(*,*) 'nat:',nat
  ! write(*,*) 'dlanc',dlanc
  ! write(*,*) 'v1',v1
  ! write(*,*) 'force', force
  ! write(*,*) 'istepperp',istepperp
  ! write(*,*) 'push',push
  ! write(*,*) 'file exists', file_exists
  !

  !
  ! do things depending on mode of the move
  !
  ! CALL seqopn( 4, 'fire', 'FORMATTED', file_exists )

  filnam = trim(tmpdir) // '/' // trim(prfx) // '.' //'fire'
  INQUIRE( file = filnam, exist = file_exists )
  OPEN( unit = 4, file = filnam, form = 'formatted', status = 'unknown', iostat = ios)
  ! write(999,*) 'move_mode:', trim(mode)
  ! write(999,*) 'open:', trim(filnam), ios, file_exists
  ! flush(999)

  !
  IF (file_exists ) THEN
     !
     READ( UNIT = 4, FMT = * ) etot, nsteppos, dt_curr, alpha
     CLOSE( UNIT = 4, STATUS = 'KEEP' )
     !
     SELECT CASE( TRIM(mode) )

     CASE( 'perp' )
        ! write(*,*) 'MOVE: perp'
        !
        IF( istepperp .eq. 0 ) THEN
           ! for the first step forget previous acc and velocity (prevent P < 0)
           acc(:,:) = 0.D0
           vel(:,:) = 0.D0
           alpha = alpha_init
           dt_curr = dt
        ELSE
           ! subtract the components that are parallel
           acc(:,:) = acc(:,:) - ddot(3*nat,acc, 1, push, 1 )*push(:,:)/ddot(3*nat,push(:,:),1, push(:,:),1)
           vel(:,:) = vel(:,:) - ddot(3*nat,vel, 1, push, 1 )*push(:,:)/ddot(3*nat,push(:,:),1, push(:,:),1)
        ENDIF
        !
     CASE( 'lanc' )
        ! write(*,*) 'MOVE: lanc'
        !
        ! set the velocity and acceleration and alpha of previous step to move correctly
        !
        vel(:,:) = 0.D0
        acc(:,:) = 0.D0
        dt_curr = dt
        alpha = 0.D0
        nsteppos = 0
        ! the step performed should be like this now translate it into the correct force
        force(:,:) = v1(:,:)*dlanc*alat*amu_ry/dt_curr**2
        ! WRITE (*,*) "Expected lanczos step:", v1(:,:)*dlanc
        !
     CASE( 'eign' )
        ! write(*,*) 'MOVE: eign'
        !
        ! for the first step forget previous acc and velocity (prevent P < 0)
        acc(:,:) = 0.D0
        vel(:,:) = 0.D0
        alpha = 0.0_DP
        force(:,:) = force(:,:)*alat*amu_ry/dt_curr**2
        !
     CASE default
        write(*,*) 'Problem with move_mode!'
     END SELECT
     !
     !
     ! write the FIRE parameters to its scratch file
     ! CALL seqopn( 4, 'fire', 'FORMATTED', file_exists )
     OPEN( unit = 4, file = filnam, form = 'formatted', status = 'unknown', iostat = ios)
     WRITE( UNIT = 4, FMT = * ) etot, nsteppos, dt_curr, alpha
     !
     CLOSE( UNIT = 4, STATUS = 'KEEP' )
     !
     !
  ELSE
     CLOSE( UNIT = 4, STATUS = 'DELETE')
  ENDIF

END SUBROUTINE move_mode




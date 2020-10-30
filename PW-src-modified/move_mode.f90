SUBROUTINE move_mode(nat, dlanc, v1, force, &
                     vel, acc, alpha_init, dt, & 
                     istepperp, push, &
                     mode)
  !
  ! unify the move routines
  !
  USE kinds, ONLY: DP
  USE cell_base, ONLY: alat
  USE constants, ONLY: amu_ry
  USE io_files, ONLY: prefix, seqopn, tmp_dir
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nat
  REAL(DP), INTENT(IN) :: dlanc
  REAL(DP), DIMENSION(3,nat), INTENT(IN) :: v1
  REAL(DP), DIMENSION(3,nat), INTENT(INOUT) :: force
  REAL(DP), DIMENSION(3,nat), INTENT(INOUT) :: vel, acc
  REAL(DP), INTENT(IN) :: alpha_init, dt  
  INTEGER, INTENT(IN) :: istepperp
  REAL(DP), DIMENSION(3,nat), INTENT(IN) :: push
  CHARACTER(LEN=4), INTENT(IN) :: mode
  REAL(DP), EXTERNAL :: ddot,dnrm2
  ! variables read from the FIRE minimization algorithm
  INTEGER :: nsteppos
  REAL(DP) :: dt_curr, alpha, etot
  LOGICAL :: file_exists
  
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
  CALL seqopn( 4, 'fire', 'FORMATTED', file_exists )
  IF (file_exists ) THEN
     READ( UNIT = 4, FMT = * ) etot, nsteppos, dt_curr, alpha
     CLOSE( UNIT = 4, STATUS = 'KEEP' )
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
     CASE( 'hess' )
        ! write(*,*) 'MOVE: hess'
        !
        acc(:,:) = 0.D0
        vel(:,:) = 0.D0
        alpha = 0.D0
        dt_curr = 0.D0
     !
     CASE default
        write(*,*) 'Problem with move_mode!'
     END SELECT
     !
     !
     ! write the FIRE parameters to its scratch file 
     CALL seqopn( 4, 'fire', 'FORMATTED', file_exists )
     WRITE( UNIT = 4, FMT = * ) etot, nsteppos, dt_curr, alpha  
     !
     CLOSE( UNIT = 4, STATUS = 'KEEP' )
     !
     !
  ELSE
     CLOSE( UNIT = 4, STATUS = 'DELETE') 
  ENDIF

END SUBROUTINE move_mode




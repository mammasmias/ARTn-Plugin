SUBROUTINE move_mode(nat, dlanc, v1, force, &
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
  INTEGER, INTENT(IN) :: istepperp
  REAL(DP), DIMENSION(3,nat), INTENT(IN) :: push
  CHARACTER(LEN=4), INTENT(IN) :: mode
  ! variables read from the FIRE minimization algorithm
  INTEGER :: istep, nsteppos
  REAL(DP) :: etot, dt_init, dt_curr, alpha, alpha_init
  REAL(DP), EXTERNAL :: ddot,dnrm2
  REAL(DP), DIMENSION(3,nat) :: tau, vel, acc
  REAL(DP), DIMENSION(3,nat) :: vel_old, acc_old
  LOGICAL :: file_exists
  !
  ! Open file written by fire
  !
  CALL seqopn( 4, 'fire', 'FORMATTED', file_exists )
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
  IF( file_exists ) THEN
     !
     ! read the file
     !
     READ( UNIT = 4, FMT = * ) etot, istep, nsteppos, dt_init, dt_curr, &
          alpha, alpha_init, tau(:,:), vel_old(:,:), acc_old(:,:)
     CLOSE( UNIT = 4, STATUS = 'KEEP' )
     !
     ! do things depending on mode of the move
     !
     SELECT CASE( TRIM(mode) )

     CASE( 'perp' )
        ! write(*,*) 'MOVE: perp'
        !
        IF( istepperp .eq. 0 ) THEN
           ! for the first step forget previous acc and velocity (prevent P < 0)
           acc_old(:,:) = 0.D0
           vel_old(:,:) = 0.D0
           alpha = alpha_init
           dt_curr = dt_init
        ELSE
           ! subtract the components that are parallel
           acc_old(:,:) = acc_old(:,:) - ddot(3*nat,acc_old, 1, push, 1 )*push(:,:)/ddot(3*nat,push(:,:),1, push(:,:),1)
           vel_old(:,:) = vel_old(:,:) - ddot(3*nat,vel_old, 1, push, 1 )*push(:,:)/ddot(3*nat,push(:,:),1, push(:,:),1)
        ENDIF
        !
     CASE( 'lanc' )
        ! write(*,*) 'MOVE: lanc'
        !
        ! set the velocity and acceleration and alpha of previous step to move correctly
        !
        vel_old(:,:) = 0.D0
        acc_old(:,:) = 0.D0
        dt_curr = dt_init
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
        acc_old(:,:) = 0.D0
        vel_old(:,:) = 0.D0
        alpha = 0.1_DP
        force(:,:) = force(:,:)*alat*amu_ry/dt_curr**2
        !
     CASE( 'hess' )
        ! write(*,*) 'MOVE: hess'
        !
        acc_old(:,:) = 0.D0
        vel_old(:,:) = 0.D0
        alpha = 0.D0
        dt_curr = 0.D0
        !
     CASE default
        write(*,*) 'Problem with move_mode!'

     END SELECT
     !
     CALL seqopn( 4, 'fire', 'FORMATTED', file_exists )
     WRITE( UNIT = 4, FMT = * ) etot, istep, nsteppos, dt_init, dt_curr, &
          alpha, alpha_init, tau(:,:), vel_old(:,:), acc_old(:,:)
     !
     CLOSE( UNIT = 4, STATUS = 'KEEP' )
     !
  ELSE
     !
     ! if file did not exist before, close and delete
     !
     CLOSE( UNIT = 4, STATUS = 'DELETE' )
     !
  END IF
  !

END SUBROUTINE move_mode




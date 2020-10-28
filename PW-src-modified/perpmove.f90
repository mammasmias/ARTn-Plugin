
SUBROUTINE perpmove(nat, istepperp, push)
  !
  ! subroutine that updates the fire acceleration and velocity so that we don't move along the push direction during perpendicular relax
  !
  USE kinds, ONLY : DP
  USE cell_base, ONLY : alat
  USE constants, ONLY : amu_ry
  USE io_files,  ONLY: prefix,seqopn,tmp_dir
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nat,istepperp
  REAL(DP),  INTENT(IN) :: push(3,nat)
  ! variables read from the FIRE minimization algorithm
  INTEGER :: istep, nsteppos
  REAL(DP) :: etot, P, dt_init, dt_curr, alpha, alpha_init
  REAL(DP), EXTERNAL :: ddot,dnrm2
  REAL(DP) :: tau(3,nat), vel(3,nat), acc(3,nat)
  REAL(DP), DIMENSION(3,nat) :: vel_old(3,nat), acc_old(3,nat)
  LOGICAL :: file_exists
  !
  ! Open the file written by fire
  !
  CALL seqopn( 4, 'fire', 'FORMATTED', file_exists )
  !
  write (*,*) "ARTn Called perpmove, number of perp steps:", istepperp
  IF ( file_exists ) THEN
     !
     ! Everything is OK! read the file
     !
     READ( UNIT = 4, FMT = * ) etot, istep, nsteppos, dt_init, dt_curr, alpha_init, &
          alpha, tau(:,:), vel_old(:,:), acc_old(:,:)
     !
     ! set the velocity and acceleration and alpha of previous step to move correctly
     !
     IF ( istepperp == 0 ) THEN
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
     WRITE (*,*) "ARTn: Perp relax out: dot(vel,push):", ddot(3*nat,vel_old, 1, push, 1 )
     WRITE (*,*) "ARTn: Perp relax out: dot(acc,push):", ddot(3*nat,acc_old, 1, push, 1 )
     ! overwrite the FIRE file
     CLOSE( UNIT = 4, STATUS = 'KEEP' )
     CALL seqopn( 4, 'fire', 'FORMATTED', file_exists )

     WRITE( UNIT = 4, FMT = * ) etot, istep, nsteppos, dt_init, dt_curr, alpha, alpha_init, &
          tau(:,:), vel_old(:,:), acc_old(:,:)
     CLOSE( UNIT = 4, STATUS = 'KEEP' )
  ELSE
     WRITE (*,*) "ARTn only works with the FIRE minimization algorithm! STOPPING "
     CLOSE( UNIT = 4, STATUS = 'DELETE' )
  END IF

END SUBROUTINE perpmove

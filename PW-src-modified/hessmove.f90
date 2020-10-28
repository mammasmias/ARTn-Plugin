
SUBROUTINE hessmove(force, nat, moves_finished )
  !
  ! subroutine that rescales the force so it corresponds to a specific move in tau ...
  !
  USE kinds, ONLY : DP
  USE ions_base, ONLY : tau
  USE cell_base, ONLY : alat
  USE control_flags,               ONLY : istep
  USE ener,               ONLY : etot
  USE constants, ONLY : amu_ry
  USE io_files,  ONLY: prefix,seqopn,tmp_dir
  !
  IMPLICIT NONE
  REAL(DP), INTENT(INOUT) :: force(3,nat)
  LOGICAL, INTENT(INOUT) :: moves_finished
  INTEGER, INTENT(IN) :: nat
  ! variables read from the FIRE minimization algorithm
  INTEGER :: nsteppos
  REAL(DP) :: P, dt_init, dt_curr, alpha, alpha_init
  REAL(DP), EXTERNAL :: ddot,dnrm2
  REAL(DP) :: vel(3,nat), acc(3,nat)
  REAL(DP), DIMENSION(3,nat) :: vel_old(3,nat), acc_old(3,nat)
  LOGICAL :: file_exists
  nsteppos = 1
  !
  ! Open the file written by fire
  !
  CALL seqopn( 4, 'fire', 'FORMATTED', file_exists )
  !
  write (*,*) "ARTn Called hess move, resetting FIRE"
  ! IF ( file_exists ) THEN
     !
     ! Everything is OK! read the file
     !
     ! READ( UNIT = 4, FMT = * ) etot, istep, nsteppos, dt_init, dt_curr, alpha, alpha_init, &
          ! tau(:,:), vel_old(:,:), acc_old(:,:)

     ! CLOSE( UNIT = 4, STATUS = 'KEEP' )
     !
     ! set the velocity and acceleration and alpha of previous step to move correctly
     !
     ! for the first step forget previous acc and velocity (prevent P < 0)
     acc_old(:,:) = 0.D0
     vel_old(:,:) = 0.D0
     alpha = 0.D0
     dt_curr = 0.D0
     !IF ( moves_finished ) THEN
     !   ! set dt to 0 so as not to move when loop restarts
     !   dt_curr = 0.0
     !   moves_finished = .false.
     !ELSE
     !   dt_curr = dt_init
     !   force(:,:) = force(:,:)*alat*amu_ry/dt_curr**2
     !ENDIF
     ! overwrite the FIRE file
     ! CALL seqopn( 4, 'fire', 'FORMATTED', file_exists )
     WRITE( UNIT = 4, FMT = * ) etot, istep, nsteppos, dt_init, dt_curr, alpha, alpha_init, &
          tau(:,:), vel_old(:,:), acc_old(:,:)
     CLOSE( UNIT = 4, STATUS = 'KEEP' )
  ! ELSE
  !    CLOSE( UNIT = 4, STATUS = 'DELETE' )
  !    dt_curr = 20.D0
  !    force(:,:) = force(:,:)*alat*amu_ry/dt_curr**2
  ! END IF

END SUBROUTINE hessmove


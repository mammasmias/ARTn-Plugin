
SUBROUTINE lancmove(nat, v1, dlanc, force )
  USE kinds, ONLY : DP
   USE cell_base, ONLY : alat
   USE constants, ONLY : amu_ry
   USE io_files,  ONLY: prefix,seqopn,tmp_dir
   !
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: nat
   REAL(DP), INTENT(IN) :: dlanc
   REAL(DP),  INTENT(IN) :: v1(3,nat)
   REAL(DP),  INTENT(INOUT) :: force(3,nat)
   REAL(DP) :: step(3,nat)
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

   IF ( file_exists ) THEN
      !
      ! Everything is OK! read the file
      !
      READ( UNIT = 4, FMT = * ) etot, istep, nsteppos, dt_init, dt_curr, alpha_init, alpha, &
           tau(:,:), vel_old(:,:), acc_old(:,:)
      !
      ! set the velocity and acceleration and alpha of previous step to move correctly
      !
      vel_old(:,:) = 0.D0
      acc_old(:,:) = 0.D0
      dt_curr = dt_init
      alpha = 0.D0
      nsteppos = 0
      ! the step performed should be like this now translate it into the correct force
      WRITE (*,*) "Expected lanczos step:", v1(:,:)*dlanc
      force(:,:) = v1(:,:)*dlanc*alat*amu_ry/dt_curr**2
      ! overwrite  fire parameters
      CLOSE( UNIT = 4, STATUS = 'KEEP' )
      CALL seqopn( 4, 'fire', 'FORMATTED', file_exists )

      WRITE( UNIT = 4, FMT = * ) etot, istep, nsteppos, dt_init, dt_curr, alpha_init, alpha, &
           tau(:,:), vel_old(:,:), acc_old(:,:)
      CLOSE( UNIT = 4, STATUS = 'KEEP' )
   ELSE
      WRITE (*,*) "ARTn only works with the FIRE minimization algorithm! STOPPING "
      CLOSE( UNIT = 4, STATUS = 'DELETE' )
   END IF

END SUBROUTINE lancmove

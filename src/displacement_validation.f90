
SUBROUTINE displacement_validation( atom_id, atom_const, push, lvalid)
  !
  ! subroutine that checks if the initial_displacement is within given parameters
  !
  USE artn_params, ONLY : DP, PI 
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: atom_id
  REAL(DP), INTENT(IN) :: atom_const(4)
  REAL(DP), INTENT(INOUT) :: push(3)
  REAL(DP), EXTERNAL :: ddot, dnrm2
  LOGICAL,         INTENT(INOUT) :: lvalid
  !
  ! Local variables
  REAL(DP)               :: cone_angle, displacement_angle
  REAL(DP)               :: dot_prod, displacement_norm, cone_dir_norm
  REAL(DP), DIMENSION(3) :: cone_dir,displacement
  !
  !
  !
  write (*,*) "ARTn: Called displacement validation: with cone_dir:",atom_const(1:3), "current push:", push(:)
  cone_dir = atom_const(1:3)
  cone_angle = atom_const(4)
  !
  displacement(:)         = push(:)
  !
  displacement_norm       = dnrm2(3,displacement,1)
  cone_dir_norm           = dnrm2(3,cone_dir,1)
  !
  dot_prod                = ddot( 3, cone_dir, 1, displacement, 1 ) / ( cone_dir_norm * displacement_norm )
  displacement_angle      = ACOS( dot_prod ) *180.0_DP / PI
  lvalid                  = ( displacement_angle < cone_angle )
  write (*,*) "Finished displacement validation",lvalid
  !
  IF ( cone_angle == 0.0_DP) THEN
     lvalid = .TRUE.
     !
     ! TODO: why is the direction multiplied by 0.1? seems kind of random ...
     !
     push(:) = cone_dir(:)
  ENDIF
  !
END SUBROUTINE displacement_validation

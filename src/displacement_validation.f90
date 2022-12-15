
!> @author
!!  Matic Poberznik,
!!  Miha Gunde
!!  Nicolas Salles
!
SUBROUTINE displacement_validation( atom_const, push, lvalid)
  !
  !> @brief
  !!   subroutine that checks if the initial_displacement is within given parameters
  !
  !> @param [in]      atom_id     id of atoms
  !> @param [in]      atom_const  Constrain Applied on the atom
  !> @param [inout]   push	  Direction of the push applied on the atoms
  !> @param [inout]   lvalid	  Flag to know if the random displacement correspond to the constrain
  !
  USE units, only : DP, PI
  !
  IMPLICIT NONE
  !INTEGER, INTENT(IN) :: atom_id
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
  !write (*,*) " * ARTn: Called displacement validation:", atom_id 
  !write(*,*) " with cone_dir:",atom_const(1:3), "current push:", push(:)
  cone_dir = atom_const(1:3)
  cone_angle = atom_const(4)
  !print*, 
  !
  displacement(:)         = push(:)
  !
  displacement_norm       = dnrm2(3,displacement,1)
  cone_dir_norm           = dnrm2(3,cone_dir,1)
  !
  dot_prod                = ddot( 3, cone_dir, 1, displacement, 1 ) / ( cone_dir_norm * displacement_norm )
  displacement_angle      = ACOS( dot_prod ) *180.0_DP / PI
  lvalid                  = ( displacement_angle < cone_angle )
  !write (*,*) "Finished displacement validation",lvalid  !&
          !, displacement_norm, cone_dir_norm, dot_prod, displacement_angle, atom_const
  !
  IF ( cone_angle == 0.0_DP) THEN
     lvalid = .TRUE.
     !
     ! TODO: why is the direction multiplied by 0.1? seems kind of random ...
     !
     push(:) = cone_dir(:)
     !print*, "**disp_valid",push, norm2(push)
  ENDIF
  !
  ! When the atom is not pushed, constrain is useless
  IF (all(push(:) .EQ. 0,1)) lvalid = .TRUE.
  !
END SUBROUTINE displacement_validation



!...........................................................................................
!subroutine apply_constrain( atom_const, push )

!  use units, only : DP
!  implicit none

!  integer, intent(in) :: atom_const(4)
!  REAL(DP), INTENT(INOUT) :: push(3)

  ! Local variables
!  REAL(DP)               :: cone_angle, displacement_angle
!  REAL(DP)               :: dot_prod, displacement_norm, cone_dir_norm
!  REAL(DP), DIMENSION(3) :: cone_dir,displacement
!  REAL(DP) :: r1, r2, r3

  
  

!end subroutine apply_constrain









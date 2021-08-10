
!> @author
!!   Matic Poberznik,
!!   Miha Gunde
SUBROUTINE perpforce( force, if_pos, push, fpara, nat )
  !
  !> @brief subroutine that subtracts parallel components to push from force
  !
  !> @param[inout]  force	Force field and become Fperp field following Push field
  !! @param[in]	    if_pos	Constrain field
  !! @param[in]	    push	Parallel field reference
  !! @param[out]    fpara	Parallel force field following Push field
  !! @param[in]	    nat		number of point in the field 
  !
  USE artn_params, ONLY : DP
  IMPLICIT NONE
  ! -- ARGUMENTS
  INTEGER,  INTENT(IN)     :: nat
  REAL(DP), INTENT(IN)     :: push(3,nat)
  REAL(DP), INTENT(INOUT)  :: force(3,nat)
  REAL(DP), INTENT(OUT)    :: fpara(3,nat)
  INTEGER,  INTENT(IN)     :: if_pos(3,nat)
  ! -- LOCAL VARIABLE
  REAL(DP) :: push_norm(3,nat)
  REAL(DP), EXTERNAL :: ddot,dnrm2
  INTEGER  :: na

  ! calculate components parallel to the push
  fpara(:,:) = ddot(3*nat,force(:,:),1,push(:,:),1) / ddot(3*nat,push(:,:),1,push(:,:),1) * push(:,:)

  ! subtract them
  force(:,:) = force(:,:) - fpara(:,:)

  ! apply constraints
  IF ( ANY(if_pos(:,:) == 0)  ) force(:,:) = force(:,:)*if_pos(:,:) 

END SUBROUTINE perpforce


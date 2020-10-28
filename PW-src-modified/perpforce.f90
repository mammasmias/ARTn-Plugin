
SUBROUTINE perpforce(force,push,fpara,nat)
  !
  ! subroutine that subtracts parallel components to push from force
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  REAL(DP), INTENT(IN)  :: push(3,nat)
  REAL(DP), INTENT(INOUT)  :: force(3,nat)
  REAL(DP), INTENT(OUT)  :: fpara(3,nat)
  REAL(DP) :: push_norm(3,nat)
  REAL(DP), EXTERNAL :: ddot,dnrm2
  INTEGER, INTENT(IN) :: nat
  INTEGER  :: na
  ! calculate components parallel to the push
  fpara(:,:) = ddot(3*nat,force(:,:),1,push(:,:),1) / ddot(3*nat,push(:,:),1,push(:,:),1) * push(:,:)
  ! subtract them
  force(:,:) = force(:,:) - fpara(:,:)

END SUBROUTINE perpforce

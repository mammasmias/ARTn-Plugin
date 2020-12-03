
SUBROUTINE sum_force(force,nat,force_tot)
  !
  ! subroutine that sums the forces on all atoms and returns the total force
  !
  USE artn_params, ONLY : DP
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nat
  REAL(DP), INTENT(IN) :: force(3,nat)
  REAL(DP), INTENT(OUT) :: force_tot
  INTEGER :: na
  force_tot = 0.0
  DO na = 1, nat
     force_tot = force_tot + force(1,na)**2 + force(2,na)**2 + force(3,na)**2
  ENDDO
  force_tot = SQRT(force_tot)
END SUBROUTINE sum_force

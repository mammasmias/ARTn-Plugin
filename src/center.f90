
SUBROUTINE center ( vec, nat)
  USE artn_params, ONLY: DP
  !
  ! takes as input a vector of size (3,nat) and centers it
  !
  IMPLICIT none
  INTEGER, INTENT(IN) :: nat
  REAL(DP), INTENT(INOUT) :: vec(3,nat)
  INTEGER :: na
  REAL(DP) :: delta(3)
  !
  delta(:) = 0.D0
  DO na = 1,nat
     delta(:) = delta(:) + vec(:,na)
  ENDDO
  !
  delta(:) = delta(:)/dble(nat)
  !
  FORALL ( na = 1:nat) vec(:,na) = vec(:,na) - delta(:)


END SUBROUTINE center

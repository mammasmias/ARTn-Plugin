!
!> @author
!!  Matic Poberznik
!!  Miha Gunde
!!  Nicolas Salles
!
!
SUBROUTINE center ( vec, nat )
  !> @brief
  !!   takes as input a vector of size (3,nat) and centers it
  !
  !> @param[in]     nat    number of atom
  !! @param[inout]  vec    output vector
  !
  USE units, ONLY: DP
  !
  !
  IMPLICIT none
  INTEGER,  INTENT(IN) :: nat             !> Size of vec array: Number of atom
  REAL(DP), INTENT(INOUT) :: vec(3,nat)   !> Vector will be centered
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


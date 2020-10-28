
SUBROUTINE push_init (nat, idum, init_step_size, push)
  !
  ! subroutine for generating a random push on all atoms (global_move scenario)
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT none
  INTEGER, INTENT(IN) :: nat, idum
  ! INTEGER, INTENT(IN) :: if_pos(3,nat)
  REAL(DP), INTENT(IN) :: init_step_size
  INTEGER :: na
  REAL(DP), INTENT(OUT) :: push(3,nat)
  REAL(DP), EXTERNAL :: ran3
  REAL(DP) :: dr2
  !
  ! generate a random push
  !
  DO na=1,nat
     DO
        push(1,na) = (0.5_DP - ran3(idum))
        push(2,na) = (0.5_DP - ran3(idum))
        push(3,na) = (0.5_DP - ran3(idum))
        dr2 = push(1,na)**2 + push(2,na)**2 + push(3,na)**2
        IF ( dr2 < 0.25_DP ) EXIT
     ENDDO
  ENDDO

  !
  ! scale the push according to initial step size
  !
  push(:,:) = init_step_size*push(:,:)

END SUBROUTINE push_init

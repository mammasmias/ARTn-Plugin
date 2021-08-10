
!> @author
!!   Matic Poberznik,
!!   Miha Gunde


SUBROUTINE sum_force( force, nat, force_tot )
  !
  !> @brief
  !!   subroutine that sums the forces on all atoms and returns the total force
  !
  !> @param [in]  nat	    Size of list: number of atoms
  !> @param [in]  force	    List of atomic forces
  !> @param [out] force_tot Sum of the list of forces
  !
  USE artn_params, ONLY : DP
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nat
  REAL(DP), INTENT(IN) :: force(3,nat)
  REAL(DP), INTENT(OUT) :: force_tot
  INTEGER :: na

  force_tot = 0.0
  DO na = 1, nat
     !force_tot = force_tot + force(1,na)**2 + force(2,na)**2 + force(3,na)**2
     force_tot = force_tot + dot_product( force(:,na), force(:,na) )
  ENDDO

  force_tot = SQRT(force_tot)

END SUBROUTINE sum_force

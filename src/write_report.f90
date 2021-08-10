
!> @author
!!   Matic Poberznik,
!!   Miha Gunde


SUBROUTINE write_report( etot, force, lowest_eigval, disp, if_pos, istep, nat, iunartout )
  !
  !> @brief
  !!   a subroutine that writes a report of the current step to the output file  
  !
  !> @param [in]  etot		energy of the system
  !> @param [in]  force		List of atomic forces
  !> @param [in]  lowest_eigval	Lowest eigenvalue obtained by lanczos
  !> @param [in]  disp		Kind of actual displacement 
  !> @param [in]  if_pos	Fix the atom or not	
  !> @param [in]  istep		actual step of ARTn 
  !> @param [in]  iunartout	Channel of output
  !
  USE artn_params, ONLY: push, MOVE 
  USE UNITS
  IMPLICIT NONE
  ! -- Arguments
  INTEGER,  INTENT(IN) :: nat, istep, iunartout
  INTEGER,  INTENT(IN) :: if_pos(3,nat)
  REAL(DP), INTENT(IN) :: force(3,nat), etot, lowest_eigval
  INTEGER,  INTENT(IN) :: disp
  ! -- Local Variables
  REAL(DP) :: fpara(3,nat), fperp(3,nat)
  REAL(DP) :: force_tot, fperp_tot, fpara_tot
  REAL(DP), EXTERNAL :: ddot
  !
  CALL sum_force(force,nat,force_tot)

  fperp(:,:) = force(:,:)

  CALL perpforce(fperp,if_pos,push,fpara,nat)

  CALL sum_force(fperp,nat,fperp_tot)

  fpara_tot = ddot(3*nat,force,1,push,1)

  !  write report 
  WRITE (iunartout,'(5X,I4,7X,A4,9X, F12.6, *(5X, F10.4))')  & !,5X, F7.4, 5X, F7.4, 5X, F7.4)') &
       & istep, MOVE(disp), etot, force_tot,fperp_tot,fpara_tot, lowest_eigval !, MAXVAL(fperp)
       !& istep, MOVE(disp), unconvert_energy(etot), force_tot, fperp_tot, fpara_tot, lowest_eigval, MAXVAL(fperp)


END SUBROUTINE write_report

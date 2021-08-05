
!SUBROUTINE write_report(etot, force, lowest_eigval, step, if_pos, istep, nat,iunartout)
SUBROUTINE write_report(etot, force, lowest_eigval, step, if_pos, istep, nat,iunartout)
  !
  ! a subroutine that writes a report of the current step to the output file  
  !
  USE artn_params, ONLY: DP, push, MOVE 
  USE UNITS
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nat, istep, iunartout
  INTEGER, INTENT(IN) :: if_pos(3,nat)
  REAL(DP), INTENT(IN) :: force(3,nat), etot, lowest_eigval
  REAL(DP) :: fpara(3,nat), fperp(3,nat)
  REAL(DP) :: force_tot, fperp_tot, fpara_tot
  !CHARACTER(LEN=4), INTENT(IN) :: step 
  INTEGER, INTENT(IN) :: step 
  REAL(DP), EXTERNAL :: ddot
  !
  CALL sum_force(force,nat,force_tot)
  fperp(:,:) = force(:,:)
  CALL perpforce(fperp,if_pos,push,fpara,nat)
  CALL sum_force(fperp,nat,fperp_tot)
  fpara_tot = ddot(3*nat,force,1,push,1)
  !  write report 
  WRITE (iunartout,'(5X,I4,7X,A4,9X, F12.6, *(5X, F10.4))')  & !,5X, F7.4, 5X, F7.4, 5X, F7.4)') &
       & istep, MOVE(step), etot, force_tot,fperp_tot,fpara_tot, lowest_eigval, MAXVAL(fperp)
       !& istep, MOVE(step), unconvert_energy(etot), force_tot, fperp_tot, fpara_tot, lowest_eigval, MAXVAL(fperp)
END SUBROUTINE write_report

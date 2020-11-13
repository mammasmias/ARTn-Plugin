
SUBROUTINE report_force(force, if_pos, forc_conv_thr_qe, eigenvec,nat,force_tot,fperp_tot,fpara_tot,iunartout)
  !
  ! a subroutine that reports the non-manipulated forces of QE (Ftot, Fpara, and Fperp)
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nat, iunartout
  INTEGER, INTENT(IN) :: if_pos(3,nat)
  REAL(DP), INTENT(IN) :: force(3,nat), eigenvec(3,nat), forc_conv_thr_qe
  REAL(DP) :: fpara(3,nat), fperp(3,nat)
  REAL(DP), INTENT(OUT) :: force_tot, fperp_tot, fpara_tot
  REAL(DP), EXTERNAL :: ddot
  !
  CALL sum_force(force,nat,force_tot)
  fperp(:,:) = force(:,:)
  CALL perpforce(fperp,if_pos,eigenvec,fpara,nat)
  CALL sum_force(fperp,nat,fperp_tot)
  fpara_tot = ddot(3*nat,force,1,eigenvec,1)
  ! CALL sum_force(fpara,nat,fpara_tot)
  WRITE (iunartout,'("ARTn Forces (Ry/a.u.):",5X, "Ftot = ",F6.3,5X,"Fperp = ",F6.3, 5X,"Fpara = ",F6.3)') &
       & force_tot,fperp_tot,fpara_tot
END SUBROUTINE report_force

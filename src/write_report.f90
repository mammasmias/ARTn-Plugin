
!> @author
!!   Matic Poberznik,
!!   Miha Gunde



SUBROUTINE write_initial_report(iunartout, filout)
  use artn_params, ONLY: engine_units, ninit, nperp, neigen, nsmooth,  &
                         init_forc_thr, forc_thr, fpara_thr, eigval_thr, &
                         push_step_size, eigen_step_size, lanc_mat_size, dlanc, &
                         push_mode
  use units, only : strg_units
  INTEGER,             INTENT(IN) :: iunartout
  CHARACTER (LEN=255), INTENT(IN) :: filout
  ! -- Local Variables
  INTEGER :: ios
  !
  ! Writes the header to the artn output file
  !
  OPEN ( UNIT = iunartout, FILE = filout, FORM = 'formatted', STATUS = 'unknown', IOSTAT = ios )
  WRITE (iunartout,'(5X, "--------------------------------------------------")')
  WRITE (iunartout,'(5X, "                ARTn plugin                       ")')
  WRITE (iunartout,'(5X, "--------------------------------------------------")')
  WRITE (iunartout,'(5X, " "                                                 )')
  WRITE (iunartout,'(5X, "               INPUT PARAMETERS                   ")')
  WRITE (iunartout,'(5X, "--------------------------------------------------")')
  WRITE (iunartout,'(5x, "engine_units:", *(x,A))') TRIM(engine_units)
  WRITE (iunartout,'(5X, "--------------------------------------------------")')
  WRITE (iunartout,'(5X, "Push and perpendicular relax:")')
  WRITE (iunartout,'(5X, "--------------------------------------------------")')
  WRITE (iunartout,'(15X,"ninit           = ", I6)') ninit
  WRITE (iunartout,'(15X,"nperp           = ", I6)') nperp
  WRITE (iunartout,'(15X,"neigen          = ", I6)') neigen
  WRITE (iunartout,'(15X,"nsmooth         = ", I6)') nsmooth
  WRITE (iunartout,'(15X,"Threshold Parameter: ")')
  WRITE (iunartout,'(15X,"init_forc_thr   = ", F6.3)') init_forc_thr
  WRITE (iunartout,'(15X,"forc_thr        = ", F6.3)') forc_thr
  WRITE (iunartout,'(15X,"fpara_thr       = ", F6.3)') fpara_thr
  WRITE (iunartout,'(15X,"eigval_thr      = ", F6.3)') eigval_thr
  WRITE (iunartout,'(15X,"push_step_size  = ", F6.1)') push_step_size
  WRITE (iunartout,'(15X,"eigen_step_size = ", F6.1)') eigen_step_size
  WRITE (iunartout,'(15X,"push_mode       = ", A6)') push_mode
  WRITE (iunartout,'(5X, "--------------------------------------------------")')
  WRITE (iunartout,'(5X, "Lanczos algorithm:")' )
  WRITE (iunartout,'(5X, "--------------------------------------------------")')
  WRITE (iunartout,'(15X, "lanc_mat_size     = ", I6)') lanc_mat_size
  WRITE (iunartout,'(15X, "dlanc          = ", F6.3)') dlanc
  WRITE (iunartout,'(5X, "--------------------------------------------------")')
  WRITE (iunartout,'(/,/)') 
  !WRITE (iunartout,*) " "
  !%! Condition on the engin_units..
  WRITE (iunartout,'(5X,"istep",4X,"ART_step",4X,"Etot",5x,"init/eig/ip/il","&
                    "3X," Ftot ",5X," Fperp ",5X," Fpara ",3X,"eigval", 6X, "delr", 2X, "npart", X,"evalf")')
  !WRITE (iunartout,'(27X, "[Ry]",17X,"-----------[Ry/a.u.]----------",3X,"Ry/a.u.^2")')
  WRITE (iunartout, strg_units )


  CLOSE ( UNIT = iunartout, STATUS = 'KEEP')

END SUBROUTINE write_initial_report





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
  USE artn_params, ONLY: push, MOVE  &
                        ,etot_init, iinit, iperp, ieigen, ilanc, delr
  USE UNITS
  IMPLICIT NONE

  ! -- Arguments
  INTEGER,  INTENT(IN) :: nat, istep, iunartout
  INTEGER,  INTENT(IN) :: if_pos(3,nat)
  REAL(DP), INTENT(IN) :: force(3,nat), etot, lowest_eigval
  INTEGER,  INTENT(IN) :: disp

  ! -- Local Variables
  CHARACTER(LEN=5) :: Mstep
  integer :: macrostep = -1, evalf, i, npart
  REAL(DP) :: fpara(3,nat), fperp(3,nat)
  REAL(DP) :: force_tot, fperp_tot, fpara_tot, detot, lowEig, dr, rc2
  REAL(DP), EXTERNAL :: ddot
  !


  CALL perpforce(force,if_pos,push,fperp,fpara,nat)

  ! ...Force processing
  CALL sum_force( force, nat, force_tot )
  CALL perpforce(force,if_pos,push,fperp,fpara,nat)
  CALL sum_force(fperp,nat,fperp_tot)
  fpara_tot = ddot(3*nat,force,1,push,1)
  ! ...Displacement processing
  npart = 0
  rc2 = 0.1*0.1
  do i = 1, nat
     if( norm2(delr(:,i)) > rc2 ) npart = npart + 1
  enddo
  call sum_force( delr, nat, dr )

  ! .. Convertion Units
  fperp_tot = unconvert_force( fperp_tot )
  fpara_tot = unconvert_force( fpara_tot )
  dEtot = unconvert_energy(etot - etot_init)
  lowEig = unconvert_hessian( lowest_eigval )
  !lowEig = lowest_eigval

  !  write report 
  !WRITE (iunartout,'(5X,I4,7X,A4,9X, F12.6, *(5X, F10.4))')  & !,5X, F7.4, 5X, F7.4, 5X, F7.4)') &
  !     & istep, MOVE(disp), etot, force_tot,fperp_tot,fpara_tot, lowest_eigval !, MAXVAL(fperp)
       !& istep, MOVE(disp), unconvert_energy(etot), force_tot, fperp_tot, fpara_tot, lowest_eigval, MAXVAL(fperp)

  !%! More Complete Output
  Mstep = "macrostep"
  !if( lbasin ) Mstep = 'Bstep'
  !if( .not.lbasin ) Mstep = 'Sstep'
  !delr = sum()
  evalf = istep
  WRITE(iunartout,5) istep, Mstep, MOVE(disp), detot, iinit, ieigen, iperp, ilanc,   &
                     force_tot, fperp_tot, fpara_tot, lowEig,     &
                     dr, npart, evalf !, a1
  !5 format(5x,i4,3x,a,x,a,F10.4,3x,4(x,i2),4(x,f10.4))
  5 format(5x,i4,3x,a,x,a,F10.4,3x,4(x,i2),5(x,f10.4),2(x,i4))

END SUBROUTINE write_report




SUBROUTINE write_end_report( iunartout, lsaddle, lpush_final, de )
  
  use units, only : DP, unconvert_energy
  implicit none

  integer, intent( in ) :: iunartout
  logical, intent( in ) :: lsaddle, lpush_final
  REAL(DP), intent( in ), value :: de
  
  if( lsaddle )then
    WRITE (iunartout,'(5X, "--------------------------------------------------")')
    WRITE (iunartout,'(5X, "    *** ARTn found a potential saddle point ***   ")')
    WRITE (iunartout,'(5X, "--------------------------------------------------")')
    WRITE (iunartout,'(15X,"E_final - E_initial =", F12.5," eV")') unconvert_energy(de)
    WRITE (iunartout,'(5X, "--------------------------------------------------")')

    IF ( lpush_final ) THEN
       WRITE (iunartout, '(5X,"       *** Pushing to adjacent minima  ***      ")')
       WRITE (iunartout,'(5X, "------------------------------------------------")')
    ENDIF

  else
    WRITE (iunartout,'(5X, "--------------------------------------------------")')
    WRITE (iunartout,'(5X, "        *** ARTn saddle search failed  ***        ")')
    WRITE (iunartout,'(5X, "--------------------------------------------------")')
  endif

END SUBROUTINE write_end_report




subroutine compute_delr( nat, pos, lat )
  use artn_params, only : delr, tau_step
   use units, only : DP
  implicit none

  INTEGER, intent( in ) :: nat
  REAL(DP), intent( in ) :: pos(3,nat), lat(3,3)

  integer :: i
  REAL(DP) :: dr(3,nat), r(3)
  REAL(DP), external :: fpbc
  
  !dr = pos - tau_step
  do i = 1, nat
     !delr(:,i) = delr(:,i) + fpbc( dr(:,i), lat )
     !call pbc( dr(:,i), lat )
     !delr(:,i) = delr(:,i) + dr(:,i)
     r = pos(:,i) - tau_step(:,i)
     call pbc( r, lat )
     delr(:,i) = delr(:,i) + r(:)
  enddo
end subroutine compute_delr





!> @author
!!   Matic Poberznik,
!!   Miha Gunde
!!   Nicolas Salles


!------------------------------------------------------------
SUBROUTINE write_initial_report(iunartout, filout)
  use artn_params, ONLY: engine_units, ninit, nperp, neigen, nsmooth,  &
                         init_forc_thr, forc_thr, fpara_thr, eigval_thr, &
                         push_step_size, eigen_step_size, lanc_mat_size, dlanc, &
                         push_mode, verbose, push_over, frelax_ene_thr, zseed
  use units, only : strg_units, unconvert_force, &
                    unconvert_energy, unconvert_hessian, unconvert_length, unit_char
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
  WRITE (iunartout,'(5x, "Verbosity Level:", *(x,i2))') verbose
  WRITE (iunartout,'(5X, "--------------------------------------------------")')
  WRITE (iunartout,'(5X, "Simulation Parameters:")')
  WRITE (iunartout,'(5X, "--------------------------------------------------")')
  WRITE (iunartout,'(13X,"* Iterators Parameter: ")')
  !WRITE (iunartout,'(15X,"Zseed           = ", I6)') zseed
  WRITE (iunartout,'(15X,"ninit           = ", I6)') ninit
  WRITE (iunartout,'(15X,"nperp           = ", I6)') nperp
  WRITE (iunartout,'(15X,"neigen          = ", I6)') neigen
  WRITE (iunartout,'(15X,"nsmooth         = ", I6)') nsmooth
  WRITE (iunartout,'(13X,"* Threshold Parameter: ")')
  WRITE (iunartout,'(15X,"init_forc_thr   = ", F6.3,2x,A)') unconvert_force( init_forc_thr ), unit_char('force')
  WRITE (iunartout,'(15X,"forc_thr        = ", F6.3,2x,A)') unconvert_force( forc_thr ), unit_char('force')
  WRITE (iunartout,'(15X,"fpara_thr       = ", F6.3,2x,A)') unconvert_force( fpara_thr ), unit_char('force')
  WRITE (iunartout,'(15X,"eigval_thr      = ", F6.3,2x,A)') unconvert_hessian( eigval_thr ), unit_char('hessian')
  WRITE (iunartout,'(15X,"frelax_ene_thr  = ", F6.3,2x,A)') unconvert_energy( frelax_ene_thr ), unit_char('energy')
  WRITE (iunartout,'(13X,"* Step size Parameter: ")')
  WRITE (iunartout,'(15X,"push_step_size  = ", F6.2,2x,A)') unconvert_length( push_step_size ), unit_char('length')
  WRITE (iunartout,'(15X,"eigen_step_size = ", F6.2,2x,A)') unconvert_length( eigen_step_size ), unit_char('length')
  WRITE (iunartout,'(15X,"push_over       = ", F6.3,2x,A)') push_over, "fraction of eigen_step_size"
  WRITE (iunartout,'(15X,"push_mode       = ", A6)') push_mode
  WRITE (iunartout,'(5X, "--------------------------------------------------")')
  WRITE (iunartout,'(5X, "Lanczos algorithm:")' )
  WRITE (iunartout,'(5X, "--------------------------------------------------")')
  WRITE (iunartout,'(15X,"lanc_mat_size   = ", I6)') lanc_mat_size
  WRITE (iunartout,'(15X,"dlanc           = ", F6.3,2x,A)') unconvert_length( dlanc ), unit_char('length')
  WRITE (iunartout,'(5X, "--------------------------------------------------")')
  WRITE (iunartout,'(/,/)') 
  !WRITE (iunartout,*) " "

  !%! Condition on the engin_units..
  select case( verbose )
    case( 0 )
    WRITE (iunartout,'(5X,"istep",4X,"ART_step",4X,"Etot",5x,"init/eign/perp/lanc/relx","&
                      "4X," Ftot ",5X," Fperp ",4X," Fpara ",4X,"eigval", 6X, "delr", 2X, "npart", X,"evalf",2X,"a1")')

    case( 1: )
    WRITE (iunartout,'(5X,"istep",4X,"ART_step",4X,"Etot",5x,"init/eign/perp/lanc/relx","&
                      "4X," Ftot ",5X," Fperp ",4X," Fpara ",4X,"eigval", 6X, "delr", 2X, "npart", X,"evalf","&
                      "2X,"B/S/R|I/P/L/E|P/B/R",4X,"a1")')
  end select

  ! -- Units
  WRITE (iunartout, strg_units )


  CLOSE ( UNIT = iunartout, STATUS = 'KEEP')

END SUBROUTINE write_initial_report





!------------------------------------------------------------
SUBROUTINE write_report( etot, force, fperp, fpara, lowest_eigval, disp, if_pos, istep, nat, iunartout, ArtnStep )
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
  !> @param [in]  ARTnStep	Flag to print at ARTn step
  !
  USE artn_params, ONLY: push, MOVE, verbose  &
                        ,etot_init, iinit, iperp, ieigen, ilanc, irelax, delr, verbose, iartn, a1 &
                        ,tau_init, lat, tau_step, delr &
                        ,lrelax, linit, lbasin, lperp, llanczos, leigen, lsaddle, lpush_final, lbackward, lrestart 
  USE UNITS
  IMPLICIT NONE

  ! -- Arguments
  INTEGER,  INTENT(IN) :: nat, istep, iunartout
  INTEGER,  INTENT(IN) :: if_pos(3,nat)
  REAL(DP), INTENT(IN) :: force(3,nat), fpara(3,nat), fperp(3,nat), etot, lowest_eigval
  INTEGER,  INTENT(IN) :: disp
  LOGICAL, intent(IN)  :: ArtnStep

  ! -- Local Variables
  CHARACTER(LEN=5) :: Mstep
  integer :: macrostep = -1, evalf, i, npart
  REAL(DP) :: force_tot, fperp_tot, fpara_tot, detot, lowEig, dr, rc2
  REAL(DP), EXTERNAL :: ddot
  LOGICAL :: C1, C2
  !

  ! ...Print only ARTn-Step
  C1 = ( .NOT.ARTnStep .AND. verbose < 2) ! When verbose = {0,1}: print only ARTnStep
  IF( C1 )RETURN
  !C2 = ( verbose == 2  .AND. ARTnStep )   ! When verbose = 2: print only noARTnStep
  !IF( C1 .OR. C2 )RETURN



  ! ...Force processing
  force_tot = MAXVAL( ABS(force) )
  fperp_tot = MAXVAL( ABS(fperp) )
  fpara_tot = MAXVAL( ABS(fpara) )



  ! .. Conversion Units
  force_tot = unconvert_force( force_tot )
  fperp_tot = unconvert_force( fperp_tot )
  fpara_tot = unconvert_force( fpara_tot )

  !write( iunartout,*) "write_report::", etot,etot_init, etot - etot_init
  dEtot = unconvert_energy(etot - etot_init)
  lowEig = unconvert_hessian( lowest_eigval )
  !lowEig = lowest_eigval


  !%! More Complete Output
  Mstep = "macrostep"
  if( lbasin ) Mstep = 'Bstep'
  if( .not.lbasin ) Mstep = 'Sstep'
  if( lrelax ) Mstep = 'Rstep'


  !delr = sum()
  evalf = istep


  dr = 0.
  npart = 0
  IF( ARTnStep )THEN

    ! ...Displacement processing
    if( iartn == 0 )then
      if( .not.allocated(tau_init) )then
        allocate( tau_init, source = tau_step )
      else
        tau_init = tau_step
      endif
    endif
    call compute_delr( nat, tau_step, tau_init, lat, delr )
    npart = 0
    rc2 = 0.1!*0.1  !! Miha: Why square? NS: Why not!
    do i = 1, nat
       if( norm2(delr(:,i)) > rc2 ) npart = npart + 1
    enddo
    !! routine sum_force is equivalent to implicit: norm2( delr )
    call sum_force( delr, nat, dr )

  ENDIF


  select case( verbose )
    case( 0 )
      WRITE(iunartout,6) iartn, Mstep, MOVE(disp), detot, iinit, ieigen, iperp, ilanc, irelax,  &
                         force_tot, fperp_tot, fpara_tot, lowEig,     &
                         dr, npart, evalf, a1
      6 format(5x,i4,3x,a,x,a,F10.4,x,5(x,i4),5(x,f10.4),2(x,i4),3X,f4.2)

    case( 1: )
      WRITE(iunartout,5) iartn, Mstep, MOVE(disp), detot, iinit, ieigen, iperp, ilanc, irelax,  &
                         force_tot, fperp_tot, fpara_tot, lowEig,     &
                         dr, npart, evalf,   &
          lbasin, lsaddle, lrelax, linit, lperp, llanczos, leigen,  lpush_final, lbackward, lrestart , a1
      5 format(5x,i4,3x,a,x,a,F10.4,x,5(x,i4),5(x,f10.4),2(x,i4),3X,10(L2),3X,f4.2)

  end select

  IF( ARTnStep )iartn = iartn + 1

END SUBROUTINE write_report



!------------------------------------------------------------
SUBROUTINE write_inter_report( u, pushfactor, de )
  use units, only : DP, unconvert_energy
  use artn_params, only : artn_resume
  implicit none

  integer, intent( in ) :: u             !> Ouput Unit 
  integer, intent( in ) :: pushfactor
  real(DP), intent( in ) :: de(*)            !> list of energie 

  SELECT CASE( pushfactor )

    CASE( 1 )
      ! de(1) = de_back
      WRITE( u, '(5X, "--------------------------------------------------")')
      WRITE( u, '(5X, "|> ARTn found adjacent minimum | backward E_act =", F12.5," eV")') unconvert_energy( de(1) )
      WRITE( u, '(5X, "--------------------------------------------------")')

    CASE( -1 )
      ! de(1) = de_back
      ! de(2) = de_fwd 
      ! de(3) = etot_init 
      ! de(4) = etot_final
      ! de(5) = etot   
      WRITE( u,'(5X, "--------------------------------------------------")')
      WRITE( u,'(5X, "    *** ARTn converged to initial minimum ***   ")')
      WRITE( u,'(5X, "--------------------------------------------------")')
      WRITE( u,'(15X,"forward  E_act =", F12.5," eV")') unconvert_energy(de(2)) 
      WRITE( u,'(15X,"backward E_act =", F12.5," eV")') unconvert_energy(de(1)) 
      WRITE( u,'(15X,"reaction dE    =", F12.5," eV")') unconvert_energy((de(5)-de(4))) 
      WRITE( u,'(15X,"dEinit - dEfinal    =", F12.5," eV")') unconvert_energy((de(3)-de(5))) 
      WRITE( u,'(5X, "--------------------------------------------------")')
      WRITE( u,'(5X, "Cofiguration Files:", X,A)') trim(artn_resume)
      WRITE( u,'(5X, "--------------------------------------------------")')


    CASE DEFAULT
      WRITE(*,*) "********* ERROR write_inter_report:: pushfactor", pushfactor, " **************"
      

  END SELECT

END SUBROUTINE write_inter_report



!------------------------------------------------------------
SUBROUTINE write_end_report( iunartout, lsaddle, lpush_final, de )
 
  use units, only : DP, unconvert_energy
  implicit none

  integer, intent( in ) :: iunartout
  logical, intent( in ) :: lsaddle, lpush_final
  REAL(DP), intent( in ), value :: de
  
  if( lsaddle )then

    WRITE (iunartout,'(5X, "--------------------------------------------------")')
    WRITE (iunartout,'(5X, "|> ARTn found a potential saddle point | E_saddle - E_initial =", F12.5," eV")') unconvert_energy(de)
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


!------------------------------------------------------------
SUBROUTINE write_fail_report( iunartout, disp, estep )

  use units, only : DP, unconvert_energy
  use artn_params, only : MOVE
  implicit none

  integer, intent( in ) :: iunartout, disp
  REAL(DP), intent( in ):: estep

  WRITE (iunartout,'(5X, "--------------------------------------------------")')
  WRITE (iunartout,'(5X, "        *** ARTn search failed at ",a," ***        ")') MOVE(DISP)
  WRITE (iunartout,'(5X, "Step Params:",f10.4)') unconvert_energy(estep)
  WRITE (iunartout,'(5X, "--------------------------------------------------")')

END SUBROUTINE write_fail_report


!------------------------------------------------------------
subroutine compute_delr( nat, pos, old_pos, lat, delr )
  use units, only : DP
  implicit none

  INTEGER, intent( in ) :: nat
  REAL(DP), intent( in ) :: pos(3,nat), lat(3,3)
  REAL(DP), intent( in ) :: old_pos(3,nat)
  REAL(DP), intent( out ) :: delr(3,nat)

  integer :: i
  REAL(DP) :: r(3)
  
  delr = 0.0
  do i = 1, nat
     r = pos(:,i) - old_pos(:,i)
     call pbc( r, lat )
     delr(:,i) = r(:)
  enddo

end subroutine compute_delr




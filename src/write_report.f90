
!> @author
!!   Matic Poberznik,
!!   Miha Gunde
!!   Nicolas Salles


!------------------------------------------------------------
SUBROUTINE write_initial_report(iunartout, filout)
  use artn_params, ONLY: engine_units, ninit, nperp, neigen, nsmooth,  &
                         init_forc_thr, forc_thr, fpara_thr, eigval_thr, &
                         push_step_size, eigen_step_size, lanczos_max_size, lanczos_disp, &
                         push_mode, verbose, push_over, frelax_ene_thr, zseed, &
                         converge_property, lanczos_eval_conv_thr
  use units, only : unconvert_force, &
                    unconvert_energy, unconvert_hessian, unconvert_length, unit_char
  INTEGER,             INTENT(IN) :: iunartout
  CHARACTER (LEN=255), INTENT(IN) :: filout
  ! -- Local Variables
  INTEGER :: ios
  !
  ! Writes the header to the artn output file
  !
  OPEN( UNIT = iunartout, FILE = filout, FORM = 'formatted', STATUS = 'unknown', POSITION='rewind', IOSTAT = ios )
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
  WRITE (iunartout,'(15X,"converge_property = ", A)') converge_property
  WRITE (iunartout,'(15X,"init_forc_thr     = ", F6.3,2x,A)') unconvert_force( init_forc_thr ), unit_char('force')
  WRITE (iunartout,'(15X,"forc_thr          = ", F6.3,2x,A)') unconvert_force( forc_thr ), unit_char('force')
  WRITE (iunartout,'(15X,"fpara_thr         = ", F6.3,2x,A)') unconvert_force( fpara_thr ), unit_char('force')
  WRITE (iunartout,'(15X,"eigval_thr        = ", F6.3,2x,A)') unconvert_hessian( eigval_thr ), unit_char('hessian')
  WRITE (iunartout,'(15X,"frelax_ene_thr    = ", F6.3,2x,A)') unconvert_energy( frelax_ene_thr ), unit_char('energy')
  WRITE (iunartout,'(13X,"* Step size Parameter: ")')
  WRITE (iunartout,'(15X,"push_step_size  = ", F6.2,2x,A)') unconvert_length( push_step_size ), unit_char('length')
  WRITE (iunartout,'(15X,"eigen_step_size = ", F6.2,2x,A)') unconvert_length( eigen_step_size ), unit_char('length')
  WRITE (iunartout,'(15X,"push_over       = ", F6.3,2x,A)') push_over, "fraction of eigen_step_size"
  WRITE (iunartout,'(15X,"push_mode       = ", A6)') push_mode
  WRITE (iunartout,'(5X, "--------------------------------------------------")')
  WRITE (iunartout,'(5X, "Lanczos algorithm:")' )
  WRITE (iunartout,'(5X, "--------------------------------------------------")')
  WRITE (iunartout,'(15X,"lanczos_max_size   = ", I6)') lanczos_max_size
  WRITE (iunartout,'(15X,"lanczos_disp           = ", G11.4,2x,A)') unconvert_length( lanczos_disp ), unit_char('length')
  WRITE (iunartout,'(15X,"lanczos_eval_conv_thr   = ", G11.4)') lanczos_eval_conv_thr
  WRITE (iunartout,'(5X, "--------------------------------------------------")')
  WRITE (iunartout,'(/,/)') 
  !WRITE (iunartout,*) " "

  CLOSE ( UNIT = iunartout, STATUS = 'KEEP')

END SUBROUTINE write_initial_report




!------------------------------------------------------------
SUBROUTINE write_header_report( iunartout )
  use artn_params, only : verbose, isearch, ifound, filout
  use units, only :  strg_units
  implicit none

  integer, intent( in ) :: iunartout 
  INTEGER               :: ios
  !integer :: ios


  OPEN( UNIT = iunartout, FILE = filout, FORM = 'formatted', STATUS = 'unknown', POSITION='append', IOSTAT = ios )
  WRITE(iunartout,'(5x,"|> ARTn research :",2(x,i0)/,5x,*(a))') isearch, ifound, repeat("-",50)

  !%! Condition on the engin_units..
  select case( verbose )
    case( 0 )
    WRITE( iunartout,'(5X,"istep",4X,"ART_step",4X,"Etot",5x,"init/eign/perp/lanc/relx","&
               &"4X," Ftot ",5X," Fperp ",4X," Fpara ",4X,"eigval", 6X, "delr", 2X, "npart", X,"evalf",2X,"a1")')

    case( 1: )
    WRITE( iunartout,'(5X,"istep",4X,"ART_step",4X,"Etot",5x,"init/eign/perp/lanc/relx","&
               &"4X," Ftot ",5X," Fperp ",4X," Fpara ",4X,"eigval", 6X, "delr", 2X, "npart", X,"evalf","&
               &"4X,"B/O/R|I/P/L/E|P/B/R",4X,"a1")')
  end select

  ! -- Units
  WRITE( iunartout, strg_units )
  CLOSE (iunartout)
END SUBROUTINE write_header_report



!------------------------------------------------------------
SUBROUTINE write_report( etot, force, fperp, fpara, lowest_eigval, if_pos, istep, nat, iunartout)
  !
  !> @brief
  !!   a subroutine that writes a report of the current step to the output file  
  !
  !> @param [in]  etot		energy of the system
  !> @param [in]  force		List of atomic forces
  !> @param [in]  fpara		List of parallel atomic forces
  !> @param [in]  fperp		List of perpendicular atomic forces
  !> @param [in]  lowest_eigval	Lowest eigenvalue obtained by lanczos
  !> @param [in]  if_pos	Fix the atom or not	
  !> @param [in]  istep		actual step of ARTn 
  !> @param [in]  iunartout	Channel of output
  !
  USE artn_params, ONLY: MOVE, verbose, rcurv, bilan, filout, ismooth, nsmooth  &
                        ,etot_init, iinit, iperp, ieigen, ilanc, irelax, delr, verbose, iartn, a1 &
                        ,tau_init, lat, tau_step, delr, converge_property &
                        ,lrelax, linit, lbasin, lperp, llanczos, leigen, lpush_over, lpush_final, lbackward, lrestart 
  USE UNITS
  IMPLICIT NONE

  ! -- Arguments
  INTEGER,  INTENT(IN) :: nat, istep, iunartout
  INTEGER,  INTENT(IN) :: if_pos(3,nat)
  REAL(DP), INTENT(IN) :: force(3,nat),   &
                          fpara(3,nat),   &
                          fperp(3,nat)
  REAL(DP), INTENT(IN) :: etot, lowest_eigval

  ! -- Local Variables
  CHARACTER(LEN=5)     :: Mstep
  INTEGER              :: evalf, i, npart
  REAL(DP)             :: force_tot, fperp_tot, fpara_tot, detot, lowEig, dr, rc2
  REAL(DP)             :: ctot, cmax
  REAL(DP), EXTERNAL   :: ddot, dsum
  INTEGER              :: disp
  INTEGER              :: ios
  !
  !...Define the displacement type 
  disp =5
  IF ( lrelax     )    disp=6
  IF ( lpush_over )    disp=7
  IF ( lperp      )    disp=3  
  IF ( lperp    .AND. iperp==0 ) THEN  
                       disp=2
      IF (.NOT.lbasin) disp=4
  ENDIF    
  IF ( llanczos .AND. ilanc==0 ) THEN  
                       disp=2
     IF (.NOT.lbasin)  disp=4
  ENDIF    
  IF ( istep==0   )    disp=1
  IF ( disp == 4 .AND. ismooth <= nsmooth .AND. nsmooth>0)&
                       disp=8
  !
  ! ...Define when to print
  IF (( .NOT.(disp==1) ).AND. &
      ( .NOT.(disp==2) ).AND. &
      ( .NOT.(disp==4) ).AND. &
      ( .NOT.(disp==8) ).AND. &
      ( .NOT.((mod(irelax,5)==0) .AND. disp==6)) .AND.&  ! can be changed to print more during relax
      ( verbose<3 ) )  RETURN
  !
  IF( disp==2 .OR. disp==4 .OR. disp==8 .OR. (disp=5.AND. istep=1) )&
      iartn = iartn + 1
  !
  ! ...Force processing
  IF( trim(converge_property) == 'norm' )THEN
    call sum_force( force*if_pos, nat, force_tot )
    call sum_force( fpara, nat, fpara_tot )
    call sum_force( fperp, nat, fperp_tot )
    !force_tot = sqrt( dsum( 3*nat, force*if_pos ) )
    !fpara_tot = sqrt( dsum( 3*nat, fpara ) )
    !fperp_tot = sqrt( dsum( 3*nat, fperp ) )
  ELSE
    force_tot = MAXVAL( ABS(force*if_pos) )
    fperp_tot = MAXVAL( ABS(fperp) )
    fpara_tot = MAXVAL( ABS(fpara) )
  ENDIF
  !
  ! .. Conversion Units
  force_tot = unconvert_force( force_tot )
  fperp_tot = unconvert_force( fperp_tot )
  fpara_tot = unconvert_force( fpara_tot )
  dEtot     = unconvert_energy(etot - etot_init)
  lowEig    = unconvert_hessian( lowest_eigval )
  !
  !%! More Complete Output
  Mstep              = "Mstep"
  IF( lbasin ) Mstep = 'Bstep'
  IF( .NOT.lbasin ) Mstep = 'Sstep'
  IF( lrelax ) Mstep = 'Rstep'
  !
  !delr = sum()
  evalf = istep+1
  dr    = 0.
  npart = 0
  IF( disp==1 ) THEN
    !
    ! ...Initialize Displacement processing
    IF( .NOT.ALLOCATED(tau_init) ) THEN
        ALLOCATE( tau_init, source = tau_step )
    ELSE
        tau_init = tau_step    
    ENDIF
  ENDIF 
  !
  IF( disp==2 .OR. disp==4 .OR. disp==6 .OR. disp==8) THEN
    !  
    ! ...Displacement processing
    call compute_delr( nat, tau_step, tau_init, lat, delr )
    npart = 0
    rc2   = 0.1!*0.1  !! Miha: Why square? NS: Why not! 
    DO i = 1, nat
      IF( norm2(delr(:,i)) > rc2 ) npart = npart + 1
    enddo
    !! routine sum_force is equivalent to implicit: norm2( delr )
    call sum_force( delr, nat, dr )
    !
  ENDIF
  !
  ! ...Fill bilan variable for the inter report
  IF ( disp==2 .OR. disp==4 .OR. disp==8) THEN 
    ctot = dr
    cmax = real(npart,DP)
  ELSE
    ctot = bilan(7)
    cmax = bilan(6)
  ENDIF
  bilan = [ detot, force_tot, fpara_tot, fperp_tot, lowEig, cmax, ctot, real(evalf,DP) ]
  !
  !
  OPEN( UNIT = iunartout, FILE = filout, FORM = 'formatted', STATUS = 'unknown', POSITION='append', IOSTAT = ios )
  SELECT CASE( verbose )
    !
    CASE( 0 )
      WRITE(iunartout,6) iartn, Mstep, MOVE(disp), detot, iinit, ieigen, iperp, ilanc, irelax,  &
                         force_tot, fperp_tot, fpara_tot, lowEig, dr, npart, evalf, a1
      6 FORMAT(5x,i4,3x,a,x,a,F10.4,x,5(x,i4),5(x,f10.4),2(x,i5),3X,f4.2)
    !
    CASE( 1: )
      WRITE(iunartout,5) iartn, Mstep, MOVE(disp), detot, iinit, ieigen, iperp, ilanc, irelax,  &
                         force_tot, fperp_tot, fpara_tot, lowEig, dr, npart, evalf, lbasin,     &
                         lpush_over, lrelax, linit, lperp, llanczos, leigen,  lpush_final,      &
                         lbackward, lrestart , a1
      5 format(5x,i4,3x,a,x,a,F10.4,x,5(x,i4),5(x,f10.4),2(x,i5),3X,10(L2),3X,f4.2)
    ! 
  END SELECT
  CLOSE(iunartout)
  !
END SUBROUTINE write_report



!------------------------------------------------------------
SUBROUTINE write_inter_report( iunartout, pushfactor, de )
  use units, only : DP, unconvert_energy, unit_char
  use artn_params, only : artn_resume, istep, bilan,filout
  implicit none

  integer, intent( in )     :: iunartout             !> Ouput Unit 
  integer, intent( in )     :: pushfactor
  real(DP), intent( in )    :: de(*)        !> list of energies 
  character(:), allocatable :: Cbilan
  INTEGER                   :: ios

  OPEN( UNIT = iunartout, FILE = filout, FORM = 'formatted', STATUS = 'unknown', POSITION='append', IOSTAT = ios )
  SELECT CASE( pushfactor )

    CASE( 1 )
      ! de(1) = de_back
      WRITE( iunartout, '(5X, "--------------------------------------------------")')
      WRITE( iunartout, '(5X, "|> ARTn converged to a forward minimum | backward E_act =", F12.5,x,a)') &
          unconvert_energy( de(1) ), unit_char('energy')
      WRITE( iunartout, '(5X, "--------------------------------------------------")')
      WRITE (iunartout,'(5X, "|> number of steps:",x, i0)') istep

    CASE( -1 )
      ! de(1) = de_back
      ! de(2) = de_fwd 
      ! de(3) = etot_init 
      ! de(4) = etot_final
      ! de(5) = etot   
      WRITE( iunartout,'(5X, "--------------------------------------------------")')
      WRITE( iunartout,'(5X, "    *** ARTn converged to a backward minimum ***   ")')
      WRITE( iunartout,'(5X, "--------------------------------------------------")')
      WRITE( iunartout,'(15X,"forward  E_act =", F12.5,x,a)') unconvert_energy(de(2)), unit_char('energy')
      WRITE( iunartout,'(15X,"backward E_act =", F12.5,x,a)') unconvert_energy(de(1)), unit_char('energy') 
      WRITE( iunartout,'(15X,"reaction dE    =", F12.5,x,a)') unconvert_energy((de(5)-de(4))), unit_char('energy')
      WRITE( iunartout,'(15X,"dEinit - dEfinal    =", F12.5,x,a)') unconvert_energy((de(3)-de(5))), unit_char('energy') 
      WRITE( iunartout,'(5X, "--------------------------------------------------")')
      WRITE( iunartout,'(5X, "|> Configuration Files:", X,A)') trim(artn_resume)
      WRITE( iunartout,'(5x, "|> Configuration Files:", X,A)') trim(artn_resume)
      WRITE( iunartout,'(5X, "--------------------------------------------------")')
      !WRITE( u,'(/)')


    CASE DEFAULT
      WRITE(*,*) "********* ERROR write_inter_report:: pushfactor", pushfactor, " **************"
      

  END SELECT
    Cbilan = '(5x,"|> DEBRIEF | dE= ",f12.5,x,"'//unit_char('energy')//' | F_{tot,para,perp}= ",3(f12.5,x),"' &
     //unit_char('force')//' | EigenVal= ", f12.5,x,"'//unit_char('hessian')//' | npart= ",f4.0,x," | delr= ",f12.5,x,"' &
     //unit_char('length')//' | evalf= ",f5.0,x,"|")'
    Write(iunartout,Cbilan) Bilan
  CLOSE(iunartout)

END SUBROUTINE write_inter_report



!------------------------------------------------------------
SUBROUTINE write_end_report( iunartout, lsaddle, lpush_final, de )
 
  use units, only : DP, unconvert_energy, unit_char
  use artn_params, only : artn_resume, istep, bilan,filout
  implicit none

  integer, intent( in ) :: iunartout
  logical, intent( in ) :: lsaddle, lpush_final
  REAL(DP), intent( in ), value :: de
  character(:), allocatable :: Cbilan
  INTEGER                   :: ios
  
  OPEN( UNIT = iunartout, FILE = filout, FORM = 'formatted', STATUS = 'unknown', POSITION='append', IOSTAT = ios )
  if( lsaddle )then

    WRITE (iunartout,'(5X, "--------------------------------------------------")')
    WRITE (iunartout,'(5X, "|> ARTn found a potential saddle point | E_saddle - E_initial =", F12.5,x,a)') &
      unconvert_energy(de), unit_char('energy')
    WRITE(iunartout,'(5X, "|> Stored in Configuration Files:", X,A)') trim(artn_resume)
    WRITE (iunartout,'(5X, "--------------------------------------------------")')

    Cbilan = '(5x,"|> DEBRIEF | dE= ",f12.5,x,"'//unit_char('energy')//' | F_{tot,para,perp}= ",3(f12.5,x),"' &
        //unit_char('force')// &
        ' | EigenVal= ", f12.5,x,"'//unit_char('hessian')//' | npart= ",f4.0,x," | delr= ",f12.5,x,"'//unit_char('length')// &
        ' | evalf= ",f5.0,x,"|")'
    Write(iunartout,Cbilan) Bilan


    IF( lpush_final ) THEN
      WRITE(iunartout,'(5X,"       *** Pushing forward to a minimum  ***      ")')
      WRITE(iunartout,'(5X, "-------------------------------------------------")')
    ELSE
      WRITE(iunartout,'(5X,"|> No push_final to Minimum :: ARTn search finished "/5x,*(a))') repeat("-",50)
      !WRITE(iunartout,'(5X,"       *** ARTn search finished ***")')
      !WRITE(iunartout,'(5X,"       *** no push_final minimal ***")')
      WRITE(iunartout,'(5X, "-------------------------------------------------"/)')
    ENDIF

  else
    WRITE (iunartout,'(5X, "--------------------------------------------------")')
    WRITE (iunartout,'(5X, "        *** ARTn saddle search failed  ***        ")')
    WRITE (iunartout,'(5X, "--------------------------------------------------")')
  endif

  WRITE (iunartout,'(5X, "|> number of steps:",x, i0)') istep
  CLOSE(iunartout)

END SUBROUTINE write_end_report


!------------------------------------------------------------
SUBROUTINE write_fail_report( iunartout, disp, estep )

  use units, only : DP, unconvert_energy, unit_char
  use artn_params, only : MOVE, ifails, error_message,filout
  implicit none

  integer, intent( in ) :: iunartout, disp
  REAL(DP), intent( in ):: estep
  INTEGER               :: ios

  ifails = ifails + 1

  OPEN  (UNIT = iunartout, FILE = filout, FORM = 'formatted', STATUS = 'unknown', POSITION='append', IOSTAT = ios )
  WRITE (iunartout,'(5X, "--------------------------------------------------")')
  WRITE (iunartout,'(5X, "        *** ARTn search failed ( ",i0," ) at ",a," *** ")') ifails, MOVE(DISP)
  WRITE (iunartout,'(5X, "Step Params: Etot = ",f10.4,x,a)') unconvert_energy(estep), unit_char('energy')
  WRITE (iunartout, '(5X, "Failure message: ",a)') trim(adjustl(error_message))
  WRITE (iunartout,'(5X, "--------------------------------------------------"//)')
  CLOSE (iunartout)
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


SUBROUTINE check_force_convergence( nat, force, if_pos, fperp, fpara, lforc_conv, lsaddle_conv ) 
  !
  !> @breif 
  !!   A subroutine that checks the force convergence of a particular step in the artn algorithm 
  !
  !> @param [in]   nat             Size of list: number of atoms
  !> @param [in]   force           Force field
  !> @param [in]   if_pos          List of atom move or not
  !> @param [in]   fperp           Perpendicular Force Field
  !> @param [in]   fpara           Parallel Force Field
  !> @param [out]  lforc_conv      Force Convergence Flag
  !> @param [out]  lsaddle_conv    Saddle-point Convergence Flag
  !
  USE units
  USE artn_params, ONLY : linit, leigen, llanczos, lperp, lrelax, &
                          ilanc, iperp, nperp, nperp_step, noperp, istep, INIT, EIGN, RELX, &
                          init_forc_thr, forc_thr, fpara_thr, tau_step, rcurv,  &
                          lowest_eigval, iunartout, restartfname, etot_step, write_restart, warning, converge_property
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: force(3,nat)
  REAL(DP), INTENT(IN) :: fperp(3,nat)
  REAL(DP), INTENT(IN) :: fpara(3,nat)
  INTEGER, INTENT(IN) :: if_pos(3,nat)
  INTEGER, INTENT(IN) :: nat
  !INTEGER, INTENT(IN) :: order(nat)
  REAL(DP) :: fperp_thr
  LOGICAL, INTENT(OUT) :: lforc_conv, lsaddle_conv
  !
  LOGICAL :: C1, C2, C3, C4
  LOGICAL, parameter :: ARTnStep = .true.
  integer :: ios
  REAL(DP) :: maxforce, maxfperp, maxfpara

  real(DP), external :: dsum

  lforc_conv = .false.
  lsaddle_conv = .false.
  !

  ! ...Open the output file
  OPEN( UNIT = iunartout, FILE = 'artn.out', FORM = 'formatted', ACCESS = 'append', STATUS = 'unknown', IOSTAT = ios )


  ! ...Compute the variable
  IF( trim(converge_property) == 'norm' )THEN
    call sum_force( force*if_pos, nat, maxforce )
    call sum_force( fpara, nat, maxfpara )
    call sum_force( fperp, nat, maxfperp )
    !maxforce = sqrt( dsum( 3*nat, force*if_pos ) )
    !maxfpara = sqrt( dsum( 3*nat, fpara ) )
    !maxfperp = sqrt( dsum( 3*nat, fperp ) )
 
  ELSE
    maxforce = MAXVAL(ABS(force*if_pos))
    maxfpara = MAXVAL(ABS(fpara))
    maxfperp = MAXVAL(ABS(fperp))
  ENDIF

  

  !
  IF ( lperp ) THEN 
     !
     ! ...Compute Force evolution
     !call compute_curve( iperp, 3*nat, tau_step, fperp )

     IF ( leigen ) THEN 

        !
        ! We are outside of the basin check both perpendicular and total force convergence
        ! 
        IF( maxforce < forc_thr  ) THEN
           lsaddle_conv = .true.

           !CALL write_restart( restartfname, nat )
           CALL write_restart( restartfname )
           CALL write_report( etot_step, force, fperp, fpara, lowest_eigval, EIGN, if_pos, istep, nat,  iunartout, ArtnStep )
           RETURN
        ENDIF



        !
        ! check whether the fperp criterion should be tightened
        ! 
        IF( maxfpara <= fpara_thr )THEN !> We are close to the saddle point 
           fperp_thr = forc_thr

           ! ...Perp-relax managment
           CALL nperp_limitation_step( 0 )

        ELSE

           fperp_thr = init_forc_thr

           ! ...Perp-relax managment
           CALL nperp_limitation_step( 0 )
        ENDIF

        !print*, " CHECK_FORCE():Nperp ", nperp, nperp_step, nperp_list

        ! 
        ! -- check perpendicular force convergence for the perp-relax 
        ! 
        C1 = ( maxfperp < fperp_thr ) ! check on the fperp field
        C2 = ( nperp > 0.AND.iperp >= nperp )    ! check on the perp-relax iteration
        C3 = ( MAXfperp < MAXfpara ) ! check wheter fperp is lower than fpara
        !C4 = ( rcurv > 0.5_DP )

        IF( C1 .and. iperp == 0 )C1 = .false.  !! Force to do at least one prep-relax
        !IF( C1.and. ABS(maxfperp - maxfpara) < maxfpara*1.20 ) C1 = .false.  !! if fperp is to o far from fpara too many perp-relax can relax to much the system

        !IF( C1 .OR. C2 .OR. C3 .OR. C4 )THEN
        IF( C1 .OR. C2 .OR. C3 )THEN
           lperp = .false.
           llanczos = .true. 
           leigen = .false. 
           IF (C1) write(iunartout,*)"Stop perp relax because fperp<fperp_thr ",&
               unconvert_force( maxfperp ), unconvert_force(fperp_thr), TRIM(converge_property) 
           IF (C2) write(iunartout,*)"Stop perp relax because iperp>nperp_thr ",&
               iperp, nperp 
           IF (C3) write(iunartout,*)"Stop perp relax because fperp<fpara ",&
               unconvert_force( maxfperp ), unconvert_force( maxfpara ), TRIM(converge_property)

           CALL write_restart( restartfname )
           CALL write_report( etot_step, force, fperp, fpara, lowest_eigval, EIGN, if_pos, istep, nat,  iunartout, ArtnStep )
           ilanc = 0

           ! ...Perp-Relax is finshed, 
           !   we count the nperp_step
           CALL nperp_limitation_step( 1 )
           !WRITE( iunartout,* ) "* NEXT NPERP ",nperp, nperp_step

        ENDIF


        ! ...Count if the fperp is always to small
        IF( C1.AND.iperp == 0 )THEN
        !  noperp = noperp + 1
        !  ! ** WARNING **
        !  if( noperp > 2 ) &
            CALL WARNING( iunartout, "No PERP-RELAX",  &
                 "Fperp is aready lower than fperp_thr", [fperp_thr])
        ENDIF




     ELSE ! ...IN  BASIN
        !
        ! In the Basin after perp-relax we always return to init mode
        !

        fperp_thr = init_forc_thr 

        ! ...Do INIT until fperp is > fperp_thr
        ! && Max Iteration of Perp-Relax
        CALL nperp_limitation_step( -1 )

        C1 = ( MAXfperp < fperp_thr ) ! check on the fperp field
        C2 = ( nperp > -1 .AND.iperp >= nperp )    ! check on the perp-relax iteration

        IF( C1 .OR. C2 )THEN
           lperp = .false.
           linit = .true.
           !write(iunartout,*)"Check_force::basin->Fperp", C1, C2
           CALL write_restart( restartfname )
           CALL write_report( etot_step, force, fperp, fpara, lowest_eigval, INIT, if_pos, istep, nat,  iunartout, ArtnStep )
        ENDIF

        ! ...Count if the fperp is always to small
        IF( C1 )THEN
          noperp = noperp + 1
          ! ** WARNING **
          if( noperp > 2 ) &
            CALL WARNING( iunartout, "Tansition Push-Init->Perp-Relax",  &
                 "The Fperp is too small after Push-INIT- change push_step_size ", [noperp])
        ENDIF


     ENDIF



     ! 
  ELSE IF ( lrelax ) THEN  
     !
     IF( MAXforce < forc_thr  )then
       lforc_conv = .true.
       !write(iunartout,*)"Check_force::Relax->Force"
       CALL write_report( etot_step, force, fperp, fpara, lowest_eigval, RELX, if_pos, istep, nat,  iunartout, ArtnStep )
     ENDIF
     !
  ENDIF


  close( iunartout )
     
ENDSUBROUTINE check_force_convergence





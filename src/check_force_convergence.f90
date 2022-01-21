
SUBROUTINE check_force_convergence( nat, force, if_pos, fperp, fpara, lforc_conv, lsaddle_conv) 
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
  USE artn_params, ONLY : linit, lbasin, leigen, llanczos, lperp, lrelax, &
                          ilanc, iperp, nperp, nperp_list, nperp_step, noperp, istep, INIT, PERP, EIGN, LANC, RELX, &
                          init_forc_thr, forc_thr, fpara_thr, push, &
                          lowest_eigval, iunartout, restartfname, etot_step, write_restart, warning
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
  LOGICAL :: C1, C2
  LOGICAL, parameter :: ARTnStep = .true.
  integer :: ios

  lforc_conv = .false.
  lsaddle_conv = .false.
  !

  OPEN ( UNIT = iunartout, FILE = 'artn.out', FORM = 'formatted', ACCESS = 'append', STATUS = 'unknown', IOSTAT = ios )

  !
  IF ( lperp ) THEN 
     !
     IF ( leigen ) THEN 

        !
        ! We are outside of the basin check both perpendicular and total force convergence
        ! 
        IF (MAXVAL( ABS(force*if_pos)) < forc_thr  ) THEN
           lsaddle_conv = .true.

           CALL write_restart( restartfname, nat )
           CALL write_report( etot_step, force, fperp, fpara, lowest_eigval, EIGN, if_pos, istep, nat,  iunartout, ArtnStep )
           RETURN
        ENDIF

        !
        ! check whether the fperp criterion should be tightened
        ! 
        IF ( MAXVAL(ABS(fpara)) <= fpara_thr ) THEN !> We are close to the saddle point 
           fperp_thr = forc_thr

           ! ...Perp-relax managment
           !if( nperp_step == 1 )nperp_list(1) = nperp
           !nperp_step = nperp_step + 1
           select case( nperp_step )
             case(:4); nperp = nperp_list( nperp_step )
             case(5:); nperp = nperp_list( 5 )
           end select
           !nperp = 0    ! Remove the Perp-Relax Iteration constrain

        ELSE

           fperp_thr = init_forc_thr

           ! ...Perp-relax managment
           !if( nperp_step == 1 )nperp_list(1) = nperp
           !nperp_step = nperp_step + 1
           select case( nperp_step )
             case(:4); nperp = nperp_list( nperp_step )
             case(5:); nperp = nperp_list( 5 )
           end select
           !nperp =  nperp_list(1)
        ENDIF

        ! 
        ! -- check perpendicular force convergence 
        ! 
        C1 = ( MAXVAL( ABS(fperp)) < fperp_thr ) ! check on the fperp field
        C2 = ( nperp > 0.AND.iperp >= nperp )    ! check on the perp-relax iteration

        IF( C1 .OR. C2  ) THEN
           lperp = .false.
           llanczos = .true. 
           leigen = .false. 

           CALL write_restart( restartfname, nat )
           CALL write_report( etot_step, force, fperp, fpara, lowest_eigval, EIGN, if_pos, istep, nat,  iunartout, ArtnStep )
           ilanc = 0

           ! ...Perp-Relax is finshed, 
           !   we count the nperp_step
           if( nperp_step == 1 )nperp_list(1) = nperp
           nperp_step = nperp_step + 1
           select case( nperp_step )
             case(:4); nperp = nperp_list( nperp_step )
             case(5:); nperp = nperp_list( 5 )
           end select
           WRITE( iunartout,* ) "* NEXT NPERP ",nperp, nperp_step

        ENDIF




     ELSE ! ...IN  BASIN
        !
        ! In the Basin after perp-relax we always return to init mode
        !

        fperp_thr = init_forc_thr 

        ! ...Do INIT until fperp is > fperp_thr
        ! && Max Iteration of Perp-Relax
        !    nperp = 0 means no nperp constrain

        C1 = ( MAXVAL( ABS(fperp)) < fperp_thr ) ! check on the fperp field
        C2 = ( nperp > 0.AND.iperp >= nperp )    ! check on the perp-relax iteration

        IF( C1 .OR. C2 ) THEN
           lperp = .false.
           linit = .true.
           !write(iunartout,*)"Check_force::basin->Fperp"
           CALL write_restart( restartfname, nat )
           CALL write_report( etot_step, force, fperp, fpara, lowest_eigval, INIT, if_pos, istep, nat,  iunartout, ArtnStep )
        ENDIF

        ! ...Count if the fperp is always to slow
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
     IF( MAXVAL( ABS(force*if_pos)) < forc_thr  )then
       lforc_conv = .true.
       !write(iunartout,*)"Check_force::Relax->Force"
       CALL write_report( etot_step, force, fperp, fpara, lowest_eigval, RELX, if_pos, istep, nat,  iunartout, ArtnStep )
     ENDIF
     !
  ENDIF


  close( iunartout )
     
ENDSUBROUTINE check_force_convergence

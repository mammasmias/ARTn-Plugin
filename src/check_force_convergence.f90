
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
                          ilanc, iperp, nperp, nperp_step, noperp, istep, INIT, EIGN, RELX, iperp_save, &
                          init_forc_thr, forc_thr, fpara_thr, tau_step, rcurv, verbose, iinit, ninit,&
                          lowest_eigval, iunartout, restartfname, etot_step, write_restart, warning, converge_property
  IMPLICIT NONE
  REAL(DP), INTENT(IN)  :: force(3,nat)
  REAL(DP), INTENT(IN)  :: fperp(3,nat)
  REAL(DP), INTENT(IN)  :: fpara(3,nat)
  INTEGER,  INTENT(IN)  :: if_pos(3,nat)
  INTEGER,  INTENT(IN)  :: nat
  !INTEGER, INTENT(IN)  :: order(nat)
  REAL(DP) :: fperp_thr
  LOGICAL,  INTENT(OUT) :: lforc_conv, lsaddle_conv
  !
  ! Local Variables
  LOGICAL               :: C0,C1, C2, C3, C4
  integer               :: ios
  REAL(DP)              :: maxforce, maxfperp, maxfpara
  real(DP), external    :: dsum
  !
  C0           = .false.
  C1           = .false.
  C2           = .false.
  C3           = .false.
  C4           = .false.
  lforc_conv   = .false.
  lsaddle_conv = .false.
  !
  ! ...Compute the variable
  IF( trim(converge_property) == 'norm' )THEN
    call sum_force( force*if_pos, nat, maxforce )
    call sum_force( fpara, nat, maxfpara )
    call sum_force( fperp, nat, maxfperp )
  ELSE
    maxforce = MAXVAL(ABS(force*if_pos))
    maxfpara = MAXVAL(ABS(fpara))
    maxfperp = MAXVAL(ABS(fperp))
  ENDIF
  !
  IF ( lperp ) THEN
     !
     ! ...Compute Force evolution
     IF ( leigen ) THEN ! ... NOT IN BASIN
        !
        ! We are outside of the basin
        !
        ! ... Is the system converged to saddle?
        C0 = ( maxforce < forc_thr )
        IF( C0  ) THEN
           lsaddle_conv = .true.
           CALL write_restart( restartfname )
           !RETURN
        ENDIF
        !
        ! ... Ccheck whether the fperp criterion should be tightened
        IF( maxfpara <= fpara_thr )THEN !> We are close to the saddle point 
           fperp_thr = forc_thr
           ! ...Perp-relax managment
           CALL nperp_limitation_step( 0 )
        ELSE
           fperp_thr = init_forc_thr
           ! ...Perp-relax managment
           CALL nperp_limitation_step( 0 )
        ENDIF
        ! 
        ! ... Check perpendicular force convergence for the perp-relax 
        C1 = ( maxfperp < fperp_thr )         ! check on the fperp field
        C2 = ( nperp > 0.AND.iperp >= nperp ) ! check on the perp-relax iteration
        C3 = ( MAXfperp < MAXfpara )          ! check wheter fperp is lower than fpara
        !C4 = ( rcurv > 0.5_DP )
        ! 
        IF( C1 .and. iperp == 0 ) C1 = .false. ! Force to do at least one prep-relax. AJ Why?
        !
        ! ... If fperp is to o far from fpara too many perp-relax can relax to much the system
        !IF( C1.and. ABS(maxfperp - maxfpara) < maxfpara*1.20 ) C1 = .false.
        !
        IF( C1 .OR. C2 .OR. C3 ) THEN
           lperp    = .false.
           llanczos = .true. 
           leigen   = .false. 
           ilanc    = 0
           !
           CALL write_restart( restartfname )
           !
           ! ...Perp-Relax is finshed,  we count the nperp_step
           CALL nperp_limitation_step( 1 )
           !
        ENDIF
        !
        !
     ELSE ! ... IN  BASIN
        !
        ! In the Basin after perp-relax we always return to init mode
        !
        fperp_thr = init_forc_thr 
        !
        ! ...Do INIT until fperp is > fperp_thr
        ! && Max Iteration of Perp-Relax
        CALL nperp_limitation_step( -1 )
        ! 
        C1 = ( MAXfperp < fperp_thr )           ! check on the fperp field
        C2 = ( nperp > -1 .AND.iperp >= nperp ) ! check on the perp-relax iteration
        !
        ! ... After perp-relax return to lanczos or init if we are still close to the minimum
        IF( C1 .OR. C2 )THEN
          IF( iinit<ninit ) THEN 
            lperp    = .false.
            linit    = .true.
          ELSE
            lperp    = .false.
            llanczos = .true.
            linit    = .false.
            ilanc    = 0
          ENDIF    
           CALL write_restart( restartfname )
        ENDIF
        !
        ! ...Count if fperp is always to small after each init push
        IF ( C1 .AND. iperp == 0) noperp = noperp + 1
        !
     ENDIF
     !
     IF (verbose >1) THEN
        !    
        ! ... Show Stop perp message
        OPEN( UNIT = iunartout, FILE = 'artn.out', FORM = 'formatted', ACCESS = 'append', STATUS = 'unknown', IOSTAT = ios )
        IF ( C0 ) WRITE(iunartout,'(5x,a46,x,f10.4,x,a1,x,f10.4,a20)') &
            "|> Stop perp relax because force < forc_thr  :",&
            unconvert_force( maxforce ),"<", unconvert_force(forc_thr), TRIM(converge_property)

        IF ( C1 .AND. iperp>0) WRITE(iunartout,'(5x,a46,x,f10.4,x,a1,x,f10.4,a20)') &
            "|> Stop perp relax because fperp < fperp_thr :",&
            unconvert_force( maxfperp ),"<", unconvert_force(fperp_thr), TRIM(converge_property)
        !
        IF ( C2 ) WRITE(iunartout,'(5x,a46,x,i3,a1,i3)') &
            "|> Stop perp relax because iperp = nperp max :",&
            iperp,"=",nperp 
        !
        IF ( C3 ) WRITE(iunartout,'(5x,a46,x,f10.4,x,a1,x,f10.4,a20)') &
            "|> Stop perp relax because fperp < fpara     :",&
            unconvert_force( maxfperp ),"<", unconvert_force( maxfpara ), TRIM(converge_property)
        !
        IF ( C1 .AND. iperp == 0) WRITE(iunartout,'(5x,a46,x,f10.4,x,a1,x,f10.4,a20)') &
            "|> No perp relax because fperp < fpara       :",&
            unconvert_force( maxfperp ),"<", unconvert_force( maxfpara ), TRIM(converge_property)
        !
        IF ( noperp>2 ) WRITE(iunartout,'(5x,a90)') &
            "|> WARNING -The Fperp is too small after each Push-INIT- You should increase push_step_size"
        CLOSE( iunartout )
        !
     ENDIF   
     !
     !... Initialize couter if perp relax is finished
     IF (.NOT. lperp ) THEN
         iperp_save = iperp  
         iperp      = 0
     ENDIF    
     ! 
  ELSE IF ( lrelax ) THEN  
     !
     ! ... Check if Minimum has been reached
     C0 = ( maxforce < forc_thr )
     IF ( C0  ) THEN
        lforc_conv = .true.
        IF ( verbose >1) THEN
           !  
           ! ... Show Stop relax message
           OPEN( UNIT = iunartout, FILE = 'artn.out', FORM = 'formatted', ACCESS = 'append', STATUS = 'unknown', IOSTAT = ios )
           WRITE(iunartout,'(5x,a46,x,f10.4,x,a1,x,f10.4,a20)') &
           "|> Stop relax because force < forc_thr       :",&
           unconvert_force( MAXforce ),"<", unconvert_force( forc_thr ), TRIM(converge_property)
           CLOSE( iunartout )
        ENDIF  
     ENDIF
     !
  ENDIF
  !
ENDSUBROUTINE check_force_convergence

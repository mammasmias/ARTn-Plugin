
!> @author
!!   Matic Poberznik
!!   Miha Gunde
!!   Nicolas Salles
!!   Antoine Jay


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
  USE units, ONLY : DP, unconvert_force
  USE artn_params, ONLY : linit, leigen, llanczos, lperp, lrelax, lbasin, nperp_step, nperp_limitation,&
                          ilanc, iperp, nperp, nperp_step, noperp, istep, iperp_save, &
                          init_forc_thr, forc_thr, fpara_thr, verbose, iinit, ninit,&
                          lowest_eigval, iunartout, restartfname, etot_step, write_restart, warning,   &
                          converge_property, ismooth, nsmooth
  IMPLICIT NONE
  REAL(DP), INTENT(IN)  :: force(3,nat)
  REAL(DP), INTENT(IN)  :: fperp(3,nat)
  REAL(DP), INTENT(IN)  :: fpara(3,nat)
  INTEGER,  INTENT(IN)  :: if_pos(3,nat)
  INTEGER,  INTENT(IN)  :: nat
  !INTEGER, INTENT(IN)  :: order(nat)
  LOGICAL,  INTENT(OUT) :: lforc_conv, lsaddle_conv
  !
  ! Local Variables
  LOGICAL               :: C0,C1, C2, C3, C4
  integer               :: ios
  REAL(DP)              :: fperp_thr
  REAL(DP)              :: maxforce, maxfperp, maxfpara
  !REAL(DP)              :: min_dir(3,nat)
  real(DP), external    :: dsum, ddot
  logical, external     :: fperp_min_alignment
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
        ! ... Is the system converged to saddle?
        C0 = ( maxforce < forc_thr )
        IF( C0 ) THEN
           lsaddle_conv = .true.
           iperp_save = iperp  !! save iperp before the write_report()
           CALL write_restart( restartfname )
           CALL write_ARTn_step_report( etot_step, force, fperp, fpara, lowest_eigval, if_pos, istep, nat,  iunartout )
           RETURN  !! ANTOINE you removed this line !!!
        ENDIF

        
        ! ... Check whether the fperp criterion should be tightened
        fperp_thr = forc_thr  !! Should be removed

        ! 
        ! ... Conditions for stopping perp_relax
        C1 = ( maxfperp < fperp_thr )          ! check on the fperp field
        C2 = ( nperp > 0.AND.iperp >= nperp )  ! check on the perp-relax iteration
        C3 = ( MAXfperp < MAXfpara )           ! check wheter fperp is lower than fpara

        IF( C3 .and. iperp == 0 ) C1 = .false. ! Force to do at least one prep-relax.
        IF( nsmooth > 0 .AND. ismooth <= nsmooth )C3 = .False.  ! Force to do a perp relax during the smooth step

        !
        ! ...Alignment between fperp and direction of minimum
        C4 = fperp_min_alignment( 0.8_DP, 0.1_DP )
        !min_dir = tau_step - tau_init 
        !min_dir = min_dir / NORM2( min_dir )
        !dtmp = ddot(3*nat,min_dir,1,push,1)
        !! IF eigenVec change suddenlly AND direction of minimum is perp to the last push
        !C4 = ( a1 < 0.8 .AND. ABS(dtmp) < 0.1 )   

        !
        ! ...Stopping condition is filled, switch to lanczos
        IF( C1 .OR. C2 .OR. C3 .OR. C4 )THEN
           lperp    = .false.
           llanczos = .true. 
           leigen   = .false. 
           ilanc    = 0
           iperp_save = iperp  !! save iperp before the write_report()
           !
           CALL write_restart( restartfname )
           CALL write_ARTn_step_report( etot_step, force, fperp, fpara, lowest_eigval, if_pos, istep, nat,  iunartout )
           !
        ENDIF
        !
        !
     ELSE ! ... IN  BASIN
        !
        !fperp_thr = init_forc_thr 
        !
        ! ... Conditions for stopping perp_relax
        !C1 = ( MAXfperp < fperp_thr )           ! check on the fperp field   !! NS: NO C1 In the BASIN
        C2 = ( nperp > -1 .AND.iperp >= nperp ) ! check on the perp-relax iteration
        !
        ! ... Stopping condition is filled, switch to lanczos or to init if we are still close to the minimum
        !IF( C1 .OR. C2 )THEN
        IF( C2 )THEN
          IF( iinit < ninit ) THEN 
            lperp    = .false.
            linit    = .true.
          ELSE
            lperp    = .false.
            llanczos = .true.
            linit    = .false.
            ilanc    = 0
          ENDIF    
          iperp_save = iperp  !! save iperp before the write_report()
          !
          CALL write_restart( restartfname )
          CALL write_ARTn_step_report( etot_step, force, fperp, fpara, lowest_eigval, if_pos, istep, nat,  iunartout )
        ENDIF
        !
        ! ...Count if fperp is always to small after each init push
        IF ( C1 .AND. iperp == 0) noperp = noperp + 1
        !
     ENDIF

     !
     !    
     ! ... Show Stop perp message
     IF( verbose >1 )THEN
        OPEN( UNIT = iunartout, FILE = 'artn.out', FORM = 'formatted', ACCESS = 'append', STATUS = 'unknown', IOSTAT = ios )
        IF ( C0 ) WRITE(iunartout,'(5x,a46,x,f10.4,x,a1,x,f10.4,a20)') &
            "|> Stop perp relax because force < forc_thr  :",&
            unconvert_force( maxforce ),"<", unconvert_force(forc_thr), TRIM(converge_property)

        !IF ( C1 .AND. iperp>0) WRITE(iunartout,'(5x,a46,x,f10.4,x,a1,x,f10.4,a20)') &
        !! NS: Remove the condition on iperp 
        IF ( C1 ) WRITE(iunartout,'(5x,a46,x,f10.4,x,a1,x,f10.4,a20)') &
            "|> Stop perp relax because fperp < fperp_thr :",&
            unconvert_force( maxfperp ),"<", unconvert_force(fperp_thr), TRIM(converge_property)
        
        IF ( C2 ) WRITE(iunartout,'(5x,a46,x,i3,a1,i3)') &
            "|> Stop perp relax because iperp = nperp max :",&
            iperp,"=",nperp 
        !
        IF ( C3 ) WRITE(iunartout,'(5x,a46,x,f10.4,x,a1,x,f10.4,a20)') &
            "|> Stop perp relax because fperp < fpara     :",&
            unconvert_force( maxfperp ),"<", unconvert_force( maxfpara ), TRIM(converge_property)
        !
        IF ( C3 .AND. iperp == 0) WRITE(iunartout,'(5x,a46,x,f10.4,x,a1,x,f10.4,a20)') &
            "|> No perp relax because fperp < fpara       :",&
            unconvert_force( maxfperp ),"<", unconvert_force( maxfpara ), TRIM(converge_property)
        !
        IF ( C4 ) WRITE(iunartout,'(5x,a46)') &
            "|> No perp relax because C4      :"
        !
        IF ( noperp > 2 ) WRITE(iunartout,'(5x,a90)') &
            "|> WARNING -The Fperp is too small after each Push-INIT- You should increase push_step_size"
        CLOSE( iunartout )
        !
     ENDIF   

     !
     !... If perp relax is finished: update counter and update number of allowed perp_relax steps
     IF (.NOT. lperp ) THEN
         !iperp_save = iperp  
         iperp      = 0
         IF ( .NOT. lbasin) THEN
            nperp_step = nperp_step + 1
            nperp = nperp_limitation(MIN(SIZE(nperp_limitation), nperp_step))
         ELSE   
            nperp = nperp_limitation(1) 
         ENDIF   
     ENDIF    
     ! 


  ELSE IF ( lrelax ) THEN  
     !
     ! ... Check if Minimum has been reached
     C0 = ( maxforce < forc_thr )
     IF ( C0 ) THEN
        lforc_conv = .true.
        CALL write_ARTn_step_report( etot_step, force, fperp, fpara, lowest_eigval, if_pos, istep, nat,  iunartout )
        !  
        ! ... Show Stop relax message
        IF( verbose > 1 )THEN
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
END SUBROUTINE check_force_convergence








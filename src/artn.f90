!
!> @author
!!   Matic Poberznik
!!   Miha Gunde
!!   Nicolas Salles
!
!> @brief
!!   Main ARTn plugin subroutine:
!
!> @details
!!   Modifies the input force to perform the ARTn algorithm
!------------------------------------------------------------------------------
SUBROUTINE artn( force, etot_eng, nat, ityp, atm, tau, order, at, if_pos, disp, displ_vec, lconv )
  !----------------------------------------------------------------------------
  !> @param[in]     force       force calculated by the engine
  !> @param[inout]  etot_eng    total energy of the engine
  !> @param[in]     nat         number of atoms
  !> @param[in]     ityp        list of type of atoms
  !> @param[in]     atm         list of the element's name relative to the atomic type
  !> @param[inout]  tau         atomic position
  !> @param[in]     order       order of atomic index in the list: force, tau, ityp
  !> @param[in]     at          lattice parameter
  !> @param[in]     if_pos      list of fixed atomic dof (0 or 1)
  !> @param[out]    disp        stage for move_mode
  !> @param[out]    displ_vec   displacement vector communicated to move_mode
  !> @param[out]    lconv       flag for controlling convergence
  !
  ! artn_params for variables and counters that need to be stored in each step
  ! DEFINED IN: artn_params_mod.f90
  !
  USE units
  USE artn_params, ONLY: iunartin, iunartout, iunstruct, &
       lrelax, linit, lperp, leigen, llanczos, lrestart, lbasin, lpush_over, lpush_final, lbackward, lmove_nextmin,  &
       irelax, istep, iperp, ieigen, iinit, ilanc, ismooth, iover, isearch, ifound, nlanc, nperp, noperp, nperp_step,  &
       if_pos_ct, lowest_eigval, etot_init, etot_step, etot_saddle, etot_final, de_back, de_fwd, &
       ninit, neigen, lanczos_max_size, nsmooth, push_mode, dist_thr, init_forc_thr, forc_thr, &
       fpara_thr, eigval_thr, frelax_ene_thr, push_step_size, current_step_size, eigen_step_size, fpush_factor, &
       push_ids, add_const, push, eigenvec, tau_step, force_step, tau_init, tau_saddle, eigen_saddle, v_in, &
       VOID, INIT, PERP, EIGN, LANC, RELX, OVER, zseed, &
       engine_units, struc_format_out, elements, ilanc_save, &
       setup_artn, read_restart, write_restart, inewchance, nnewchance,&
       push_over, ran3, a1, old_lanczos_vec, lend, fill_param_step, &
       filin, filout, sadfname, initpfname, eigenfname, restartfname, warning, flag_false,  &
       prefix_min, nmin, prefix_sad, nsaddle, artn_resume, natoms, old_lowest_eigval, &
       lanczos_always_random, etot_diff_limit, error_message, prev_push, SMTH
  !
  IMPLICIT NONE

  ! -- ARGUMENTS
  INTEGER, value,   INTENT(IN)    :: nat              !> number of atoms
  REAL(DP),         INTENT(INOUT) :: etot_eng         !> total energy in current step
  INTEGER,          INTENT(IN)    :: order(nat)       !> Engine order of atom
  REAL(DP),         INTENT(IN)    :: at(3,3)          !> lattice parameters in alat units
  INTEGER,          INTENT(IN)    :: ityp(nat)        !> atom types
  INTEGER,          INTENT(IN)    :: if_pos(3,nat)    !> coordinates fixed by engine
  CHARACTER(LEN=3), INTENT(IN)    :: atm(*)           !> name of atom corresponding to ityp
  REAL(DP),         INTENT(IN)    :: force(3,nat)     !> force calculated by the engine
  REAL(DP),         INTENT(INOUT) :: tau(3,nat)       !> atomic positions (needed for output only)
  REAL(DP),         INTENT(OUT)   :: displ_vec(3,nat) !> displacement vector communicated to move mode
  INTEGER,          INTENT(OUT)   :: disp             !> Stage for move_mode
  LOGICAL,          INTENT(OUT)   :: lconv            !> flag for controlling convergence

  ! -- LOCAL VARIABLES
  REAL(DP), EXTERNAL              :: dnrm2, ddot      ! lapack functions
  INTEGER                         :: na, icoor        ! integers for loops
  INTEGER                         :: iidum            ! integers for loops
  REAL(DP)                        :: fpara(3,nat)     ! force parallel to push/eigenvec
  REAL(DP)                        :: fperp(3,nat)     ! force parallel to push/eigenvec
  REAL(DP)                        :: fpara_tot        ! total force in parallel direction
  INTEGER                         :: ios ,i           ! file IOSTAT
  LOGICAL                         :: lforc_conv       ! flag true when forces are converged 
  LOGICAL                         :: lsaddle_conv     ! flag true when saddle is reached
  LOGICAL                         :: ArtnStep         ! Is it an artn step?
  LOGICAL                         :: lerror           ! flag for an error from the engine
  character(len=256)              :: outfile          ! file where are written the steps
  REAL(DP)                        :: z

  !
  ! The ARTn algorithm proceeds as follows:
  ! (1) push atoms in the direction specified by user & relax in the perpendicular direction;
  ! (2) use the lanczos algorithm calculate the lowest eigenvalue/eigenvec
  ! (3) a negative eigenvalue, update push direction otherwise push again
  ! (4) follow the lanczos direction twoard the saddle point
  ! (5) push twoards adjacent minimum & initial minimum
  !
  ! ... Flags that controls convergence
  lconv        = .false.
  lforc_conv   = .false.
  lsaddle_conv = .false.
  !
  ! ... fpara_tot is used to scale the magnitude of the eigenvector
  fpara_tot = 0.D0




  !
  ! ... Initialize artn
  IF( istep == 0 )THEN !! ---------------------------------------------------------------------------------------------  ISTEP = 0
    !
    ! ...Initialize if it is the first search
    IF( isearch == 0 )CALL setup_artn( nat, iunartin, filin )
    !
    ! ...Fill the *_step Arrays and parameters
    CALL Fill_param_step( nat, at, order, tau, etot_eng, force, lerror )
    IF ( lerror ) THEN
       ! 
       ! ... Something went wrong in filling the arrays!
       disp =void 
       call write_fail_report( iunartout, disp, etot_eng )
       lconv = .true.
       !
    ENDIF
    !
    ! ... Initialize pushvect and eigenvec accoriding to user's choice
    call start_guess( zseed, nat, order, force, push, eigenvec )
    !
    IF( lrestart ) THEN
      !
      ! ...Signal that it is a restart
      OPEN ( UNIT = iunartout, FILE = filout, FORM = 'formatted', STATUS = 'old', POSITION = 'append', IOSTAT = ios )
      WRITE (iunartout, '(5x,a/)') "|> Restarted previous ARTn calculation"
      CLOSE ( UNIT = iunartout, STATUS = 'KEEP')
      !
      ! ...Read the FLAGS, FORCES, POSITIONS, ENERGY, ...
      CALL read_restart( restartfname, nat, order, ityp, lerror )
      IF( lerror )THEN
        error_message = 'RESTART FILE DOES NOT EXIST'
        call write_fail_report( iunartout, disp, etot_step )
        lconv = .true.
      ENDIF
      !
      ! ...Overwirte the engine Arrays
      tau(:,:) = tau_step(:,order(:))
      !
    ELSE
      ! 
      ! ...Create The input
      !!    To be able to do multiple research in the same run we keep 
      !!    in memory "isearch" how many time we pass here and open an output 
      !!    file only once
      IF( isearch == 0 ) CALL write_initial_report( iunartout, filout )
      isearch = isearch + 1

      ! ...Initial parameter
      etot_init = etot_step
      !
    ENDIF
    !
    ! ...Split the force field in para/perp field following the push field
    !CALL perpforce( force_step, if_pos, push, fperp, fpara, nat)
    CALL splitfield( 3*nat, force_step, if_pos, push, fperp, fpara )

    ! ...Start to write the output
    CALL write_header_report( iunartout )
    CALL write_report( etot_step, force_step, fperp, fpara, lowest_eigval, if_pos, istep, nat,  iunartout )

    ! ...Write the structure
    CALL write_struct( at, nat, tau, order, elements, ityp, push, etot_eng, 1.0_DP, iunstruct, struc_format_out, initpfname )
    artn_resume = '* Start: '//trim(initpfname)//'.'//trim(struc_format_out)
    !



  ELSE !! ---------------------------------------------------------------------------------------------------  ISTEP > 0
    !! receive variables from the engine, split force into perp and para, and check if it is converged
    !
    ! ...Fill the *_step Arrays
    CALL Fill_param_step( nat, at, order, tau, etot_eng, force, lerror )
    IF( lerror ) THEN
       !! somehing went wrong
       call write_fail_report( iunartout, void, etot_step )
       !! finish current search
       displ_vec = 0.0_DP
       lconv = .true.
    ENDIF
    !
    ! ...Split the force field in para/perp field following the push field
    CALL splitfield( 3*nat, force_step, if_pos, push, fperp, fpara )

    ! ...Write Output
    CALL write_report( etot_step, force_step, fperp, fpara, lowest_eigval, if_pos, istep, nat,  iunartout)

    ! ...Check the convergence forces  
    CALL check_force_convergence( nat, force_step, if_pos, fperp, fpara, lforc_conv, lsaddle_conv )
    !
  ENDIF





  !
  disp = VOID
  !
  ! initial displacement , then switch off linit, and pass to lperp
  !
  IF ( linit ) THEN
     !
     !=============================
     ! Send a push with initial push vector and decide what to do next: perp_relax, or lanczos
     !=============================
     ! linit flag is touched by:
     !   - initialize_artn(),
     !   - check_force_convergence()
     !   - here
     !.............................

     ! ...User cancel the INIT push
     IF ( istep ==0 .AND. ninit==0 ) THEN
        !
        ! Pass to lanczos 
        llanczos = .true.
        linit    = .false.
        lperp    = .false.
        !
     ELSE
        !
        ! Do init push, and switch to perp relax for next step
        iinit = iinit + 1
        disp  = INIT
        prev_push = disp !! save previous push
        !
        ! displacement equal to the push
        displ_vec = push
        !
        !call info_field( iunartout, nat, displ_vec, "init::displ_vec" )
        ! ...set up the flags for next step (we do an initial push, then we need to relax perpendiculary)
        linit = .false.
        lperp = .true.
        !
     ENDIF
     ilanc = 0
     !
  ELSE IF ( lperp ) THEN
     !
     !===============================================
     ! Relax forces perpendicular to eigenvector/push
     !===============================================
     ! lperp is touched by:
     !   - initialize_artn(),
     !   - check_force_convergence()
     !   - here
     !.............................
     !
     ! ...Some pre-processing on the fperp
     !fperp_tot = ddot(3*nat, fperp(:,:), 1, fperp(:,:), 1)
     !current_step_size = MIN(eigen_step_size,ABS(fpara_tot)/MAX( ABS(lowest_eigval), 0.01_DP ))
     !
     !call compute_curve( iperp, 3*nat, tau_step, fperp )
     !
     !
     disp = PERP
     displ_vec(:,:) = fperp(:,:)
     !
     iperp = iperp + 1
     !
     !> Here we do a last verification on displ_vec to detect
     !! the box explosion
     !! -> Stop the search if one of displacement has 5 number
     z = 0.0_DP
     do i = 1,nat
        z = max( z, norm2(displ_vec(:,i)) )
     enddo
     IF( nat /= natoms .OR. z > 1.0e4 )THEN
       error_message = "Box explosion"
       lconv = .true.  !! Stop the research
     ENDIF
     ! 
     !call info_field( iunartout, nat, displ_vec, "perp::displ_vec" )
     !
  ELSE IF ( leigen  )THEN
     !================================================
     ! Push in the direction of the lowest eigenvector
     !================================================
     !
     ! leigen is .true. after we obtain a good eigenvector
     ! if we have a good lanczos eigenvector use it as push vector
     !
     !
     disp    = EIGN
     ismooth = ismooth + 1
     ieigen  = ieigen  + 1

     ! ...reset the iterator of previous step
     ilanc   = 0

     !
     IF( nsmooth > 0 .AND. ismooth <= nsmooth ) THEN
       CALL smooth_interpol( ismooth, nsmooth, nat, force_step, push, eigenvec )
       prev_push = SMTH !! save previous push
     ELSE 
       push(:,:) = eigenvec(:,:)
       prev_push = disp !! save previous push
     ENDIF
     !
     ! rescale the eigenvector according to the current force in the parallel direction
     ! see Cances_JCP130: some improvements of the ART technique doi:10.1063/1.3088532
     ! 0.13 is taken from ARTn, 0.5 eV/Angs^2 corresponds roughly to 0.01 Ry/Bohr^2
     !
     ! ...Recompute the norm of fpara because eigenvec change a bit
     fpara_tot = ddot(3*nat, force_step(:,:), 1, eigenvec(:,:), 1)
     current_step_size = -SIGN(-1.0_DP,fpara_tot)*MIN(eigen_step_size,ABS(fpara_tot)/MAX( ABS(lowest_eigval), 0.01_DP ))
     !
     ! Put some test on current_step_size
     !
     displ_vec(:,:) = eigenvec(:,:)*current_step_size
     ! 
     IF ( ieigen >= neigen  ) THEN
        ! do a perpendicular relax
        lperp = .true.
     ENDIF
     !
     CALL write_struct( at, nat, tau, order, elements, ityp, force_step, &
          etot_eng, 1.0_DP, iunstruct, struc_format_out, eigenfname )
     !
  END IF




  !
  ! The saddle point is reached -> confirmed by check_force_convergence()
  !
  !> SHOULD BE A ROUTINE but not :: it's because we call write_struct() that needs
  !!  arguments exist only in artn()
  IF( lsaddle_conv )THEN
     !
     !> store the saddle point energy
     etot_saddle = etot_step
     tau_saddle = tau_step
     eigen_saddle = eigenvec
     !
     !lsaddle = .true.
     lpush_over = .true.
     ifound = ifound + 1
     !
     !CALL write_struct( at, nat, tau, order, elements, ityp, force_step, 1.0_DP, iunstruct, struc_format_out, sadfname )
     call make_filename( outfile, prefix_sad, nsaddle )
     CALL write_struct( at, nat, tau, order, elements, ityp, force_step, &
          etot_eng, 1.0_DP, iunstruct, struc_format_out, outfile )
     artn_resume = trim(artn_resume)//" | "//trim(outfile)//'.'//trim(struc_format_out)
     !
     CALL write_end_report( iunartout, lpush_over, lpush_final, etot_step - etot_init )
     !
     !> If the saddle point is lower in energy
     !!  than the initial point: Mode refine
     IF ( etot_step < etot_init ) THEN
        ! ...HERE Warning to says we should be in refine saddle mode
        OPEN ( UNIT = iunartout, FILE = filout, FORM = 'formatted', STATUS = 'old', POSITION = 'append', IOSTAT = ios )
        WRITE( iunartout, '(5x,a)' ) "|> WARNING::E_Saddle < E_init => Should be a saddle refine mode"
        CLOSE(iunartout)
     ENDIF
     !
  ENDIF



  !
  ! ...If saddle point is reached
  ! This block do only Push to adjacent minima after the saddle point
  !
  !IF ( lsaddle ) THEN
  IF ( lpush_over ) THEN
     !
     ! do we do a final push ?
     !
     IF ( lpush_final ) THEN
        ! set convergence and other flags to false
        lconv    = .false.
        lperp    = .false.
        leigen   = .false.
        llanczos = .false.
        !
        ! normalize eigenvector
        IF( lbackward ) THEN
          eigenvec(:,:) = eigen_saddle(:,:)
          lbackward     = .false.
          etot_step     = etot_saddle
        ELSE
          !! Normalize it to be sure
          eigenvec(:,:) = eigenvec(:,:)/dnrm2(3*nat,eigenvec,1)
        ENDIF
        !
        !
        ! ...PUSH_OVER works => If diff Energy is negative
        IF( etot_step - etot_saddle < frelax_ene_thr ) THEN
           ! we started going downhill ...
           if( .NOT.lrelax )irelax = 0
           lrelax = .true.
           lpush_over = .false.
           !
        ELSE  !< It is a PUSH_OVER the saddle point
           disp =EIGN 
           CALL PUSH_OVER_PROCEDURE( iover, nat, tau, eigenvec, fpush_factor, order, displ_vec, lconv )
           !
        END IF
        !
     ELSE  ! --- NO FINAL_PUSH
        !
        !> At this point the saddle point is already wrote by write_struct before
        !! Here we finish the ARTn search.
        !! Preparation of the possible new ARTn search following this step.
        !! - Cleaning the flag/parameter
        !! - write in output saying no more research
        !! - return a configuration in which a new ARTn search can start
        !
        !call write_end_report( iunartout, lsaddle, lpush_final, 0.0_DP )
        OPEN ( UNIT = iunartout, FILE = filout, FORM = 'formatted', STATUS = 'old', POSITION = 'append', IOSTAT = ios )
        WRITE(iunartout,'(5x,a/)') "|> NO FINAL_PUSH :: Return to the start configuration "
        CLOSE(iunartout)

        !call clean_artn()  !! No final_push

        ! ...Return to the initial comfiguration
        tau(:,:) = tau_init(:,order(:))

        ! ...Tell to the engine it is finished
        call flag_false()

        lconv = .true.

        ! ...Set the force to zero
        displ_vec = 0.0_DP

        ! ...Here we dont load the next minimum because it does not exist

     ENDIF
  ENDIF




  !
  ! perform a FIRE relaxation (only for reaching adjacent minima)
  !
  RELAX: IF ( lrelax ) THEN
     !
     ! reset force
     !
     disp      = RELX
     displ_vec = force_step
     irelax    = irelax + 1

     prev_push = disp !! Save the previous displacement 

     !
     ! The convergence is reached:
     !  - Switch the push_over or
     !  - Finish the ARTn search
     !
     IF ( lforc_conv ) THEN
        !
        IF ( fpush_factor == 1.0 ) THEN
           !
           ! ...It found the adjacent minimum!
           !   We save it and return to the saddle point
           CALL make_filename( outfile, prefix_min, nmin )
           CALL write_struct( at, nat, tau, order, elements, ityp, force_step, &
                etot_eng, 1.0_DP, iunstruct, struc_format_out, outfile )
           artn_resume = trim(artn_resume)//" | "//trim(outfile)//'.'//trim(struc_format_out)
           !
           ! ...Save the minimum if it is new
           call save_min( nat, tau_step )
           disp = RELX
           !
           ! ...restart from saddle point
           tau(:,:)      = tau_saddle(:,order(:))
           eigenvec(:,:) = eigen_saddle(:,:)
           lbackward     = .true.
           !
           ! ...Return to Push_Over Step in opposit direction
           !lsaddle = .true.
           lpush_over = .true.
           lrelax     = .false.
           !
           etot_final = etot_step
           de_back = etot_saddle - etot_final
           !
           call write_inter_report( iunartout, int(fpush_factor), [de_back] )
           !
           ! ...reverse direction for the push_over
           fpush_factor = -1.0
           irelax = 0
           !
        !ELSEIF( .NOT.lend )THEN  !< If already pass before no need to rewrite again
        ELSE  !< If already pass before no need to rewrite again
           !
           ! ...It found the starting minimum! (should be the initial configuration)
           CALL make_filename( outfile, prefix_min, nmin )
           CALL write_struct( at, nat, tau, order, elements, ityp, &
                force_step, etot_eng, 1.0_DP, iunstruct, struc_format_out, outfile )
           ! ...Save the structure name file to print it
           artn_resume = trim(artn_resume)//" | "//trim(outfile)//'.'//trim(struc_format_out)
           !
           ! ...Communicate to the engine it is finished
           CALL flag_false()
           !
           lconv = .true.
           lend = lconv  !! Maybe don't need anymore
           ! 
           ! ...Save the Energy difference
           de_fwd = etot_saddle - etot_step
           !
           call write_inter_report( iunartout, int(fpush_factor), &
                [de_back, de_fwd, etot_init, etot_final, etot_step] )
           ! 
        END IF
        !
     END IF
     !
  END IF RELAX




  !> WHAT FOR THIS BLOCK??? ANTOINE??
  !!  This should be in check_force()
  IF( etot_step - etot_init > etot_diff_limit ) then
     error_message = 'ENERGY EXCEEDS THE LIMIT'
     call write_fail_report( iunartout, disp, etot_step )
     lconv = .true.
  ENDIF





  !
  ! check if we should perform the lanczos algorithm
  !
  ! Lanczos at the end. Reason: check convergence at saddle before going into
  ! un-needed lanczos near saddle.
  !
  LANCZOS_: IF ( llanczos ) THEN
     !
     !==========================================
     ! Perform Lanczos algo, one step at a time
     !==========================================
     !
     disp =LANC
     IF (ilanc == 0 ) THEN
        !
        ! first iternation of current lanczos call
        !
        v_in(:,:) = eigenvec(:,:)
        !
        IF( lanczos_always_random ) THEN
           ! generate random initial vector
           CALL random_number(z)
           z = z *1e8
           iidum = INT(z)
           DO na = 1, nat
              v_in(:,na) = (/0.5_DP - ran3(iidum),0.5_DP - ran3(iidum),0.5_DP - ran3(iidum)/)
           ENDDO
           ! normalize
           v_in = v_in / norm2(v_in)
        ENDIF
        !
        ! reset the eigenvalue flag
        !
        leigen = .false.
        !
        ! allocate memory for previous lanczos vec
        !
        if( .not. allocated( old_lanczos_vec ) ) allocate( old_lanczos_vec, source = v_in )
        a1 = 0.0
     ENDIF
     !
     ! apply constraints from the engine. Works only with engines which fill if_pos!! (not lammps)
     !
     IF ( ANY(if_pos(:,:) == 0) ) THEN
        DO na=1,nat
           DO icoor=1,3
              IF (if_pos(icoor,na) == 1 ) if_pos_ct = if_pos_ct + 1
           ENDDO
        END DO
        IF ( if_pos_ct < nlanc .and. if_pos_ct /= 0 ) nlanc = if_pos_ct
        v_in(:,:) = v_in(:,:)*if_pos(:,:)
        force_step(:,:) = force_step(:,:)*if_pos(:,:)
     ENDIF
     !
     !
     CALL lanczos( nat, v_in, push, force_step, &
          ilanc, nlanc, lowest_eigval, eigenvec, displ_vec)
     !
     ilanc = ilanc + 1
     !
     ! Lanczos has converged:
     ! nlanc = number of steps it took to converge,
     !
     IF ( ilanc > nlanc ) THEN
        !
        ! check lowest eigenvalue, decide what to do in next step
        !
        ilanc_save =ilanc
        IF ( lowest_eigval < eigval_thr ) THEN
           ! structure is out of the basin (above inflection),
           ! in next step make a push with the eigenvector
           !! Next Mstep outside the basin
           lbasin = .false.
           ! ...push in eigenvector direction
           leigen = .true.
           ! ...Save the eigenvector
           ! ...No yet perp relax
           lperp  = .false.
           old_lowest_eigval = lowest_eigval
           !
        ELSE
           !
           ! ...If we lose the eigval
           IF ( .NOT. lbasin .AND. lowest_eigval > 0.0) THEN
              ! 
              IF (inewchance < nnewchance) THEN
                 ! Hope by continue pushing along init we find something 
                 inewchance = inewchance +1
                 ismooth      = 0
              ELSE
                 error_message = 'EIGENVALUE LOST'
                 call write_fail_report( iunartout, disp, lowest_eigval )
                 lconv = .true.
              ENDIF
              !
           ENDIF
           !
           ! structure is still in basin (under unflection),
           ! in next step it move following push vetor (can be a previous eigenvec)
           !! Next Mstep inside the Basin
           !lowest_eigval = 0.D0
           leigen = .false.
           linit  = .true.
           lbasin = .true.
           noperp = 0      !> count the init-perp fail
           nperp_step = 1  !> count the out-basin perp relax step
           !
        ENDIF
        !
        ! ...Compare the eigenvec with the previous one
        !
        a1 = ddot( 3*nat, eigenvec, 1, old_lanczos_vec, 1 )
        a1 = abs( a1 )
        ! set current eignevec for comparison in next step
        old_lanczos_vec = eigenvec
        ! deallocate( old_lanczos_vec )
        !
        ! finish lanczos for now
        !
        llanczos = .false.
        !
        ! reset lanczos size for next call
        !
        nlanc = lanczos_max_size
        !
     ENDIF
     ! 
  ENDIF LANCZOS_



  !
  !! --- Finalization Block
  !
  IF( lconv )THEN
    !
    OPEN ( UNIT = iunartout, FILE = filout, FORM = 'formatted', STATUS = 'old', POSITION = 'append', IOSTAT = ios )
    WRITE (iunartout,'(5x, "|> BLOCK FINALIZE..")')
    WRITE (iunartout,'(5X, "|> number of steps:",x, i0)') istep
    CLOSE (iunartout)
    !> SCHEMA FINILIZATION
    lend = lconv
    !
    IF( lerror ) THEN
       ! STOP the search
       error_message = 'STOPPING DUE TO ERROR'
       call write_fail_report( iunartout, void, etot_step )
       STOP
    ENDIF
    !
    ! ...Here we should load the next minimum if the user ask
    IF( lmove_nextmin )THEN
      CALL move_nextmin( nat, tau )
    ELSE
      tau(:,:) = tau_init(:,order(:))
    ENDIF
    !
    ! ...Force = 0.0
    displ_vec = 0.0_DP
    !disp = VOID
    disp = RELX
    !
    ! ...The search IS FINISHED
    RETURN
    !
  ENDIF
  !
  ! ...Increment the ARTn-step
  istep = istep + 1
  !
END SUBROUTINE artn







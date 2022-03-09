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
       lrelax, linit, lperp, leigen, llanczos, lrestart, lbasin, lpush_over, lpush_final, lbackward, lmove_nextmin, lread_param,  &
       irelax, istep, iperp, ieigen, iinit, ilanc, ismooth, iover, ifails, isearch, ifound, nlanc, nperp, noperp, nperp_step,  &
       if_pos_ct, lowest_eigval, etot_init, etot_step, etot_saddle, etot_final, de_saddle, de_back, de_fwd, &
       ninit, neigen, lanc_mat_size, nsmooth, push_mode, dist_thr, init_forc_thr, forc_thr, &
       fpara_thr, eigval_thr, frelax_ene_thr, push_step_size, current_step_size, dlanc, eigen_step_size, fpush_factor, &
       push_ids, add_const, push, eigenvec, tau_step, force_step, tau_init, tau_saddle, eigen_saddle, v_in, &
       VOID, INIT, PERP, EIGN, LANC, RELX, OVER, zseed, &
       engine_units, struc_format_out, elements, &
       initialize_artn, read_restart, write_restart, &
       push_over, ran3, a1, old_lanczos_vec, lend, lat, fill_param_step, &
       filin, filout, sadfname, initpfname, eigenfname, restartfname, warning, flag_false,  &
       prefix_min, nmin, prefix_sad, nsaddle, artn_resume
  !
  IMPLICIT NONE
  ! -- ARGUMENTS
  INTEGER,  INTENT(IN), value ::    nat       !> number of atoms

  REAL(DP), INTENT(INOUT) ::    etot_eng         !> total energy in current step
  INTEGER,  INTENT(IN) ::    order(nat)       !> Engine order of atom
  REAL(DP), INTENT(IN) ::    at(3,3)          !> lattice parameters in alat units
  INTEGER,  INTENT(IN) ::    ityp(nat)        !> atom types
  INTEGER,  INTENT(IN) ::    if_pos(3,nat)    !> coordinates fixed by engine
  CHARACTER(LEN=3),   INTENT(IN) :: atm(*)    !> name of atom corresponding to ityp

  REAL(DP), INTENT(IN) :: force(3,nat)     !> force calculated by the engine
  REAL(DP), INTENT(INOUT) :: tau(3,nat)       !> atomic positions (needed for output only)

  REAL(DP), INTENT(OUT)  :: displ_vec(3,nat)  !> displacement vector communicated to move mode
  INTEGER,          INTENT(OUT) :: disp       !> Stage for move_mode
  LOGICAL,          INTENT(OUT) :: lconv      !> flag for controlling convergence

  ! -- LOCAL VARIABLES
  REAL(DP), EXTERNAL :: dnrm2, ddot           ! lapack functions
  INTEGER :: na, icoor, idum                  ! integers for loops
  !
  !
  REAL(DP)  :: fpara(3,nat)                   ! force parallel to push/eigenvec
  REAL(DP)  :: fperp(3,nat)                   ! force parallel to push/eigenvec
  REAL(DP)  :: fpara_tot                      ! total force in parallel direction
  REAL(DP)  :: smoothing_factor               ! mixing factor for smooth transition between eigenvec and push
  REAL(DP)  :: etot!, lat(3,3)
  INTEGER   :: ios ,i                         ! file IOSTAT
  !CHARACTER( LEN=255) :: filin !, filout, sadfname, initpfname, eigenfname, restartfname
  LOGICAL :: lforc_conv, lsaddle_conv, ArtnStep
  !character(:), allocatable :: outfile
  character(len=256) :: outfile

  integer :: natom
  LOGICAL, PARAMETER :: noARTnStep = .false.
  REAL(DP) :: z

  !
  ! The ARTn algorithm proceeds as follows:
  ! (1) push atoms in the direction specified by user & relax in the perpendicular direction;
  ! (2) use the lanczos algorithm calculate the lowest eigenvalue/eigenvec
  ! (3) a negative eigenvalue, update push direction otherwise push again
  ! (4) follow the lanczos direction twoard the saddle point
  ! (5) push twoards adjacent minimum & initial minimum
  !
  ! flag that controls convergence
  !
  lconv = .false.
  lforc_conv = .false.
  lsaddle_conv = .false.
  !
  ! fpara_tot is used to scale the magnitude of the eigenvector
  !
  fpara_tot = 0.D0
  !
  !filin = 'artn.in'
  !> Move output file name in module : Can be customize




  !
  ! initialize artn
  !


  IF( istep == 0 )THEN

    !write(*,*) "ISTEP = 0: lread_param; ifails", lread_param, ifails, ( .NOT.lread_param.AND.ifails == 0 )
    !> @brief lread_param and ifail are there to not read input multiple time
    !!  we need 2 parameter to be able to write the output header once 

    !> Don't initialize ARTn if it fails before
    !NOFAILS: if( .NOT.lread_param.AND.ifails == 0 )then  
    ONCE: IF( isearch == 0 )THEN

      ! ...Read the input parameters
      CALL initialize_artn( nat, iunartin, filin )


      ! set initial random seed from input (could be moved to initialize_artn)
      ! value zseed = 0 means generate random seed
      idum = zseed
      IF( idum .EQ. 0) THEN
        !
        ! generate random seed
        CALL random_number(z)
        z = z *1e8
        idum = INT(z)
      ENDIF
      !> Save the seed for DEBUG
      open( newunit=zseed, file="random_seed.dat" )
      write( zseed, * )" zseed = ", int(z)
      close( zseed )

    endif ONCE


    ! ...Fill the *_step Arrays
    !CALL Fill_param_step( nat, at, order, tau, etot_eng, force )
    lat = at
    tau_step(:,order(:)) = tau(:,:)
    etot_step = convert_energy( etot_eng )
    force_step(:,order(:)) = convert_force( force(:,:) )



    IF( lrestart ) THEN

      ! ...Signal that it is a restart
      OPEN ( UNIT = iunartout, FILE = filout, FORM = 'formatted', STATUS = 'old', POSITION = 'append', IOSTAT = ios )
      WRITE (iunartout, *) "Restarted previous ARTn calculation"
      CLOSE ( UNIT = iunartout, STATUS = 'KEEP')

      ! ...Read the FLAGS, FORCES, POSITIONS, ENERGY, ...
      CALL read_restart( restartfname, nat, order, ityp )

      ! ...Overwirte the engine Arrays
      tau(:,:) = tau_step(:,order(:))

      isearch = isearch + 1

    ELSE

      if( isearch == 0 )CALL write_initial_report( iunartout, filout )
      isearch = isearch + 1
      ! ...Initial parameter
      etot_init = etot_step

    ENDIF

    !
    ! --------------------------
    !

    !------------------------------------------------------------------------------    
    !> @brief 
    !!   Here we have to initialize the push and eigenvec thanks to different way
    !!   depending the user choice.
    !!   - push_init() works for random 
    ! ...Define the Initial Push
    CALL push_init(nat, tau, order, lat, idum, push_ids, dist_thr, add_const, push_step_size, push , push_mode)


    ! ...Generate first lanczos eigenvector 
    add_const(:,:) = 0.0
    CALL push_init(nat, tau, order, lat, idum, push_ids, dist_thr, add_const, eigen_step_size, eigenvec , 'all ')
    !------------------------------------------------------------------------------    



    ! ...Split the force field in para/perp field following the push field
    !CALL perpforce( force_step, if_pos, push, fperp, fpara, nat)
    CALL splitfield( 3*nat, force_step, if_pos, push, fperp, fpara )



  ELSE !>     ISTEP > 0

    ! ...Fill the *_step Arrays
    !CALL Fill_param_step( nat, at, order, tau, etot_eng, force )
    lat = at
    tau_step(:,order(:)) = tau(:,:)
    etot_step = convert_energy( etot_eng )
    force_step(:,order(:)) = convert_force( force(:,:) )


    CALL perpforce( force_step, if_pos, push, fperp, fpara, nat)
    CALL check_force_convergence( nat, force_step, if_pos, fperp, fpara, lforc_conv, lsaddle_conv )

  ENDIF



  ! ...Initialize the displacement
  disp = VOID


  !
  ! Open the output file for writing
  !
  OPEN ( UNIT = iunartout, FILE = filout, FORM = 'formatted', ACCESS = 'append', STATUS = 'unknown', IOSTAT = ios )


  ! -- Start a ARTn search
  !write( iunartout,*)"|> test print: istep", istep
  if( istep == 0 )then
     ! ...Write Zero step
     CALL write_header_report( iunartout )
     CALL write_report( etot_step, force_step, fperp, fpara, lowest_eigval, VOID, if_pos, istep, nat, iunartout, .true. )

     CALL write_struct( at, nat, tau, order, elements, ityp, push, etot_eng, &
          1.0_DP, iunstruct, struc_format_out, initpfname )
     artn_resume = '* Start: '//trim(initpfname)
  endif


  !
  ! initial displacement , then switch off linit, and pass to lperp
  !
  IF ( linit ) THEN
     !
     !=============================
     ! generate initial push vector
     !=============================
     ! linit is touched by:
     !   - initialize_artn(),
     !   - check_force_convergence()
     !   - here
     !.............................



     ! ...Start lanczos when number of init steps is reached
     IF ( iinit >= ninit ) THEN

        ! ...For actual step
        llanczos = .true.
        ilanc = 0
        ! ...For next Step
        linit = .false.
        lperp = .false.


     ELSE  ! ...Init Push

        iinit = iinit + 1
        disp = INIT


        ! ...modify the force to be equal to the push
        displ_vec = push
 

        CALL write_report( etot_step, force_step, fperp, fpara, lowest_eigval, &
             disp, if_pos, istep, nat, iunartout, noARTnStep )

        !call info_field( iunartout, nat, displ_vec, "init::displ_vec" )


        ! ...set up the flags (we did an initial push, now we need to relax perpendiculary)
        linit = .false.
        lperp = .true.
        iperp = 0

     ENDIF
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
     ! If eigenvalue is good, overwrite push with eigenvec (NS: Why?!)
     !  If we come back in bassin we use this eigenvec and not push_init
     !
     !
     ! Subtract parrallel components to push from force.
     !
     disp = PERP
     !
     !
     ! If the force_perp component are small we continue
     ! to push
     !
     displ_vec(:,:) = fperp(:,:)
     CALL write_report( etot_step, force_step, fperp, fpara, lowest_eigval, disp, &
          if_pos, istep, nat, iunartout, noARTnStep )
     !
     iperp = iperp + 1
     !

     !> Here we do a last verification on displ_vec to prevent 
     !! the box explosion
     !! -> Stop the search if one of displacement has 5 number
     z = 0.0_DP
     do i = 1,nat
        z = max( z, norm2(displ_vec(:,i)) )
     enddo
     IF( z > 1.0e4 )THEN

       write(*,'(5x,"|> BOX EXPLOSION - STEP",x,i0," - STOP ARTn RESEARCH ")') istep
       write(iunartout,'(5x,"|> BOX EXPLOSION - STEP",x,i0," - STOP ARTn RESEARCH ")') istep
       write(iunartout,'(5x,150("-")/)') 
      
       !ifails = ifails + 1
       

       !> SCHEMA FINILIZATION
       ! ...Laod the start configuration
       tau(:,:) = tau_init(:,order(:))

       ! ...Force = 0.0
       displ_vec = 0.0_DP

       ! ...The research IS FINISHED
       !CALL clean_artn()  ! Explosion Box

       !disp = RELX
       lconv = .true.  !! Stop the research 
       !return

       
     ENDIF

     !call info_field( iunartout, nat, displ_vec, "perp::displ_vec" )

     !
  ELSE IF ( leigen  )THEN
     !================================================
     ! Push in the direction of the lowest eigenvector
     !================================================
     !
     ! leigen is .true. after we obtain a good eigenvector
     !
     ! if we have a good lanczos eigenvector use it as push vector
     !
     !
     disp = EIGN


     IF( nsmooth > 0 ) &
       CALL Smooth_interpol( ismooth, nsmooth, nat, force_step, push, eigenvec )

     ! ...Overwrite the initial Push with Eigenvector
     IF( nsmooth == 0.OR.ismooth > nsmooth )then
       !write(iunartout,'(x,"DEBUG::EIGEN::Overwrite push = eigenvec")')
       push(:,:) = eigenvec(:,:)  
     ENDIF


     !
     ! rescale the eigenvector according to the current force in the parallel direction
     ! see Cances_JCP130: some improvements of the ART technique doi:10.1063/1.3088532
     !
     ! 0.13 is taken from ARTn, 0.5 eV/Angs^2 corresponds roughly to 0.01 Ry/Bohr^2
     !
     ! ...Recompute the norm of fpara because eigenvec change a bit
     fpara_tot = ddot(3*nat, force_step(:,:), 1, eigenvec(:,:), 1)
     current_step_size = MIN(eigen_step_size,ABS(fpara_tot)/MAX( ABS(lowest_eigval), 0.01_DP ))
     !current_step_size = MIN(eigen_step_size,ABS(MAXVAL(fpara))/MAX( ABS(lowest_eigval), 0.01_DP ))
     !
     !%! Put some test on current_step_size
     !
     displ_vec(:,:) = eigenvec(:,:)*current_step_size


     ! count the number of steps made with the eigenvector
     ieigen = ieigen + 1

     IF ( ieigen == neigen  ) THEN
        ! do a perpendicular relax
        lperp = .true.
        iperp = 0
     ENDIF

     CALL write_struct( at, nat, tau, order, elements, ityp, force_step, &
          etot_eng, 1.0_DP, iunstruct, struc_format_out, eigenfname )
     CALL write_report( etot_step, force_step, fperp, fpara, lowest_eigval, disp, if_pos, istep, nat,  iunartout, noARTnStep )

  END IF



  !
  ! The saddle point is reached -> confirmed by check_force_convergence()
  !
  !> SHOULD BE A ROUTINE but not :: it's because we call write_struct() that needs 
  !!  arguments exist only in artn() 
  IF( lsaddle_conv )THEN

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
     artn_resume = trim(artn_resume)//" | "//trim(outfile)
     !


     !> If the saddle point is lower in energy
     !!  than the initial point: Mode refine
     IF ( etot_step < etot_init ) THEN
        ! ...HERE Warning to says we should be in refine saddle mode
        write( iunartout, '(5x,a)' ) "!> WARNING::E_Saddle < E_init => Should be a saddle refine mode"
     ENDIF

     !CALL write_end_report( iunartout, lsaddle, lpush_final, etot_step - etot_init )
     CALL write_end_report( iunartout, lpush_over, lpush_final, etot_step - etot_init )

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
        lconv = .false.
        lperp = .false.
        leigen = .false.
        llanczos = .false.
        !
        ! normalize eigenvector
        if( lbackward )then
          eigenvec(:,:) = eigen_saddle(:,:)
          lbackward = .false.
        else
          !! Normalize it to be sure
          eigenvec(:,:) = eigenvec(:,:)/dnrm2(3*nat,eigenvec,1)
        endif
        !
        !
        ! ...PUSH_OVER works => If diff Energy is negative
        IF( etot_step - etot_saddle < frelax_ene_thr ) THEN
        !IF ( etot_step - etot_saddle < 0.0_DP ) THEN

           ! we started going downhill ...
           if( .NOT.lrelax )irelax = 0
           lrelax = .true.
           !lsaddle = .false.
           lpush_over = .false.



        ELSE  !< It is a PUSH_OVER the saddle point
           disp = EIGN

           CALL PUSH_OVER_PROCEDURE( iover, nat, tau, eigenvec, fpush_factor, order, displ_vec, lconv )

           !
           CALL write_report( etot_step, force_step, fperp, fpara, lowest_eigval, &
                OVER, if_pos, istep, nat,  iunartout, noARTnStep )

        END IF


     ELSE  ! --- NO FINAL_PUSH

        !> At this point the saddle point is already wrote by write_struct before
        !! Here we finish the ARTn search.
        !! Preparation of the possible new ARTn search following this step.
        !! - Cleaning the flag/parameter
        !! - write in output saying no more research
        !! - return a configuration in which a new ARTn search can start

        !call write_end_report( iunartout, lsaddle, lpush_final, 0.0_DP )
        write(iunartout,'(5x,a/)') "|> NO FINAL_PUSH :: Return to the start configuration "
        write(*,'(5x,a)') "|> NO FINAL_PUSH :: Return to the start configuration "

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
     disp = RELX
     displ_vec = force_step

     ! ...Print the relaxation state
     ArtnStep = noArtnStep
     if( mod(irelax,5) == 0 ) ArtnStep = .true.  !> The 5 can be custom parameter : nrprint
     CALL write_report( etot_step, force_step, fperp, fpara, lowest_eigval, disp, &
          if_pos, istep, nat, iunartout, ARTnStep )

     irelax = irelax + 1


     !
     ! The convergence is reached: 
     !  - Switch the push_over or
     !  - Finish the ARTn search
     !
     IF ( lforc_conv ) THEN

        IF ( fpush_factor == 1.0 ) THEN


           ! ...It found the adjacent minimum!
           !   We save it and return to the saddle point
           CALL make_filename( outfile, prefix_min, nmin )
           CALL write_struct( at, nat, tau, order, elements, ityp, force_step, &
                etot_eng, 1.0_DP, iunstruct, struc_format_out, outfile )
           artn_resume = trim(artn_resume)//" | "//trim(outfile)


           ! ...Save the minimum if it is new
           call save_min( nat, tau_step )

           disp = RELX

           ! ...restart from saddle point
           tau(:,:) = tau_saddle(:,order(:))
           eigenvec(:,:) = eigen_saddle(:,:)
           lbackward = .true.

           ! ...Return to Push_Over Step in opposit direction
           !lsaddle = .true.
           lpush_over = .true.
           lrelax = .false.

           etot_final = etot_step
           de_back = etot_saddle - etot_final

           CALL write_report( etot_step, force_step, fperp, fpara, lowest_eigval, &
                disp, if_pos, istep, nat, iunartout, .true. )

           call write_inter_report( iunartout, int(fpush_factor), [de_back] )

           ! ...reverse direction for the push_over
           fpush_factor = -1.0
           irelax = 0


        !ELSEIF( .NOT.lend )THEN  !< If already pass before no need to rewrite again
        ELSE  !< If already pass before no need to rewrite again

           ! ...It found the starting minimum! (should be the initial configuration)
           CALL make_filename( outfile, prefix_min, nmin )
           CALL write_struct( at, nat, tau, order, elements, ityp, &
                force_step, etot_eng, 1.0_DP, iunstruct, struc_format_out, outfile )
           ! ...Save the structure name file to print it
           artn_resume = trim(artn_resume)//" | "//trim(outfile)

           ! ...Tell to the machin it is finished
           call flag_false()

           lconv = .true.
           lend = lconv  !! Maybe don't need anymore

           ! ...Save the Energy difference
           de_fwd = etot_saddle - etot_step


           CALL write_report( etot_step, force_step, fperp, fpara, lowest_eigval, &
                disp, if_pos, istep, nat, iunartout, .true. )

           call write_inter_report( iunartout, int(fpush_factor), &
                [de_back, de_fwd, etot_init, etot_final, etot_step] )

        END IF
        !

     END IF
     !
  END IF RELAX




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
     disp = LANC
     CALL write_report( etot_step, force_step, fperp, fpara, lowest_eigval, disp, &
          if_pos, istep, nat,  iunartout, noARTnStep )

     IF (ilanc == 0 ) THEN
        !
        ! first iternation of current lanczos call
        !
        v_in(:,:) = eigenvec(:,:)
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
     ! apply constraints from the engine (QE)
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
     iperp = 0

     !
     ! Lanczos has converged:
     ! nlanc = number of steps it took to converge,
     ! and ilanc = ilanc + 1
     !
     IF ( ilanc > nlanc ) THEN
        !
        ! check lowest eigenvalue, decide what to do in next step
        !
        IF ( lowest_eigval < eigval_thr ) THEN
           ! structure is out of the basin (above inflection),
           ! in next step make a push with the eigenvector
           !! Next Mstep outside the basin
           lbasin = .false.
           ! ...push in eigenvector direction
           leigen = .true.
           ieigen = 0
           ! ...Save the eigenvector
           !push(:,:) = eigenvec(:,:)  !! Move it in leigen block
           ! ...No yet perp relax
           lperp = .false.
           iperp = 0
           !
        ELSE
           ! structure is still in basin (under unflection),
           ! in next step it move following push vetor (can be a previous eigenvec)
           !! Next Mstep inside the Basin
           !lowest_eigval = 0.D0
           leigen = .false.
           linit  = .true.
           lbasin = .true.
           iperp =  0
           noperp = 0      !> count the init-perp fail
           nperp_step = 1  !> count the out-basin perp relax step
           ismooth = 1     !> Initialise the smoothy step
           iinit = iinit - 1
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
        nlanc = lanc_mat_size

     ENDIF

  ENDIF LANCZOS_



  !
  !! --- Finalization Block
  !
  IF( lconv )THEN

    write(iunartout,'(5x,"|> BLOCK FINALIZE..")')
    !> SCHEMA FINILIZATION
    lend = lconv

    ! ...Here we should load the next minimum if the user ask
    IF( lmove_nextmin )THEN
      CALL move_nextmin( nat, tau )
    ELSE
      tau(:,:) = tau_init(:,order(:))
    ENDIF

    ! ...Force = 0.0
    displ_vec = 0.0_DP

    ! ...The research IS FINISHED
    CLOSE( UNIT = iunartout, STATUS = 'KEEP' )
    return

  ENDIF




  ! ...Increment the ARTn-step
  istep = istep + 1

  ! CALL write_restart(restartfname,nat)
  ! ...Close the output file
  CLOSE (UNIT = iunartout, STATUS = 'KEEP')





END SUBROUTINE artn




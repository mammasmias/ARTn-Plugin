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
!!   modifies the input force to perform the ARTn algorithm
!------------------------------------------------------------------------------
SUBROUTINE artn( force, etot_eng, nat, ityp, atm, tau, order, at, if_pos, disp, displ_vec, lconv )
  !----------------------------------------------------------------------------
  !
  ! artn_params for variables and counters that need to be stored in each step
  ! DEFINED IN: artn_params_mod.f90
  !
  USE units
  USE artn_params, ONLY: iunartin, iunartout, iunstruct, &
       lrelax, linit, lperp, leigen, llanczos, lrestart, lbasin, lsaddle, lpush_final, lbackward, &
       irelax, istep, iperp, ieigen, iinit, ilanc, ismooth, iover, nlanc, nperp, noperp, nperp_step, if_pos_ct, &
       lowest_eigval, etot_init, etot_step, etot_saddle, etot_final, de_saddle, de_back, de_fwd, &
       ninit, neigen, lanc_mat_size, nsmooth, push_mode, dist_thr, init_forc_thr, forc_thr, &
       fpara_thr, eigval_thr, frelax_ene_thr, push_step_size, current_step_size, dlanc, eigen_step_size, fpush_factor, &
       push_ids, add_const, push, eigenvec, tau_step, force_step, tau_saddle, eigen_saddle, v_in, &
       VOID, INIT, PERP, EIGN, LANC, RELX, OVER, zseed, &
       engine_units, struc_format_out, elements, &
       initialize_artn, read_restart, write_restart, &
       push_over, ran3, a1, old_lanczos_vec, lend, lat, fill_param_step, &
       filout, sadfname, initpfname, eigenfname, restartfname, warning
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

  REAL(DP), INTENT(INOUT) :: force(3,nat)     !> force calculated by the engine
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
  CHARACTER( LEN=255) :: filin !, filout, sadfname, initpfname, eigenfname, restartfname
  LOGICAL :: lforc_conv, lsaddle_conv, ArtnStep

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
  filin = 'artn.in'
  !> Moved in module: Can be customize
  !filout = 'artn.out'
  !sadfname = 'saddle'
  !initpfname = 'initp'
  !eigenfname = 'latest_eigenvec'
  !restartfname = 'artn.restart'





  !
  ! initialize artn
  !


  IF( istep == 0 )THEN

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

    ! ...Fill the *_step Arrays
    !CALL Fill_param_step( nat, order, tau, etot_eng, force )
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
      !force(:,:) = force_step(:,order(:))

    ELSE 

      CALL write_initial_report( iunartout, filout )
      ! ...Initial parameter
      etot_init = etot_step

    ENDIF

    !
    ! --------------------------
    !

    ! ...Define the Initial Push
    CALL push_init(nat, tau, order, lat, idum, push_ids, dist_thr, add_const, push_step_size, push , push_mode)


    ! ...Generate first lanczos eigenvector (NS: Maybe put in push_init)
    DO i = 1, nat
       eigenvec(:,i) = (/0.5_DP - ran3(idum),0.5_DP - ran3(idum),0.5_DP - ran3(idum)/) * if_pos(:,i)
    END DO
    eigenvec(:,:) = eigenvec(:,:)/dnrm2(3*nat,eigenvec,1) * dnrm2(3*nat,push,1)

    !> Test
    !eigenvec = push


    CALL perpforce( force_step, if_pos, push, fperp, fpara, nat)


    
  ELSE !>     ISTEP > 0

    ! ...Fill the *_step Arrays
    !CALL Fill_param_step( nat, order, tau, etot_eng, force )
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

  if( istep == 0 )then
     ! ...Write Zero step  
     CALL write_report( etot_step, force_step, fperp, fpara, lowest_eigval, VOID, if_pos, istep, nat, iunartout, .true. )
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
     

     ! ...The push counter (controls if we should call lanczos or keep pushing)
     IF( iinit == 1 )THEN
       CALL write_struct( at, nat, tau, order, elements, ityp, push, 1.0_DP, iunstruct, struc_format_out, initpfname )
     ENDIF
     
     
     ! ...Start lanczos when number of init steps is reached
     IF ( iinit >= ninit ) THEN

        llanczos = .true.
        ilanc = 0
        !write(iunartout,*)"ARTn(linit)::initialize ilanc"
        linit = .false.
        lperp = .false.


     ELSE  ! ...Init Push

        iinit = iinit + 1
        disp = INIT


        ! ...modify the force to be equal to the push
        displ_vec = push

        CALL write_report( etot_step, force_step, fperp, fpara, lowest_eigval, disp, if_pos, istep, nat, iunartout, noARTnStep )


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
     !
     IF( leigen ) THEN
        push(:,:) = eigenvec(:,:)
     ENDIF
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
     CALL write_report( etot_step, force_step, fperp, fpara, lowest_eigval, disp, if_pos, istep, nat, iunartout, noARTnStep )
     !
     iperp = iperp + 1
     !
     !
  ELSE IF ( leigen  ) THEN
     !================================================
     ! Push in the direction of the lowest eigenvector
     !================================================
     !
     ! leigen is always .true. after we obtain a good eigenvector
     ! except during lanczos iterations
     !
     ! if we have a good lanczos eigenvector use it as push vector
     !
     !
     disp = EIGN

     !CALL Smooth_interpol( ismooth, nat, force_step, push, eigenvec, fpara_tot )
     !>>>>
     smoothing_factor = 1.0_DP*ismooth/nsmooth
     !!
     fpara_tot = ddot(3*nat, force_step(:,:), 1, eigenvec(:,:), 1)

     eigenvec(:,:) = (1.0_DP - smoothing_factor)*push(:,:) &
          -SIGN(1.0_DP,fpara_tot)*smoothing_factor*eigenvec(:,:)
     !
     IF ( ismooth < nsmooth) ismooth = ismooth + 1
     !<<<<<
     !
     ! rescale the eigenvector according to the current force in the parallel direction
     ! see Cances_JCP130: some improvements of the ART technique doi:10.1063/1.3088532
     !
     ! 0.13 is taken from ARTn, 0.5 eV/Angs^2 corresponds roughly to 0.01 Ry/Bohr^2
     !%! Should be a parameter: push_over?
     !
     current_step_size = MIN(eigen_step_size,ABS(fpara_tot)/MAX( ABS(lowest_eigval), 0.01_DP ))
     !
     eigenvec(:,:) = eigenvec(:,:)*current_step_size

     displ_vec(:,:) = eigenvec(:,:)
     ! count the number of steps made with the eigenvector
     !
     ieigen = ieigen + 1
     IF ( ieigen == neigen  ) THEN
        ! do a perpendicular relax
        lperp = .true.
        iperp = 0

        ! return to initial number of lanczos steps
        !ilanc = 0  !< initialize it when turn llanczos  = T
        ! 
     ENDIF
     CALL write_struct( at, nat, tau, order, elements, ityp, force, 1.0_DP, iunstruct, struc_format_out, eigenfname )
     CALL write_report( etot_step, force_step, fperp, fpara, lowest_eigval, disp, if_pos, istep, nat,  iunartout, noARTnStep )
      
  END IF



  !
  ! check for convergence of total forces (only after eigevec was obtained)
  !
  IF ( lsaddle_conv ) THEN

     IF ( etot_step > etot_init ) THEN

        ! store the saddle point energy
        etot_saddle = etot_step
        tau_saddle = tau_step
        eigen_saddle = eigenvec
        !
        lsaddle = .true.
        !
        CALL write_struct( at, nat, tau, order, elements, ityp, force, 1.0_DP, iunstruct, struc_format_out, sadfname )
        !
        CALL write_end_report( iunartout, lsaddle, lpush_final, etot_step - etot_init )
        ! 
     ELSE
        CALL write_end_report( iunartout, lsaddle, lpush_final, etot_step )
     ENDIF
  ENDIF
  !
  ! Push to adjacent minima after the saddle point
  !
  IF ( lsaddle ) THEN
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
          eigenvec(:,:) = eigenvec(:,:)/dnrm2(3*nat,eigenvec,1)
        endif
        !
        !
        ! ...If diff Energy is negative
        IF ( etot_step - etot_saddle < frelax_ene_thr ) THEN
           ! we started going downhill ...
           if( .NOT.lrelax )irelax = 0
           lrelax = .true.

        ELSE  !< It is a PUSH_OVER the saddle point  
           disp = EIGN
           iover = iover + 1

           ! ** WARNING **
           if( iover > 2 ) &
             CALL WARNING( iunartout, "PUSH-OVER","Too many push over the saddle point-> PARAM: Push_Over ", [iover])
           !
           CALL write_report( etot_step, force_step, fperp, fpara, lowest_eigval, OVER, if_pos, istep, nat,  iunartout, noARTnStep )
           ! 
           displ_vec(:,:) = fpush_factor*eigenvec(:,:)*eigen_step_size * push_over
        END IF

     ELSE

        lconv = .true.

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

     !
     ArtnStep = noArtnStep
     if( mod(irelax,5) == 0 ) ArtnStep = .true.  !> The 5 can be custom parameter : nrprint
     CALL write_report( etot_step, force_step, fperp, fpara, lowest_eigval, disp, if_pos, istep, nat, iunartout, ARTnStep )

     irelax = irelax + 1
     !write(iunartout, *) "ARTn(lrelax)::increment IRELAX"

     !
     ! check for convergence
     !
     IF ( lforc_conv ) THEN
        IF ( fpush_factor == 1.0 ) THEN

           disp = RELX

           ! restart from saddle point
           tau(:,:) = tau_saddle(:,order(:))
           eigenvec(:,:) = eigen_saddle(:,:)
           !DO i = 1,nat
           !   tau(:,i) = tau_saddle(:,order(i))
           !   eigenvec(:,i) = eigen_saddle(:,i)
           !ENDDO
           lbackward = .true.

           lrelax = .false.
           etot_final = etot_step
           de_back = etot_saddle - etot_final

           CALL write_report( etot_step, force_step, fperp, fpara, lowest_eigval, disp, if_pos, istep, nat, iunartout, .true. )
           call write_inter_report( iunartout, int(fpush_factor), [de_back] )

           ! reverse direction of push
           fpush_factor = -1.0
           irelax = 0

        ELSEIF( .NOT.lend )THEN  !< If already pass before no need to rewrite again

           lconv = .true.
           lend = lconv
           de_fwd = etot_saddle - etot_step

           CALL write_report( etot_step, force_step, fperp, fpara, lowest_eigval, disp, if_pos, istep, nat, iunartout, .true. )
           call write_inter_report( iunartout, int(fpush_factor), [de_back, de_fwd, etot_init, etot_final, etot_step] )

        END IF
        !
     END IF
     !
  END IF RELAX

  !
  ! write restart before the lanczos. If restart happens during lanczos,
  ! the first step will repeat the last lanczos step from before restart.
  ! Reason for this: if we put restart after lanczos, then the restarted lanczos
  ! does not come back to initial point properly.
  !
  !CALL write_restart(restartfname,nat)
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
     CALL write_report( etot_step, force_step, fperp, fpara, lowest_eigval, disp, if_pos, istep, nat,  iunartout, noARTnStep )

     IF (ilanc == 0 ) THEN
        IF ( .not. leigen ) THEN
           !
           ! add_const is set to zero so no constraints are put on the atoms
           !
           ! add_const(:,:) = 0.D0
           ! CALL push_init(nat,idum,push_ids,add_const,1.D0,v_in,'all ')
           ! the push vector is already centered within the push_init subroutine
           ! note push vector is normalized inside lanczos
           v_in(:,:) = eigenvec(:,:)
           ! input a different vector for testing purposes
        ELSE
           ! rescale the lanczos eigenvec back to original size
           v_in(:,:) = eigenvec(:,:)/current_step_size
           ! reset the eigenvalue flag to continue lanczos
           leigen = .false.
        ENDIF
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
        !%! force(:,:) = force(:,:)*if_pos(:,:)
        force_step(:,:) = force_step(:,:)*if_pos(:,:)
     ENDIF
     !
     !
     CALL lanczos( nat, force_step, displ_vec, v_in, dlanc, nlanc, ilanc, lowest_eigval,  eigenvec, push)
     ilanc = ilanc + 1
     iperp = 0

     !
     ! when lanczos converges, nlanc = number of steps it took to converge,
     ! and ilanc = ilanc + 1
     !
     IF ( ilanc > nlanc ) THEN
        !
        ! max number of lanczos steps exceeded ; reset counters
        !
        nlanc = lanc_mat_size
        !ilanc = 0
        !
        llanczos = .false.
        !
        ! check the eigenvalue, if it's lower than threshold use the eigenvector, otherwise push again ...
        !
        IF ( lowest_eigval < eigval_thr ) THEN
           !* make a push with the eigenvector
           !! Next Mstep outside the basin
           leigen = .true.
           lbasin = .false.
           ieigen = 0
           lperp = .false.
           iperp = 0

           !
        ELSE
           !* make an initial push
           !! Next Mstep inside the Basin
           leigen = .false.
           linit  = .true.
           lbasin = .true.
           iperp =  0
           noperp = 0      !> count the init-perp fail
           nperp_step = 0  !> count the out-basin perp relax step
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

     ENDIF

  ENDIF LANCZOS_



  ! ...Increment the ARTn-step
  istep = istep + 1
  ! CALL write_restart(restartfname,nat)
  ! ...Close the output file
  CLOSE (UNIT = iunartout, STATUS = 'KEEP')



  ! ...Unconvert the Force and position
  force = unconvert_force( force )
  !tau = unconvert_length( tau )
  !at = unconvert_length( at )


END SUBROUTINE artn




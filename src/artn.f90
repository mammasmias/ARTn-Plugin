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
SUBROUTINE artn( force, etot_eng, nat, ityp, atm, tau, order, at, if_pos, disp, lconv )
  !----------------------------------------------------------------------------
  !
  ! artn_params for variables and counters that need to be stored in each step
  ! DEFINED IN: artn_params_mod.f90
  !
  USE units
  USE artn_params, ONLY: iunartin, iunartout, iunstruct, &
       lartn, lrelax,lpush_init,lperp,leigen,llanczos, lrestart, lsaddle, lpush_final, lbackward, &
       istep, iperp, ieigen, ipush, ilanc, ismooth, nlanc, if_pos_ct, &
       lowest_eigval, etot_init, etot_step, etot_saddle, etot_final, de_saddle, de_back, de_fwd, &
       npush, neigen, nlanc_init, nsmooth, push_mode, dist_thr, convcrit_init, convcrit_final, &
       fpara_convcrit, eigval_thr, relax_thr, push_step_size, current_step_size, dlanc, eigen_step_size, fpush_factor, &
       push_ids,add_const, push, eigenvec, tau_step, force_step, tau_saddle, eigen_saddle, v_in, &
       VOID, INIT, PERP, EIGN, LANC, RELX, zseed, &
       engine_units, struc_format_out, elements, &
       initialize_artn, write_initial_report, read_restart, write_restart, &
       push_over, ran3
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

  INTEGER,          INTENT(OUT) :: disp       !> Stage for move_mode
  LOGICAL,          INTENT(OUT) :: lconv      !> flag for controlling convergence

  ! -- LOCAL VARIABLES
  !REAL(DP), EXTERNAL :: ran3, dnrm2, ddot     ! lapack functions
  REAL(DP), EXTERNAL :: dnrm2, ddot     ! lapack functions
  INTEGER :: na, icoor, idum                  ! integers for loops
  !
  REAL(DP)  :: force_in(3,nat)                ! stores non-modified force
  REAL(DP)  :: fpara(3,nat)                   ! force parallel to push/eigenvec
  REAL(DP)  :: fpara_tot                      ! total force in parallel direction
  REAL(DP)  :: smoothing_factor               ! mixing factor for smooth transition between eigenvec and push
  REAL(DP)  :: lat(3,3), etot
  INTEGER   :: ios ,i                         ! file IOSTAT
  CHARACTER( LEN=255) :: filin, filout, sadfname, initpfname, eigenfname, restartfname

  integer :: natom

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
  !
  ! fpara_tot is used to scale the magnitude of the eigenvector
  !
  fpara_tot = 0.D0
  !
  filin = 'artn.in'
  filout = 'artn.out'
  sadfname = 'saddle'
  initpfname = 'initp'
  eigenfname = 'latest_eigenvec'
  restartfname = 'artn.restart'
  !
  ! initialize artn
  !
  IF( istep == 0 ) THEN
     ! read the input parameters
     CALL initialize_artn( nat, iunartin, filin )
     !
     ! check if a restart was requested
     !
     IF ( lrestart ) THEN
        OPEN ( UNIT = iunartout, FILE = filout, FORM = 'formatted', STATUS = 'old', POSITION = 'append', IOSTAT = ios )
        WRITE (iunartout, *) "Restarted previous ARTn calculation"
        CLOSE ( UNIT = iunartout, STATUS = 'KEEP')
        CALL read_restart(restartfname,nat)
        !
        ! ...Unconvert Energy/Forces because it will be convert just after
        tau = tau_step
        force = unconvert_force( force_step )
        etot_eng = unconvert_energy( etot_step )
     ELSE
        CALL write_initial_report(iunartout, filout)
        ! store energy of initial state
        etot_init = convert_energy( etot_eng )
     ENDIF
  ENDIF

  ! set initial random seed
  IF( zseed .ne. 0 ) idum = zseed


  ! ...Store & convert original force in a.u.
  force_in = convert_force( force )
  force = force_in
  force_step = force

  ! ...Convert the Energy
  etot = convert_energy( etot_eng )
  etot_step = etot
  ! ...Initialize the displacement
  disp = VOID
  ! store positions of current step
  tau_step = tau

  !
  ! Open the output file for writing
  !
  OPEN ( UNIT = iunartout, FILE = filout, FORM = 'formatted', ACCESS = 'append', STATUS = 'unknown', IOSTAT = ios )

  !
  ! do initial push, then switch off lpush_init, and pass to lperp
  !
  IF ( lpush_init ) THEN
     !
     !=============================
     ! generate initial push vector
     !=============================
     !
     CALL push_init(nat, tau, order, at, idum, push_ids, dist_thr, add_const, push_step_size, push , push_mode)
     !
     ! set up the flags (we did an initial push, now we need to relax perpendiculary)
     lpush_init = .false.
     lperp = .true.
     ! the push counter (controls if we should call lanczos or keep pushing)
     ipush = ipush + 1
     !
     ! modify the force to be equal to the push
     !
     !force(:,:) =  push(:,:) ! Use push in move mode
     disp = INIT
     !
     CALL write_report( etot, force_in, lowest_eigval, disp, if_pos, istep, nat, iunartout )
     !
     CALL write_struct( at, nat, tau, order, elements, ityp, push, 1.0_DP, iunstruct, struc_format_out, initpfname )
     !
     ! Set first eigenvec as push, this is for sure step one of ARTn algo.
     ! eigenvec(:,:) = push(:,:)
     !
     ! generate random initial eigenvec (used as input for Lanczos)
     DO i = 1, nat
        eigenvec(:,i) = (/0.5_DP - ran3(idum),0.5_DP - ran3(idum),0.5_DP - ran3(idum)/) * if_pos(:,i)
     END DO
     eigenvec(:,:) = eigenvec(:,:)/dnrm2(3*nat,eigenvec,1) * dnrm2(3*nat,push,1)

     !
  ELSE IF ( lperp .and. .not. llanczos) THEN
     !
     !=============================
     ! Relax forces perpendicular to eigenvector
     !=============================
     !
     ! If we don't have eigenvector yet, overwrite it with push direction
     ! IF ( .not. leigen ) THEN
     !    ! before lanczos relax perpendiculary with respect to the push
     !    eigenvec(:,:) = push(:,:)
     ! ENDIF
     !
     ! If eigenvalue is good, overwrite push with eigenvec
     !
     IF( leigen ) THEN
        push(:,:) = eigenvec(:,:)
     ENDIF
     !
     ! Subtract parrallel components to push from force.
     ! Modify and overwrite force such that force = force - fpara = fperp
     !
     ! CALL perpforce(force,if_pos,eigenvec,fpara,nat)
     CALL perpforce( force, if_pos, push, fpara, nat )
     !
     IF (MAXVAL(ABS(fpara)) <= fpara_convcrit) THEN
        ! tighten perpendicular convergence criterion
        convcrit_init = convcrit_final
     END IF
     !
     disp = PERP
     !
     iperp = iperp + 1

     !
     ! If the force_perp component are small we continue
     ! to push
     !
     IF (( MAXVAL( ABS(force) )) < convcrit_init ) THEN  ! Here force = force_perp
        !
        ! we reached convergence in the perpendicular direction, go to pushing
        !
        lperp = .false.
        iperp = 0
        !
        IF ( ipush < npush ) THEN
           !
           ! The number of desired pushed still not reached:
           ! push in the specified direction
           !
           ! force(:,:) =  eigenvec(:,:)
           force(:,:) =  push(:,:)
           disp = INIT
           !
           ipush = ipush + 1
           !
           lperp = .true.
           !
           CALL write_report(etot,force_in, lowest_eigval, disp, if_pos, istep, nat,  iunartout)
           !
        ELSE IF ( ipush >= npush  ) THEN
           !
           ! Number of desired pushes has been reached:
           ! use true force & start lanczos
           !
           force(:,:) = force_in(:,:)
           llanczos = .true.
        END IF
     ELSE
        !
        ! keep doing lperp
        !
        CALL write_report(etot,force_in, lowest_eigval, disp, if_pos, istep, nat,  iunartout)
        !
     END IF
     !
  ELSE IF ( leigen .and. .not. lperp ) THEN
     !=============================
     ! Push in the direction of the lowest eigenvector
     !=============================
     !
     ! leigen is always .true. after we obtain a good eigenvector
     ! except during lanczos iterations
     !
     ! if we have a good lanczos eigenvector use it as push vector
     !
     force(:,:) = eigenvec(:,:)
     disp = EIGN
     ieigen = ieigen + 1
     !
     CALL write_report(etot,force_in, lowest_eigval, disp, if_pos, istep, nat,  iunartout)

     ! count the number of steps made with the eigenvector
     CALL write_struct( at, nat, tau, order, elements, ityp, force, 1.0_DP, iunstruct, struc_format_out, eigenfname)
     !
     IF ( ieigen == neigen  ) THEN
        ! do a perpendicular relax
        lperp = .true.
        iperp = 0
        ! return to initial number of lanczos steps
        ilanc = 0
     ENDIF
  END IF
  !
  ! check for convergence of total forces (only after eigevec was obtained)
  !
  IF (MAXVAL(ABS(force_in*if_pos)) < convcrit_final .AND. leigen .AND. .NOT. lsaddle ) THEN
     force(:,:) = force_in(:,:)
     IF ( etot > etot_init ) THEN
        ! store the saddle point energy
        etot_saddle = etot
        !tau_saddle = tau
        do i = 1, nat
           tau_saddle(:,order(i)) = tau(:,i) !> The list follows the atomic order
           eigen_saddle(:,order(i)) = eigenvec(:,i)
        enddo
        !do i = 1, nat
        !   print*, i, tau_saddle(:,i)
        !enddo
        !
        lsaddle = .true.
        !
        !CALL write_struct( at, nat, tau, order, atm, ityp, force, 1.0_DP, iunstruct, struc_format_out, sadfname)
        CALL write_struct( at, nat, tau, order, elements, ityp, force, 1.0_DP, iunstruct, struc_format_out, sadfname)
        !
        WRITE (iunartout,'(5X, "--------------------------------------------------")')
        WRITE (iunartout,'(5X, "    *** ARTn found a potential saddle point ***   ")')
        WRITE (iunartout,'(5X, "--------------------------------------------------")')
        !WRITE (iunartout,'(15X,"E_final - E_initial =", F12.5," eV")') (etot - etot_init)*RY2EV
        WRITE (iunartout,'(15X,"E_final - E_initial =", F12.5," eV")') unconvert_energy((etot - etot_init)) !*RY2EV
        WRITE (iunartout,'(5X, "--------------------------------------------------")')

        !write( iunartout, * ) " ARTn::", MAXVAL(ABS(force_in*if_pos)), convcrit_final, leigen, lsaddle

        IF ( lpush_final ) THEN
           WRITE (iunartout, '(5X,"       *** Pushing to adjacent minima  ***      ")')
           WRITE (iunartout,'(5X, "------------------------------------------------")')
        ENDIF
     ELSE
        WRITE (iunartout,'(5X, "--------------------------------------------------")')
        WRITE (iunartout,'(5X, "        *** ARTn saddle search failed  ***        ")')
        WRITE (iunartout,'(5X, "--------------------------------------------------")')
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
          eigenvec(:,:) = eigen_saddle(:,order(:))
          lbackward = .false.
        else
          eigenvec(:,:) = eigenvec(:,:)/dnrm2(3*nat,eigenvec,1)
        endif
        !
        force(:,:) = fpush_factor*eigenvec(:,:)*eigen_step_size
        !
        ! ...If diff Energy is negative
        IF ( etot - etot_saddle < relax_thr ) THEN
           ! we started going downhill ...
           lrelax = .true.
        ELSE
           !
           disp = EIGN
           !
           CALL write_report(etot,force_in, lowest_eigval, disp, if_pos, istep, nat,  iunartout)
        END IF
     ELSE
        lconv = .true.
     ENDIF
  ENDIF
  !
  ! perform a FIRE relaxation (only for reaching adjacent minima)
  !
  IF ( lrelax ) THEN
     !
     ! reset force
     !
     disp = RELX
     force(:,:) = force_in(:,:)
     !
     CALL write_report(etot, force_in, lowest_eigval, disp, if_pos, istep, nat, iunartout)
     !
     ! check for convergence
     !
     IF ( MAXVAL(ABS(force_in(:,:)*if_pos(:,:))) <= convcrit_final  ) THEN
        IF ( fpush_factor == 1.0 ) THEN

           disp = RELX
           ! reverse direction of push
           fpush_factor = -1.0

           ! restart from saddle point
           do i = 1,nat
           tau(:,i) = tau_saddle(:,order(i))
           eigenvec(:,i) = eigen_saddle(:,order(i))
           enddo
           lbackward = .true.

           !do i = 1,nat
           ! print*, i, order(i), tau(:,i)
           !enddo

           lrelax = .false.
           etot_final = etot
           de_back = etot_saddle - etot_final

           WRITE (iunartout,'(5X, "--------------------------------------------------")')
           WRITE (iunartout,'(5X, "    *** ARTn found adjacent minimum ***   ")')
           WRITE (iunartout,'(5X, "--------------------------------------------------")')
           !WRITE (iunartout,'(15X,"backward E_act =", F12.5," eV")') de_back*RY2EV
           WRITE (iunartout,'(15X,"backward E_act =", F12.5," eV")') unconvert_energy(de_back) !*RY2EV
           WRITE (iunartout,'(5X, "--------------------------------------------------")')
        ELSE
           lconv = .true.
           de_fwd = etot_saddle - etot
           WRITE (iunartout,'(5X, "--------------------------------------------------")')
           WRITE (iunartout,'(5X, "    *** ARTn converged to initial minimum ***   ")')
           WRITE (iunartout,'(5X, "--------------------------------------------------")')
           !WRITE (iunartout,'(15X,"forward  E_act =", F12.5," eV")') de_fwd*RY2EV
           !WRITE (iunartout,'(15X,"backward E_act =", F12.5," eV")') de_back*RY2EV
           !WRITE (iunartout,'(15X,"reaction dE    =", F12.5," eV")') (etot-etot_final) *RY2EV
           !WRITE (iunartout,'(15X,"dEinit - dEfinal    =", F12.5," eV")') (etot_init-etot) *RY2EV
           WRITE (iunartout,'(15X,"forward  E_act =", F12.5," eV")') unconvert_energy(de_fwd) !*RY2EV
           WRITE (iunartout,'(15X,"backward E_act =", F12.5," eV")') unconvert_energy(de_back) !*RY2EV
           WRITE (iunartout,'(15X,"reaction dE    =", F12.5," eV")') unconvert_energy((etot-etot_final)) ! *RY2EV
           WRITE (iunartout,'(15X,"dEinit - dEfinal    =", F12.5," eV")') unconvert_energy((etot_init-etot)) ! *RY2EV
           WRITE (iunartout,'(5X, "--------------------------------------------------")')
        END IF
        !
     END IF
     !
  END IF

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
  IF ( llanczos ) THEN
     !
     !=============================
     ! Perform Lanczos algo, one step at a time
     !=============================
     !
     disp = LANC
     CALL write_report(etot,force_in, lowest_eigval, disp, if_pos, istep, nat,  iunartout)
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
        force(:,:) = force(:,:)*if_pos(:,:)
     ENDIF
     !
     !
     CALL lanczos( nat, force, v_in, dlanc, nlanc, ilanc, lowest_eigval,  eigenvec, push)

     iperp = 0

     !
     ! when lanczos converges, nlanc = number of steps it took to converge,
     ! and ilanc = ilanc + 1
     !
     IF ( ilanc > nlanc ) THEN
        !
        ! max number of lanczos steps exceeded ; reset counters
        !
        nlanc = nlanc_init
        ilanc = 0
        !
        llanczos = .false.
        !
        ! check the eigenvalue, if it's lower than threshold use the eigenvector, otherwise push again ...
        !
        IF ( lowest_eigval < eigval_thr ) THEN
           !
           smoothing_factor = 1.0_DP*ismooth/nsmooth
           !
           !
           fpara_tot = ddot(3*nat,force_in(:,:),1,eigenvec(:,:),1)
           !
           eigenvec(:,:) = (1.0_DP - smoothing_factor)*push(:,:) &
                -SIGN(1.0_DP,fpara_tot)*smoothing_factor*eigenvec(:,:)
           !
           IF ( ismooth < nsmooth) ismooth = ismooth + 1
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
           !
           leigen = .true.
           !
           ieigen = 0
           lperp = .false.
           iperp = 0
           !
        ELSE
           ! make an initial push
           leigen = .false.
           lperp = .true.
           iperp =  0
           ipush = ipush - 1
           !
        ENDIF
     ENDIF
  ENDIF

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




! 
! 
! Main ARTn plugin subroutine:
!        modifies the input force to perform the ARTn algorithm 
!------------------------------------------------------------------------------
SUBROUTINE artn(force,etot,nat,ityp,atm,tau,at,alat,istep,if_pos,vel,dt,fire_alpha_init,lconv,prefix,tmp_dir)
  !----------------------------------------------------------------------------
  !
  ! artn_params for variables and counters that need to be stored in each step   
  ! DEFINED IN: artn_params_mod.f90
  ! 
  USE artn_params, ONLY: DP, RY2EV,B2A, iunartin, iunartout, iunstruct, &
       lrelax,lpush_init,lperp,leigen,llanczos, lsaddle, lpush_final, &
       iperp, ieigen, ipush, ilanc, ismooth, nlanc, if_pos_ct, &
       lowest_eigval, etot_init, etot_saddle, etot_final, de_saddle, de_back, de_fwd, &
       npush, neigen, nlanc_init, nsmooth, push_mode, dist_thr, convcrit_init, convcrit_final, &
       fpara_convcrit, eigval_thr, relax_thr, push_step_size, current_step_size, dlanc, eigen_step_size, fpush_factor, &
       push_ids,add_const, push, eigenvec, tau_saddle, initialize_artn
  ! 
  IMPLICIT NONE
  REAL(DP), INTENT(INOUT) :: force(3,nat)     ! force calculated by the engine
  REAL(DP), INTENT(INOUT) :: vel(3,nat)       ! velocity of previous FIRE step
  REAL(DP), INTENT(INOUT) :: tau(3,nat)       ! atomic positions (needed for output only)
  REAL(DP), INTENT(IN) ::    etot             ! total energy in current step
  REAL(DP), INTENT(IN) ::    dt               ! default time step in FIRE  
  REAL(DP), INTENT(IN) ::    fire_alpha_init  ! initial value of alpha in FIRE 
  REAL(DP), INTENT(IN) ::    alat             ! lattice parameter of QE
  REAL(DP), INTENT(IN) ::    at(3,nat)        ! lattice parameters in alat units 
  INTEGER,  INTENT(IN) ::    nat              ! number of atoms
  INTEGER,  INTENT(IN) ::    ityp(nat)        ! atom types
  INTEGER,  INTENT(IN) ::    istep            ! current step
  INTEGER,  INTENT(IN) ::    if_pos(3,nat)    ! coordinates fixed by engine 
  CHARACTER(LEN=3),   INTENT(IN) :: atm(*)    ! name of atom corresponding to ityp
  CHARACTER(LEN=255), INTENT(IN) :: tmp_dir   ! scratch directory of engine 
  CHARACTER(LEN=255), INTENT(IN) :: prefix    ! prefix for scratch files of engine 
  LOGICAL, INTENT(OUT) :: lconv               ! flag for controlling convergence 
  REAL(DP), EXTERNAL :: ran3, dnrm2, ddot     ! lapack functions 
  INTEGER :: na, icoor, idum                  ! integers for loops 
  !
  REAL(DP)  :: force_in(3,nat)                ! stores non-modified force 
  REAL(DP)  :: fpara(3,nat)                   ! force parallel to push/eigenvec   
  REAL(DP)  :: v_in(3,nat)                    ! input vector for lanczos 
  REAL(DP)  :: fpara_tot                      ! total force in parallel direction
  REAL(DP) ::  smoothing_factor               ! mixing factor for smooth transition between eigenvec and push
  INTEGER   :: ios                            ! file IOSTAT  
  CHARACTER( LEN=255) :: filin, filout, sadfname, initpfname, eigenfname, restartfname
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
  ! store original force
  force_in(:,:) = force(:,:)
  !
  ! fpara_tot is used to scale the magnitude of the eigenvector 
  ! 
  fpara_tot = 0.D0
  !
  filin = 'artn.in'
  filout = 'artn.out'
  sadfname = 'saddle.xsf'
  initpfname = 'initp.xsf'
  eigenfname = 'latest_eigenvec.xsf'
  !
  ! initialize artn 
  !  
  IF (istep == 0 ) THEN
     ! read the input parameters 
     CALL initialize_artn(nat,iunartin,iunartout,filin,filout)
     ! store the total energy of the initial state
     etot_init = etot 
  ENDIF
  ! 
  ! Open the output file for writing   
  ! 
  OPEN ( UNIT = iunartout, FILE = filout, FORM = 'formatted', ACCESS = 'append', STATUS = 'unknown', IOSTAT = ios )
  
  ! 
  ! start initial push
  !
  IF ( lpush_init ) THEN
     !
     ! initial push 
     !
     CALL push_init(nat, tau, at, alat, idum, push_ids, dist_thr, add_const, push_step_size, push ,push_mode)
     ! set up the flags (we did an initial push, now we need to relax perpendiculary) 
     lpush_init = .false.
     lperp = .true.
     ! the push counter (controls if we should call lanczos or keep pushing)
     ipush = ipush + 1
     ! 
     ! modify the force to be equal to the push
     !
     force(:,:) =  push(:,:)
     CALL move_mode( nat, dlanc, v_in, force, &
          vel, fire_alpha_init, dt,  &
          iperp, push, 'eign', prefix, tmp_dir)
     !
     CALL write_report(etot,force_in, lowest_eigval, 'push' , if_pos, istep, nat,  iunartout)
     !
     CALL write_struct(alat, at, nat, tau, atm, ityp, force, 1.0_DP, iunstruct, 'xsf', initpfname)
     ! 
  ELSE IF ( lperp ) THEN
     !
     !  subtract parrallel components to push from force
     !
     IF ( .not. leigen) THEN
        ! before lanczos relax perpendiculary with respect to the push otherwise 
        eigenvec(:,:) = push(:,:) 
     ENDIF
     ! 
     CALL perpforce(force,if_pos,eigenvec,fpara,nat)
     ! 
     IF (MAXVAL(ABS(fpara)) <= fpara_convcrit) THEN
        ! tighten perpendicular convergance criterion
        convcrit_init = convcrit_final  
     END IF
     !
     CALL move_mode( nat,  dlanc, v_in, force, &
          vel, fire_alpha_init, dt,  &
          iperp, eigenvec, 'perp', prefix, tmp_dir)
     ! 
     iperp = iperp + 1
     !
     IF (( MAXVAL( ABS(force) )) < convcrit_init ) THEN
        !
        ! we reached convergence in the perpendicular direction of push
        !
        lperp = .false.
        iperp = 0
        ! 
        IF ( ipush < npush ) THEN
           ! continue pushing in the specified direction
           force(:,:) =  eigenvec(:,:)
           CALL move_mode( nat, dlanc, v_in, force, &
                vel, fire_alpha_init, dt,  &
                iperp, eigenvec, 'eign', prefix, tmp_dir)
           ! 
           ipush = ipush + 1
           ! 
           lperp = .true.
           ! 
           CALL write_report(etot,force_in, lowest_eigval, 'push' , if_pos, istep, nat,  iunartout)
           ! 
        ELSE IF ( ipush >= npush  ) THEN
           ! regenerate force & start lanczos
           force(:,:) = force_in(:,:)
           llanczos = .true.
        END IF
     ELSE
        CALL write_report(etot,force_in, lowest_eigval, 'perp' , if_pos, istep, nat,  iunartout)       
     END IF
     ! leigen is always .true. after we obtain a good eigenvector
     ! except during lanczos iterations 
  ELSE IF ( leigen .and. .not. lperp ) THEN
     !
     ! if we have a good lanczos eigenvector use it
     !  
     force(:,:) = eigenvec(:,:)
     CALL move_mode( nat, dlanc, v_in, force, &
          vel, fire_alpha_init, dt,  &
          iperp, eigenvec, 'eign', prefix, tmp_dir)
     ! update eigenstep counter 
     ieigen = ieigen + 1
     !      
     CALL write_report(etot,force_in, lowest_eigval, 'eign' , if_pos, istep, nat,  iunartout)
     ! count the number of steps made with the eigenvector
     CALL write_struct(alat, at, nat, tau, atm, ityp, force, 1.0_DP, iunstruct, 'xsf', eigenfname)
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
  ! check if we should perform the lanczos algorithm
  !
  IF ( llanczos ) THEN 
     !
     CALL write_report(etot,force_in, lowest_eigval, 'lanc' , if_pos, istep, nat,  iunartout)
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
     CALL lanczos( nat, force, vel, fire_alpha_init, dt, &
          v_in, dlanc, nlanc, ilanc, lowest_eigval,  eigenvec, push, prefix, tmp_dir)
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
           ! 
           current_step_size = MIN(eigen_step_size,ABS(fpara_tot)/MAX(ABS(lowest_eigval),0.01_DP))
           ! 
           eigenvec(:,:) = eigenvec(:,:)*current_step_size
           !           
           leigen = .true.
           !
           ieigen = 0
           lperp = .false.
           iperp = 0
             
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
  !
  ! check for convergence of total forces (only after eigevec was obtained)  
  !
  IF (MAXVAL(ABS(force_in*if_pos)) < convcrit_final .AND. leigen .AND. .NOT. lsaddle ) THEN
     force(:,:) = force_in(:,:)
     IF ( etot > etot_init ) THEN
        ! store the saddle point energy  
        etot_saddle = etot
        tau_saddle = tau
        ! 
        lsaddle = .true.
        CALL write_struct(alat, at, nat, tau, atm, ityp, force, 1.0_DP, iunstruct, 'xsf', sadfname)
        !
        WRITE (iunartout,'(5X, "--------------------------------------------------")')
        WRITE (iunartout,'(5X, "    *** ARTn found a potential saddle point ***   ")')
        WRITE (iunartout,'(5X, "--------------------------------------------------")')
        WRITE (iunartout,'(15X,"E_final - E_initial =", F12.5," eV")') (etot - etot_init)*RY2EV
        WRITE (iunartout,'(5X, "--------------------------------------------------")')
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
        ! set convergence and other flags to false ...
        lconv = .false.
        lperp = .false.
        leigen = .false. 
        ! 
        ! normalize eigenvector 
        eigenvec(:,:) = eigenvec(:,:)/dnrm2(3*nat,eigenvec,1)
        !
        force(:,:) = fpush_factor*eigenvec(:,:)*eigen_step_size
        !
        IF ( etot - etot_saddle < relax_thr ) THEN
           ! we started going downhill ... 
           lrelax = .true.
        ELSE
           ! 
           CALL move_mode( nat, dlanc, v_in, force, &
          vel, fire_alpha_init, dt,  &
          iperp, eigenvec, 'eign', prefix, tmp_dir)
           !
           CALL write_report(etot,force_in, lowest_eigval, 'eign' , if_pos, istep, nat,  iunartout)
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
     force(:,:) = force_in(:,:)
     ! 
     CALL write_report(etot, force_in, lowest_eigval, 'relx', if_pos, istep, nat, iunartout)
     !
     ! check for convergence 
     ! 
     IF ( MAXVAL(ABS(force_in(:,:)*if_pos(:,:))) <= convcrit_final  ) THEN
        IF ( fpush_factor == 1.0 ) THEN
           ! make sure QE doesn't converge
           forc_conv_thr_qe = 1.0D-8
           ! reverse direction of push 
           fpush_factor = -1.0
           ! restart from saddle point
           tau(:,:) = tau_saddle(:,:)
           lrelax = .false.
           etot_final = etot 
           de_back = etot_saddle - etot_final 
           WRITE (iunartout,'(5X, "--------------------------------------------------")')
           WRITE (iunartout,'(5X, "    *** ARTn found adjacent minimum ***   ")')
           WRITE (iunartout,'(5X, "--------------------------------------------------")')
           WRITE (iunartout,'(15X,"backward E_act =", F12.5," eV")') de_back*RY2EV
           WRITE (iunartout,'(5X, "--------------------------------------------------")')
        ELSE
           lconv = .true.
           de_fwd = etot_saddle - etot
           WRITE (iunartout,'(5X, "--------------------------------------------------")')
           WRITE (iunartout,'(5X, "    *** ARTn converged to initial minimum ***   ")')
           WRITE (iunartout,'(5X, "--------------------------------------------------")')
           WRITE (iunartout,'(15X,"forward  E_act =", F12.5," eV")') de_fwd*RY2EV
           WRITE (iunartout,'(15X,"backward E_act =", F12.5," eV")') de_back*RY2EV
           WRITE (iunartout,'(15X,"reaction dE    =", F12.5," eV")') (etot-etot_final) *RY2EV
           WRITE (iunartout,'(15X,"dEinit - dEfinal    =", F12.5," eV")') (etot_init-etot) *RY2EV
           WRITE (iunartout,'(5X, "--------------------------------------------------")')
        END IF
        ! 
     END IF
     ! 
  END IF

  CLOSE (UNIT = iunartout, STATUS = 'KEEP') 
END SUBROUTINE artn 

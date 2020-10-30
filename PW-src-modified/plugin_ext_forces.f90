!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE plugin_ext_forces(force)
  !----------------------------------------------------------------------------
  !
  !
  USE mp,               ONLY : mp_bcast
  USE mp_images,        ONLY : intra_image_comm
  USE io_global,        ONLY : stdout, ionode, ionode_id
  USE kinds,            ONLY : DP
  !
  USE plugin_flags
  ! modifications start here
  ! we take some stuff from QE
  USE ions_base,     ONLY : nat, tau, if_pos, ityp, amass
  USE cell_base,     ONLY : alat, at
  USE control_flags, ONLY : istep
  USE dynamics_module, ONLY : vel, acc, dt, fire_alpha_init 
  USE io_files,      ONLY : prefix,seqopn,tmp_dir
  !
  IMPLICIT NONE
  !
  ! local varibles
  !
  REAL(DP), INTENT(INOUT)  :: force(3,nat)
  INTEGER :: na, idum, nlanc, nlanccalls, npush, neigenstep, io
  INTEGER, PARAMETER :: iunart = 50, iunlanc = 51
  INTEGER, PARAMETER :: npushmin = 3, nlanccallsmax = 1, neigenstepmax = 1
  INTEGER ::  nlanciter, nlanciter_init, natpush, istepperp, if_pos_ct, i
  INTEGER ::  ihess,jhess
  REAL(DP) :: dlanc
  ! flags for controlling ARTn
  LOGICAL :: lpush_init, leigen, llanczos, lperp, file_exists, lhess, lpush_list
  !
  REAL(DP)  :: force_in(3,nat), push(3,nat), fpara(3,nat), pushdir(3,nat), eigenvec(3,nat), v_in(3,nat)
  REAL(DP) :: hess(3*nat,3*nat), hess_symm(3*nat,3*nat), hess_eigvals(3*nat), vibmode(3,nat)
  REAL(DP) :: fpara_tot, fperp_tot, force_tot
  REAL(DP) :: convcrit,init_step_size, step_size, current_step_size
  REAL(DP), ALLOCATABLE :: add_const(:,:)
  INTEGER :: if_pos_init(3,nat)
  INTEGER, ALLOCATABLE :: push_ids(:)
  REAL(DP) :: ran3,dnrm2, ddot
  REAL(DP) :: lowest_eigval
  convcrit = 1.0d-2
  init_step_size = 1.5
  step_size = 0.5
  !
  ! The basic idea is:
  ! (1) push atoms in the direction specified by user & relax in the perpendicular direction;
  ! (2) use the lanczos algorithm calculate the lowest eigenvalue/eigenvec
  ! (3) a negative eigenvalue, update push direction otherwise push again
  ! (4) follow the lanczos direction twoard the saddle point
  !
  ! Set up control flags for the first step
  !
  ! lpush_init = .true.
  lhess = .false.
  lpush_init = .true.
  lperp = .false.
  llanczos = .false.
  leigen = .false.
  !
  ! set up counters
  !
  eigenvec(:,:) = 0.D0
  push(:,:) = 0.D0
  pushdir(:,:) = 0.D0
  lowest_eigval  = 0.D0
  neigenstep = 0
  nlanc = 0
  istepperp = 0
  nlanccalls = 1
  nlanciter = 16
  nlanciter_init = 16
  natpush = 1
  if_pos_ct = 0
  ! counters for the calculation of the hessian
  ihess = 1
  jhess = 1
  
  WRITE (*,*) "ARTn read FIRE Parms:",dt,fire_alpha_init 
  ! store original force
  force_in(:,:) = force(:,:)
  ALLOCATE(add_const(4,nat),source=0.D0)
  ALLOCATE(push_ids(natpush))
  ! define which atoms are to be pushed and the constraints ...
  push_ids = (/1/)
  add_const(:,1) = (/1.0, 0.0, 0.0, 0.0/)
  lpush_list = .true.
  !
  fperp_tot = 0.D0
  fpara_tot = 0.D0
  force_tot = 0.D0
  !
  ! First read the scratch file and update flags
  !
  CALL seqopn( iunart, 'artn', 'FORMATTED', file_exists )
  IF ( file_exists ) THEN
     !
     ! ... the scratch file is read
     !
     READ( UNIT =  iunart, FMT = *)  lpush_init, lperp, npush, nlanc, nlanciter, &
           nlanccalls, neigenstep, istepperp, ihess, jhess, llanczos, leigen, lhess,  &
           push(:,:), eigenvec(:,:), lowest_eigval, &
           current_step_size, convcrit, pushdir(:,:)
     !
     CLOSE( UNIT =  iunart, STATUS = 'KEEP' )
  ELSE
     CLOSE( UNIT = iunart, STATUS = 'DELETE')
  END IF
  ! for now set the hess calc flag here
  IF ( lhess ) THEN
     CALL calc_hessian(force,ihess,jhess,hess)
     write (*,*) "Hessian calc:",ihess,jhess
     IF ( ihess == 3*nat .and. jhess > 3*nat ) THEN
        lhess = .false.
        write (*,*) "Hessian calc finished:"
        ! symmetrize hess
        DO ihess=1,3*nat
           DO jhess=1,3*nat
              hess_symm(ihess,jhess) = 0.5_DP*(hess(ihess,jhess) + hess(jhess,ihess))
           ENDDO
        ENDDO
        DO na=1,3*nat
           write (*,*) hess_symm(:,na)
        ENDDO
        CALL diag(3*nat, hess_symm, hess_eigvals, 1 )
        write (*,*) "Hessian eigenvalues:",hess_eigvals(:)
        write (*,*) "Hessian eigenvecs:", hess_symm(:,:)
        OPEN (UNIT = 405)
        DO na=1,3*nat
           vibmode(1,:) = hess_symm(1:nat,na)
           vibmode(2,:) = hess_symm(nat+1:2*nat,na)
           vibmode(3,:) = hess_symm(2*nat+1:3*nat,na)
           write (405,*) "PRIMCOORD",na
           write (405,*) nat, "1"
           DO i=1,nat
              write (405,*) tau(:,i)*alat*0.52917720859_DP,vibmode(:,i)
           ENDDO
        ENDDO

        CLOSE (UNIT = 405, STATUS = 'KEEP')
        STOP 1
     ENDIF
  ENDIF

  WRITE (*,*) "ARTn: initial READ at step:",istep, "Initial push",&
       & lpush_init, "Perp relax:", lperp, "Lanczos:", llanczos, &
       & "Good eigenvalue:", leigen
  !
  ! start initial push
  !
  IF ( lpush_init ) THEN
     !
     ! check type of initial push
     !
     WRITE (*,*) "ARTn initial push at step:",istep
     IF ( lpush_list ) THEN
        CALL push_init_list(nat, natpush, idum, push_ids, add_const, init_step_size, push )
     ELSE
        CALL push_init(nat,idum, init_step_size, push)
     ENDIF
     ! set up the flags (we did an initial push, now we need to relax perpendiculary
     lpush_init = .false.
     lperp = .true.
     istepperp = 0
     ! the push counter (controls if we should call lanczos or keep pushing)
     npush = 1
     force(:,:) =  push(:,:)
     WRITE (*,*) "Starting ARTn with initial push: ", push(:,:)
  ELSE IF ( lperp .and. .not. llanczos ) THEN
     !
     !  subtract parrallel components to push from force
     !
     WRITE (*,*) "ARTn: perpedicular relax at step:", istep
     IF ( leigen ) THEN
        WRITE (*,*) "Doing a perpendicular relax with obtained eigenvec:", eigenvec(:,:)
        push(:,:) = eigenvec(:,:)
     ENDIF
     CALL perpforce(force,push,fpara,nat)
     ! CALL perpmove(nat,istepperp,push)
     CALL move_mode( nat, dlanc, v_in, force, &
          vel, acc, fire_alpha_init, dt, & 
          istepperp, push, 'perp')
     istepperp = istepperp + 1
     !
     !
     IF (( MAXVAL( ABS(force) )) < convcrit ) THEN
        !
        ! we reached convergence in the perpendicular direction of push
        !
        WRITE (*,*) "ARTn : reached convergence in perpendicular direction at step:", istep, &
             "maxval of force was:",( MAXVAL( ABS(force) ))
        lperp = .false.
        istepperp = 0
        IF ( npush < npushmin ) THEN
           ! continue pushing in the specified direction
           WRITE (*,*) "ARTn: Continuing in push direction"
           force(:,:) =  push(:,:)
           npush = npush + 1
           lperp = .true.

        ELSE IF ( npush >= npushmin  ) THEN
           ! start lanczos algorithm
           force(:,:) = force(:,:) + fpara(:,:)
           write (*,*) "ARTn: start lanzcos in step:",istep
           llanczos = .true.
        END IF
     END IF
  ELSE IF ( leigen .and. .not. lperp ) THEN
     !
     ! if we have a good lanczos eigenvector use it
     !
     ! if_pos(:,:) = if_pos_init(:,:)
     fpara_tot = ddot(3*nat,force(:,:),1,eigenvec(:,:),1)
     WRITE (*,*) "ARTn: Magnitude of parallel force at push:",fpara_tot
     IF (ABS(fpara_tot) <= 0.005_DP) THEN
        WRITE (*,*) "ARTn: Tightening perpendicular convergence criterion"
        ! decrease perpendicular convergance criterion
        convcrit = 0.1_DP*convcrit
     END IF
     ! rescale the eigenvector according to the current force in the parallel direction
     current_step_size = MIN(step_size,ABS(fpara_tot)/MAX(ABS(lowest_eigval),0.5_DP))
     eigenvec(:,:) = -SIGN(1.0D0,fpara_tot)*eigenvec(:,:)*current_step_size
     write (*,*) "ARTn: eigenvec rescaled by:",current_step_size
     write (*,*) "ARTn: used lanczos eigenvector: ", eigenvec(:,:)
     force(:,:) = eigenvec(:,:)
     ! call eigenmove(force,nat)
     CALL move_mode( nat, dlanc, v_in, force, &
          vel, acc, fire_alpha_init, dt,  & 
          istepperp, push, 'eign')
     neigenstep = neigenstep + 1
     ! count the number of steps made with the eigenvector
     IF ( neigenstep == neigenstepmax  ) THEN
        write (*,*) "ARTn: Made the maximum number of steps in the eigenvec direction"
        ! do a perpendicular relax
        ! return to initial number of lanczos steps
        lperp = .true.
        istepperp = 0
        nlanc = 0
     ENDIF
  END IF
  !
  ! check if we should perform the lanczos algorithm
  !
  IF ( llanczos ) THEN
     !
     WRITE (*,*) "ARTn: Called lanczos, at step:",istep,"iteration:", nlanc, "of max:",nlanciter
     IF (nlanc == 0 ) THEN
        write(555,*) 'calling lanczos with size:',nlanciter
        IF ( .not. leigen ) THEN
           ! generate a random lanczos vector for input
           pushdir(:,:) = push(:,:)
           CALL push_init(nat, idum, 1.D0, v_in )
           CALL center (v_in(:,:), nat)
        ELSE
           ! rescale the lanczos eigenvec back to original size
           pushdir(:,:)  = eigenvec(:,:)/current_step_size
           v_in(:,:) = pushdir(:,:)
           ! CALL center (v_in(:,:), nat)
           ! reset the eigenvalue flag to continue lanczos
           leigen = .false.
        ENDIF
     ENDIF
     ! write (*,*) "Coordinates at lanczos step:", nlanc, tau(:,:)
     ! modify the force according to the constraint
     IF ( ANY(if_pos(:,:) == 0) ) THEN
        DO na=1,nat
           DO i=1,3
              IF (if_pos(i,na) == 1 ) if_pos_ct = if_pos_ct + 1
           ENDDO
        END DO
        IF ( if_pos_ct < nlanciter .and. if_pos_ct /= 0 ) nlanciter = if_pos_ct

        WRITE (*,*) "ARTn: applying constraint"
        v_in(:,:) = v_in(:,:)*if_pos(:,:)
        force(:,:) = force(:,:)*if_pos(:,:)
     ENDIF
     CALL lanczos( nat, force, vel, acc, fire_alpha_init, dt, &
          v_in, nlanciter, nlanc, lowest_eigval,  eigenvec, pushdir )

     !
     ! when lanczos converges, nlanciter = number of steps it took to converge,
     ! and nlanc = nlanciter + 1
     !
     IF ( nlanc > nlanciter ) THEN
        !
        ! max number of lanczos steps exceeded ; reset counters
        !
        nlanciter = nlanciter_init
        !
        llanczos = .false.
        ! reset constraints
        ! write(*,*) "Reseting if_pos:",if_pos_init(:,:)
        ! set lanczos counter to 0
        nlanc = 0
        !
        ! check the eigenvalue, if it's negative continue take the projection, otherwise push again ...
        !
        IF ( lowest_eigval < 0.0_DP ) THEN
          !  WRITE (*,*) "Coordinates of lanczos end:", tau(:,:)
           WRITE (*,*) "ARTn: found negative eigenvalue, using eigvec:", norm2(eigenvec)
           write(*,'(3f13.6)') eigenvec(:,:)
           ! CALL center(eigenvec,nat)
           WRITE (406,*) "PRIMCOORD"
           WRITE (406,*) nat, "1"
           DO na=1,nat
              write (406,*) tau(:,na)*alat*0.52917720859_DP,eigenvec(:,na)
           ENDDO
           CLOSE (UNIT = 406,STATUS = 'KEEP')
           leigen = .true.
           ! lhess = .true.
           neigenstep = 0
           lperp = .false.
           istepperp = 0
        ELSE
           ! npush = 0
           leigen = .false.
           ! lhess = .true.
           lperp = .true.
           istepperp =  0
           ! force(:,:) = push(:,:)
           npush = npush - 1
        ENDIF
     ENDIF
  ENDIF
  !
  ! report forces during perpendicular relaxation
  !
  IF (  lperp  ) THEN
     IF ( leigen ) THEN
        CALL report_force(force_in,eigenvec,nat,force_tot,fperp_tot,fpara_tot)
     ELSE
        CALL report_force(force_in,push,nat,force_tot,fperp_tot,fpara_tot)
     ENDIF
  ENDIF

  CALL seqopn( iunart, 'artn', 'FORMATTED', file_exists )


  WRITE (UNIT = iunart, FMT = * ) lpush_init, lperp, npush, nlanc, nlanciter, nlanccalls, &
       & neigenstep, istepperp, ihess, jhess, llanczos, leigen, lhess, &
       & push(:,:), eigenvec(:,:), lowest_eigval, &
       & current_step_size, convcrit, pushdir(:,:)
  CLOSE (UNIT = iunart, STATUS = 'KEEP')

  DEALLOCATE(add_const,push_ids)

  ! write(456,*) nat
  ! write(456,*) 'Lattice="',at(1,:)*alat*0.529177_DP, at(2,:)*alat*0.529177_DP, at(3,:)*alat*0.529177_DP,'"'
  ! do i = 1, nat
  !    write(456,*) ityp(i),tau(:,i)*alat*0.529177_DP, force(:,i)
  ! end do
  ! flush(456)

END SUBROUTINE plugin_ext_forces

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

SUBROUTINE hessmove(force, nat, moves_finished )
  !
  ! subroutine that rescales the force so it corresponds to a specific move in tau ...
  !
  USE kinds, ONLY : DP
  USE ions_base, ONLY : tau
  USE cell_base, ONLY : alat
  USE control_flags,               ONLY : istep
  USE ener,               ONLY : etot
  USE constants, ONLY : amu_ry
  USE io_files,  ONLY: prefix,seqopn,tmp_dir
  !
  IMPLICIT NONE
  REAL(DP), INTENT(INOUT) :: force(3,nat)
  LOGICAL, INTENT(INOUT) :: moves_finished
  INTEGER, INTENT(IN) :: nat
  ! variables read from the FIRE minimization algorithm
  INTEGER :: nsteppos
  REAL(DP) :: P, dt_init, dt_curr, alpha, alpha_init
  REAL(DP), EXTERNAL :: ddot,dnrm2
  REAL(DP) :: vel(3,nat), acc(3,nat)
  REAL(DP), DIMENSION(3,nat) :: vel_old(3,nat), acc_old(3,nat)
  LOGICAL :: file_exists
  nsteppos = 1
  !
  ! Open the file written by fire
  !
  CALL seqopn( 4, 'fire', 'FORMATTED', file_exists )
  !
  write (*,*) "ARTn Called hess move, resetting FIRE"
  ! IF ( file_exists ) THEN
     !
     ! Everything is OK! read the file
     !
     ! READ( UNIT = 4, FMT = * ) etot, istep, nsteppos, dt_init, dt_curr, alpha, alpha_init, &
          ! tau(:,:), vel_old(:,:), acc_old(:,:)

     ! CLOSE( UNIT = 4, STATUS = 'KEEP' )
     !
     ! set the velocity and acceleration and alpha of previous step to move correctly
     !
     ! for the first step forget previous acc and velocity (prevent P < 0)
     acc_old(:,:) = 0.D0
     vel_old(:,:) = 0.D0
     alpha = 0.D0
     dt_curr = 0.D0
     !IF ( moves_finished ) THEN
     !   ! set dt to 0 so as not to move when loop restarts
     !   dt_curr = 0.0
     !   moves_finished = .false.
     !ELSE
     !   dt_curr = dt_init
     !   force(:,:) = force(:,:)*alat*amu_ry/dt_curr**2
     !ENDIF
     ! overwrite the FIRE file
     ! CALL seqopn( 4, 'fire', 'FORMATTED', file_exists )
     WRITE( UNIT = 4, FMT = * ) etot, istep, nsteppos, dt_init, dt_curr, alpha, alpha_init, &
          tau(:,:), vel_old(:,:), acc_old(:,:)
     CLOSE( UNIT = 4, STATUS = 'KEEP' )
  ! ELSE
  !    CLOSE( UNIT = 4, STATUS = 'DELETE' )
  !    dt_curr = 20.D0
  !    force(:,:) = force(:,:)*alat*amu_ry/dt_curr**2
  ! END IF

END SUBROUTINE hessmove


SUBROUTINE calc_hessian(force,ihess,jhess,hmat)
  !
  !  a subroutine that should calculate the hessian at the current position
  USE kinds, ONLY:DP
  USE cell_base, ONLY : alat
  USE ions_base,     ONLY : nat, tau, if_pos
  USE control_flags, ONLY : istep
  USE ener, ONLY: etot
  USE io_files,  ONLY: prefix,seqopn,tmp_dir
  IMPLICIT none
  REAL(DP), INTENT(INOUT) :: force(3,nat)
  INTEGER, INTENT(INOUT) :: ihess, jhess
  INTEGER, PARAMETER :: iunhess = 533
  REAL(DP), INTENT(OUT) :: hmat(3*nat,3*nat)
  REAL(DP) :: force_in(3,nat), ftmp(3*nat), tau_tmp(3*nat), tau_init(3,nat)
  REAL(DP) :: dr, ddE, dF, etot_init
  LOGICAL :: posmove,negmove,neutmovei, neutmovej, moves_finished, file_exists
  dr = 0.1D0/alat
  ! store the input gradient at current position
  hmat(:,:) = 0.D0
  ! the force should be set to 0 so that we move only the specific coordinate we are looking at ...
  tau_tmp(:) = 0.D0
  ! each element of the hessian matrix should be calculated as H_{i,j} = (grad(E)_dRij+ - grad(E)_dRij-)/dRij
  posmove = .true.
  negmove = .false.
  moves_finished = .false.
  ! write (*,*) "Hessian init:", posmove, file_exists

  dF = 0.D0
  CALL seqopn( iunhess, 'artnhess', 'FORMATTED', file_exists )
  IF ( file_exists ) THEN
     !
     ! read hessian data from the previous iteration  ...
     !
     READ( UNIT = iunhess, FMT = * ) hmat(:,:), posmove,negmove, moves_finished, tau_init(:,:), &
          etot_init, dF
     !
     CLOSE(UNIT = iunhess, STATUS = 'KEEP')
  ELSE
     CLOSE(UNIT = iunhess, STATUS = 'DELETE')
  END IF
  IF ( ihess == 1 .and. jhess == 1 .and. posmove ) THEN
     ! store the initial position and etot
     tau_init(:,:) = tau(:,:)
     etot_init = etot
  ENDIF
  !
  ! due to the way plugin_force_ext is called we need to update the force after the step is made ( in the second call to calc_hess; )
  !
  tau_tmp(1:nat) = tau_init(1,:)
  tau_tmp(nat+1:2*nat) = tau_init(2,:)
  tau_tmp(2*nat+1:3*nat) = tau_init(3,:)
  ftmp(1:nat) = force(1,:)
  ftmp(nat+1:2*nat) = force(2,:)
  ftmp(2*nat+1:3*nat) = force(3,:)
  write (*,*) "Hess parms:", posmove,negmove, etot, etot_init
  write (*,*) "Hess step:", ihess,jhess,istep
  write (*,*) "Hess force in:", force(:,:)
  IF ( ihess <= 3*nat  ) THEN
     IF (jhess <= 3*nat ) THEN
        IF ( posmove ) THEN
           ! increment both indices by step dr
           tau_tmp(jhess) = tau_tmp(jhess) + 0.5*dr
           ! write (*,*) "tau shifts: (j)", tau_tmp(jhess)
           posmove = .false.
           negmove = .true.
        ELSEIF (negmove) THEN
           write (*,*) "Hess added force:", ftmp(ihess)
           dF = dF - ftmp(ihess)
           ! increment by -dr
           ! NOTE: because we incremented by +dr in previous step do -2*dr
           tau_tmp(jhess) = tau_tmp(jhess) - 0.5*dr
           negmove = .false.
           moves_finished = .true.
           ! go to next j-atom
        ELSEIF (moves_finished) THEN
           write (*,*) "Hess subtracted force:", ftmp(ihess)
           ! write (*,*) "Subtracted etot:",etot
           dF = dF + ftmp(ihess)
           write (*,*) "Hess dF:", dF
           ! move to next step, reset tau
           hmat(ihess,jhess) = dF/(dr*alat)
           dF = 0.D0
           write (*,*) "Hess wrote:", hmat(ihess,jhess), "at pos:",ihess,jhess
           jhess = jhess + 1
           ! tau(:,:) = tau_init(:,:)
           moves_finished = .false.
           posmove = .true.
        ENDIF
        tau(1,:) = tau_tmp(1:nat)
        tau(2,:) = tau_tmp(nat+1:2*nat)
        tau(3,:) = tau_tmp(2*nat+1:3*nat)
        ! force(:,:) = RESHAPE(ftmp,(/3, nat/))
        ! force(:,:) = 0.01D0
        CALL hessmove(force,nat,moves_finished)
     ELSE
        ihess = ihess + 1
        ! reset jhess when loop is finished
        jhess = 1
     END IF
  END IF
  WRITE (*,*) "Hessian calculation at:",ihess,jhess
  WRITE (*,*) "Current hessian:", hmat(:,:)
  CALL seqopn( iunhess, 'artnhess', 'FORMATTED', file_exists )
  WRITE (UNIT = iunhess, FMT = * ) hmat(:,:), posmove,negmove, moves_finished, tau_init(:,:), &
       etot_init, dF
  CLOSE (UNIT = iunhess, STATUS = 'KEEP')
END SUBROUTINE calc_hessian

SUBROUTINE center ( vec, nat)
  USE kinds, ONLY: DP
  !
  ! takes as input a vector of size (3,nat) and centers it
  !
  IMPLICIT none
  INTEGER, INTENT(IN) :: nat
  REAL(DP), INTENT(INOUT) :: vec(3,nat)
  INTEGER :: na
  REAL(DP) :: delta(3)
  !
  delta(:) = 0.D0
  DO na = 1,nat
     delta(:) = delta(:) + vec(:,na)
  ENDDO
  !
  delta(:) = delta(:)/dble(nat)
  !
  FORALL ( na = 1:nat) vec(:,na) = vec(:,na) - delta(:)


END SUBROUTINE center
SUBROUTINE diag(n, A, eigvals, vec)
  USE kinds,            ONLY : DP
  !! assuming a general square matrix (can be nonsymmetric).
  !! On output A is overwritten by eigenvectors in rows, if vec=0, then
  !! A is just 0.0 on output.
  !!
  !! n       : dimension of matrix A
  !! A       : matrix to be diagonalised, overwritten by eigenvectors on output
  !! eigvals : output vector of eigenvalues, not sorted!
  !! vec     : 0 if don't want to compute eigenvectors, 1 otherwise
  !!
  IMPLICIT NONE
  INTEGER,              intent(in) :: n
  REAL(DP), DIMENSION(n,n), intent(inout) :: A
  REAL(DP), DIMENSION(n),   intent(out) :: eigvals
  INTEGER,              intent(in) :: vec
  REAL(DP), DIMENSION(n) :: eigvals_i !! imaginary part of the eigenvalues
  REAL(DP), DIMENSION(n,n) :: eigvec
  INTEGER :: lda
  INTEGER :: lwork
  REAL(DP) :: Dummy(1000)
  INTEGER :: info
  CHARACTER(len=1) :: getvec
  getvec = 'N'
  if( vec == 1 ) getvec='V'
  lda = n
  eigvals_i(:) = 0.0
  eigvec(:,:) = 0.0
  !! test workspace
  lwork = -1
  !call sgeev('N',getvec, n, A, lda, eigvals, eigvals_i, &
  !     dummy, 1, eigvec, n, dummy, lwork, info)
  call dgeev('N',getvec, n, A, lda, eigvals, eigvals_i, &
       dummy, 1, eigvec, n, dummy, lwork, info)
  !! choose optimal size of workspace (as in example from intel website)
  lwork = min( 1000, nint(dummy(1)) )
  !! compute stuffs
  !call sgeev('N',getvec, n, A, lda, eigvals, eigvals_i, &
  !     dummy, 1, eigvec, n, dummy, lwork, info)
  call dgeev('N',getvec, n, A, lda, eigvals, eigvals_i, &
       dummy, 1, eigvec, n, dummy, lwork, info)
  !! overwrite a on output
  A(:,:) = eigvec(:,:)
END SUBROUTINE diag

SUBROUTINE eigenmove(force,nat)
  !
  ! subroutine that updates the fire acceleration and velocity so that we don't move along the push direction during perpendicular relax
  !
  USE kinds, ONLY : DP
  USE cell_base, ONLY : alat
  USE constants, ONLY : amu_ry
  USE io_files,  ONLY: prefix,seqopn,tmp_dir
  !
  IMPLICIT NONE
  REAL(DP), INTENT(INOUT) :: force(3,nat)
  INTEGER, INTENT(IN) :: nat
  ! variables read from the FIRE minimization algorithm
  INTEGER :: istep, nsteppos
  REAL(DP) :: etot, P, dt_init, dt_curr, alpha, alpha_init
  REAL(DP), EXTERNAL :: ddot,dnrm2
  REAL(DP) :: tau(3,nat), vel(3,nat), acc(3,nat)
  REAL(DP), DIMENSION(3,nat) :: vel_old(3,nat), acc_old(3,nat)
  LOGICAL :: file_exists
  !
  ! Open the file written by fire
  !
  CALL seqopn( 4, 'fire', 'FORMATTED', file_exists )
  !
  write (*,*) "ARTn Called eigen move, resetting FIRE"
  IF ( file_exists ) THEN
     !
     ! Everything is OK! read the file
     !
     READ( UNIT = 4, FMT = * ) etot, istep, nsteppos, dt_init, dt_curr, alpha, alpha_init, &
          tau(:,:), vel_old(:,:), acc_old(:,:)
     !
     ! set the velocity and acceleration and alpha of previous step to move correctly
     !
     ! for the first step forget previous acc and velocity (prevent P < 0)
     acc_old(:,:) = 0.D0
     vel_old(:,:) = 0.D0
     alpha = 0.1_DP
     force(:,:) = force(:,:)*alat*amu_ry/dt_curr**2
     ! overwrite the FIRE file
     CLOSE( UNIT = 4, STATUS = 'KEEP' )
     CALL seqopn( 4, 'fire', 'FORMATTED', file_exists )

     WRITE( UNIT = 4, FMT = * ) etot, istep, nsteppos, dt_init, dt_curr, alpha, alpha_init, &
          tau(:,:), vel_old(:,:), acc_old(:,:)
     CLOSE( UNIT = 4, STATUS = 'KEEP' )
  ELSE
     WRITE (*,*) "ARTn only works with the FIRE minimization algorithm! STOPPING "
     CLOSE( UNIT = 4, STATUS = 'DELETE' )
  END IF
END SUBROUTINE eigenmove


SUBROUTINE perpmove(nat, istepperp, push)
  !
  ! subroutine that updates the fire acceleration and velocity so that we don't move along the push direction during perpendicular relax
  !
  USE kinds, ONLY : DP
  USE cell_base, ONLY : alat
  USE constants, ONLY : amu_ry
  USE io_files,  ONLY: prefix,seqopn,tmp_dir
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nat,istepperp
  REAL(DP),  INTENT(IN) :: push(3,nat)
  ! variables read from the FIRE minimization algorithm
  INTEGER :: istep, nsteppos
  REAL(DP) :: etot, P, dt_init, dt_curr, alpha, alpha_init
  REAL(DP), EXTERNAL :: ddot,dnrm2
  REAL(DP) :: tau(3,nat), vel(3,nat), acc(3,nat)
  REAL(DP), DIMENSION(3,nat) :: vel_old(3,nat), acc_old(3,nat)
  LOGICAL :: file_exists
  !
  ! Open the file written by fire
  !
  CALL seqopn( 4, 'fire', 'FORMATTED', file_exists )
  !
  write (*,*) "ARTn Called perpmove, number of perp steps:", istepperp
  IF ( file_exists ) THEN
     !
     ! Everything is OK! read the file
     !
     READ( UNIT = 4, FMT = * ) etot, istep, nsteppos, dt_init, dt_curr, alpha_init, &
          alpha, tau(:,:), vel_old(:,:), acc_old(:,:)
     !
     ! set the velocity and acceleration and alpha of previous step to move correctly
     !
     IF ( istepperp == 0 ) THEN
        ! for the first step forget previous acc and velocity (prevent P < 0)
        acc_old(:,:) = 0.D0
        vel_old(:,:) = 0.D0
        alpha = alpha_init
        dt_curr = dt_init
     ELSE
        ! subtract the components that are parallel
        acc_old(:,:) = acc_old(:,:) - ddot(3*nat,acc_old, 1, push, 1 )*push(:,:)/ddot(3*nat,push(:,:),1, push(:,:),1)
        vel_old(:,:) = vel_old(:,:) - ddot(3*nat,vel_old, 1, push, 1 )*push(:,:)/ddot(3*nat,push(:,:),1, push(:,:),1)
     ENDIF
     WRITE (*,*) "ARTn: Perp relax out: dot(vel,push):", ddot(3*nat,vel_old, 1, push, 1 )
     WRITE (*,*) "ARTn: Perp relax out: dot(acc,push):", ddot(3*nat,acc_old, 1, push, 1 )
     ! overwrite the FIRE file
     CLOSE( UNIT = 4, STATUS = 'KEEP' )
     CALL seqopn( 4, 'fire', 'FORMATTED', file_exists )

     WRITE( UNIT = 4, FMT = * ) etot, istep, nsteppos, dt_init, dt_curr, alpha, alpha_init, &
          tau(:,:), vel_old(:,:), acc_old(:,:)
     CLOSE( UNIT = 4, STATUS = 'KEEP' )
  ELSE
     WRITE (*,*) "ARTn only works with the FIRE minimization algorithm! STOPPING "
     CLOSE( UNIT = 4, STATUS = 'DELETE' )
  END IF

END SUBROUTINE perpmove

SUBROUTINE lancmove(nat, v1, dlanc, force )
  USE kinds, ONLY : DP
   USE cell_base, ONLY : alat
   USE constants, ONLY : amu_ry
   USE io_files,  ONLY: prefix,seqopn,tmp_dir
   !
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: nat
   REAL(DP), INTENT(IN) :: dlanc
   REAL(DP),  INTENT(IN) :: v1(3,nat)
   REAL(DP),  INTENT(INOUT) :: force(3,nat)
   REAL(DP) :: step(3,nat)
   ! variables read from the FIRE minimization algorithm
   INTEGER :: istep, nsteppos
   REAL(DP) :: etot, P, dt_init, dt_curr, alpha, alpha_init
   REAL(DP), EXTERNAL :: ddot,dnrm2
   REAL(DP) :: tau(3,nat), vel(3,nat), acc(3,nat)
   REAL(DP), DIMENSION(3,nat) :: vel_old(3,nat), acc_old(3,nat)
   LOGICAL :: file_exists
   !
   ! Open the file written by fire
   !

   CALL seqopn( 4, 'fire', 'FORMATTED', file_exists )

   IF ( file_exists ) THEN
      !
      ! Everything is OK! read the file
      !
      READ( UNIT = 4, FMT = * ) etot, istep, nsteppos, dt_init, dt_curr, alpha_init, alpha, &
           tau(:,:), vel_old(:,:), acc_old(:,:)
      !
      ! set the velocity and acceleration and alpha of previous step to move correctly
      !
      vel_old(:,:) = 0.D0
      acc_old(:,:) = 0.D0
      dt_curr = dt_init
      alpha = 0.D0
      nsteppos = 0
      ! the step performed should be like this now translate it into the correct force
      WRITE (*,*) "Expected lanczos step:", v1(:,:)*dlanc
      force(:,:) = v1(:,:)*dlanc*alat*amu_ry/dt_curr**2
      ! overwrite  fire parameters
      CLOSE( UNIT = 4, STATUS = 'KEEP' )
      CALL seqopn( 4, 'fire', 'FORMATTED', file_exists )

      WRITE( UNIT = 4, FMT = * ) etot, istep, nsteppos, dt_init, dt_curr, alpha_init, alpha, &
           tau(:,:), vel_old(:,:), acc_old(:,:)
      CLOSE( UNIT = 4, STATUS = 'KEEP' )
   ELSE
      WRITE (*,*) "ARTn only works with the FIRE minimization algorithm! STOPPING "
      CLOSE( UNIT = 4, STATUS = 'DELETE' )
   END IF

END SUBROUTINE lancmove

SUBROUTINE lanczos( nat, force, vel, acc, alpha_init, dt, &
     v_in, nlanciter, nlanc, lowest_eigval, lowest_eigvec, pushdir)
  USE kinds,            ONLY : DP
  USE io_files, ONLY: prefix,seqopn,tmp_dir
  !
  ! Lanczos subroutine for the ARTn algorithm; based on the lanczos subroutine as written by M. Gunde
  !
  IMPLICIT NONE
  INTEGER,                INTENT(IN) :: nat
  REAL(DP), DIMENSION(3,nat), INTENT(IN) :: v_in
  REAL(DP), DIMENSION(3,nat), INTENT(INOUT) :: vel, acc
  REAL(DP), INTENT(IN) :: alpha_init, dt   
  REAL(DP), DIMENSION(3,nat), INTENT(INOUT) :: force
  REAL(DP), DIMENSION(3,nat), INTENT(INOUT) :: lowest_eigvec
  REAL(DP), INTENT(INOUT) :: lowest_eigval
  REAL(DP), DIMENSION(3,nat), INTENT(IN) :: pushdir
  !
  INTEGER :: i, j, io, id_min
  INTEGER, PARAMETER ::  iunlanc = 51
  INTEGER, INTENT(INOUT) :: nlanciter
  INTEGER, INTENT(INOUT) :: nlanc
  REAL(DP) :: dlanc
  REAL(DP), PARAMETER :: eigvec_thr = 1.0D-3, eigval_thr = 1.0D-2
  REAL(DP), DIMENSION(3,nat) :: force_old, lowest_eigvec_old
  REAL(DP), ALLOCATABLE :: v0(:,:), v1(:,:), v2(:,:), q(:,:), eigvals(:)
  REAL(DP), ALLOCATABLE :: Vmat(:,:,:), Vmat_mul(:,:), H(:,:), Hstep(:,:)
  REAL(DP) :: lowest_eigvec_tmp(3*nat)
  REAL(DP) :: dir
  REAL(DP), EXTERNAL :: ran3,dnrm2,ddot
  REAL(DP) :: alpha, beta, lowest_eigval_old, eigvec_diff, largest_eigvec_diff, eigval_diff
  LOGICAL :: file_exists

  !! allocate vectors
  ALLOCATE( q(3,nat) )
  ALLOCATE( v0(3,nat))
  ALLOCATE( v1(3,nat))
  ALLOCATE( v2(3,nat))

  !! allocate matrices
  ALLOCATE( Vmat_mul(3*nat,nlanciter), source=0.D0)
  ALLOCATE( Vmat(3,nat,1:nlanciter), source=0.D0 )
  ALLOCATE( H(1:nlanciter,1:nlanciter), source=0.D0 )
  ALLOCATE( Hstep(1:nlanciter,1:nlanciter), source=0.D0 )
  CALL seqopn( iunlanc, 'artnlanc', 'FORMATTED', file_exists )
  !
  ! initialize lanczos counter and variables
  !
  ! store the eignvalue and eigenvec of previous iteration
  write (*,*) "Read lowest Eigenvalue:",lowest_eigval
  lowest_eigvec_old(:,:) = lowest_eigvec(:,:)
  lowest_eigval_old = lowest_eigval
  write (*,*) "Stored Eigenvalue:",lowest_eigval_old
  ! parameter for lanczos moves
  dlanc = 1.0D-2
  IF ( file_exists ) THEN
     !
     ! read lanczos data from the previous iteration  ...
     !
     READ( UNIT = iunlanc, FMT = * )  Vmat(:,:,1:nlanciter), H(1:nlanciter,1:nlanciter),  &
         force_old(:,:)
     !
     IF (nlanc > nlanciter) THEN
        ! if we reached the final lanczos step then delete the lanczos file
        CLOSE( UNIT = iunlanc, STATUS = 'DELETE' )
     ELSE
        CLOSE( UNIT = iunlanc, STATUS = 'KEEP' )
     ENDIF
  ELSE
     CLOSE(UNIT = iunlanc, STATUS = 'DELETE')
  END IF
  !
  ! in the first lanczos step we should give a random push to initiate the algorithm and then calc the force of the new pos with qe.
  !
  IF ( nlanc  == 0 ) THEN
     !!
     !! normalize initial vector
     !!
     ! NOTE: the inital vector to lanczos is random
     v0(:,:) = v_in(:,:)
     v0(:,:) = v0(:,:) / dnrm2( 3*nat, v0, 1 )
     ! store v0
     Vmat(:,:,1) = v0(:,:)
     write (*,*) "ARTn Lanczos: initial vec:", v0(:,:)
     ! write lanczos data to file for future cycles
     nlanc = nlanc + 1
     CALL seqopn( iunlanc, 'artnlanc', 'FORMATTED', file_exists )

     WRITE (UNIT = iunlanc, FMT = * ) Vmat(:,:,1:nlanciter), H(1:nlanciter,1:nlanciter), &
          force(:,:)
     CLOSE (UNIT = iunlanc, STATUS = 'KEEP')

     ! CALL lancmove(nat, v0, dlanc, force )
     CALL move_mode( nat, dlanc, v0, force, &
          vel, acc, alpha_init, dt, & 
          0, pushdir, 'lanc')
     ! GO BACK (make move, get new force)
     RETURN
  ELSEIF (nlanc == 1 ) THEN
     !! q(:,:) now represents {Force} = -[Hessian]{R_1}, where {R_1} = {R} + {dR} = {R} + {v0}*d_lanc,
     !! and {v0} is the Lanczos vector v0(:). We need q(:) to represent [Hessian]{v0}, thus
     !! first do q(:) = -[Hessian]( {R_1} - {R} ) = -[Hessian]{dR}
     q(:,:) = force(:,:) - force_old(:,:)
     !! now do q(:) = [Hessian]{dR}/d_lanc = [Hessian]{v0}
     q(:,:) = -q(:,:) / dlanc
     alpha = ddot(3*nat,Vmat(:,:,1),1,q(:,:),1)
     ! check if the eigenvalue
     !IF ( lowest_eigval_old /= 0.D0) THEN
     !   IF ( (alpha - lowest_eigval_old)/lowest_eigval_old .lt. eigval_thr ) THEN
     !      write (*,*) "ARTn Lanczos: new eigenvalue close to the old eigenvalue, STOPPING"
     !      ! set the lanczos counter to current number of lanczos iterations
     !      nlanciter = nlanc
     !   ENDIF
     !ENDIF
     v1(:,:) = q(:,:) - alpha*Vmat(:,:,1)
     beta = dnrm2(3*nat,v1,1)
     ! store the vecs for future cycles ...
     Vmat(:,:,2) = v1(:,:)/beta
     write (*,*) "ARTn Lanczos: second vec:", v1(:,:)/beta
     H(1,1) = alpha
     H(2,1) = beta
     H(1,2) = beta

     write(555,'(i4,2f9.4)') nlanc, alpha, (alpha-lowest_eigval_old)/lowest_eigval_old
     flush(555)
     nlanc = nlanc + 1
     CALL seqopn( iunlanc, 'artnlanc', 'FORMATTED', file_exists )
     WRITE (UNIT = iunlanc, FMT = * ) Vmat(:,:,1:nlanciter), H(1:nlanciter,1:nlanciter), &
          force(:,:)
     CLOSE (UNIT = iunlanc, STATUS = 'KEEP')
     ! CALL lancmove(nat, v1/beta, dlanc, force )
     CALL move_mode( nat, dlanc, v1/beta, force, &
          vel, acc, alpha_init, dt,  & 
          0, pushdir, 'lanc')
     ! CALL move_mode( nat, dlanc, v1/beta, force, 0, pushdir, 'lanc')
     ! GO BACK  (make move, get new force)
     lowest_eigval = alpha
     RETURN
  ELSEIF (nlanc > 1 .and. nlanc <= nlanciter ) THEN
     !
     v0(:,:) = Vmat(:,:,nlanc-1)
     v1(:,:) = Vmat(:,:,nlanc)
     !
     q(:,:) = force(:,:) - force_old(:,:)
     q(:,:) = -q(:,:)/dlanc
     !
     alpha = ddot(3*nat,v1(:,:),1,q(:,:),1)
     H(nlanc,nlanc) = alpha
     !
     ! Do a diagonalization here ; check if eigenvalues are converged
     !
     ALLOCATE( eigvals(nlanc) )
     ! store the H matrix, because its overwritten by eigvecs on diagonalization
     Hstep(:,:) = H(:,:)

     CALL diag(nlanc, H(1:nlanc,1:nlanc), eigvals, 1 )

     lowest_eigval = eigvals(1)
     id_min = 1
     !
     ! get the lowest eigenvalue
     !
     DO i = 1, nlanc
        IF (eigvals(i) < lowest_eigval ) THEN
           lowest_eigval = eigvals(i)
           id_min = i
        ENDIF
     ENDDO
     !
     ! reshape the V mat for now (FIX this so it uses the DGEMM from lapack)
     !
     Vmat_mul(1:nat,:) = Vmat(1,:,:)
     Vmat_mul(nat+1:2*nat,:) = Vmat(2,:,:)
     Vmat_mul(2*nat+1:3*nat,:) = Vmat(3,:,:)

     lowest_eigvec_tmp(:) = matmul(Vmat_mul(:,:),H(:,id_min))
     lowest_eigvec(1,:) = lowest_eigvec_tmp(1:nat)
     lowest_eigvec(2,:) = lowest_eigvec_tmp(nat+1:2*nat)
     lowest_eigvec(3,:) = lowest_eigvec_tmp(2*nat+1:3*nat)
     !
     ! check if the obtained eigenvec points in the same direction as the input eigvec
     !
     dir = ddot(3*nat,lowest_eigvec,1, pushdir, 1)
     IF ( dir < 0.D0 ) THEN
        WRITE (*,*) "Dir is less than 0, flipping eigvec:", lowest_eigvec(:,:)
        lowest_eigvec(:,:) = -1.D0*lowest_eigvec(:,:)
     ENDIF
     WRITE (*,*) "ARTn Lanczos: Eigenvalues:", eigvals(:), "at step:", nlanc
     !
     ! calculate the largest difference between current eigenvec vs previous
     !
     !largest_eigvec_diff = 0.D0
     !DO i=1,nat
     !   DO j=1,3
     !      eigvec_diff = ABS(lowest_eigvec(j,i)) - ABS(lowest_eigvec_old(j,i))
     !      IF ( ABS(eigvec_diff) > largest_eigvec_diff ) THEN
     !         largest_eigvec_diff = eigvec_diff
     !      ENDIF
     !   ENDDO
     !ENDDO
     !
     ! Check for the convergence of the lanczos eigenvalue
     !
     eigval_diff = (lowest_eigval - lowest_eigval_old)/lowest_eigval_old
     write(555,'(i4,2f9.4)') nlanc, lowest_eigval, eigval_diff
     flush(555)
     write (*,*) "Eigenvalue difference:", eigval_diff
     IF ( ABS(eigval_diff) <= eigval_thr ) THEN
        !
        WRITE (*,*) "ARTn Lanczos: eigenvalue converged at step:", nlanc, "eigenthr:", eigval_thr
        !
        nlanciter = nlanc
        CALL seqopn( iunlanc, 'artnlanc', 'FORMATTED', file_exists )
        WRITE (UNIT = iunlanc, FMT = * ) Vmat(:,:,1:nlanciter), Hstep(1:nlanciter,1:nlanciter), &
          force(:,:)
        CLOSE (UNIT = iunlanc, STATUS = 'KEEP')
        ! write(555,*) 'eigvec:'
        ! write(555,'(3f13.6)') lowest_eigvec
        ! flush(555)
     ENDIF
     !
     ! if lanczos is not yet converged generate new matrix elements
     !
     IF ( nlanc < nlanciter ) THEN
        beta = Hstep(nlanc,nlanc-1)
        !
        v2(:,:) = q(:,:) - alpha*v1(:,:) - beta*v0(:,:)
        !
        ! orthogonalize vectors in accordance with previous ones ...
        DO j = 1, nlanc - 1
           v2(:,:) = v2(:,:) - ddot(3*nat, v2 ,1, Vmat(:,:,j),1)*Vmat(:,:,j)
        ENDDO
        !
        !  stop the lanczos algorithm when v2 is very small ...
        !
        IF ( dnrm2(3*nat, v2(:,:), 1) < 1.0D-15 ) THEN
           !
           nlanciter = nlanc
           write (*,*) "ARTn Lanczos: Stopped because Lanczos vector very small!", dnrm2(3*nat, v2(:,:), 1)
           !
           ! Backtrack to initial position
           !
           nlanc = nlanc + 1
           v1(:,:) = 0.D0
           DO i=1,nlanciter
              v1(:,:) = v1(:,:) - Vmat(:,:,i)
           ENDDO
           ! make the move back & delete the lanczos file
           ! CALL lancmove(nat, v1(:,:), dlanc, force(:,:) )
           ! CALL move_mode( nat, dlanc, v1, force, 0, pushdir, 'lanc')
           CALL move_mode( nat, dlanc, v1, force, &
                vel, acc, alpha_init, dt,  & 
                0, pushdir, 'lanc')
           CALL seqopn( iunlanc, 'artnlanc', 'FORMATTED', file_exists )
           CLOSE (UNIT = iunlanc, STATUS = 'DELETE')
           RETURN
        ELSE
           !
           v2(:,:) = v2(:,:)/dnrm2(3*nat,v2(:,:),1)
           write (*,*) "ARTn Lanczos: next vec:", v2(:,:)
           Vmat(:,:,nlanc+1) = v2(:,:)
           ! store values to H
           beta = ddot(3*nat, q, 1, v2, 1)
           Hstep(nlanc+1,nlanc ) = beta
           Hstep(nlanc, nlanc+1) = beta
           nlanc = nlanc +1

           CALL seqopn( iunlanc, 'artnlanc', 'FORMATTED', file_exists )

           WRITE (UNIT = iunlanc, FMT = * ) Vmat(:,:,1:nlanciter), Hstep(1:nlanciter,1:nlanciter), &
                force(:,:)
           CLOSE (UNIT = iunlanc, STATUS = 'KEEP')

           DEALLOCATE( eigvals )
           !
           ! for the next move v2 will be v1 ...; call move with v2
           !
           ! CALL lancmove(nat, v2(:,:), dlanc, force(:,:) )
           CALL move_mode( nat, dlanc, v2, force, &
                vel, acc, alpha_init, dt,  & 
                0, pushdir, 'lanc')
           ! CALL move_mode( nat, dlanc, v2, force, 0, pushdir, 'lanc')
           !  GO BACK (make move, get new force)
           RETURN
        ENDIF
     ELSE
        !
        ! backtrack to the initial position
        !
        nlanc = nlanc + 1
        v1(:,:) = 0.D0
        DO i=1,nlanciter
           v1(:,:) = v1(:,:) - Vmat(:,:,i)
        ENDDO
        ! make the final move & delete lanczos file
        ! CALL lancmove(nat, v1(:,:), dlanc, force(:,:) )
        CALL move_mode( nat, dlanc, v1, force, &
             vel, acc, alpha_init, dt, & 
             0, pushdir, 'lanc')
        ! CALL move_mode( nat, dlanc, v1, force, 0, pushdir, 'lanc')
        CALL seqopn( iunlanc, 'artnlanc', 'FORMATTED', file_exists )
        CLOSE (UNIT = iunlanc, STATUS = 'DELETE')
        RETURN
     END IF
  ENDIF
  ! deallocate vectors
  DEALLOCATE( q, v0, v1, v2 )
  ! deallocate matrices
  DEALLOCATE(Vmat, Vmat_mul, H, Hstep)

END SUBROUTINE lanczos

SUBROUTINE displacement_validation( atom_id, atom_const, push, lvalid)
  !
  ! subroutine that checks if the initial_displacement is within given parameters
  !
  USE constants, ONLY : PI
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: atom_id
  REAL(DP), INTENT(IN) :: atom_const(4)
  REAL(DP), INTENT(INOUT) :: push(3)
  REAL(DP), EXTERNAL :: ddot, dnrm2
  LOGICAL,         INTENT(INOUT) :: lvalid
  !
  ! Local variables
  REAL(DP)               :: cone_angle, displacement_angle
  REAL(DP)               :: dot_prod, displacement_norm, cone_dir_norm
  REAL(DP), DIMENSION(3) :: cone_dir,displacement
  !
  !
  !
  write (*,*) "ARTn: Called displacement validation: with cone_dir:",atom_const(1:3), "current push:", push(:)
  cone_dir = atom_const(1:3)
  cone_angle = atom_const(4)
  !
  displacement(:)         = push(:)
  !
  displacement_norm       = dnrm2(3,displacement,1)
  cone_dir_norm           = dnrm2(3,cone_dir,1)
  !
  dot_prod                = ddot( 3, cone_dir, 1, displacement, 1 ) / ( cone_dir_norm * displacement_norm )
  displacement_angle      = ACOS( dot_prod ) *180.0_DP / PI
  lvalid                  = ( displacement_angle < cone_angle )
  write (*,*) "Finished displacement validation",lvalid
  !
  IF ( cone_angle == 0.0_DP) THEN
     lvalid = .TRUE.
     !
     ! TODO: why is the direction multiplied by 0.1? seems kind of random ...
     !
     push(:) = cone_dir(:)
  ENDIF
  !
END SUBROUTINE displacement_validation


SUBROUTINE push_init (nat, idum, init_step_size, push)
  !
  ! subroutine for generating a random push on all atoms (global_move scenario)
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT none
  INTEGER, INTENT(IN) :: nat, idum
  ! INTEGER, INTENT(IN) :: if_pos(3,nat)
  REAL(DP), INTENT(IN) :: init_step_size
  INTEGER :: na
  REAL(DP), INTENT(OUT) :: push(3,nat)
  REAL(DP), EXTERNAL :: ran3
  REAL(DP) :: dr2
  !
  ! generate a random push
  !
  DO na=1,nat
     DO
        push(1,na) = (0.5_DP - ran3(idum))
        push(2,na) = (0.5_DP - ran3(idum))
        push(3,na) = (0.5_DP - ran3(idum))
        dr2 = push(1,na)**2 + push(2,na)**2 + push(3,na)**2
        IF ( dr2 < 0.25_DP ) EXIT
     ENDDO
  ENDDO

  !
  ! scale the push according to initial step size
  !
  push(:,:) = init_step_size*push(:,:)

END SUBROUTINE push_init

SUBROUTINE push_init_list (nat, natpush, idum, push_ids, add_const, init_step_size, push)
  !
  ! subroutine that generates a random push to a list of atoms specified by user;
  !
  !
  ! the user should supply: number and list of atoms to push; and add_constraints on these atoms
  !
  USE kinds, ONLY : DP
  IMPLICIT none
  INTEGER, INTENT(IN) :: nat,natpush, idum
  INTEGER :: na
  INTEGER, INTENT(IN) :: push_ids(natpush)
  REAL(DP), INTENT(IN) :: init_step_size
  REAL(DP), INTENT(IN) :: add_const(4,nat)
  REAL(DP), INTENT(OUT) :: push(3,nat)
  REAL(DP), EXTERNAL :: ran3
  REAL(DP) :: dr2, pushat(3)
  LOGICAL :: lvalid
  INTEGER :: atom_displaced(nat)
  !
  !  read the list of pushed atoms
  !
  push(:,:) = 0.d0
  lvalid = .false.
  DO na=1,nat
     IF (ANY(push_ids == na)) THEN
        atom_displaced(na) = 1
     ELSE
        atom_displaced(na) = 0
     ENDIF
  ENDDO
  !
  !
  DO na=1,nat
     IF (atom_displaced(na) == 1 ) THEN
        DO
           push(:,na) = (/0.5_DP - ran3(idum),0.5_DP - ran3(idum),0.5_DP - ran3(idum)/)
           dr2 = push(1,na)**2 + push(2,na)**2 + push(3,na)**2
           ! check if the atom is constrained
           IF (ANY(ABS(add_const(:,na)) > 0.D0)) THEN
              ! write (*,*) "Atom",na, "Constrained to:", add_const(:,na)
              ! check if the displacement is within the chosen constraint
              CALL displacement_validation(na,add_const(:,na),push(:,na),lvalid )
              ! write (*,*) push(:)
              IF ( .not. lvalid ) CYCLE
           ENDIF
           IF ( dr2 < 0.25_DP ) EXIT
        ENDDO
     ENDIF
  ENDDO
  !
  ! scale the push according to initial step size
  !
  push(:,:) = init_step_size*push(:,:)

END SUBROUTINE push_init_list

SUBROUTINE report_force(force,eigenvec,nat,force_tot,fperp_tot,fpara_tot)
  !
  ! a subroutine that reports the non-manipulated forces of QE (Ftot, Fpara, and Fperp)
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nat
  REAL(DP), INTENT(IN) :: force(3,nat), eigenvec(3,nat)
  REAL(DP) :: fpara(3,nat), fperp(3,nat)
  REAL(DP), INTENT(OUT) :: force_tot, fperp_tot, fpara_tot
  REAL(DP), EXTERNAL :: ddot
  !
  CALL sum_force(force,nat,force_tot)
  fperp(:,:) = force(:,:)
  CALL perpforce(fperp,eigenvec,fpara,nat)
  CALL sum_force(fperp,nat,fperp_tot)
  fpara_tot = ddot(3*nat,force,1,eigenvec,1)
  ! CALL sum_force(fpara,nat,fpara_tot)
  WRITE (*,'("ARTn Forces (Ry/a.u.):",5X, "Ftot = ",F6.3,5X,"Fperp = ",F6.3, 5X,"Fpara = ",F6.3)') &
       & force_tot,fperp_tot,fpara_tot
END SUBROUTINE report_force

SUBROUTINE sum_force(force,nat,force_tot)
  !
  ! subroutine that sums the forces on all atoms and returns the total force
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nat
  REAL(DP), INTENT(IN) :: force(3,nat)
  REAL(DP), INTENT(OUT) :: force_tot
  INTEGER :: na
  force_tot = 0.0
  DO na = 1, nat
     force_tot = force_tot + force(1,na)**2 + force(2,na)**2 + force(3,na)**2
  ENDDO
  force_tot = SQRT(force_tot)
END SUBROUTINE sum_force

SUBROUTINE perpforce(force,push,fpara,nat)
  !
  ! subroutine that subtracts parallel components to push from force
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  REAL(DP), INTENT(IN)  :: push(3,nat)
  REAL(DP), INTENT(INOUT)  :: force(3,nat)
  REAL(DP), INTENT(OUT)  :: fpara(3,nat)
  REAL(DP) :: push_norm(3,nat)
  REAL(DP), EXTERNAL :: ddot,dnrm2
  INTEGER, INTENT(IN) :: nat
  INTEGER  :: na
  ! calculate components parallel to the push
  fpara(:,:) = ddot(3*nat,force(:,:),1,push(:,:),1) / ddot(3*nat,push(:,:),1,push(:,:),1) * push(:,:)
  ! subtract them
  force(:,:) = force(:,:) - fpara(:,:)

END SUBROUTINE perpforce
SUBROUTINE move_mode(nat, dlanc, v1, force, &
                     vel, acc, alpha_init, dt, & 
                     istepperp, push, &
                     mode)
  !
  ! unify the move routines
  !
  USE kinds, ONLY: DP
  USE cell_base, ONLY: alat
  USE constants, ONLY: amu_ry
  USE io_files, ONLY: prefix, seqopn, tmp_dir
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nat
  REAL(DP), INTENT(IN) :: dlanc
  REAL(DP), DIMENSION(3,nat), INTENT(IN) :: v1
  REAL(DP), DIMENSION(3,nat), INTENT(INOUT) :: force
  REAL(DP), DIMENSION(3,nat), INTENT(INOUT) :: vel, acc
  REAL(DP), INTENT(IN) :: alpha_init, dt  
  INTEGER, INTENT(IN) :: istepperp
  REAL(DP), DIMENSION(3,nat), INTENT(IN) :: push
  CHARACTER(LEN=4), INTENT(IN) :: mode
  REAL(DP), EXTERNAL :: ddot,dnrm2
  ! variables read from the FIRE minimization algorithm
  INTEGER :: nsteppos
  REAL(DP) :: dt_curr, alpha, etot
  LOGICAL :: file_exists
  
  !
  ! write(*,*) 'in MOVE:'
  ! write(*,*) 'mode:',trim(mode)
  ! write(*,*) 'nat:',nat
  ! write(*,*) 'dlanc',dlanc
  ! write(*,*) 'v1',v1
  ! write(*,*) 'force', force
  ! write(*,*) 'istepperp',istepperp
  ! write(*,*) 'push',push
  ! write(*,*) 'file exists', file_exists
  !
     
  !
  ! do things depending on mode of the move
  !
  CALL seqopn( 4, 'fire', 'FORMATTED', file_exists )
  IF (file_exists ) THEN
     READ( UNIT = 4, FMT = * ) etot, nsteppos, dt_curr, alpha
     CLOSE( UNIT = 4, STATUS = 'KEEP' )
     SELECT CASE( TRIM(mode) )
        
     CASE( 'perp' )
        ! write(*,*) 'MOVE: perp'
        !
        IF( istepperp .eq. 0 ) THEN
           ! for the first step forget previous acc and velocity (prevent P < 0)
           acc(:,:) = 0.D0
           vel(:,:) = 0.D0
           alpha = alpha_init
           dt_curr = dt
        ELSE
           ! subtract the components that are parallel
           acc(:,:) = acc(:,:) - ddot(3*nat,acc, 1, push, 1 )*push(:,:)/ddot(3*nat,push(:,:),1, push(:,:),1)
           vel(:,:) = vel(:,:) - ddot(3*nat,vel, 1, push, 1 )*push(:,:)/ddot(3*nat,push(:,:),1, push(:,:),1)
        ENDIF
        !
     CASE( 'lanc' )
        ! write(*,*) 'MOVE: lanc'
        !
        ! set the velocity and acceleration and alpha of previous step to move correctly
        !
        vel(:,:) = 0.D0
        acc(:,:) = 0.D0
        dt_curr = dt
        alpha = 0.D0
        nsteppos = 0
        ! the step performed should be like this now translate it into the correct force
        force(:,:) = v1(:,:)*dlanc*alat*amu_ry/dt_curr**2
        ! WRITE (*,*) "Expected lanczos step:", v1(:,:)*dlanc
        !
     CASE( 'eign' )
        ! write(*,*) 'MOVE: eign'
        !
        ! for the first step forget previous acc and velocity (prevent P < 0)
        acc(:,:) = 0.D0
        vel(:,:) = 0.D0
        alpha = 0.0_DP
        force(:,:) = force(:,:)*alat*amu_ry/dt_curr**2
        !
     CASE( 'hess' )
        ! write(*,*) 'MOVE: hess'
        !
        acc(:,:) = 0.D0
        vel(:,:) = 0.D0
        alpha = 0.D0
        dt_curr = 0.D0
     !
     CASE default
        write(*,*) 'Problem with move_mode!'
     END SELECT
     !
     !
     ! write the FIRE parameters to its scratch file 
     CALL seqopn( 4, 'fire', 'FORMATTED', file_exists )
     WRITE( UNIT = 4, FMT = * ) etot, nsteppos, dt_curr, alpha  
     !
     CLOSE( UNIT = 4, STATUS = 'KEEP' )
     !
     !
  ELSE
     CLOSE( UNIT = 4, STATUS = 'DELETE') 
  ENDIF

END SUBROUTINE move_mode




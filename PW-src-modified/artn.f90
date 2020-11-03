!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE artn(force,nat,alat,istep,if_pos,vel,acc,dt,fire_alpha_init,prefix,tmp_dir)
  !----------------------------------------------------------------------------
  !
  !
  USE kinds,            ONLY : DP
  IMPLICIT NONE
  !
  ! local varibles
  !
  REAL(DP), INTENT(INOUT)  :: force(3,nat),vel(3,nat),acc(3,nat)
  REAL(DP), INTENT(IN) :: dt, alat, fire_alpha_init
  INTEGER, INTENT(IN) :: nat, istep, if_pos(3,nat)
  CHARACTER(LEN=255), INTENT(IN) :: prefix, tmp_dir 
  INTEGER, PARAMETER :: iunart = 50, iunlanc = 51, iunartin=52
  INTEGER :: na, idum, nlanc, npush, neigenstep,neigenstepmax,npushmin
  INTEGER ::  nlanciter, nlanciter_init, natpush, istepperp, if_pos_ct, i
  REAL(DP) :: dlanc
  ! flags for controlling ARTn
  LOGICAL :: lpush_init, leigen, llanczos, lperp, file_exists, lpush_list
  !
  REAL(DP)  :: force_in(3,nat), push(3,nat), fpara(3,nat), pushdir(3,nat), eigenvec(3,nat), v_in(3,nat)
  REAL(DP) :: fpara_tot, fperp_tot, force_tot
  REAL(DP) :: convcrit,init_step_size, step_size, current_step_size
  REAL(DP) :: add_const(4,nat)
  INTEGER :: push_ids(nat)
  REAL(DP), EXTERNAL :: ran3,dnrm2, ddot
  REAL(DP) :: lowest_eigval, eigval_thr 
  INTEGER :: ios
  CHARACTER( LEN=255) :: filnam
  !
  ! The basic idea is:
  ! (1) push atoms in the direction specified by user & relax in the perpendicular direction;
  ! (2) use the lanczos algorithm calculate the lowest eigenvalue/eigenvec
  ! (3) a negative eigenvalue, update push direction otherwise push again
  ! (4) follow the lanczos direction twoard the saddle point
  
  ! input namelist 
  NAMELIST/artn_parameters/ convcrit,init_step_size,lpush_list,npushmin,push_ids,add_const,dlanc,nlanciter_init,step_size
  !
  ! Set up control flags for the first step
  !
  ! lhess = .false.
  lpush_init = .true.
  lperp = .false.
  llanczos = .false.
  leigen = .false.
  !
  add_const(:,:) = 0.D0
  push_ids(:) = 0
  push(:,:) = 0.D0
  ! 
  istepperp = 0
  !
  nlanc = 0
  nlanciter = 0 
  if_pos_ct = 0
  !
  neigenstep = 0
  neigenstepmax = 1
  ! initialize push vector
  ! initialize vecs used in the lanczos algorithm 
  eigenvec(:,:) = 0.D0
  pushdir(:,:) = 0.D0
  ! set up lowest eigenvalue 
  lowest_eigval  = 0.D0
  ! store original force
  force_in(:,:) = force(:,:)
  !
  fperp_tot = 0.D0
  fpara_tot = 0.D0
  force_tot = 0.D0
  !
  ! input parameters 
  ! 
  ! npushmin = 3 
  ! convcrit = 1.0d-2
  ! init_step_size = 1.5
  ! parameters for the push
  ! lpush_list = .true.
  ! define which atoms are to be pushed and the constraints ...
  ! push_ids = (/1/)
  ! add_const(:,1) = (/1.0, 0.0, 0.0, 0.0/)
  ! lanczos parameters
  ! dlanc = 1.D-2 
  ! nlanciter_init = 16
  ! eigenpush parms
  ! step_size = 0.5
  

  !
  filnam = 'artn.in'
  OPEN( UNIT = iunartin, FILE = filnam, FORM = 'formatted', STATUS = 'unknown', IOSTAT = ios)
  READ( NML = artn_parameters, UNIT = iunartin)
  CLOSE ( UNIT = iunartin, STATUS = 'KEEP')  
  ! set nlanciter to the one from input ;
  nlanciter = nlanciter_init
  !
  ! First read the scratch file and update flags
  !
  filnam = trim(tmp_dir) // '/' // trim(prefix) // '.' // 'artn'
  INQUIRE( file = filnam, exist = file_exists )
  OPEN( unit = iunart, file = filnam, form = 'formatted', status = 'unknown', iostat = ios )
  !
  IF ( file_exists ) THEN
     !
     ! ... the scratch file is read
     !
     READ( UNIT =  iunart, FMT = *)  lpush_init, lperp, npush, nlanc, nlanciter, &
           neigenstep, istepperp, llanczos, leigen,  &
           push(:,:), eigenvec(:,:), lowest_eigval, &
           current_step_size, convcrit, pushdir(:,:)
     !
     CLOSE( UNIT =  iunart, STATUS = 'KEEP' )
  ELSE
     CLOSE( UNIT = iunart, STATUS = 'DELETE')
  END IF
  !
  ! start initial push
  !
  IF ( lpush_init ) THEN
     !
     ! check type of initial push
     !
     IF ( lpush_list ) THEN
        CALL push_init_list(nat, idum, push_ids, add_const, init_step_size, push )
     ELSE
        CALL push_init(nat,idum, init_step_size, push)
     ENDIF
     ! set up the flags (we did an initial push, now we need to relax perpendiculary
     lpush_init = .false.
     lperp = .true.
     istepperp = 0
     ! the push counter (controls if we should call lanczos or keep pushing)
     npush = 1
     ! set the force equal to the push 
     force(:,:) =  push(:,:)
  ELSE IF ( lperp ) THEN
     !
     !  subtract parrallel components to push from force
     !
     IF ( leigen ) THEN
        push(:,:) = eigenvec(:,:)
     ENDIF
     CALL perpforce(force,push,fpara,nat)
     CALL move_mode( nat, alat, dlanc, v_in, force, &
          vel, acc, fire_alpha_init, dt, &
          istepperp, push, 'perp', prefix, tmp_dir)
     istepperp = istepperp + 1
     !
     !
     IF (( MAXVAL( ABS(force) )) < convcrit ) THEN
        !
        ! we reached convergence in the perpendicular direction of push
        !
        lperp = .false.
        istepperp = 0
        IF ( npush < npushmin ) THEN
           ! continue pushing in the specified direction
           force(:,:) =  push(:,:)
           npush = npush + 1
           lperp = .true.
        ELSE IF ( npush >= npushmin  ) THEN
           ! start lanczos algorithm
           force(:,:) = force_in(:,:) 
           llanczos = .true.
        END IF
     END IF
  ELSE IF ( leigen .and. .not. lperp ) THEN
     !
     ! if we have a good lanczos eigenvector use it
     !
     fpara_tot = ddot(3*nat,force(:,:),1,eigenvec(:,:),1)
     ! check magnitude of parallel force 
     IF (ABS(fpara_tot) <= 0.005_DP) THEN
        ! decrease perpendicular convergance criterion
        convcrit = 0.1_DP*convcrit
     END IF
     ! rescale the eigenvector according to the current force in the parallel direction
     current_step_size = MIN(step_size,ABS(fpara_tot)/MAX(ABS(lowest_eigval),0.5_DP))
     eigenvec(:,:) = -SIGN(1.0D0,fpara_tot)*eigenvec(:,:)*current_step_size
     ! 
     ! 
     force(:,:) = eigenvec(:,:)
     CALL move_mode( nat, alat, dlanc, v_in, force, &
          vel, acc, fire_alpha_init, dt,  &
          istepperp, push, 'eign', prefix, tmp_dir)
     ! update eigenstep counter 
     neigenstep = neigenstep + 1
     ! count the number of steps made with the eigenvector
     IF ( neigenstep == neigenstepmax  ) THEN
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
     IF (nlanc == 0 ) THEN
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

     ! apply constraints from QE 
     IF ( ANY(if_pos(:,:) == 0) ) THEN
        DO na=1,nat
           DO i=1,3
              IF (if_pos(i,na) == 1 ) if_pos_ct = if_pos_ct + 1
           ENDDO
        END DO
        IF ( if_pos_ct < nlanciter .and. if_pos_ct /= 0 ) nlanciter = if_pos_ct
        v_in(:,:) = v_in(:,:)*if_pos(:,:)
        force(:,:) = force(:,:)*if_pos(:,:)
     ENDIF
     CALL lanczos( nat, alat, force, vel, acc, fire_alpha_init, dt, &
          v_in, dlanc, nlanciter, nlanc, lowest_eigval,  eigenvec, pushdir, prefix, tmp_dir )

     !
     ! when lanczos converges, nlanciter = number of steps it took to converge,
     ! and nlanc = nlanciter + 1
     !
     IF ( nlanc > nlanciter ) THEN
        !
        ! max number of lanczos steps exceeded ; reset counters
        !
        nlanciter = nlanciter_init
        nlanc = 0
        !
        llanczos = .false.
        !
        ! check the eigenvalue, if it's negative continue take the projection, otherwise push again ...
        !
        IF ( lowest_eigval < 0.0_DP ) THEN
           WRITE (*,*) "ARTn: found negative eigenvalue, using eigvec:", norm2(eigenvec)
           ! CALL center(eigenvec,nat)
           !WRITE (406,*) "PRIMCOORD"
           !WRITE (406,*) nat, "1"
           !DO na=1,nat
           !   write (406,*) tau(:,na)*alat*0.52917720859_DP,eigenvec(:,na)
           !ENDDO
           !CLOSE (UNIT = 406,STATUS = 'KEEP')
           leigen = .true.
           neigenstep = 0
           lperp = .false.
           istepperp = 0
        ELSE
           leigen = .false.
           lperp = .true.
           istepperp =  0
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

  OPEN( unit = iunart, file = filnam, form = 'formatted', status = 'unknown', iostat = ios )

  WRITE (UNIT = iunart, FMT = * ) lpush_init, lperp, npush, nlanc, nlanciter, &
       & neigenstep, istepperp, llanczos, leigen, &
       & push(:,:), eigenvec(:,:), lowest_eigval, &
       & current_step_size, convcrit, pushdir(:,:)
  CLOSE (UNIT = iunart, STATUS = 'KEEP')

END SUBROUTINE artn 


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

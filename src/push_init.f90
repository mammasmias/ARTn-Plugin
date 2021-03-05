
SUBROUTINE push_init (nat, tau, at, alat, idum, push_ids, dist_thr, add_const, init_step_size, push, mode)
  !
  ! subroutine that generates the initial push; options are specified by mode: 
  !           (1) 'all' generates a push on all atoms 
  !           (2) 'list' generates a push on a list of atoms
  !           (3) 'rad' generates a push on a list of atoms and all atoms within dist_thr 
  ! the user should supply: number and list of atoms to push; and add_constraints on these atoms
  !
  USE artn_params, ONLY : DP
  IMPLICIT none
  INTEGER, INTENT(IN) :: nat,idum
  INTEGER :: na, ia 
  INTEGER, INTENT(IN) :: push_ids(nat)
  REAL(DP), INTENT(IN) :: alat, dist_thr, init_step_size
  REAL(DP), INTENT(IN) :: tau(3,nat), at(3,3), add_const(4,nat)
  CHARACTER(LEN=4), INTENT(IN) :: mode
  REAL(DP), INTENT(OUT) :: push(3,nat)
  REAL(DP), EXTERNAL :: ran3, dnrm2 
  REAL(DP) :: dr2, pushat(3)
  REAL(DP) :: dist(3), tau0(3) 
  LOGICAL :: lvalid
  INTEGER :: atom_displaced(nat)
  !
  push(:,:) = 0.d0
  atom_displaced(:) = 0
  lvalid = .false.
  !
  !  read the list of pushed atoms
  !
  IF ( mode == 'all' ) THEN
     ! displace all atoms 
     atom_displaced(:) = 1
  ELSE IF ( mode == 'list' ) THEN 
     ! displace only atoms in list 
     DO na=1,nat
        IF (ANY(push_ids == na)) THEN
           atom_displaced(na) = 1
        ENDIF
     ENDDO
  ELSE IF ( mode == 'rad' ) THEN 
     ! displace only atoms in list and all atoms within chosen a cutoff radius ...
     DO na=1,nat
        IF (ANY(push_ids == na)) THEN
           atom_displaced(na) = 1           
           !
           tau0 = tau(:,na)*alat
           DO ia = 1,nat
              IF (ia /= na) THEN 
                 dist(:) = tau(:,ia)*alat - tau0(:)
                 
                 CALL pbc( dist, at, alat)
                 IF ( dnrm2(3,dist,1) <= dist_thr ) THEN
                    ! found an atom within dist_thr 
                    atom_displaced(ia) = 1
                 ENDIF
              ENDIF
           ENDDO
        ENDIF
     ENDDO
              
  ENDIF
  !
  !
  DO na=1,nat
     IF (atom_displaced(na) == 1 ) THEN
        DO
           push(:,na) = (/0.5_DP - ran3(idum),0.5_DP - ran3(idum),0.5_DP - ran3(idum)/)
           dr2 = push(1,na)**2 + push(2,na)**2 + push(3,na)**2
           ! check if the atom is constrained
           IF (ANY(ABS(add_const(:,na)) > 0.D0)) THEN
              ! check if the displacement is within the chosen constraint
              CALL displacement_validation(na,add_const(:,na),push(:,na),lvalid )
              IF ( .not. lvalid ) CYCLE
           ENDIF
           IF ( dr2 < 0.25_DP ) EXIT
        ENDDO
     ENDIF
  ENDDO
  ! if all atoms are pushed center the push vector to avoid translational motion 
  IF ( mode == 'all')  CALL center(push(:,:), nat)
  ! 
  !
  ! normalize so that maxval of the vector is 1.0
  !
  push(:,:) = push(:,:)/MAXVAL(ABS(push(:,:)))
  ! scale initial push vector according to step size 
  push(:,:) = init_step_size*push(:,:)

END SUBROUTINE push_init

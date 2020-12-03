
SUBROUTINE push_init (nat, idum, push_ids, add_const, init_step_size, push,mode)
  !
  ! subroutine that generates a random push to a list of atoms specified by user;
  !
  !
  ! the user should supply: number and list of atoms to push; and add_constraints on these atoms
  !
  USE kinds, ONLY : DP
  IMPLICIT none
  INTEGER, INTENT(IN) :: nat,idum
  INTEGER :: na
  INTEGER, INTENT(IN) :: push_ids(nat)
  REAL(DP), INTENT(IN) :: init_step_size
  REAL(DP), INTENT(IN) :: add_const(4,nat)
  CHARACTER(LEN=4), INTENT(IN) :: mode
  REAL(DP), INTENT(OUT) :: push(3,nat)
  REAL(DP), EXTERNAL :: ran3
  REAL(DP) :: dr2, pushat(3)
  LOGICAL :: lvalid
  INTEGER :: atom_displaced(nat)
  !
  push(:,:) = 0.d0
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
        ELSE
           atom_displaced(na) = 0
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
  !
  ! scale the push according to initial step size
  !
  push(:,:) = init_step_size*push(:,:)

END SUBROUTINE push_init

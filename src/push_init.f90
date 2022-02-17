
!> @author
!!   Matic Poberznik,
!!   Miha Gunde


SUBROUTINE push_init( nat, tau, order, at, idum, push_ids, dist_thr, add_const, init_step_size, push, mode)
  !
  !> @brief
  !!   subroutine that generates the initial push; options are specified by mode: 
  !!           (1) 'all' generates a push on all atoms 
  !!           (2) 'list' generates a push on a list of atoms
  !!           (3) 'rad' generates a push on a list of atoms and all atoms within dist_thr 
  !!   the user should supply: number and list of atoms to push; and add_constraints on these atoms
  !
  !> @param [in]    nat		    Size of list: number of atoms 
  !> @param [in]    idum	    looks like it is the seed for random gen
  !> @param [in]    push_ids	    List of atoms on which apply a push
  !> @param [in]    order	    order of atom in the list
  !> @param [in]    dist_thr	    Threshold on the distance interatomic
  !> @param [in]    init_step_size  length of initial step
  !> @param [in]    tau		    atomic position
  !> @param [in]    at		    Box length
  !> @param [inout] add_const	    list of atomic constrain
  !> @param [in]    mode	    Actual kind displacement 
  !> @param [out]   push	    list of push applied on the atoms (ORDERED)
  !
  USE artn_params, ONLY : DP, ran3, istep, iunartout, warning
  IMPLICIT none
  ! -- ARGUMENTS
  INTEGER,          INTENT(IN)  :: nat,idum
  INTEGER,          INTENT(IN)  :: push_ids(nat)
  INTEGER,          INTENT(IN)  :: order(nat)           !%! f: i --> id
  REAL(DP),         INTENT(IN)  :: dist_thr,    &
                                   init_step_size
  REAL(DP),         INTENT(IN)  :: tau(3,nat),  &
                                   at(3,3)
  REAL(DP),         INTENT(INOUT) ::  add_const(4,nat)
  CHARACTER(LEN=4), INTENT(IN)  :: mode
  REAL(DP),         INTENT(OUT) :: push(3,nat)
  !
  ! -- LOCAL VARIABLE
  INTEGER :: na, ia , iglob
  REAL(DP) :: dr2, pushat(3)
  REAL(DP) :: dist(3), tau0(3) 
  LOGICAL :: lvalid
  INTEGER :: atom_displaced(nat)
  REAL(DP), EXTERNAL :: dnrm2, fpbc
  !
  push(:,:) = 0.d0
  atom_displaced(:) = 0
  lvalid = .false.

  !write(iunartout,*)" PUSH_INIT"
  !write(*,*)" PUSH_INIT"


  !
  !  read the list of pushed atoms
  !
  IF ( mode == 'all' ) THEN
     ! displace all atoms 
     atom_displaced(:) = 1

  ELSE IF ( mode == 'list' ) THEN 
     ! displace only atoms in list 
     DO na=1,nat
        iglob = order(na)
        IF( ANY(push_ids == iglob) )THEN
           atom_displaced(na) = 1
           !print*, " * PUSH_INIT::Atom_displeced", na, atom_displaced(na), order(na), push_ids
        ENDIF
     ENDDO

  ELSE IF ( mode == 'rad' ) THEN 
     if( sum(push_ids) == 0) &
       call warning( iunartout, "PUSH_INIT()",&
            "push_mode = 'rad' need a list of atoms: define push_ids keyword ", push_ids )
     ! displace only atoms in list and all atoms within chosen a cutoff radius ...
     DO na=1,nat
        !IF (ANY(push_ids == na)) THEN
        iglob = order(na)
        IF( ANY(push_ids == iglob) )THEN
           atom_displaced(na) = 1   !%! Array based on local index i           
           !
           tau0 = tau(:,na)
           DO ia = 1,nat
              IF( ia /= na ) THEN 
                 dist(:) = tau(:,ia) - tau0(:)
                 
                 CALL pbc( dist, at)
                 IF ( dnrm2(3,dist,1) <= dist_thr ) THEN
                    ! found an atom within dist_thr 
                    atom_displaced(ia) = 1
                 ENDIF
              ENDIF
           ENDDO
        ENDIF
     ENDDO

  ENDIF

  !%! Order the ADD_CONST Array: i = order(i)
  add_const(:,:) = add_const(:,order(:))



  !
  !%! Now All the information are converted in local index
  INDEX:DO na=1,nat
     IF (atom_displaced(na) == 1 ) THEN
        RDM:DO
           push(:,na) = (/0.5_DP - ran3(idum),0.5_DP - ran3(idum),0.5_DP - ran3(idum)/)
           dr2 = push(1,na)**2 + push(2,na)**2 + push(3,na)**2

           ! check if the atom is constrained
           IF(ANY(ABS(add_const(:,na)) > 0.D0)) THEN
              ! check if the displacement is within the chosen constraint
              CALL displacement_validation(na,add_const(:,na),push(:,na),lvalid )
              IF( .not. lvalid )THEN; CYCLE RDM      ! draw another random vector
              ELSE;                   CYCLE INDEX    ! go to the next atom index
              ENDIF
           ENDIF
           IF ( dr2 < 0.25_DP ) EXIT

        ENDDO RDM 
     ENDIF
  ENDDO INDEX
  ! if all atoms are pushed center the push vector to avoid translational motion 
  IF ( mode == 'all')  CALL center(push(:,:), nat)
  ! 
  !
  ! normalize so that maxval of the vector is 1.0
  !
  push(:,:) = push(:,:)/MAXVAL(ABS(push(:,:)))

  ! scale initial push vector according to step size (ORDERED) 
  push(:,order(:)) = init_step_size*push(:,:)


! do na =1,nat
!    !if( ANY(ABS(add_const(:,na)) > 0.D0) )print*, "PUSH_INIT:", na, order(na), push(:,na)
!    if( atom_displaced(na) == 1 )print*, "PUSH_INIT:", na, order(na), push(:,order(na))
!    !print*, "PUSH_INIT:", na, order(na), push(:,na)
! enddo




END SUBROUTINE push_init

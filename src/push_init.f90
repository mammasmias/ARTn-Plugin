!> @author
!!   Matic Poberznik,
!!   Miha Gunde
!!
SUBROUTINE push_init( nat, tau, order, lat, idum, push_ids, dist_thr, add_const, init_step_size, push, mode)
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
  USE units, only : DP
  USE artn_params, ONLY : ran3, iunartout, warning, force_step, random_array
  IMPLICIT none
  ! -- ARGUMENTS
  INTEGER,          INTENT(IN)  :: nat,idum
  INTEGER,          INTENT(IN)  :: push_ids(nat)
  INTEGER,          INTENT(IN)  :: order(nat)           !%! f: i --> id
  REAL(DP),         INTENT(IN)  :: dist_thr,    &
                                   init_step_size
  REAL(DP),         INTENT(IN)  :: tau(3,nat),  &
                                   lat(3,3)
  REAL(DP),         INTENT(INOUT) ::  add_const(4,nat)
  CHARACTER(*),     INTENT(IN)  :: mode
  REAL(DP),         INTENT(OUT) :: push(3,nat)
  !
  ! -- LOCAL VARIABLE
  INTEGER :: na, ia , iglob
  REAL(DP) :: dr2, bias(3,nat)
  REAL(DP) :: dist(3), tau0(3), vmax
  LOGICAL :: lvalid, lcenter
  INTEGER :: atom_displaced(nat)
  REAL(DP), EXTERNAL :: dnrm2, fpbc
  !
  push(:,:) = 0.0_DP
  atom_displaced(:) = 0
  lvalid = .false.
  lcenter = .false.
  bias = 0.0_DP

  !
  !  read the list of pushed atoms
  !
  SELECT CASE( trim(mode) )

    CASE( 'all' )  !! displace all atoms 

      !atom_displaced(:) = 1
      bias = 1.0_DP
      lcenter = .true.
 

    CASE( 'list' ) !! displace only atoms in list 

      DO na=1,nat
         !iglob = order(na)
         !IF( ANY(push_ids == iglob) )THEN
         IF( ANY(push_ids == na) )THEN
            atom_displaced(na) = 1
            bias(:,na) = 1.0_DP
         ENDIF
      ENDDO

     
    CASE( 'rad' ) !! radius around atoms

      IF( sum(push_ids) == 0 ) &
       call warning( iunartout, "PUSH_INIT()",&
            "push_mode = 'rad' need a list of atoms: define push_ids keyword ", push_ids )

      ! displace only atoms in list and all atoms within chosen a cutoff radius ...
      DO na=1,nat
         IF( ANY(push_ids == na) )THEN
           !iglob = order(na)
           !IF( ANY(push_ids == iglob) )THEN
           !atom_displaced(na) = 1   !%! Array based on local index i           
           bias(:,na) = 1.0_DP
           !
           tau0 = tau(:,na)
           DO ia = 1,nat
              IF( ia /= na ) THEN 
                 dist(:) = tau(:,ia) - tau0(:)
                 
                 CALL pbc( dist, lat)
                 IF ( dnrm2(3,dist,1) <= dist_thr ) THEN
                    ! found an atom within dist_thr 
                    !atom_displaced(ia) = 1
                    bias(:,ia) = 1.0_DP
                 ENDIF
              ENDIF
           ENDDO
        ENDIF
      ENDDO



    ! ...Generalize the bias 
    CASE( 'bias_force' )
      !! define an array bias = force_step that is used in random_array()
      !!  to bias the 
      bias = force_step / dnrm2( 3*nat, force_step, 1)
      lcenter = .true.



    CASE( 'list_force' )
      !! Equivalent to Miha list on the force
      !bias = merge( 1.0_DP, 0.0_DP, force_step > 1e-16 )  !! Component by component
      do na=1,nat
         bias(:,na) = merge( 1.0_DP, 0.0_DP, norm2(force_step(:,na)) > 1e-16 )  !! On the norm(force) as Miha did
         !print*, "push_init", na, bias(:,na), push_ids(na)
      enddo


  END SELECT


  !%! Order the ADD_CONST Array: i = order(i)
  !add_const(:,:) = add_const(:,order(:))  !! Now all array are ordered


  !
  ! ...Now All the information are converted in local index

  INDEX:DO na=1,nat


        ia = 0
        RDM:DO
           ia = ia + 1
 
           push(:,na) = (/ (0.5_DP - ran3(idum)) * bias(1,na),   &
                           (0.5_DP - ran3(idum)) * bias(2,na),   &
                           (0.5_DP - ran3(idum)) * bias(3,na) /)
           dr2 = push(1,na)**2 + push(2,na)**2 + push(3,na)**2

           !if( atom_displaced(na) == 1 ) &
           !  print'("PUSH_INIT::DRAW ",i0,4(x,g10.3),x,i0)', na, push(:,na), dr2, ia


           ! check if the atom is constrained
           IF( ANY(ABS(add_const(:,na)) > 0.0_DP) ) THEN

              ! check if the displacement is within the chosen constraint
              CALL displacement_validation( add_const(:,na), push(:,na), lvalid )

              IF( .not. lvalid )THEN;      CYCLE RDM      ! draw another random vector
              ELSEIF( dr2 < 0.25_DP )THEN; CYCLE INDEX    ! go to the next atom index
              ENDIF

           ENDIF

           ! ...Isotrop Random Condition
           IF ( dr2 < 0.25_DP ) CYCLE INDEX  !! next atom

        ENDDO RDM 


  ENDDO INDEX
  

  !
  ! ...if all atoms are pushed center the push vector to avoid translational motion 
  !IF ( mode == 'all')  CALL center(push(:,:), nat)
  IF( lcenter )CALL center(push(:,:), nat)


  !
  ! ...normalize so that the norm of the largest displacement of an atom is 1.0
  vmax = 0.0_DP
  do na = 1,nat
     vmax = max( vmax, norm2(push(:,na)) )
  enddo
  push(:,:) = push(:,:)/ vmax
  
  !
  ! ...scale initial push vector according to step size (ORDERED) 
  push = init_step_size*push


END SUBROUTINE push_init






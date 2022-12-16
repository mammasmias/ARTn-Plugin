
!......................................................
!> @author
!!   Matic Poberznik
!!   Miha Gunde
!!   Nicolas Salles

!> @brief
!!   provide a 3 random number \in [-0.5:0.5] with norm < 0.25
! 
!> @param[in]    idum    seed for rng
!! @param[inout] vec     output vector
!
!> @note
!!   the random vector is inside a cercle of radius 0.5
!!   because x, y, z \in [-.5:.5]
!
subroutine random_displacement( idum, vec )
  !
  use units, only : DP
  use artn_params, only : ran3
  implicit none

  integer, intent(in) :: idum
  real(DP), intent(inout ) :: vec(3)

  real(DP) :: dr
  real(DP), external :: dnrm2

  RDM:DO
     vec(:) = (/ 0.5_DP - ran3(idum), 0.5_DP - ran3(idum), 0.5_DP - ran3(idum) /)
     dr = dnrm2( 3, vec, 1 )
     IF ( dr < 0.25_DP ) RETURN
  ENDDO RDM

end subroutine random_displacement


!......................................................
!> @author
!!   Matic Poberznik
!!   Miha Gunde
!!   Nicolas Salles

!> @brief 
!!   provide a random displacement to the atom's ID neighbors relative to the 
!!   threshold distance Rcut
!
!> @param[in]    idum    seed 
!> @param[in]    nat     number of atom
!> @param[in]    id      atom's Id
!> @param[in]    rcut    distance threshold
!> @param[out]   vec     output displacement
!
subroutine neigh_random_displacement( idum, nat, id, rcut, vec )
  !
  use units, only : DP, unconvert_length
  use artn_params, only : lat, tau_step, push_ids
  implicit none

  integer, intent( in ) :: idum, id, nat
  real(DP), intent( in ) :: rcut
  real(DP), intent( out ) :: vec(3,nat)

  integer :: na
  real(DP) :: x0(3), dr(3), d, rc
  real(DP), external :: dnrm2

  !
  ! -- WARNING : The position and lattice are stil in engine units
  !       we unconvert the rcut which should be in 
  !
  rc = unconvert_length( rcut )

  x0 = tau_step(:,id)
  DO na = 1,nat
     IF( id == na)cycle
     !IF( ANY(push_ids == na) )cycle
     dr(:) = tau_step(:,na) - x0(:)

     CALL pbc( dr, lat)
     d = dnrm2(3,dr,1) 
     IF( d <= rc )THEN
       ! found an atom within dist_thr 
       !call random_displacement( idum, na, vec(:,na)) 
       call random_displacement( idum, vec(:,na)) 
       !print*, id, na, d, "neigh random disp:", vec(:,na)
     ENDIF
  ENDDO


end subroutine neigh_random_displacement



!.....................................................................................................
!> @authors
!!   Matic Poberznik
!!   Miha Gunde
!!   Nicolas Salles
!
!> @brief Read field direction from file
!
!> @par Purpose
!  ============
!
!> @verbatim
!>   Read the configuration from a file formatted xyz but as we want to customise 
!>   the push the position are the push, no position means random displacement
!>   Can list only a part of particle in the system.
!> @endverbatim
!
!> @ingroup Control
!
!> @param[in]     nat       number of atoms  
!> @param[in]     idum      seed for random number generator
!> @param[out]    vec       initial push
!> @param[in]     filename  input file name
!
!> @snippet read_guess.f90 read_guess
!>
SUBROUTINE READ_GUESS( idum, nat, vec, filename )
  !
  !> [read_guess]
  use units,       only : DP, unconvert_length, read_line, parser
  use artn_params, only : warning, iunartout, dist_thr, push_ids, push_step_size
  ! use tools
  implicit none

  integer,      intent( in ) :: nat, idum
  REAL(DP),     intent( out ) :: vec(3,nat)
  character(*), intent( in ) :: filename

  character(len=256) :: line
  character(:), allocatable :: words(:)
  integer :: i, n, u0, nwords, idx, j
  logical :: ok, neiglist


  !PRINT*, "   ** ENTER IN READ_GUESS()"

  ! ...Look at the file

  inquire( file=filename, exist=ok )
  if( .not.ok )CALL warning( iunartout, 'READ_GUESS','File does not find')



  !  ...Initialization

  if( allocated(push_ids) )deallocate(push_ids)
  neiglist = .false.
  if( dist_thr > 0.0e-8 ) neiglist = .true.
  !print*, "DIST_THR", dist_thr, unconvert_length( dist_thr )



  ! ...Read file

  OPEN( newunit=u0, file=filename, ACTION="READ" )

  READ(u0,*) n
  allocate(push_ids(n))
  READ(u0,*)

  do i = 1, n

     idx = 0
     call read_line( u0, line )
     nwords = parser( trim(line), " ", words )

     select case( nwords )

       !! Only the atom index
       case( 1 )
         IF( is_numeric(words(1)) )read(words(1),*) idx
         push_ids(i) = idx
         call random_displacement( idum, vec(:,idx) )
         vec(:,idx) = vec(:,idx) * push_step_size
         !print*, idx, "random disp:", vec(:,idx)


       !! Atom index and push direction constrain
       case( 2: )
         IF( is_numeric(words(1)) )then
           read(words(1),*) idx
         else
           call warning( iunartout, 'READ_GUESS', 'index  proposed are not valid', words )
         endif
         push_ids(i) = idx

         !print*, "   ** push_ids", idx
         do j = 2,4
            IF( is_numeric(trim(words(j))) )then
              read(words(j),*) vec(j-1,idx)
            !ELSEIF( words(j) == "*" )THEN         !! Idea for more flexibility 
            !  mask(j-1,idx)
            ELSE
              call warning( iunartout, 'READ_GUESS', 'Displacement propose are not valid', words )
            ENDIF
         enddo
         !!> @warning :: Maybe put a test the norm of the user vector to compare to the 
         !!   push_step_size parameters. 
         !print*, idx, "constrain disp:", vec(:,idx)


       case default
         call warning( iunartout, 'READ_GUESS', 'Empty line' )
         exit

     end select

     ! ...Add the neigbors
     if( neiglist )call neigh_random_displacement( idum, nat, idx, dist_thr, vec )

  enddo


  CLOSE( u0 )
  !> [read_guess]

CONTAINS
  !......................................................
  !> @brief
  !!   test if the string represent a number or not
  !
  !> @param[in]    string   input string
  !! @return       logical  
  !
  elemental FUNCTION is_numeric(string)
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: string
    LOGICAL :: is_numeric
    REAL :: x
    INTEGER :: e,n
    CHARACTER(len=12) :: fmt

    n = LEN_TRIM(string)
    WRITE(fmt,'("(F",I0,".0)")') n
    READ(string,fmt,IOSTAT=e) x
    is_numeric = (e == 0)
  END FUNCTION is_numeric

END SUBROUTINE READ_GUESS






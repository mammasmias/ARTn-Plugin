


SUBROUTINE move_mode(nat, force, vel, etot, nsteppos, dt_curr, alpha, alpha_init, dt_init, cmode )
  !
  ! translate specified move to appropriate force and set FIRE parameters accordingly  
  !
  use iso_c_binding, only : c_char
  USE artn_params, ONLY: DP, AMU_RY, iperp, push0 => push, push=>eigenvec, dlanc !,  &
  !
  IMPLICIT NONE

  ! -- Arguments
  INTEGER, INTENT(IN)                       :: nat

  REAL(DP), DIMENSION(3,nat), INTENT(INOUT) :: force
  REAL(DP), DIMENSION(3,nat), INTENT(INOUT) :: vel

  REAL(DP), INTENT(IN)                      :: alpha_init, dt_init
  REAL(DP), INTENT(INOUT)                   :: etot, alpha, dt_curr
  INTEGER,  INTENT(INOUT)                   :: nsteppos
  
  CHARACTER(LEN=1,KIND=c_char), INTENT(IN)  :: cmode(:)
  ! 
  ! -- Local Variable
  CHARACTER(LEN=:), allocatable             :: mode
  !CHARACTER(LEN=:), allocatable, external    :: ctrim
  REAL(DP), EXTERNAL               :: ddot,dnrm2

  !
  ! do things depending on mode of the move
  ! NOTE force units of Ry/a.u. are assumed ... 
  !

  ! .. Convert C_Char to Fortran String
  !print*, " * ARTn::MOVE_MODE::", cmode
  mode = ctrim( cmode )


  print*, " * ARTn::MOVE_MODE::mode ", trim(mode), len_trim(mode)
  print*, " * ARTn::MOVE_MODE::iperp ", iperp



  SELECT CASE( TRIM(mode) )

  CASE( 'init' )
     !
     etot = 0.D0
     vel(:,:) = 0.D0
     alpha = 0.0_DP
     dt_curr = dt_init
     nsteppos = 0
     force = push0*amu_ry/dt_curr**2

  CASE( 'perp' )
     !
     !IF( iperp .eq. 0 ) THEN
     IF( iperp - 1 .eq. 0 ) THEN  !%! Because I invrement iperp before to enter in move_mode
        ! for the first step forget previous velocity (prevent P < 0)
        etot = 0.D0
        vel(:,:) = 0.D0
        alpha = alpha_init
        dt_curr = dt_init
     ELSE
        ! subtract the components that are parallel
        vel(:,:) = vel(:,:) - ddot(3*nat,vel, 1, push, 1 )*push(:,:)/ddot(3*nat,push(:,:),1, push(:,:),1)
     ENDIF
        !
  CASE( 'lanc' )
     !
     ! set the velocity and acceleration and alpha of previous step to move correctly
     !
     etot = 0.D0
     vel(:,:) = 0.D0
     dt_curr = dt_init
     alpha = 0.D0
     nsteppos = 0
     ! the step performed should be like this now translate it into the correct force
     force(:,:) = force(:,:)*dlanc*amu_ry/dt_curr**2
     !force(:,:) = eigenvec(:,:)*dlanc*amu_ry/dt_curr**2   ! Should be like that
     !
  CASE( 'eign' )
     !
     etot = 0.D0
     vel(:,:) = 0.D0
     alpha = 0.0_DP
     dt_curr = dt_init
     nsteppos = 0
     force(:,:) = force(:,:)*amu_ry/dt_curr**2
     !force(:,:) = eigenvec(:,:)*amu_ry/dt_curr**2   ! Should be like that
     !

  CASE( 'relx' )
     !forc_thr = 10D-8    !! QE dependent
     alpha = alpha_init
     dt_curr = dt_init

  CASE default
     write(*,*) 'Problem with move_mode!'

  END SELECT
 

 CONTAINS

!................................................................................
!> Convert C_CHAR (array of letter) in string char (fortran style)
!
!> @param[in] char_in   One length char
!! @return    String Fortran Style

function ctrim( char_in ) !BIND( C )
  use iso_c_binding, only : c_null_char
  implicit none
  !integer, intent( in ) :: n
  !character(len=1), dimension(n), intent( in ) :: char_in
  character(len=1), dimension(:), intent( in ) :: char_in
  character(len=:), allocatable :: ctrim
  character(len=256) :: ctmp
  !
  integer :: i, idx, n
  !
  n = size( char_in )
  !write(*,*) " * CTRIM::",n, char_in
  idx=1
  do i = 1,n
     if( char_in(i) == " " ) cycle
     if( char_in(i) == c_null_char) exit
     ctmp(idx:idx) = char_in( i )
     idx = idx + 1
  enddo
  allocate(character(idx-1) :: ctrim )
  !print*, ' size ',idx, len(ctrim), len(trim(ctmp)), trim(ctmp)
  do i = 1,idx-1
     ctrim(i:i) = ctmp(i:i)
  enddo
  !ctrim(idx:idx) = c_null_char
  !print*, 'leave ctrim', len(ctrim), ctrim
  !
end function ctrim

END SUBROUTINE move_mode




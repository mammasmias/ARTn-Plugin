
!> @author
!!   Nicolas Salles
!!   Matic Poberznik
!!   Miha Gunde

!> @brief
!!   Allow to load the new minimum found at the end of the research
!
!> @param[in]      nat     number of atom
!! @param[inout]   pos     atomic position
!
SUBROUTINE move_nextmin( nat, pos )
  !
  USE UNITS, only : DP
  USE artn_params, only : tau_nextmin, etot_init, etot_final
  implicit none

  ! -- arguments
  INTEGER, INTENT(in) :: nat
  REAL(DP), INTENT(inout) :: pos(3,nat)

  !LOGICAL, intent(in) :: lnext


  if( .not.allocated(tau_nextmin) )then

    write(*,'(5x,"*** ARTn:: Not other minimum saved ")')
    return

  else

    !pos = unconvert_length( tau_nextmin )
    pos = tau_nextmin
    etot_init = etot_final
    write(*,'(5x,"*** ARTn:: Next minimum loaded ")')

  endif

END SUBROUTINE move_nextmin


!> @author
!!   Nicolas Salles
!!   Matic Poberznik
!!   Miha Gunde

!> @brief 
!!   Save the configuration of new minimum if the distance from the 
!!   initial minimum is greater than a threshold Rc
!
!> @param[in]      nat     number of atom
!! @param[inout]   pos     atomic position
!
SUBROUTINE save_min( nat, pos )
  USE UNITS, only : DP, unconvert_length
  USE artn_params, only : tau_init, lat, tau_nextmin
  implicit none

  ! -- Arguments
  INTEGER, INTENT(in) :: nat
  REAL(DP), INTENT(inout) :: pos(3,nat)  ! it is in ARTn units (bohr)

  ! -- local variables
  REAL(DP) :: delr(3,nat), dr2, dr1, Rc
  REAL(DP), external :: dsum

  Rc = unconvert_length( 0.5_DP )

  call compute_delr( nat, pos, tau_init, lat, delr )
  !call sum_force( delr, nat, dr1 )
  dr1 = dsum( 3*nat, delr )

  !...Comparison in bohr
  if( dr1 > Rc )then

    ! ...We never saved another minimum
    if( .not.allocated(tau_nextmin) )then
      allocate( tau_nextmin, source=pos )
    else

     ! ...We already saved a minimum and the auestion is:
     !!   is it new/farther as the init/start position
     call compute_delr( nat, tau_nextmin, tau_init, lat, delr )
     !call sum_force( delr, nat, dr2 )
     dr2 = dsum( 3*nat, delr )

     ! ...We save the minimum farther than the initial positon
     if( dr1 > dr2 ) tau_nextmin = pos

    endif

  else
    write(*,*) " *** ARTn:: The minimum found is the initial minimum "

  endif

END SUBROUTINE save_min




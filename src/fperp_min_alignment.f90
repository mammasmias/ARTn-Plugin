
!> @author
!!   Matic Poberjnik
!!   Miha Gunde
!!   Nicolas Salles

!.......................................................................................
logical function fperp_min_alignment( thr1, thr2 )result( res )
  !> @brief 
  !!   compute the 2 condition:
  !!    - eigenVec has been suddenlly changed
  !!    - direction of minimum is perp to the last push
  !
  !> @param[in] thr1    threshold on the eigenvec alignement
  !> @param[in] thr2    threshold in the fperp - direction of minimum alignment
  !
  USE units, only : DP
  USE artn_params, only : a1, tau_step, tau_init, push, natoms
  implicit none

  real(DP), intent(in) :: thr1, thr2

  REAL(DP) :: min_dir(3,natoms), dtmp
  REAL(DP), external :: ddot

  min_dir = tau_step - tau_init
  min_dir = min_dir / NORM2( min_dir )
  dtmp = ddot(3*natoms,min_dir,1,push,1)

  !! IF eigenVec change suddenlly AND direction of minimum is perp to the last push
  res = ( a1 < thr1 .AND. ABS(dtmp) < thr2 )

end function fperp_min_alignment

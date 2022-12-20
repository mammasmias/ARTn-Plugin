
! ....................................................................................
!> @author
!!   Matic Poberznik
!!   Miha Gunde
!!   Nicolas Salles

!> @brief
!!   subroutine that sums the forces on all atoms and returns the total force
!
!> @param [in]  nat       Size of list: number of atoms
!> @param [in]  force     List of atomic forces
!> @param [out] force_tot Sum of the list of forces
!
SUBROUTINE sum_force( force, nat, force_tot )
  !
  USE artn_params, ONLY : DP
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nat
  REAL(DP), INTENT(IN) :: force(3,nat)
  REAL(DP), INTENT(OUT) :: force_tot
  INTEGER :: na

  force_tot = 0.0
  DO na = 1, nat
     !force_tot = force_tot + force(1,na)**2 + force(2,na)**2 + force(3,na)**2
     force_tot = force_tot + dot_product( force(:,na), force(:,na) )
  ENDDO

  force_tot = SQRT(force_tot)

END SUBROUTINE sum_force






! ....................................................................................
!> @author
!!   Matic Poberznik
!!   Miha Gunde
!!   Nicolas Salles

!> @brief
!!   sum the component square of the field in the mood of ddot of lib lapack
!!   The unroll loop can be faster than previous version
!
!> @param[in]   n     number of field's component 
!! @param[in]   f     field f(n)
!! @return      res   sum of square if the field f
!
FUNCTION dsum( n, f )result( res )
  !
  use units, only : DP
  implicit none

  integer, intent(in) :: n
  real(DP), intent(in) :: f(*)

  integer :: i, m, mp1
  real(DP) :: tmp, res

  res = 0.0_DP
  IF( n <= 0 )return

  ! ...Sum the unfit dimension
  tmp = 0.0_DP
  m = mod(n,5)
  if( m /= 0 )then
    do i = 1,m
       tmp = tmp + f(i)**2
    enddo
    if( n < 5 )then
      res = tmp
      return
    endif
  endif

  ! ...Unroll the loop
  mp1 = m + 1
  do i = mp1,n,5
     tmp = tmp + f(i)**2 + f(i+1)**2 + f(i+2)**2 + f(i+3)**2 + f(i+4)**2 
  enddo
  res = tmp
  return
END FUNCTION dsum






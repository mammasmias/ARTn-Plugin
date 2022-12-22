
!> @author
!!  Matic Poberznik
!!  Miha Gunde
!!  Nicolas Salles

!> @brief
!!   A function that takes into account periodic boundary conditions,
!!   based on the pbc function of the contraints_module of QE
!
!> @param [inout] vec   input vector in atomic units
!> @param [in]    at    lattice vectors
!
SUBROUTINE pbc( vec, at )

  USE artn_params, ONLY : DP 
  IMPLICIT none 
  ! -- ARGUMENTS
  REAL(DP), INTENT(INOUT) :: vec(3) !> input vector in atomic units  
  REAL(DP), INTENT(IN) :: at(3,3)   !> lattice vectors
  ! -- LOCAL VARIABLES
  REAL(DP) :: bg(3,3) ! inverse of at(3,3) 
  !
  ! calculate the reciprocal lattice parameters of at 
  ! 
  CALL invmat3x3(at,bg)
  !
  ! convert to crystal coords 
  vec(:) = matmul(vec(:),bg(:,:))

  ! move the vector to original box  
  vec(:) = vec(:) - anint(vec(:))

  ! convert back to cartesian coordinates 
  vec(:) = matmul(at(:,:), vec(:))
  !   
END SUBROUTINE pbc


!> @author
!!  Matic Poberznik
!!  Miha Gunde
!!  Nicolas Salles

  !> @brief 
  !!   A function that takes into account periodic boundary conditions,
  !!   based on the pbc function of the contraints_module of QE
  !

function fpbc( vec, at )RESULT( res )
  ! 
  USE artn_params, ONLY : DP
  IMPLICIT none

  REAL(DP), INTENT(IN)  :: vec(3) ! input vector in atomic units  
  REAL(DP), INTENT(IN)  :: at(3,3)   ! lattice vectors
  REAL(DP)              :: res(3)   ! lattice vectors

  REAL(DP)                :: bg(3,3) ! inverse of at(3,3) 
  !
  ! calculate the reciprocal lattice parameters of at 
  ! 
  CALL invmat3x3(at,bg)
  !
  ! convert to crystal coords: U = R.B^-1 
  res(:) = matmul(vec(:),bg(:,:))
  ! move the vector to original box  
  res(:) = res(:) - anint(res(:))
  ! convert back to cartesian coordinates: R = B.U ??? 
  res(:) = matmul(at(:,:), res(:))
  !   
END function fpbc






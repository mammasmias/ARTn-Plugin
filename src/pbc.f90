!SUBROUTINE pbc(vec, at, alat)
SUBROUTINE pbc( vec, at )
  USE artn_params, ONLY : DP 
  IMPLICIT none 
  ! 
  ! A function that takes into account periodic boundary conditions,
  ! based on the pbc function of the contraints_module of QE
  ! 
  REAL(DP), INTENT(INOUT) :: vec(3) ! input vector in atomic units  
  REAL(DP), INTENT(IN) :: at(3,3)   ! lattice vectors
  !REAL(DP), INTENT(IN) :: alat      ! a lattice parameter
  REAL(DP) :: bg(3,3) ! inverse of at(3,3) 
  !
  ! calculate the reciprocal lattice parameters of at 
  ! 
  CALL invmat3x3(at,bg)
  !
  ! convert to crystal coords 
  !vec(:) = matmul(vec(:),bg(:,:))/alat
  vec(:) = matmul(vec(:),bg(:,:))
  ! move the vector to original box  
  vec(:) = vec(:) - anint(vec(:))
  ! convert back to cartesian coordinates 
  !vec(:) = matmul(at(:,:), vec(:))*alat
  vec(:) = matmul(at(:,:), vec(:))
  !   
END SUBROUTINE pbc

function fpbc( vec, at )RESULT( res )
  USE artn_params, ONLY : DP
  IMPLICIT none
  ! 
  ! A function that takes into account periodic boundary conditions,
  ! based on the pbc function of the contraints_module of QE
  ! 
  REAL(DP), INTENT(IN)  :: vec(3) ! input vector in atomic units  
  REAL(DP), INTENT(IN)  :: at(3,3)   ! lattice vectors
  REAL(DP)              :: res(3)   ! lattice vectors
  !REAL(DP), INTENT(IN) :: alat      ! a lattice parameter
  REAL(DP)                :: bg(3,3) ! inverse of at(3,3) 
  !
  ! calculate the reciprocal lattice parameters of at 
  ! 
  CALL invmat3x3(at,bg)
  !
  ! convert to crystal coords 
  res(:) = matmul(vec(:),bg(:,:))
  ! move the vector to original box  
  res(:) = res(:) - anint(res(:))
  ! convert back to cartesian coordinates 
  res(:) = matmul(at(:,:), res(:))
  !   
END function fpbc






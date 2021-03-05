SUBROUTINE pbc(vec, at, alat)
  USE artn_params, ONLY : DP 
  IMPLICIT none 
  ! 
  ! A function that takes into account periodic boundary conditions,
  ! based on the pbc function of the contraints_module of QE
  ! 
  REAL(DP), INTENT(INOUT) :: vec(3) ! input vector in atomic units  
  REAL(DP), INTENT(IN) :: at(3,3)   ! lattice vectors
  REAL(DP), INTENT(IN) :: alat      ! a lattice parameter
  REAL(DP) :: bg(3,3) ! inverse of at(3,3) 
  !
  ! calculate the reciprocal lattice parameters of at 
  ! 
  CALL invmat3x3(at,bg)
  !
  ! convert to crystal coords 
  vec(:) = matmul(vec(:),bg(:,:))/alat
  ! move the vector to original box  
  vec(:) = vec(:) - anint(vec(:))
  ! convert back to cartesian coordinates 
  vec(:) = matmul(at(:,:), vec(:))*alat
  !   
END SUBROUTINE pbc


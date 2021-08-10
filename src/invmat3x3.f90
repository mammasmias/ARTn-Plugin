
!> @author
!!  Matic Poberjnik,
!!  Miha Gunde


SUBROUTINE invmat3x3(mat,inv)
  USE artn_params, ONLY : DP 
  IMPLICIT none
  !
  !> @brief
  !!   Subroutine that calculates the inverse of a 3x3 matrix 
  !
  !> @param [in]  mat	Matrix to inverse
  !> @param [out] inv	Inverse of the Matrix
  !
  REAL(DP), INTENT(IN) :: mat(3,3)
  REAL(DP), INTENT(OUT) :: inv(3,3)
  !
  REAL(DP) :: det, invdet
  ! 
  det = 0.0_DP
  !
  ! calculate the determinant ...
  ! 
  det = det + mat(1,1)*mat(2,2)*mat(3,3) &
            + mat(1,2)*mat(2,3)*mat(3,1) &
            + mat(1,3)*mat(2,1)*mat(3,2) &
            - mat(1,3)*mat(2,2)*mat(3,1) &
            - mat(1,2)*mat(2,1)*mat(3,3) &
            - mat(1,1)*mat(2,3)*mat(3,2)
  ! invert the determinant  
  invdet = 1/det
  ! calculate the inverse matrix
  inv(1,1) = invdet  * ( mat(2,2)*mat(3,3) - mat(2,3)*mat(3,2) )
  inv(2,1) = -invdet * ( mat(2,1)*mat(3,3) - mat(2,3)*mat(3,1) )
  inv(3,1) = invdet  * ( mat(2,1)*mat(3,2) - mat(2,2)*mat(3,1) )
  inv(1,2) = -invdet * ( mat(1,2)*mat(3,3) - mat(1,3)*mat(3,2) )
  inv(2,2) = invdet  * ( mat(1,1)*mat(3,3) - mat(1,3)*mat(3,1) )
  inv(3,2) = -invdet * ( mat(1,1)*mat(3,2) - mat(1,2)*mat(3,1) )
  inv(1,3) = invdet  * ( mat(1,2)*mat(2,3) - mat(1,3)*mat(2,2) )
  inv(2,3) = -invdet * ( mat(1,1)*mat(2,3) - mat(1,3)*mat(2,1) )
  inv(3,3) = invdet  * ( mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1) )
  ! 
END SUBROUTINE invmat3x3

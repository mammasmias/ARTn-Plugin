SUBROUTINE diag(n, A, eigvals, vec)
  USE artn_params,            ONLY : DP
  !! assuming a general square matrix (can be nonsymmetric).
  !! On output A is overwritten by eigenvectors in rows, if vec=0, then
  !! A is just 0.0 on output.
  !!
  !! n       : dimension of matrix A
  !! A       : matrix to be diagonalised, overwritten by eigenvectors on output
  !! eigvals : output vector of eigenvalues, not sorted!
  !! vec     : 0 if don't want to compute eigenvectors, 1 otherwise
  !!
  IMPLICIT NONE
  INTEGER,              intent(in) :: n
  REAL(DP), DIMENSION(n,n), intent(inout) :: A
  REAL(DP), DIMENSION(n),   intent(out) :: eigvals
  INTEGER,              intent(in) :: vec
  REAL(DP), DIMENSION(n) :: eigvals_i !! imaginary part of the eigenvalues
  REAL(DP), DIMENSION(n,n) :: eigvec
  INTEGER :: lda
  INTEGER :: lwork
  REAL(DP) :: Dummy(1000)
  INTEGER :: info
  CHARACTER(len=1) :: getvec
  getvec = 'N'
  if( vec == 1 ) getvec='V'
  lda = n
  eigvals_i(:) = 0.0
  eigvec(:,:) = 0.0
  !! test workspace
  lwork = -1
  !call sgeev('N',getvec, n, A, lda, eigvals, eigvals_i, &
  !     dummy, 1, eigvec, n, dummy, lwork, info)
  call dgeev('N',getvec, n, A, lda, eigvals, eigvals_i, &
       dummy, 1, eigvec, n, dummy, lwork, info)
  !! choose optimal size of workspace (as in example from intel website)
  lwork = min( 1000, nint(dummy(1)) )
  !! compute stuffs
  !call sgeev('N',getvec, n, A, lda, eigvals, eigvals_i, &
  !     dummy, 1, eigvec, n, dummy, lwork, info)
  call dgeev('N',getvec, n, A, lda, eigvals, eigvals_i, &
       dummy, 1, eigvec, n, dummy, lwork, info)
  !! overwrite a on output
  A(:,:) = eigvec(:,:)
END SUBROUTINE diag

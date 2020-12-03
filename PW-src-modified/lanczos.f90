
SUBROUTINE lanczos( nat, alat, force, vel, acc, alpha_init, dt, &
     v_in, dlanc, nlanciter, nlanc, lowest_eigval, lowest_eigvec, pushdir, prfx,tmpdir )
  USE artn_params,            ONLY: DP, Vmat, H, force_old, initialize_lanczos 
  !
  ! Lanczos subroutine for the ARTn algorithm; based on the lanczos subroutine as written by M. Gunde
  !
  IMPLICIT NONE
  INTEGER,                INTENT(IN) :: nat
  REAL(DP), DIMENSION(3,nat), INTENT(IN) :: v_in
  REAL(DP), DIMENSION(3,nat), INTENT(IN) :: pushdir
  REAL(DP), INTENT(IN) :: dlanc
  REAL(DP), INTENT(IN) :: alpha_init, dt, alat
  CHARACTER(LEN=255), INTENT(IN) :: tmpdir, prfx
  REAL(DP), DIMENSION(3,nat), INTENT(INOUT) :: vel, acc
  REAL(DP), DIMENSION(3,nat), INTENT(INOUT) :: force
  REAL(DP), DIMENSION(3,nat), INTENT(INOUT) :: lowest_eigvec
  REAL(DP), INTENT(INOUT) :: lowest_eigval
  INTEGER, INTENT(INOUT) :: nlanciter
  INTEGER, INTENT(INOUT) :: nlanc
  ! 
  INTEGER :: i, j, io, id_min
  REAL(DP), PARAMETER :: eigvec_thr = 1.0D-3, eigval_thr = 1.0D-2
  REAL(DP), ALLOCATABLE :: v1(:,:), q(:,:), eigvals(:)
  REAL(DP), ALLOCATABLE :: Vmat_mul(:,:), Hstep(:,:)
  REAL(DP) :: lowest_eigvec_tmp(3*nat)
  REAL(DP) :: dir
  REAL(DP), EXTERNAL :: ran3,dnrm2,ddot
  REAL(DP) :: alpha, beta, lowest_eigval_old, eigvec_diff, largest_eigvec_diff, eigval_diff
  ! allocate vectors and put to zero
  ALLOCATE( q(3,nat), source=0.D0 )
  ALLOCATE( v1(3,nat), source=0.D0)
  ! allocate matrices and put to zero
  ALLOCATE( Vmat_mul(3*nat,nlanciter), source=0.D0)
  ALLOCATE( Hstep(1:nlanciter,1:nlanciter), source=0.D0 )
  ! 
  ! store the eigenvalue of the previous iteration
  lowest_eigval_old = lowest_eigval
  !
  IF ( nlanc  == 0 ) THEN
     !
     ! initialize the Lanczos algorithm (allocate the matrices) 
     !
     CALL initialize_lanczos(nlanciter,nat) 
     !
     ! normalize initial vector
     ! NOTE: the inital vector to lanczos is random
     !
     v1(:,:) = v_in(:,:)
     v1(:,:) = v1(:,:) / dnrm2( 3*nat, v1, 1 )
     !
     ! store this vector
     !
     Vmat(:,:,1) = v1(:,:)
     nlanc = nlanc + 1
     !
  ELSEIF (nlanc == 1 ) THEN
     ! Generate lanczos vector, v = q - alpha*v0
     !
     ! q(:,:) now represents {Force} = -[Hessian]{R_1}, where {R_1} = {R} + {dR} = {R} + {v0}*d_lanc,
     ! and {v0} is the Lanczos vector v0(:). We need q(:) to represent [Hessian]{v0}, thus
     ! first do q(:) = -[Hessian]( {R_1} - {R} ) = -[Hessian]{dR}
     q(:,:) = force(:,:) - force_old(:,:)
     !! now do q(:) = [Hessian]{dR}/d_lanc = [Hessian]{v0}
     q(:,:) = -q(:,:) / dlanc
     !
     alpha = ddot(3*nat,Vmat(:,:,1),1,q(:,:),1)
     !
     v1(:,:) = q(:,:) - alpha*Vmat(:,:,1)
     !
     ! beta is the norm of v, used for next step
     !
     beta = dnrm2(3*nat,v1,1)
     v1(:,:) = v1(:,:) / beta
     !
     ! store the vecs for future cycles
     !
     Vmat(:,:,2) = v1(:,:)
     H(1,1) = alpha
     H(2,1) = beta
     H(1,2) = beta

     write(555,'(i4,2f9.4)') nlanc, alpha, (alpha-lowest_eigval_old)/lowest_eigval_old
     flush(555)

     lowest_eigval = alpha
     nlanc = nlanc + 1
     !
  ELSEIF (nlanc > 1 .and. nlanc <= nlanciter ) THEN
     !
     ! Generate lanczos vector, v = q - alpha*v1 - beta*v0
     !
     q(:,:) = force(:,:) - force_old(:,:)
     q(:,:) = -q(:,:)/dlanc
     !
     ! alpha = dot( v1, q )
     !
     alpha = ddot(3*nat,vmat(:,:,nlanc),1,q(:,:),1)
     H(nlanc,nlanc) = alpha
     !
     ! Do a diagonalization here ; check if eigenvalues are converged in this step
     !
     ALLOCATE( eigvals(nlanc) )
     ! store the H matrix, because its overwritten by eigvecs on diagonalization
     Hstep(:,:) = H(:,:)

     CALL diag(nlanc, Hstep(1:nlanc,1:nlanc), eigvals, 1 )

     lowest_eigval = eigvals(1)
     id_min = 1
     !
     ! get the lowest eigenvalue
     !
     DO i = 1, nlanc
        IF (eigvals(i) < lowest_eigval ) THEN
           lowest_eigval = eigvals(i)
           id_min = i
        ENDIF
     ENDDO
     !
     ! reshape the V mat for now (FIX this so it uses the DGEMM from lapack)
     !
     Vmat_mul(1:nat,:) = Vmat(1,:,:)
     Vmat_mul(nat+1:2*nat,:) = Vmat(2,:,:)
     Vmat_mul(2*nat+1:3*nat,:) = Vmat(3,:,:)
     !
     ! H now stores eigvecs of H
     ! eigvecs in coordinate space are computed as matmul(V, eigvec_H )
     !
     lowest_eigvec_tmp(:) = matmul(Vmat_mul(:,:),Hstep(:,id_min))
     lowest_eigvec(1,:) = lowest_eigvec_tmp(1:nat)
     lowest_eigvec(2,:) = lowest_eigvec_tmp(nat+1:2*nat)
     lowest_eigvec(3,:) = lowest_eigvec_tmp(2*nat+1:3*nat)
     !
     ! check if the obtained eigenvec points in the same direction as the input eigvec
     !
     dir = ddot(3*nat,lowest_eigvec,1, pushdir, 1)
     IF ( dir < 0.D0 ) THEN
        lowest_eigvec(:,:) = -1.D0*lowest_eigvec(:,:)
     ENDIF
     !
     DEALLOCATE( eigvals )
     !
     ! Check for the convergence of the lanczos eigenvalue
     !
     eigval_diff = (lowest_eigval - lowest_eigval_old)/lowest_eigval_old
     write(555,'(i4,2f9.4)') nlanc, lowest_eigval, eigval_diff
     flush(555)
     IF ( ABS(eigval_diff) <= eigval_thr ) THEN
        !
        ! lanczos has converged
        ! set max number of iternations to current iteration  
        !
        nlanciter = nlanc
        ! increase lanczos counter for last step 
        nlanc = nlanc + 1
        !
     ENDIF
     !
     ! if lanczos is not yet converged generate new matrix elements
     !
     IF ( nlanc < nlanciter ) THEN
        beta = H(nlanc,nlanc-1)
        !
        v1(:,:) = q(:,:) - alpha*vmat(:,:,nlanc) - beta*vmat(:,:,nlanc-1)
        !
        ! orthogonalize vectors in accordance with previous ones ...
        DO j = 1, nlanc - 1
           v1(:,:) = v1(:,:) - ddot(3*nat, v1 ,1, Vmat(:,:,j),1)*Vmat(:,:,j)
        ENDDO
        !
        !  do stuff with the new vector
        !
        IF ( dnrm2(3*nat, v1(:,:), 1) < 1.0D-15 ) THEN
           !
           ! new lanczos vector very small, stop (converge)
           !
           nlanciter = nlanc
           !
           ! Backtrack to initial position (sum all lanczos vectors generated)
           !
           nlanc = nlanc + 1
           !
        ELSE
           !
           ! normalize the new vector
           !
           v1(:,:) = v1(:,:)/dnrm2(3*nat,v1(:,:),1)
           Vmat(:,:,nlanc+1) = v1(:,:)
           ! store values to H
           beta = ddot(3*nat, q, 1, v1, 1)
           H(nlanc+1,nlanc ) = beta
           H(nlanc, nlanc+1) = beta
           nlanc = nlanc + 1
           !
        ENDIF
        !
     ELSE
        ! increas counter if lanczos is not converged in nlanciter
       nlanc = nlanc + 1 
     END IF
     !
  ENDIF
  ! store the force
  force_old(:,:) = force(:,:)
  !
  ! write data for next lanczos step
  !
  IF( nlanc > nlanciter ) THEN
     ! lanczos has converged, no need to write anything, delete file
     ! backtrack to initial position here
     v1(:,:) = 0.D0
     DO i = 1, nlanciter
           v1(:,:) = v1(:,:) - Vmat(:,:,i)
     ENDDO
  ENDIF
  !
  ! write data for move
  !
  CALL move_mode( nat, alat, dlanc, v1, force, &
       vel, acc, alpha_init, dt, &
       0, pushdir, 'lanc', prfx, tmpdir)
  !
  ! deallocate the matrices used only in the iteration 
  !
  DEALLOCATE( q, v1 )
  DEALLOCATE(Vmat_mul, Hstep)

END SUBROUTINE lanczos

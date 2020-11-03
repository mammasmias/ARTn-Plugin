
SUBROUTINE lanczos( nat, alat, force, vel, acc, alpha_init, dt, &
     v_in, dlanc, nlanciter, nlanc, lowest_eigval, lowest_eigvec, pushdir, prfx, tmpdir )
  USE kinds,            ONLY : DP
  !
  ! Lanczos subroutine for the ARTn algorithm; based on the lanczos subroutine as written by M. Gunde
  !
  IMPLICIT NONE
  INTEGER,                INTENT(IN) :: nat
  REAL(DP), DIMENSION(3,nat), INTENT(IN) :: v_in
  REAL(DP), DIMENSION(3,nat), INTENT(INOUT) :: vel, acc
  REAL(DP), INTENT(IN) :: alpha_init, dt, alat
  REAL(DP), DIMENSION(3,nat), INTENT(INOUT) :: force
  REAL(DP), DIMENSION(3,nat), INTENT(INOUT) :: lowest_eigvec
  REAL(DP), INTENT(INOUT) :: lowest_eigval
  REAL(DP), DIMENSION(3,nat), INTENT(IN) :: pushdir
  REAL(DP), INTENT(IN) :: dlanc
  CHARACTER(LEN=255), INTENT(IN) :: prfx, tmpdir
  !
  INTEGER :: i, j, io, id_min
  INTEGER, PARAMETER ::  iunlanc = 51
  INTEGER, INTENT(INOUT) :: nlanciter
  INTEGER, INTENT(INOUT) :: nlanc
  REAL(DP), PARAMETER :: eigvec_thr = 1.0D-3, eigval_thr = 1.0D-2
  REAL(DP), DIMENSION(3,nat) :: force_old, lowest_eigvec_old
  REAL(DP), ALLOCATABLE :: v1(:,:), q(:,:), eigvals(:)
  REAL(DP), ALLOCATABLE :: Vmat(:,:,:), Vmat_mul(:,:), H(:,:), Hstep(:,:)
  REAL(DP) :: lowest_eigvec_tmp(3*nat)
  REAL(DP) :: dir
  REAL(DP), EXTERNAL :: ran3,dnrm2,ddot
  REAL(DP) :: alpha, beta, lowest_eigval_old, eigvec_diff, largest_eigvec_diff, eigval_diff
  LOGICAL :: file_exists
  CHARACTER(LEN=255) :: filnam
  INTEGER :: ios

  ! allocate vectors and put to zero
  ALLOCATE( q(3,nat), source=0.D0 )
  ALLOCATE( v1(3,nat), source=0.D0)

  ! allocate matrices and put to zero
  ALLOCATE( Vmat_mul(3*nat,nlanciter), source=0.D0)
  ALLOCATE( Vmat(3,nat,1:nlanciter), source=0.D0 )
  ALLOCATE( H(1:nlanciter,1:nlanciter), source=0.D0 )
  ALLOCATE( Hstep(1:nlanciter,1:nlanciter), source=0.D0 )
  force_old(:,:) = 0.D0
  !
  ! open the lanczos data file
  !
  filnam = trim(tmpdir) // '/' // trim(prfx) // '.' // 'artnlanc'
  INQUIRE( file = filnam, exist = file_exists )
  OPEN( unit = iunlanc, file = filnam, form = 'formatted', status = 'unknown', iostat = ios)
  ! write(999,*) 'lanczos'
  ! write(999,*) 'open file:',trim(filnam), ios, file_exists
  ! flush(999)
  !
  ! initialize lanczos counter and variables
  !
  ! store the eignvalue and eigenvec of previous iteration
  write (*,*) "Read lowest Eigenvalue:",lowest_eigval
  lowest_eigvec_old(:,:) = lowest_eigvec(:,:)
  lowest_eigval_old = lowest_eigval
  write (*,*) "Stored Eigenvalue:",lowest_eigval_old
  !
  !
  IF ( file_exists ) THEN
     !
     ! read lanczos data from the previous iteration, if it exists
     !
     READ( UNIT = iunlanc, FMT = * )  Vmat(:,:,1:nlanciter), H(1:nlanciter,1:nlanciter),  &
         force_old(:,:)
     REWIND( UNIT = iunlanc )
     !
  END IF
  !
  ! in the first lanczos step we should give a random push to initiate the algorithm and then calc the force of the new pos with qe.
  !
  IF ( nlanc  == 0 ) THEN
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
     write (*,*) "ARTn Lanczos: initial vec:", v1(:,:)
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
     alpha = ddot(3*nat,Vmat(:,:,1),1,q(:,:),1)
     ! check if the eigenvalue
     !IF ( lowest_eigval_old /= 0.D0) THEN
     !   IF ( (alpha - lowest_eigval_old)/lowest_eigval_old .lt. eigval_thr ) THEN
     !      write (*,*) "ARTn Lanczos: new eigenvalue close to the old eigenvalue, STOPPING"
     !      ! set the lanczos counter to current number of lanczos iterations
     !      nlanciter = nlanc
     !   ENDIF
     !ENDIF
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
     write (*,*) "ARTn Lanczos: second vec:", v1(:,:)
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
     ! H now sotres eigvecs of H
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
        WRITE (*,*) "Dir is less than 0, flipping eigvec:", lowest_eigvec(:,:)
        lowest_eigvec(:,:) = -1.D0*lowest_eigvec(:,:)
     ENDIF
     WRITE (*,*) "ARTn Lanczos: Eigenvalues:", eigvals(:), "at step:", nlanc
     !
     DEALLOCATE( eigvals )
     ! calculate the largest difference between current eigenvec vs previous
     !
     !largest_eigvec_diff = 0.D0
     !DO i=1,nat
     !   DO j=1,3
     !      eigvec_diff = ABS(lowest_eigvec(j,i)) - ABS(lowest_eigvec_old(j,i))
     !      IF ( ABS(eigvec_diff) > largest_eigvec_diff ) THEN
     !         largest_eigvec_diff = eigvec_diff
     !      ENDIF
     !   ENDDO
     !ENDDO
     !
     ! Check for the convergence of the lanczos eigenvalue
     !
     eigval_diff = (lowest_eigval - lowest_eigval_old)/lowest_eigval_old
     write(555,'(i4,2f9.4)') nlanc, lowest_eigval, eigval_diff
     flush(555)
     write (*,*) "Eigenvalue difference:", eigval_diff
     IF ( ABS(eigval_diff) <= eigval_thr ) THEN
        !
        ! lanczos has converged
        !
        WRITE (*,*) "ARTn Lanczos: eigenvalue converged at step:", nlanc, "eigenthr:", eigval_thr
        !
        ! set max number of iternations to current number
        !
        nlanciter = nlanc
        !
        ! Backtrack to initial position (sum all lanczos vectors generated)
        !
        nlanc = nlanc + 1
        v1(:,:) = 0.D0
        DO i = 1, nlanciter
           v1(:,:) = v1(:,:) - Vmat(:,:,i)
        ENDDO
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
           write (*,*) "ARTn Lanczos: Stopped because Lanczos vector very small!", dnrm2(3*nat, v1(:,:), 1)
           !
           ! Backtrack to initial position (sum all lanczos vectors generated)
           !
           nlanc = nlanc + 1
           v1(:,:) = 0.D0
           DO i = 1, nlanciter
              v1(:,:) = v1(:,:) - Vmat(:,:,i)
           ENDDO
           !
        ELSE
           !
           ! normalize the new vector
           !
           v1(:,:) = v1(:,:)/dnrm2(3*nat,v1(:,:),1)
           write (*,*) "ARTn Lanczos: next vec:", v1(:,:)
           Vmat(:,:,nlanc+1) = v1(:,:)
           ! store values to H
           beta = ddot(3*nat, q, 1, v1, 1)
           H(nlanc+1,nlanc ) = beta
           H(nlanc, nlanc+1) = beta
           nlanc = nlanc + 1
           !
        ENDIF
        !
     END IF
     !
  ENDIF
  !
  ! write data for next lanczos step
  !
  IF( nlanc > nlanciter ) THEN
     ! lanczos has converged, no need to write anything, delete file
     CLOSE( UNIT = iunlanc, STATUS = 'DELETE' )
  ELSE
     ! write data and close file
     WRITE (UNIT = iunlanc, FMT = * ) Vmat(:,:,1:nlanciter), H(1:nlanciter,1:nlanciter), &
          force(:,:)
     CLOSE (UNIT = iunlanc, STATUS = 'KEEP')
     !
  ENDIF
  !
  ! write data for move
  !
  CALL move_mode( nat, alat, dlanc, v1, force, &
       vel, acc, alpha_init, dt, &
       0, pushdir, 'lanc', prfx, tmpdir)
  !
  ! deallocate
  !
  DEALLOCATE( q, v1 )
  DEALLOCATE(Vmat, Vmat_mul, H, Hstep)

END SUBROUTINE lanczos


!> @author
!!  Matic Poberznik,
!!  Miha Gunde


SUBROUTINE lanczos( nat, v_in, pushdir, force, &
     ilanc, nlanc, lowest_eigval, lowest_eigvec, displ_vec )

  USE artn_params, ONLY: DP, Vmat, H, force_old, dlanc, eval_conv_thr, &
                         lanczos_min_size
  USE units, ONLY: unconvert_length, unconvert_hessian
  !
  !> @brief
  !!   Lanczos subroutine for the ARTn algorithm;
  !! based on the lanczos subroutine as written by M. Gunde
  !!
  !! The idea is to overwrite the 'force' with the vector of desired move,
  !! according to Lanczos diagonalisation algorithm. The array 'force' on input
  !! contans the real forces, but gets overwritten with the desired move on output.
  !
  !> @param [in]      nat	       number of atoms
  !> @param [in]      v_in	      Input lanczos vector: only used in first step of each lanczos call
  !> @param [in]      pushdir	      List of Direction of push on atoms
  !> @param [inout]   force	      on input: array of Forces on the atoms
  !> @param [inout]   lowest_eigvec   Lowest eigenvector obtained by lanczos algo
  !> @param [inout]   lowest_eigval   Lowest eigenvalue obtained by lanczos algo
  !> @param [inout]   ilanc	      current step of lanczos
  !
  IMPLICIT NONE
  ! -- ARGUMENTS
  INTEGER,                    INTENT(IN) :: nat
  REAL(DP), DIMENSION(3,nat), INTENT(IN) :: v_in
  REAL(DP), DIMENSION(3,nat), INTENT(IN) :: pushdir
  REAL(DP), DIMENSION(3,nat), INTENT(IN) :: force
  INTEGER,                    INTENT(INOUT) :: ilanc
  INTEGER,                    INTENT(INOUT) :: nlanc
  REAL(DP),                   INTENT(INOUT) :: lowest_eigval
  REAL(DP), DIMENSION(3,nat), INTENT(INOUT) :: lowest_eigvec
  REAL(DP), DIMENSION(3,nat), INTENT(OUT) :: displ_vec

  ! -- LOCAL VARIABLES
  INTEGER :: i, j, io, id_min
  REAL(DP), ALLOCATABLE :: v1(:,:), q(:,:), eigvals(:)
  REAL(DP) :: dir
  REAL(DP), EXTERNAL :: ran3,dnrm2,ddot
  REAL(DP) :: alpha, beta, lowest_eigval_old, eigval_diff

  ! Try to remove a temporary array when call diag
  REAL(DP) :: Htmp(ilanc,ilanc), Hstep(nlanc,nlanc)


  ! allocate vectors and put to zero
  ALLOCATE( q(3,nat), source=0.D0 )
  ALLOCATE( v1(3,nat), source=0.D0)

  !
  ! store the eigenvalue of the previous iteration
  !
  IF( ilanc > 0 )THEN
    lowest_eigval_old = lowest_eigval
  ELSE
    lowest_eigval_old = 0.0_DP
  ENDIF


  !
  ! decide what to do based on which step ilanc we are in
  !
  IF ( ilanc  == 0 ) THEN
     !
     ! initialization of the lanczos: save the original force, and
     ! the initial lanczos vector
     !
     ! write(785,*) 'entering lanc with size:',nlanc, unconvert_length(dlanc), dlanc
     !
     ! store the force of the initial position
     !
     force_old = force(:,:)
     !
     ! normalize initial vector
     !
     v1(:,:) = v_in(:,:) / dnrm2( 3*nat, v_in, 1 )
     !
     ! store this vector in matrix of Lanczos vectors
     !
     Vmat(:,:,1) = v1(:,:)
     !
  ELSEIF (ilanc == 1 ) THEN
     !
     ! First step of the lanczos algorithm:
     !
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
     ! beta is the norm of v1, used for next ilanc step
     !
     beta = dnrm2(3*nat,v1,1)
     v1(:,:) = v1(:,:) / beta
     !
     ! store the vector for future cycles
     !
     Vmat(:,:,2) = v1(:,:)
     H(1,1) = alpha
     H(2,1) = beta
     H(1,2) = beta
     !
     ! there is only one possible eigval in this step: alpha
     ! the corresponding eigvec is unchanged
     !
     lowest_eigval = alpha
     !
     ! check for convergence in this step
     !
     eigval_diff = (lowest_eigval - lowest_eigval_old)/lowest_eigval_old
     ! write(785,*) 1, lowest_eigval_old, lowest_eigval, abs(eigval_diff)
     !
     IF ( abs(lowest_eigval_old) > 0.0_DP ) THEN
!        IF ( ABS(eigval_diff) <= eval_conv_thr ) THEN
        IF ( ilanc .ge. lanczos_min_size .and. ABS(eigval_diff) <= eval_conv_thr ) THEN
           !
           ! lanczos has converged
           ! set max number of iternations to current iteration
           !
           !write(785,*) 'converged step1'
           nlanc = ilanc
           ! increase lanczos counter for last step
           ! lowest_eigvec(:,:) = v_in(:,:)
           !
           ! the displ_vec going out should be: -Vmat(:,:,1), so
           ! put v1 to 0.0, then subtract Vmat(:,:,ilanc) few lines later
           v1(:,:) = 0.0
        ENDIF
     ENDIF
     !
     ! correct v1 so that the move is made from the initial position
     v1(:,:) = v1(:,:) - Vmat(:,:,ilanc)
     !
  ELSEIF (ilanc > 1 .and. ilanc <= nlanc ) THEN
     !
     ! first generate current alpha
     !
     q(:,:) = force(:,:) - force_old(:,:)
     q(:,:) = -q(:,:)/dlanc
     !
     ! alpha = dot( v1, q )
     !
     alpha = ddot(3*nat,vmat(:,:,ilanc),1,q(:,:),1)
     H(ilanc,ilanc) = alpha
     !
     ! then check convergence of the H matrix up to this step
     !
     ALLOCATE( eigvals(ilanc) )
     ! store the H matrix, because its overwritten by eigvecs on diagonalization
     Hstep(:,:) = H(:,:)
     Htmp = H(1:ilanc,1:ilanc)  !%! NS: add this step to remove a warning

     CALL diag(ilanc, Htmp, eigvals, 1 )
     Hstep(1:ilanc,1:ilanc) = Htmp
     !
     ! find the lowest eigenvalue in the vector
     !
     lowest_eigval = eigvals(1)
     id_min = 1
     DO i = 1, ilanc
        IF (eigvals(i) < lowest_eigval ) THEN
           lowest_eigval = eigvals(i)
           id_min = i
        ENDIF
     ENDDO
     !write (*,*) "Debug eigval", lowest_eigval
     !
     ! generate eigenvector in real space, corresponding to lowest eigenvalue
     !
     ! Hstep now stores eigvecs of H
     ! eigvecs in coordinate space are computed as matmul(V, lowest_eigvec_H )
     !
     ! Multiply matrices (V_1 | ... | V_ilanc)*H(min)=eigen(min) using dgemm of lapack
     !
     ! The call to dgemm contains:
     ! (see http://www.math.utah.edu/software/lapack/lapack-blas/dgemm.html)
     ! 'N'    ... do not transpose Vmat
     ! 'N'    ... do not transpose Hstep(:,id_min)
     ! 3*nat  ... rows of Vmat(1:3,1:nat)
     ! 1      ... columns of Hstep
     ! ilanc  ... columns of Vmat, rows of Hstep
     ! 1.0_DP ... alpha for dgemm
     ! Vmat(:,:,1:ilanc) ... Vmat of current step
     ! 3*nat             ... first dimension of Vmat
     ! Hstep(:,id_min)   ... eigenvector with lowest eigenvalue of H
     ! ilanc             ... first dimension of Hstep
     ! 0.0_DP            ... beta of dgemm
     ! lowest_eigvec     ... resulting eigenvector dimensions (1:3,1:nat)
     ! 3*nat             ... first dimension of lowest_eigvec
     !
     CALL dgemm('N','N',3*nat,1,ilanc,1.0_DP,Vmat(:,:,1:ilanc),3*nat,Hstep(:,id_min),ilanc,0.0_DP,lowest_eigvec,3*nat)
     !
     ! check if the obtained eigenvec points in the same direction as the pushing direction
     !
     dir = ddot(3*nat,lowest_eigvec,1, pushdir, 1)
     !write (*,*) "Lanczos::Debug dir:", dir
     IF ( dir < 0.D0 ) THEN
        lowest_eigvec(:,:) = -1.D0*lowest_eigvec(:,:)
     ENDIF
     !
     DEALLOCATE( eigvals )
     !
     ! Check for the convergence of the lanczos eigenvalue
     !
     eigval_diff = (lowest_eigval - lowest_eigval_old)/lowest_eigval_old
     !write (*,*) "Debug eigval:", ilanc, lowest_eigval_old, lowest_eigval, abs(eigval_diff)
     ! write(785,*) ilanc, lowest_eigval_old, lowest_eigval, abs(eigval_diff)
     !
     !IF ( ABS(eigval_diff) <= eval_conv_thr ) THEN
     IF ( ilanc.ge.lanczos_min_size .and. ABS(eigval_diff) <= eval_conv_thr ) THEN
        ! write(*,*) 'converged! in:',ilanc
        !write(785,*) 'converged! in:',ilanc
        !
        ! lanczos has converged
        ! set max number of iternations to current iteration
        nlanc = ilanc
        !
     ENDIF
     !
     ! if lanczos is not yet converged generate new lanczos vector
     !
     IF ( ilanc < nlanc ) THEN
        beta = H(ilanc,ilanc-1)
        !
        ! Generate lanczos vector, v = q - alpha*v1 - beta*v0
        v1(:,:) = q(:,:) - alpha*vmat(:,:,ilanc) - beta*vmat(:,:,ilanc-1)
        !
        ! orthogonalize vectors in accordance with previous ones ...
        DO j = 1, ilanc - 1
           v1(:,:) = v1(:,:) - ddot(3*nat, v1 ,1, Vmat(:,:,j),1)*Vmat(:,:,j)
        ENDDO
        !
        !  do stuff with the new vector
        !
        IF ( dnrm2(3*nat, v1(:,:), 1) < 1.0D-15 ) THEN
           !
           ! new lanczos vector very small, stop (converge)
           !
           ! write(785,*) 'new vector small!'
           !
        ELSE
           !
           ! normalize the new vector
           !
           v1(:,:) = v1(:,:)/dnrm2(3*nat,v1(:,:),1)
           !
           ! save it into Vmat
           Vmat(:,:,ilanc+1) = v1(:,:)
           !
           ! generate new beta
           beta = ddot(3*nat, q, 1, v1, 1)
           !
           ! store value to H
           H(ilanc+1,ilanc ) = beta
           H(ilanc, ilanc+1) = beta
           !
        ENDIF
        !
     END IF
     !
     ! correct v1 so that the move is made from the initial position
     v1(:,:) = v1(:,:) - Vmat(:,:,ilanc)
     !
  ENDIF
  !
  ! write data for next lanczos step
  !
  ! IF( ilanc > nlanc ) THEN
  !    write(785,*) 'moving back'
  !    ! final move back to initial position
  !    v1(:,:) = 0.D0
  !    v1(:,:) = v1(:,:) - Vmat(:,:,nlanc)
  ! ENDIF

  !
  ! Overwrite displ_vec by the next vector displacement
  !
  displ_vec(:,:) = v1(:,:)

  DEALLOCATE( q, v1 )

  ! flush(785)

END SUBROUTINE lanczos

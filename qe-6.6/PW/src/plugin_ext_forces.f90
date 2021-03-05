!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE plugin_ext_forces()
  !----------------------------------------------------------------------------
  !
  !
  USE mp,               ONLY : mp_bcast
  USE mp_images,        ONLY : intra_image_comm
  USE io_global,        ONLY : stdout, ionode, ionode_id
  USE kinds,            ONLY : DP
  !
  USE plugin_flags
  ! modifications start here
  ! we take some stuff from QE
  USE ions_base,     ONLY : nat, tau, if_pos, ityp, atm, amass
  USE cell_base,     ONLY : alat, at
  USE force_mod,     ONLY : force
  USE ener,          ONLY : etot 
  USE relax,         ONLY : epsf, epse
  USE control_flags, ONLY : istep
  USE dynamics_module, ONLY : vel, dt, fire_alpha_init
  USE io_files,      ONLY : prefix,tmp_dir
  !
  IMPLICIT NONE
  ! 
  LOGICAL :: lconv
  !
  ! ARTn convergence flag 
  ! 
  lconv = .false. 
  !
  IF ( ionode ) THEN
     CALL artn(force,etot,epsf,nat,ityp,atm,tau,at,alat,istep,if_pos,vel,dt,fire_alpha_init,lconv,prefix,tmp_dir) 
  ENDIF
  IF ( lconv ) THEN
     WRITE (*,*) "ARTn calculation converged, stopping" 
     STOP 1
  END IF
  
END SUBROUTINE plugin_ext_forces
MODULE artn_params
  !
  ! This module contains all global variables that are used in the ARTn plugin 
  !
  IMPLICIT none
  SAVE
  ! constants 
  INTEGER, PARAMETER ::  DP = selected_real_kind(14,200) ! double precision
  REAL(DP), PARAMETER :: PI     = 3.14159265358979323846_DP ! pi 
  REAL(DP), PARAMETER :: RY2EV =  13.605691930242388_DP ! Ry to eV conversion 
  REAL(DP), PARAMETER :: B2A =  0.529177210903_DP ! bohr to angstrom conversion
  REAL(DP), PARAMETER :: AMU_RY = 911.44424310865645_DP ! calculated from QE using DP 
  INTEGER, PARAMETER :: iunartin = 52  ! fortran file unit for ARTn input file  
  INTEGER, PARAMETER :: iunartout = 53 ! fortran file unit for ARTn output file 
  INTEGER, PARAMETER :: iunsaddle = 54 ! fortran file unit for writing the saddle coords 
  ! control flags
  LOGICAL :: lpush_init ! initial push 
  LOGICAL :: lperp      ! perpendicular relax 
  LOGICAL :: leigen     ! push with lanczos eigenvector
  LOGICAL :: llanczos   ! lanczos algorithm
  ! counters
  INTEGER :: istepperp  ! number of steps in perpendicular relaxation
  INTEGER :: neigenstep ! number of steps made with eigenvector
  INTEGER :: npush      ! number of pushes made
  INTEGER :: nlanc      ! current lanczos iteration
  INTEGER :: nlanciter  ! max number of lanczos iterations
  INTEGER :: if_pos_ct  ! counter used to determine the number of fixed coordinates 
  ! lanczos variables
  REAL(DP) :: lowest_eigval
  !                                           ! 
  ! arrays that are needed by ARTn internally ! 
  !                                           ! 
  REAL(DP), ALLOCATABLE :: push(:,:)     ! initial push vector
  REAL(DP), ALLOCATABLE :: eigenvec(:,:) ! lanczos eigenvector
  !
  REAL(DP) :: etot_init ! store the total energy of the initial state
  !                                               !
  ! arrays that are used by the Lanczos algorithm ! 
  !                                               ! 
  REAL(DP), ALLOCATABLE :: H(:,:) ! tridiagonal matrix
  REAL(DP), ALLOCATABLE :: Vmat(:,:,:) ! matrix containing the laczos vectors
  REAL(DP), ALLOCATABLE :: force_old(:,:) ! force in the previous step 
  !------------------------------------------------------------! 
  ! variables that are read from the input  start here   
  !------------------------------------------------------------!  
  !
  NAMELIST/artn_parameters/ lrelax,npushmin, neigenstepmax, nlanciter_init, push_mode, dist_thr, &
       convcrit_init,convcrit_final,fpara_convcrit, eigval_thr, &
       init_step_size, dlanc, step_size, &
       push_ids,add_const
  ! 
  LOGICAL :: lrelax     ! do we want to ignore artn and do a regular relax calculation
  ! 
  INTEGER :: npushmin              ! number of initial pushes before lanczos start
  INTEGER :: neigenstepmax         ! number of steps made with eigenvector before perp relax 
  INTEGER :: nlanciter_init        ! maximum number of lanczos iterations 
  CHARACTER (LEN = 4) :: push_mode ! type of initial push (all , list or rad) 
  ! convergence criteria
  REAL(DP) :: dist_thr       ! distance threshold for push mode "rad" 
  REAL(DP) :: convcrit_init  ! initial perp force convergence criterion for perp relax
  REAL(DP) :: convcrit_final ! tightened force convergence criterion when near the saddle point
  REAL(DP) :: fpara_convcrit ! parallel force convergence criterion, used to determine when to tighten convcrit_final
  REAL(DP) :: eigval_thr 
  ! step sizes
  REAL(DP) :: init_step_size ! step size of inital push in angstrom 
  REAL(DP) :: step_size      ! step size for a step with the lanczos eigenvector
  REAL(DP) :: current_step_size ! controls the current size of eigenvector step
  REAL(DP) :: dlanc         ! step size in the lanczos algorithm 
  ! arrays related to constraints 
  INTEGER,  ALLOCATABLE :: push_ids(:)    ! IDs of atoms to be pushed 
  REAL(DP), ALLOCATABLE :: add_const(:,:) ! constraints on initial push
  
CONTAINS
  !
  SUBROUTINE initialize_artn(nat,iunartin,iunartout,filnam,filout)
    !
    ! sets defaults, reads input and creates ARTn output file
    ! 
    IMPLICIT none
    INTEGER, INTENT(IN) :: nat,iunartin,iunartout
    CHARACTER (LEN=255), INTENT(IN) :: filnam,filout
    LOGICAL :: file_exists
    INTEGER :: ios
    !
    ! set up defaults for flags and counters 
    !
    lrelax = .false. 
    lpush_init = .true.
    lperp = .false.
    llanczos = .false.
    leigen = .false.
    !
    npush = 0 
    istepperp = 0
    nlanc = 0
    neigenstep = 0
    if_pos_ct = 0 
    !
    lowest_eigval = 0.D0
    !
    ! Defaults for input parameters
    ! 
    npushmin = 3
    neigenstepmax = 1 
    !
    dist_thr = 0.0_DP 
    ! 
    convcrit_init = 1.0d-2
    convcrit_final = 1.0d-3
    fpara_convcrit = 0.5d-2
    eigval_thr = -0.01_DP ! in Ry/bohr^2 corresponds to 0.5 eV/Angs^2  
    ! 
    init_step_size = 0.3
    step_size = 0.2
    !
    push_mode = 'all' 
    !
    dlanc = 1.D-2 
    nlanciter_init = 16
    !
    ! Allocate the arrays
    !
    IF ( .not. ALLOCATED(add_const)) ALLOCATE(add_const(4,nat), source = 0.D0)
    IF ( .not. ALLOCATED(push_ids)) ALLOCATE(push_ids(nat), source = 0)
    IF ( .not. ALLOCATED(push)) ALLOCATE(push(3,nat), source = 0.D0)
    IF ( .not. ALLOCATED(eigenvec)) ALLOCATE(eigenvec(3,nat), source = 0.D0)
    INQUIRE( file = filnam, exist = file_exists )
    ! 
    IF ( file_exists ) THEN
       ! read the ARTn input file  
       OPEN( UNIT = iunartin, FILE = filnam, FORM = 'formatted', STATUS = 'unknown', IOSTAT = ios)
       READ( NML = artn_parameters, UNIT = iunartin)
       CLOSE ( UNIT = iunartin, STATUS = 'KEEP')
       ! open the ARTn output file 
       OPEN ( UNIT = iunartout, FILE = filout, FORM = 'formatted', STATUS = 'unknown', IOSTAT = ios )
       WRITE (iunartout,'(5X, "--------------------------------------------------")')
       WRITE (iunartout,'(5X, "                ARTn plugin                       ")')
       WRITE (iunartout,'(5X, "--------------------------------------------------")')
       WRITE (iunartout,'(5X, " "                                                 )')
       WRITE (iunartout,'(5X, "               INPUT PARAMETERS                   ")')
       WRITE (iunartout,'(5X, "--------------------------------------------------")')
       WRITE (iunartout,'(5X, "Push and perpendicular relax:")')
       WRITE (iunartout,'(5X, "--------------------------------------------------")') 
       WRITE (iunartout,'(15X,"npushmin       = ", I6)') npushmin 
       WRITE (iunartout,'(15X,"convcrit_init  = ", F6.3)') convcrit_init 
       WRITE (iunartout,'(15X,"convcrit_final = ", F6.3)') convcrit_final
       WRITE (iunartout,'(15X,"fpara_convcrit = ", F6.3)') fpara_convcrit 
       WRITE (iunartout,'(15X,"eigval_thr     = ", F6.3)') eigval_thr  
       WRITE (iunartout,'(15X,"init_step_size = ", F6.1)') init_step_size 
       WRITE (iunartout,'(15X,"step_size      = ", F6.1)') step_size
       WRITE (iunartout,'(15X,"push_mode      = ", A6)') push_mode
       WRITE (iunartout,'(5X, "--------------------------------------------------")')
       WRITE (iunartout,'(5X, "Lanczos algorithm:")' )
       WRITE (iunartout,'(5X, "--------------------------------------------------")')
       WRITE (iunartout,'(15X, "nlanciter_init = ", I6)') nlanciter_init
       WRITE (iunartout,'(15X, "dlanc          = ", F6.3)') dlanc
       WRITE (iunartout,'(5X, "--------------------------------------------------")')
       WRITE (iunartout,*) " "
       WRITE (iunartout,*) " "
       WRITE (iunartout,'(5X,"istep",4X,"ART_step",12X,"Etot",9X," Ftot ",5X," Fperp ",5X," Fpara ",6X,"eigval")') 
       WRITE (iunartout,'(34X, "[Ry]",9X,"-----------[Ry/a.u.]----------",6X,"Ry/a.u.^2")')
       CLOSE ( UNIT = iunartout, STATUS = 'KEEP')
    ELSE
       WRITE(*,*) "ARTn: Input file does not exist!"
       lpush_init = .false.
       lrelax = .true.
       RETURN 
    ENDIF
    ! set initial number of lanczos iterations 
    nlanciter = nlanciter_init
    ! convert push/lanczos/eigenvec step size to bohr (because force units are in Ry/bohr) 
    step_size = step_size/B2A
    init_step_size = init_step_size/B2A
    dlanc = dlanc/B2A
    dist_thr = dist_thr/B2A
    ! 
  END SUBROUTINE initialize_artn
  !
  SUBROUTINE initialize_lanczos(nlanciter,nat)
    IMPLICIT none
    INTEGER, INTENT(IN) :: nlanciter,nat 
    ! allocate the matrices 
    IF ( .not. ALLOCATED(H)) ALLOCATE( H(1:nlanciter,1:nlanciter), source = 0.D0 )
    IF ( .not. ALLOCATED(Vmat)) ALLOCATE( Vmat(3,nat,1:nlanciter), source = 0.D0 )
    IF ( .not. ALLOCATED(force_old) ) ALLOCATE( force_old(3,nat), source = 0.D0 )
  END SUBROUTINE initialize_lanczos
  ! 
END MODULE artn_params

SUBROUTINE center ( vec, nat)
  USE artn_params, ONLY: DP
  !
  ! takes as input a vector of size (3,nat) and centers it
  !
  IMPLICIT none
  INTEGER, INTENT(IN) :: nat
  REAL(DP), INTENT(INOUT) :: vec(3,nat)
  INTEGER :: na
  REAL(DP) :: delta(3)
  !
  delta(:) = 0.D0
  DO na = 1,nat
     delta(:) = delta(:) + vec(:,na)
  ENDDO
  !
  delta(:) = delta(:)/dble(nat)
  !
  FORALL ( na = 1:nat) vec(:,na) = vec(:,na) - delta(:)


END SUBROUTINE center
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

SUBROUTINE lanczos( nat, force, vel, alpha_init, dt, &
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
  REAL(DP), INTENT(IN) :: alpha_init, dt
  CHARACTER(LEN=255), INTENT(IN) :: tmpdir, prfx
  REAL(DP), DIMENSION(3,nat), INTENT(INOUT) :: vel 
  REAL(DP), DIMENSION(3,nat), INTENT(INOUT) :: force
  REAL(DP), DIMENSION(3,nat), INTENT(INOUT) :: lowest_eigvec
  REAL(DP), INTENT(INOUT) :: lowest_eigval
  INTEGER, INTENT(INOUT) :: nlanciter
  INTEGER, INTENT(INOUT) :: nlanc
  ! 
  INTEGER :: i, j, io, id_min
  REAL(DP), PARAMETER :: eigval_thr = 1.0D-2
  REAL(DP), ALLOCATABLE :: v1(:,:), q(:,:), eigvals(:)
  REAL(DP), ALLOCATABLE ::  Hstep(:,:)
  REAL(DP) :: dir
  REAL(DP), EXTERNAL :: ran3,dnrm2,ddot
  REAL(DP) :: alpha, beta, lowest_eigval_old, eigvec_diff, largest_eigvec_diff, eigval_diff
  ! allocate vectors and put to zero
  ALLOCATE( q(3,nat), source=0.D0 )
  ALLOCATE( v1(3,nat), source=0.D0)
  ! allocate matrices and put to zero
  ALLOCATE( Hstep(1:nlanciter,1:nlanciter), source=0.D0 )
  ! 
  ! store the eigenvalue of the previous iteration
  lowest_eigval_old = lowest_eigval
  !
  IF ( nlanc  == 0 ) THEN
     ! store the force of the initial position 
     force_old = force(:,:) 
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
     ! correct v1 so that the move is made from the initial position 
     v1(:,:) = v1(:,:) - Vmat(:,:,nlanc -1) 
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
     ! Hstep now stores eigvecs of H
     ! eigvecs in coordinate space are computed as matmul(V, lowest_eigvec_H )
     !
     ! Multiply matrices (V_1 | ... | V_nlanc)*H(min)=eigen(min) using dgemm of lapack  
     !  
     !
     ! The call to dgemm contains:
     ! (see http://www.math.utah.edu/software/lapack/lapack-blas/dgemm.html)
     ! 'N'    ... do not transpose Vmat
     ! 'N'    ... do not transpose Hstep(:,id_min)
     ! 3*nat  ... rows of Vmat(1:3,1:nat)
     ! 1      ... columns of Hstep
     ! nlanc  ... columns of Vmat, rows of Hstep
     ! 1.0_DP ... alpha for dgemm
     ! Vmat(:,:,1:nlanc) ... Vmat of current step
     ! 3*nat             ... first dimension of Vmat
     ! Hstep(:,id_min)   ... eigenvector with lowest eigenvalue of H 
     ! nlanc             ... first dimension of Hstep 
     ! 0.0_DP            ... beta of dgemm
     ! lowest_eigvec     ... resulting eigenvector dimensions (1:3,1:nat)
     ! 3*nat             ... first dimension of lowest_eigvec
     ! 
     CALL dgemm('N','N',3*nat,1,nlanc,1.0_DP,Vmat(:,:,1:nlanc),3*nat,Hstep(:,id_min),nlanc,0.0_DP,lowest_eigvec,3*nat)
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
    
    ! correct v1 so that the move is made from the initial position
    v1(:,:) = v1(:,:) - Vmat(:,:,nlanc-1)
    !    
 ENDIF
  ! store the force
  ! force_old(:,:) = force(:,:)
  !
  ! write data for next lanczos step
  !
 IF( nlanc > nlanciter ) THEN
    ! final move back to initial position 
    v1(:,:) = 0.D0
    v1(:,:) = v1(:,:) - Vmat(:,:,nlanciter)  
 ENDIF
 !
 ! write data for move
 !
 CALL move_mode( nat, dlanc, v1, force, &
      vel, alpha_init, dt, &
      0, pushdir, 'lanc', prfx, tmpdir)
 !
 ! deallocate the matrices used only in the iteration 
 !
 DEALLOCATE( q, v1 )
 DEALLOCATE(Hstep)

END SUBROUTINE lanczos

SUBROUTINE displacement_validation( atom_id, atom_const, push, lvalid)
  !
  ! subroutine that checks if the initial_displacement is within given parameters
  !
  USE artn_params, ONLY : DP, PI 
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: atom_id
  REAL(DP), INTENT(IN) :: atom_const(4)
  REAL(DP), INTENT(INOUT) :: push(3)
  REAL(DP), EXTERNAL :: ddot, dnrm2
  LOGICAL,         INTENT(INOUT) :: lvalid
  !
  ! Local variables
  REAL(DP)               :: cone_angle, displacement_angle
  REAL(DP)               :: dot_prod, displacement_norm, cone_dir_norm
  REAL(DP), DIMENSION(3) :: cone_dir,displacement
  !
  !
  !
  write (*,*) "ARTn: Called displacement validation: with cone_dir:",atom_const(1:3), "current push:", push(:)
  cone_dir = atom_const(1:3)
  cone_angle = atom_const(4)
  !
  displacement(:)         = push(:)
  !
  displacement_norm       = dnrm2(3,displacement,1)
  cone_dir_norm           = dnrm2(3,cone_dir,1)
  !
  dot_prod                = ddot( 3, cone_dir, 1, displacement, 1 ) / ( cone_dir_norm * displacement_norm )
  displacement_angle      = ACOS( dot_prod ) *180.0_DP / PI
  lvalid                  = ( displacement_angle < cone_angle )
  write (*,*) "Finished displacement validation",lvalid
  !
  IF ( cone_angle == 0.0_DP) THEN
     lvalid = .TRUE.
     !
     ! TODO: why is the direction multiplied by 0.1? seems kind of random ...
     !
     push(:) = cone_dir(:)
  ENDIF
  !
END SUBROUTINE displacement_validation

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

SUBROUTINE invmat3x3(mat,inv)
  USE artn_params, ONLY : DP 
  IMPLICIT none
  !
  ! Subroutine that calculates the inverse of a 3x3 matrix 
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

SUBROUTINE push_init (nat, tau, at, alat, idum, push_ids, dist_thr, add_const, init_step_size, push, mode)
  !
  ! subroutine that generates the initial push; options are specified by mode: 
  !           (1) 'all' generates a push on all atoms 
  !           (2) 'list' generates a push on a list of atoms
  !           (3) 'rad' generates a push on a list of atoms and all atoms within dist_thr 
  ! the user should supply: number and list of atoms to push; and add_constraints on these atoms
  !
  USE artn_params, ONLY : DP
  IMPLICIT none
  INTEGER, INTENT(IN) :: nat,idum
  INTEGER :: na, ia 
  INTEGER, INTENT(IN) :: push_ids(nat)
  REAL(DP), INTENT(IN) :: alat, dist_thr, init_step_size
  REAL(DP), INTENT(IN) :: tau(3,nat), at(3,3), add_const(4,nat)
  CHARACTER(LEN=4), INTENT(IN) :: mode
  REAL(DP), INTENT(OUT) :: push(3,nat)
  REAL(DP), EXTERNAL :: ran3, dnrm2 
  REAL(DP) :: dr2, pushat(3)
  REAL(DP) :: dist(3), tau0(3) 
  LOGICAL :: lvalid
  INTEGER :: atom_displaced(nat)
  !
  push(:,:) = 0.d0
  atom_displaced(:) = 0
  lvalid = .false.
  !
  !  read the list of pushed atoms
  !
  IF ( mode == 'all' ) THEN
     ! displace all atoms 
     atom_displaced(:) = 1
  ELSE IF ( mode == 'list' ) THEN 
     ! displace only atoms in list 
     DO na=1,nat
        IF (ANY(push_ids == na)) THEN
           atom_displaced(na) = 1
        ENDIF
     ENDDO
  ELSE IF ( mode == 'rad' ) THEN 
     ! displace only atoms in list and all atoms within chosen a cutoff radius ...
     DO na=1,nat
        IF (ANY(push_ids == na)) THEN
           atom_displaced(na) = 1           
           !
           tau0 = tau(:,na)*alat
           DO ia = 1,nat
              IF (ia /= na) THEN 
                 dist(:) = tau(:,ia)*alat - tau0(:)
                 
                 CALL pbc( dist, at, alat)
                 IF ( dnrm2(3,dist,1) <= dist_thr ) THEN
                    ! found an atom within dist_thr 
                    atom_displaced(ia) = 1
                 ENDIF
              ENDIF
           ENDDO
        ENDIF
     ENDDO
              
  ENDIF
  !
  !
  DO na=1,nat
     IF (atom_displaced(na) == 1 ) THEN
        DO
           push(:,na) = (/0.5_DP - ran3(idum),0.5_DP - ran3(idum),0.5_DP - ran3(idum)/)
           dr2 = push(1,na)**2 + push(2,na)**2 + push(3,na)**2
           ! check if the atom is constrained
           IF (ANY(ABS(add_const(:,na)) > 0.D0)) THEN
              ! check if the displacement is within the chosen constraint
              CALL displacement_validation(na,add_const(:,na),push(:,na),lvalid )
              IF ( .not. lvalid ) CYCLE
           ENDIF
           IF ( dr2 < 0.25_DP ) EXIT
        ENDDO
     ENDIF
  ENDDO
  ! if all atoms are pushed center the push vector to avoid translational motion 
  IF ( mode == 'all')  CALL center(push(:,:), nat)
  ! 
  !
  ! normalize so that maxval of the vector is 1.0
  !
  push(:,:) = push(:,:)/MAXVAL(ABS(push(:,:)))
  ! scale initial push vector according to step size 
  push(:,:) = init_step_size*push(:,:)

END SUBROUTINE push_init

SUBROUTINE sum_force(force,nat,force_tot)
  !
  ! subroutine that sums the forces on all atoms and returns the total force
  !
  USE artn_params, ONLY : DP
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nat
  REAL(DP), INTENT(IN) :: force(3,nat)
  REAL(DP), INTENT(OUT) :: force_tot
  INTEGER :: na
  force_tot = 0.0
  DO na = 1, nat
     force_tot = force_tot + force(1,na)**2 + force(2,na)**2 + force(3,na)**2
  ENDDO
  force_tot = SQRT(force_tot)
END SUBROUTINE sum_force

SUBROUTINE perpforce(force,if_pos,push,fpara,nat)
  !
  ! subroutine that subtracts parallel components to push from force
  !
  USE artn_params, ONLY : DP
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nat
  REAL(DP), INTENT(IN)  :: push(3,nat)
  REAL(DP), INTENT(INOUT)  :: force(3,nat)
  REAL(DP), INTENT(OUT)  :: fpara(3,nat)
  INTEGER, INTENT(IN) :: if_pos(3,nat)
  REAL(DP) :: push_norm(3,nat)
  REAL(DP), EXTERNAL :: ddot,dnrm2
  INTEGER  :: na
  ! calculate components parallel to the push
  fpara(:,:) = ddot(3*nat,force(:,:),1,push(:,:),1) / ddot(3*nat,push(:,:),1,push(:,:),1) * push(:,:)
  ! subtract them
  force(:,:) = force(:,:) - fpara(:,:)
  ! apply constraints
  IF ( ANY(if_pos(:,:) == 0)  ) force(:,:) = force(:,:)*if_pos(:,:) 

END SUBROUTINE perpforce
SUBROUTINE move_mode(nat, dlanc, v1, force, &
                     vel, alpha_init, dt, &
                     istepperp, push, &
                     mode, prfx, tmpdir )
  !
  ! translate specified move to appropriate force and set FIRE parameters accordingly  
  !
  USE artn_params, ONLY: DP, AMU_RY 
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nat
  REAL(DP), INTENT(IN) :: dlanc  
  REAL(DP), DIMENSION(3,nat), INTENT(IN) :: v1
  REAL(DP), DIMENSION(3,nat), INTENT(INOUT) :: force
  REAL(DP), DIMENSION(3,nat), INTENT(INOUT) :: vel
  REAL(DP), INTENT(IN) :: alpha_init, dt
  INTEGER, INTENT(IN) :: istepperp
  REAL(DP), DIMENSION(3,nat), INTENT(IN) :: push
  CHARACTER(LEN=4), INTENT(IN) :: mode
  CHARACTER(LEN=255), INTENT(IN) :: tmpdir, prfx
  ! 
  REAL(DP), EXTERNAL :: ddot,dnrm2
  ! variables read from the FIRE minimization algorithm
  INTEGER :: nsteppos
  INTEGER :: ios
  REAL(DP) :: dt_curr, alpha, etot
  LOGICAL :: file_exists
  CHARACTER(len=256) :: filnam
  ! etot is set to zero to ensure that a value is written down in the initial push 
  etot = 0.D0
  !
  ! do things depending on mode of the move
  ! NOTE force units of Ry/a.u. are assumed ... 
  !
  filnam = trim(tmpdir) // '/' // trim(prfx) // '.' //'fire'
  INQUIRE( file = filnam, exist = file_exists )
  OPEN( unit = 4, file = filnam, form = 'formatted', status = 'unknown', iostat = ios)
  !
  IF (file_exists ) THEN
     ! if file exists read the data, otherwise just close it 
     READ( UNIT = 4, FMT = * ) etot, nsteppos, dt_curr, alpha
     CLOSE( UNIT = 4, STATUS = 'KEEP' )
  ELSE
     CLOSE( UNIT = 4, STATUS = 'DELETE')
  ENDIF

  SELECT CASE( TRIM(mode) )

  CASE( 'perp' )
     !
     IF( istepperp .eq. 0 ) THEN
        ! for the first step forget previous velocity (prevent P < 0)
        etot = 0.D0
        vel(:,:) = 0.D0
        alpha = alpha_init
        dt_curr = dt
     ELSE
        ! subtract the components that are parallel
        vel(:,:) = vel(:,:) - ddot(3*nat,vel, 1, push, 1 )*push(:,:)/ddot(3*nat,push(:,:),1, push(:,:),1)
     ENDIF
        !
  CASE( 'lanc' )
     !
     ! set the velocity and acceleration and alpha of previous step to move correctly
     !
     etot = 0.D0
     vel(:,:) = 0.D0
     dt_curr = dt
     alpha = 0.D0
     nsteppos = 0
     ! the step performed should be like this now translate it into the correct force
     force(:,:) = v1(:,:)*dlanc*amu_ry/dt_curr**2
     !
  CASE( 'eign' )
     !
     etot = 0.D0
     vel(:,:) = 0.D0
     alpha = 0.0_DP
     dt_curr = dt
     nsteppos = 0
     force(:,:) = force(:,:)*amu_ry/dt_curr**2
     !
  CASE default
     write(*,*) 'Problem with move_mode!'
  END SELECT
  !
  ! write the FIRE parameters to its scratch file
  ! 
  OPEN( unit = 4, file = filnam, form = 'formatted', status = 'unknown', iostat = ios)
  WRITE( UNIT = 4, FMT = * ) etot, nsteppos, dt_curr, alpha
  !
  CLOSE( UNIT = 4, STATUS = 'KEEP' )
  !
  !
 

END SUBROUTINE move_mode




SUBROUTINE write_report(etot, force, lowest_eigval, step, if_pos, istep, nat,iunartout)
  !
  ! a subroutine that writes a report of the current step to the output file  
  !
  USE artn_params, ONLY: DP, push 
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nat, istep, iunartout
  INTEGER, INTENT(IN) :: if_pos(3,nat)
  REAL(DP), INTENT(IN) :: force(3,nat), etot, lowest_eigval
  REAL(DP) :: fpara(3,nat), fperp(3,nat)
  REAL(DP) :: force_tot, fperp_tot, fpara_tot
  CHARACTER(LEN=4), INTENT(IN) :: step 
  REAL(DP), EXTERNAL :: ddot
  !
  CALL sum_force(force,nat,force_tot)
  fperp(:,:) = force(:,:)
  CALL perpforce(fperp,if_pos,push,fpara,nat)
  CALL sum_force(fperp,nat,fperp_tot)
  fpara_tot = ddot(3*nat,force,1,push,1)
  !  write report 
  WRITE (iunartout,'(5X,I4,7X,A4,9X, F12.6, 5X, F7.4,5X, F7.4, 5X, F7.4, 5X, F7.4)') &
       & istep, step, etot, force_tot,fperp_tot,fpara_tot, lowest_eigval
END SUBROUTINE write_report
SUBROUTINE write_struct(alat, at, nat, tau, atm, ityp, force, fscale, ounit, form, fname)
  !
  ! A subroutine that writes the structure to a file (based on xsf_struct of QE)  
  !
  USE artn_params, ONLY: DP,B2A 
  IMPLICIT NONE 
  INTEGER, INTENT(IN) :: nat             ! number of atoms 
  INTEGER, INTENT(IN) :: ityp(nat)       ! atom type
  CHARACTER(LEN=3), INTENT(IN) :: atm(*) ! contains information on atomic types 
  INTEGER, INTENT(IN) :: ounit           ! output fortran unit 
  REAL(DP), INTENT(IN) :: alat           ! alat of QE   
  REAL(DP), INTENT(IN) :: tau(3,nat)     ! atomic positions 
  REAL(DP), INTENT(IN) :: at(3,3)        ! lattice parameters in alat units 
  REAL(DP), INTENT(IN) :: force(3,nat)   ! forces
  REAL(DP), INTENT(IN) :: fscale         ! factor for scaling the force 
  CHARACTER(LEN=3), INTENT(IN) :: form   ! format of the structure file (default xsf)
  CHARACTER(LEN=255), INTENT(IN) :: fname  ! file name 
  ! 
  !
  INTEGER :: i, j, na, ios 
  REAL(DP) :: at_angs (3,3)
  OPEN ( UNIT = ounit, FILE = fname, FORM = 'formatted',  STATUS = 'unknown', IOSTAT = ios )
  IF ( form == 'xsf' ) THEN 
     
     DO i=1,3
        DO j=1,3
           at_angs(j,i) = at(j,i)*alat*B2A
        ENDDO
     ENDDO

     WRITE(ounit,*) 'CRYSTAL'
     WRITE(ounit,*) 'PRIMVEC'
     WRITE(ounit,'(2(3F15.9/),3f15.9)') at_angs
     WRITE(ounit,*) 'PRIMCOORD'
     WRITE(ounit,*) nat, 1
     
     DO na=1,nat
        ! positions are in Angstroms
        WRITE(ounit,'(a3,3x,6f15.9)') atm(ityp(na)), &
             tau(1,na)*alat*B2A, &
             tau(2,na)*alat*B2A, &
             tau(3,na)*alat*B2A, &
             force(1,na)*fscale, &
             force(2,na)*fscale, &
             force(3,na)*fscale
     ENDDO
  ELSE
     WRITE (ounit,*) "Specified structure format not supported"
  ENDIF
  CLOSE (UNIT = ounit , STATUS = 'KEEP')
END SUBROUTINE write_struct
! 
! 
! Main ARTn plugin subroutine:
!        modifies the input force to perform the ARTn algorithm 
!----------------------------------------------------------------------------
SUBROUTINE artn(force,etot,forc_conv_thr_qe,nat,ityp,atm,tau,at,alat,istep,if_pos,vel,dt,fire_alpha_init,lconv,prefix,tmp_dir)
  !----------------------------------------------------------------------------
  !
  ! artn_params for variables and counters that need to be stored   
  ! DEFINED IN: artn_params_mod.f90
  ! 
  USE artn_params, ONLY: DP, RY2EV,B2A, iunartin, iunartout, iunsaddle, &
       lrelax,lpush_init,lperp,leigen,llanczos, &
       istepperp, neigenstep, npush, nlanc, nlanciter, if_pos_ct, &
       lowest_eigval, etot_init, &
       npushmin, neigenstepmax, nlanciter_init, push_mode, dist_thr, convcrit_init, convcrit_final, &
       fpara_convcrit, eigval_thr, init_step_size, current_step_size, dlanc, step_size, &
       push_ids,add_const, push, eigenvec, initialize_artn
  ! 
  IMPLICIT NONE
  REAL(DP), INTENT(INOUT) :: force(3,nat)     ! force calculated by the engine
  REAL(DP), INTENT(INOUT) :: vel(3,nat)       ! velocity of previous FIRE step
  REAL(DP), INTENT(IN) ::    etot             ! total energy in current step
  REAL(DP), INTENT(IN) ::    forc_conv_thr_qe ! force convergence threshold of the engine
  REAL(DP), INTENT(IN) ::    dt               ! default time step in FIRE  
  REAL(DP), INTENT(IN) ::    fire_alpha_init  ! initial value of alpha in FIRE 
  REAL(DP), INTENT(IN) ::    alat             ! lattice parameter of QE
  REAL(DP), INTENT(IN) ::    tau(3,nat)       ! atomic positions (needed for output only)
  REAL(DP), INTENT(IN) ::    at(3,nat)        ! lattice parameters in alat units 
  INTEGER,  INTENT(IN) ::    nat              ! number of atoms
  INTEGER,  INTENT(IN) ::    ityp(nat)        ! atom types
  INTEGER,  INTENT(IN) ::    istep            ! current step
  INTEGER,  INTENT(IN) ::    if_pos(3,nat)    ! coordinates fixed by engine 
  CHARACTER(LEN=3),   INTENT(IN) :: atm(*)    ! name of atom corresponding to ityp
  CHARACTER(LEN=255), INTENT(IN) :: tmp_dir   ! scratch directory of engine 
  CHARACTER(LEN=255), INTENT(IN) :: prefix    ! prefix for scratch files of engine 
  LOGICAL, INTENT(OUT) :: lconv               ! flag for controlling convergence 
  REAL(DP), EXTERNAL :: ran3, dnrm2, ddot     ! lapack functions 
  INTEGER :: na, icoor, idum                  ! integers for loops 
  !
  REAL(DP)  :: force_in(3,nat)                ! stores non-modified force 
  REAL(DP)  :: fpara(3,nat)                   ! force parallel to push/eigenvec   
  REAL(DP)  :: v_in(3,nat)                    ! input vector for lanczos 
  REAL(DP)  :: fpara_tot                      ! total force in parallel direction 
  INTEGER   :: ios                            ! file IOSTAT  
  CHARACTER( LEN=255) :: filin, filout, sadfname, initpfname, eigenfname
  !
  ! The ARTn algorithm proceeds as follows:
  ! (1) push atoms in the direction specified by user & relax in the perpendicular direction;
  ! (2) use the lanczos algorithm calculate the lowest eigenvalue/eigenvec
  ! (3) a negative eigenvalue, update push direction otherwise push again
  ! (4) follow the lanczos direction twoard the saddle point
  ! 
  ! 
  ! flag that controls convergence
  !
  lconv = .false.
  !
  ! store original force
  force_in(:,:) = force(:,:)
  !
  ! fpara_tot is used to scale the magnitude of the eigenvector 
  ! 
  fpara_tot = 0.D0
  !
  filin = 'artn.in'
  filout = 'artn.out'
  sadfname = 'saddle.xsf'
  initpfname = 'initp.xsf'
  eigenfname = 'latest_eigenvec.xsf'
  !
  ! initialize artn 
  !  
  IF (istep == 0 ) THEN
     ! read the input parameters 
     CALL initialize_artn(nat,iunartin,iunartout,filin,filout)
     ! store the total energy of the initial state
     etot_init = etot 
  ENDIF
  IF ( lrelax ) RETURN 
  ! 
  ! Open the output file for writing   
  ! 
  OPEN ( UNIT = iunartout, FILE = filout, FORM = 'formatted', ACCESS = 'append', STATUS = 'unknown', IOSTAT = ios )
  ! 
  ! start initial push
  !
  IF ( lpush_init ) THEN
     !
     ! initial push 
     !
     CALL push_init(nat, tau, at, alat, idum, push_ids, dist_thr, add_const, init_step_size, push ,push_mode)
     ! set up the flags (we did an initial push, now we need to relax perpendiculary) 
     lpush_init = .false.
     lperp = .true.
     ! the push counter (controls if we should call lanczos or keep pushing)
     npush = npush + 1
     ! 
     ! modify the force to be equal to the push
     !
     force(:,:) =  push(:,:)
     CALL move_mode( nat, dlanc, v_in, force, &
          vel, fire_alpha_init, dt,  &
          istepperp, push, 'eign', prefix, tmp_dir)
     !
     CALL write_report(etot,force_in, lowest_eigval, 'push' , if_pos, istep, nat,  iunartout)
     !
     CALL write_struct(alat, at, nat, tau, atm, ityp, force, 1.0_DP, 556, 'xsf', initpfname)
     ! 
  ELSE IF ( lperp ) THEN
     !
     !  subtract parrallel components to push from force
     !
     CALL perpforce(force,if_pos,push,fpara,nat)
     !
     CALL move_mode( nat,  dlanc, v_in, force, &
          vel, fire_alpha_init, dt,  &
          istepperp, push, 'perp', prefix, tmp_dir)
     istepperp = istepperp + 1
     !
     IF (( MAXVAL( ABS(force) )) < convcrit_init ) THEN
        !
        ! we reached convergence in the perpendicular direction of push
        !
        lperp = .false.
        istepperp = 0
        IF ( npush < npushmin ) THEN
           ! continue pushing in the specified direction
           force(:,:) =  push(:,:)
           write (*,*) "ARTn push:", force(:,:)
           CALL move_mode( nat, dlanc, v_in, force, &
                vel, fire_alpha_init, dt,  &
                istepperp, push, 'eign', prefix, tmp_dir)
           npush = npush + 1
           lperp = .true.
           CALL write_report(etot,force_in, lowest_eigval, 'push' , if_pos, istep, nat,  iunartout)
        ELSE IF ( npush >= npushmin  ) THEN
           ! regenerate force & start lanczos
           force(:,:) = force_in(:,:)
           llanczos = .true.
        END IF
     ELSE
        CALL write_report(etot,force_in, lowest_eigval, 'perp' , if_pos, istep, nat,  iunartout)       
     END IF
     ! leigen is always .true. after we obtain a good eigenvector
     ! except during lanczos iterations 
  ELSE IF ( leigen .and. .not. lperp ) THEN
     !
     ! if we have a good lanczos eigenvector use it
     !
     fpara_tot = ddot(3*nat,force(:,:),1,eigenvec(:,:),1)
     fpara(:,:) = fpara_tot*eigenvec(:,:)
     ! check maxmimu value of parallel force  
     IF (MAXVAL(ABS(fpara)) <= fpara_convcrit) THEN
        ! tighten perpendicular convergance criterion
        convcrit_init = convcrit_final  
     END IF
     ! rescale the eigenvector according to the current force in the parallel direction
     ! see Cances_JCP130: some improvements of the ART technique doi:10.1063/1.3088532
     ! 
     ! 0.13 is taken from ARTn, 0.5 eV/Angs^2 corresponds roughly to 0.01 Ry/Bohr^2
     ! 
     current_step_size = MIN(step_size,ABS(fpara_tot)/MAX(ABS(lowest_eigval),0.01_DP))
     ! 
     eigenvec(:,:) = -SIGN(1.0D0,fpara_tot)*eigenvec(:,:)*current_step_size 
     !  
     force(:,:) = eigenvec(:,:)
     CALL move_mode( nat, dlanc, v_in, force, &
          vel, fire_alpha_init, dt,  &
          istepperp, push, 'eign', prefix, tmp_dir)
     ! update eigenstep counter 
     neigenstep = neigenstep + 1
     
     CALL write_report(etot,force_in, lowest_eigval, 'eign' , if_pos, istep, nat,  iunartout)
     ! count the number of steps made with the eigenvector
     CALL write_struct(alat, at, nat, tau, atm, ityp, force, 0.1_DP, 556, 'xsf', eigenfname)
     ! 
     IF ( neigenstep == neigenstepmax  ) THEN
        ! do a perpendicular relax
        lperp = .true.
        istepperp = 0
        ! return to initial number of lanczos steps
        nlanc = 0
     ENDIF
  END IF
  !
  ! check if we should perform the lanczos algorithm
  !
  IF ( llanczos ) THEN 
     !
     CALL write_report(etot,force_in, lowest_eigval, 'lanc' , if_pos, istep, nat,  iunartout)
     IF (nlanc == 0 ) THEN
        IF ( .not. leigen ) THEN
           !
           ! add_const is set to zero so no constraints are put on the atoms
           ! 
           ! add_const(:,:) = 0.D0 
           ! CALL push_init(nat,idum,push_ids,add_const,1.D0,v_in,'all ')
           ! the push vector is already centered within the push_init subroutine 
           ! note push vector is normalized inside lanczos
           v_in(:,:) = push(:,:)
           ! input a different vector for testing purposes
           !DO na=1,nat
           !   ! do a translation 
           !   v_in(:,na) = push(:,5) 
           !ENDDO
        ELSE
           ! rescale the lanczos eigenvec back to original size
           v_in(:,:) = push(:,:)/current_step_size 
           ! reset the eigenvalue flag to continue lanczos
           leigen = .false.
        ENDIF
     ENDIF
     !
     ! apply constraints from the engine (QE)  
     ! 
     IF ( ANY(if_pos(:,:) == 0) ) THEN
        DO na=1,nat
           DO icoor=1,3
              IF (if_pos(icoor,na) == 1 ) if_pos_ct = if_pos_ct + 1
           ENDDO
        END DO
        IF ( if_pos_ct < nlanciter .and. if_pos_ct /= 0 ) nlanciter = if_pos_ct
        v_in(:,:) = v_in(:,:)*if_pos(:,:)
        force(:,:) = force(:,:)*if_pos(:,:)
     ENDIF
     !
     CALL lanczos( nat, force, vel, fire_alpha_init, dt, &
          v_in, dlanc, nlanciter, nlanc, lowest_eigval,  eigenvec, push, prefix, tmp_dir)
     !
     ! when lanczos converges, nlanciter = number of steps it took to converge,
     ! and nlanc = nlanciter + 1
     !
     IF ( nlanc > nlanciter ) THEN
        !
        ! max number of lanczos steps exceeded ; reset counters
        !
        nlanciter = nlanciter_init
        nlanc = 0
        !
        llanczos = .false.
        !
        ! check the eigenvalue, if it's lower than threshold use the eigenvector, otherwise push again ...
        !
        IF ( lowest_eigval < eigval_thr ) THEN
           ! set push to be equal to the eigenvec check eigenvec here eventually ...
           push(:,:) = eigenvec(:,:)
           ! 
           leigen = .true.
           neigenstep = 0
           lperp = .false.
           istepperp = 0
        ELSE
           !
           leigen = .false.
           lperp = .true.
           istepperp =  0
           npush = npush - 1
        ENDIF
     ENDIF
  ENDIF
  !
  ! check for convergence of total forces 
  !
  IF (MAXVAL(ABS(force_in*if_pos)) < convcrit_final .AND. istep /= 0 ) THEN
     force(:,:) = force_in(:,:)
     lconv = .true.
     CALL write_struct(alat, at, nat, tau, atm, ityp, force, 1.0_DP, 556, 'xsf', sadfname)
     WRITE (iunartout,'(5X, "--------------------------------------------------")')
     WRITE (iunartout,'(5X, "          *** ARTn converged to saddle ***        ")')
     WRITE (iunartout,'(5X, "--------------------------------------------------")')
     WRITE (iunartout,'(15X,"E_final - E_initial =", F12.5," eV")') (etot - etot_init)*RY2EV
  ENDIF
  CLOSE (UNIT = iunartout, STATUS = 'KEEP') 
END SUBROUTINE artn 

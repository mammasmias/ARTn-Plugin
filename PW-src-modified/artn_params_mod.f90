MODULE artn_params
  !
  ! This module contains all global variables that are used in the ARTn plugin 
  !
  IMPLICIT none
  SAVE
  ! constants 
  INTEGER, PARAMETER :: DP = selected_real_kind(14,200) ! double precision 
  REAL(DP), PARAMETER :: RY2EV =  13.605691930242388_DP ! Ry to eV conversion 
  REAL(DP), PARAMETER :: B2A =  0.529177210903_DP ! bohr to angstrom conversion 
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
  NAMELIST/artn_parameters/ npushmin, neigenstepmax, nlanciter_init, push_mode,&
       convcrit_init,convcrit_final,fpara_convcrit, eigval_thr, &
       init_step_size, dlanc, step_size, &
       push_ids,add_const
  INTEGER :: npushmin              ! number of initial pushes before lanczos start
  INTEGER :: neigenstepmax         ! number of steps made with eigenvector before perp relax 
  INTEGER :: nlanciter_init        ! maximum number of lanczos iterations 
  CHARACTER (LEN = 4) :: push_mode ! type of initial push (all or list) 
  ! convergence criteria 
  REAL(DP) :: convcrit_init  ! initial perp force convergence criterion for perp relax
  REAL(DP) :: convcrit_final ! tightened force convergence criterion when near the saddle point
  REAL(DP) :: fpara_convcrit ! parallel force convergence criterion, used to determine when to tighten convcrit_final
  REAL(DP) :: eigval_thr 
  ! step sizes
  REAL(DP) :: init_step_size ! step size of inital push
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
    convcrit_init = 1.0d-2
    convcrit_final = 1.0d-3
    fpara_convcrit = 0.5d-2
    eigval_thr = -0.05_DP 
    ! 
    init_step_size = 1.5
    step_size = 0.5
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
       WRITE (iunartout,'(5X,"istep",4X,"ART_step",12X,"Etot",9X," Ftot ",5X," Fperp ",5X," Fpara ")') 
       WRITE (iunartout,'(34X, "[Ry]",9X,"-----------[Ry/a.u.]----------")')
       CLOSE ( UNIT = iunartout, STATUS = 'KEEP')
    ELSE
       WRITE(*,*) "ARTn: Input file does not exist!"
       RETURN 
ENDIF
    ! set initial number of lanczos iterations 
    nlanciter = nlanciter_init 
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

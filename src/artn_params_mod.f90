!
!> @author
!!   Matic Poberznik
!!   Miha Gunde
!!   Nicolas Salles
!
!> @brief
!!   This module contains all global variables that are used in the ARTn plugin
!
MODULE artn_params
  !
  !
  USE units, ONLY : DP
  IMPLICIT none
  SAVE
  ! constants unit pipe
  !INTEGER, PARAMETER ::  DP = selected_real_kind(14,200) ! double precision
  INTEGER, PARAMETER :: iunartin = 52    !> fortran file unit for ARTn input file
  INTEGER, PARAMETER :: iunartout = 53   !> fortran file unit for ARTn output file
  INTEGER, PARAMETER :: iunartres = 54   !> fortran file unit for ARTn restart file
  INTEGER, PARAMETER :: iunstruct = 556  !> fortran file unit for writing the structure
  INTEGER, PARAMETER :: iunrestart = 557 !> fortran file unit for writing the structure
  ! file names
  CHARACTER(LEN=255) :: filin = 'artn.in'
  CHARACTER(LEN=255) :: filout = 'artn.out'
  CHARACTER(LEN=255) :: sadfname = 'saddle'
  CHARACTER(LEN=255) :: initpfname = 'initp'
  CHARACTER(LEN=255) :: eigenfname = 'latest_eigenvec'
  CHARACTER(LEN=255) :: restartfname = 'artn.restart'
  CHARACTER(LEN=255) :: prefix_min = 'min'
  CHARACTER(LEN=255) :: prefix_sad = 'sad'
  CHARACTER(LEN=255) :: artn_resume
  ! Constante move
  INTEGER :: VOID = 1, INIT = 2, PERP = 3, EIGN = 4, LANC = 5, RELX = 6, OVER = 7
  CHARACTER(LEN=4) :: MOVE(7)
  PARAMETER( MOVE = [ 'void', 'init', 'perp', 'eign', 'lanc', 'relx', 'over' ])
  ! control flags
  LOGICAL :: linit      !> initial push OF THE MACROSTEP
  LOGICAL :: lperp      !> perpendicular relax
  LOGICAL :: leigen     !> push with lanczos eigenvector
  LOGICAL :: llanczos   !> lanczos algorithm
  LOGICAL :: lbasin     !> true while in basin 
  LOGICAL :: lsaddle    !> saddle point obtained
  LOGICAL :: lbackward  !> backward saddle point obtained
  LOGICAL :: lmove_nextmin  !> backward saddle point obtained
  LOGICAL :: lread_param  !> flag read artn params
  !
  LOGICAL :: lend
  INTEGER :: verbose    !> Verbose Level
  ! counters
  INTEGER :: istep, iartn, ifails
  INTEGER :: iperp      !> number of steps in perpendicular relaxation
  INTEGER :: nperp, nperp_step, noperp
  INTEGER :: nperp_list(5) = [ 4, 8, 12, 16, 0 ]
  !INTEGER :: nperp_list(5) = [ 4, 0, 0, 0, 0 ]
  INTEGER :: iover
  INTEGER :: irelax     !> Number of relaxation iteration
  INTEGER :: ieigen     !> number of steps made with eigenvector
  INTEGER :: iinit      !> number of pushes made
  INTEGER :: ilanc      !> current lanczos iteration
  INTEGER :: nlanc      !> max number of lanczos iterations
  INTEGER :: ismooth    !> number of smoothing steps
  INTEGER :: if_pos_ct  !> counter used to determine the number of fixed coordinates
  INTEGER :: zseed      !> random number generator seed
  INTEGER :: ifound     !> Number of saddle point found
  INTEGER :: isearch    !> Number of saddle point research

  INTEGER :: nmin       !> count the number of minimum found
  INTEGER :: nsaddle    !> count the number of saddle point found
  ! lanczos variables
  REAL(DP) :: lowest_eigval !> Lowest eigenvalues obtained by lanczos algorithm
  !                                           !
  ! arrays that are needed by ARTn internally !
  !                                           !
  REAL(DP) :: lat(3,3)
  REAL(DP), ALLOCATABLE :: tau_init(:,:)
  REAL(DP), ALLOCATABLE :: tau_nextmin(:,:)
  REAL(DP), ALLOCATABLE :: delr(:,:)
  REAL(DP), ALLOCATABLE :: push(:,:)             !> initial push vector
  REAL(DP), ALLOCATABLE, target :: eigenvec(:,:) !> lanczos eigenvector
  REAL(DP), ALLOCATABLE :: tau_step(:,:)         !> current coordinates (restart)
  REAL(DP), ALLOCATABLE :: force_step(:,:)       !> current force (restart)
  REAL(DP), ALLOCATABLE :: tau_saddle(:,:)       !> coordinates of saddle point
  REAL(DP), ALLOCATABLE :: eigen_saddle(:,:)     !> saddle point eigenvector
  !
  ! stored total energies and energy differences
  !
  REAL(DP) :: etot_init    !>  the total energy of the initial state
  REAL(DP) :: etot_step    !>  the total energy in the current step
  REAL(DP) :: etot_saddle  !>  the total energy of the saddle point
  REAL(DP) :: etot_final   !>  the total energy of the next minimum along eigenvector
  REAL(DP) :: de_saddle !> change in E from starting point
  REAL(DP) :: de_back   !> backward barrier
  REAL(DP) :: de_fwd    !> forward barrier
  !                                               !
  ! arrays that are used by the Lanczos algorithm !
  !                                               !
  REAL(DP) :: a1  !> dot product between previous and actual min lanczos vector
  REAL(DP), ALLOCATABLE :: old_lanczos_vec(:,:) !> Store the previous lanczos vec
  REAL(DP), ALLOCATABLE :: H(:,:)               !> tridiagonal matrix
  REAL(DP), ALLOCATABLE :: Vmat(:,:,:)          !> matrix containing the laczos vectors
  REAL(DP), ALLOCATABLE :: force_old(:,:)       !> force in the previous step
  REAL(DP), ALLOCATABLE :: v_in(:,:)            !> first lanczos eigenvector
  !------------------------------------------------------------!
  ! variables that are read from the input  start here
  !------------------------------------------------------------!
  !
  LOGICAL :: lrestart    !> do we want to restart the calculation
  LOGICAL :: lrelax      !> do start the relaxation to adjacent minima from the saddle point
  LOGICAL :: lpush_final !> push to adjacent minimum along eigenvector
  !
  INTEGER :: ninit                !> number of initial pushes before lanczos start
  INTEGER :: neigen               !> number of steps made with eigenvector before perp relax
  INTEGER :: lanc_mat_size        !> size of the lanczos tridiagonal matrix 
  INTEGER :: nsmooth              !> number of smoothing steps from push to eigenvec
  CHARACTER(LEN = 4) :: push_mode !> type of initial push (all , list or rad)
  ! convergence criteria
  REAL(DP) :: dist_thr       !> distance threshold for push mode "rad"
  REAL(DP) :: init_forc_thr  !> initial perp force threshold for perp relax convergence
  REAL(DP) :: forc_thr  !> tightened force convergence criterion when near the saddle point
  REAL(DP) :: fpara_thr !> parallel force convergence criterion, used to determine when to tighten convcrit_final
  REAL(DP) :: eigval_thr     !> threshold for eigenvalue
  REAL(DP) :: frelax_ene_thr      !> threshold to start relaxation to adjacent minima
  ! step sizes
  REAL(DP) :: push_step_size        !> step size of inital push in angstrom
  REAL(DP) :: eigen_step_size       !> step size for a step with the lanczos eigenvector
  REAL(DP) :: current_step_size     !> controls the current size of eigenvector step
  REAL(DP) :: fpush_factor          !> factor for the final push 
  REAL(DP), target :: dlanc         !> step size in the lanczos algorithm 
  REAL(DP) :: push_over             !> EigenVec fraction Push_over the saddle point for the relax
  ! Default Values (in Ry, au)
  REAL(DP), PARAMETER :: NAN = HUGE( dlanc )  !! Biggest number in DP representation
  REAL(DP), PARAMETER :: def_dist_thr = 0.0_DP,       def_init_forc_thr = 1.0d-2,   &
                         def_forc_thr = 1.0d-3,       def_fpara_thr = 0.5d-2,  &
                         def_eigval_thr = -0.01_DP,   def_frelax_ene_thr  = 0.00_DP,    &
                         def_push_step_size = 0.3,    def_eigen_step_size = 0.2,    &
                         def_dlanc = 1.D-2
  ! arrays related to constraints
  INTEGER,  ALLOCATABLE :: push_ids(:)    !> IDs of atoms to be pushed
  REAL(DP), ALLOCATABLE :: add_const(:,:) !> constraints on initial push
  !
  CHARACTER(LEN=256) :: engine_units
  CHARACTER(LEN=10) :: struc_format_out
  CHARACTER(LEN=3), ALLOCATABLE :: elements(:)
  CHARACTER(:), allocatable :: converge_property
  !
  NAMELIST/artn_parameters/ lrestart, lrelax, lpush_final, lmove_nextmin, &
       ninit, neigen, nperp, lanc_mat_size, nsmooth, push_mode, dist_thr,  &
       init_forc_thr,forc_thr, fpara_thr, eigval_thr, frelax_ene_thr, &
       push_step_size, dlanc, eigen_step_size, current_step_size, push_over, &
       push_ids, add_const, engine_units, zseed, struc_format_out, elements, &
       verbose, filout, sadfname, initpfname, eigenfname, restartfname, &
       converge_property

  interface warning
    module procedure :: warning_nothing, warning_int, warning_real
  end interface
  !
CONTAINS
  !
  !
  !
  SUBROUTINE initialize_artn( nat, iunartin, filnam )
    !
    !> @breif
    !!   Sets defaults, reads input and creates ARTn output file
    !
    !> @param[in] nat Number of Atoms
    !> @param[in] iunartin Channel of input
    !> @param[in] filnam Input file name
    !
    USE units
    IMPLICIT none
    ! -- Arguments
    INTEGER,             INTENT(IN) :: nat,iunartin
    CHARACTER (LEN=255), INTENT(IN) :: filnam
    ! -- Local Variables
    LOGICAL :: file_exists, verb
    INTEGER :: ios, mem
    character(len=256) :: ftmp, ctmp

    verb = .true.
    verb = .false.

    INQUIRE( file = filnam, exist = file_exists )

    write(*,'(5x,a)') "|> Initialize_ARTn()"

    IF( .not.file_exists )THEN

       WRITE(*,*) "ARTn: Input file does not exist!"
       lrelax = .true.
       RETURN

    ELSE !%! FILE EXIST
      !
      ! set up defaults for flags and counters
      !
      lrelax = .false.
      linit = .true.
      lbasin = .true.
      lperp = .false.
      llanczos = .false.
      leigen = .false.
      lsaddle = .false.
      lpush_final = .false.
      lbackward = .true.
      lrestart = .false.
      lmove_nextmin = .false.
      lread_param = .false.

      lend = .false.
      verbose = 0
      !
      ifails = 0
      iartn = 0
      istep = 0
      iinit = 0
      iperp = 0
      ilanc = 0
      ieigen = 0
      ismooth = 1
      if_pos_ct = 0
      irelax = 0
      iover = 0
      zseed = 0
      ifound = 0
      isearch = 0
      !
      lowest_eigval = 0.D0
      fpush_factor = 1.0
      push_over = 1.0_DP
      !
      ! Defaults for input parameters
      !
      ninit = 3
      nperp_step = 1
      nperp = nperp_list( nperp_step )
      noperp = 0 
      neigen = 1
      nsmooth = 1
      nmin = 0
      nsaddle = 0
      !
      dist_thr = NAN
      !
      init_forc_thr = NAN
      forc_thr = NAN
      fpara_thr = NAN
      eigval_thr = NAN ! 0.1 Ry/bohr^2 corresponds to 0.5 eV/Angs^2
      frelax_ene_thr  = NAN ! in Ry; ( etot - etot_saddle ) < frelax_ene_thr
      !
      push_step_size = NAN
      eigen_step_size = NAN
      !
      push_mode = 'all'
      struc_format_out = 'xsf'
      !
      dlanc = NAN
      lanc_mat_size = 16
      !
      engine_units = 'qe'
      !
      ! Default convergence parameter
      !
      converge_property = "maxval"
      !
      ! Allocate the arrays
      !
      IF ( .not. ALLOCATED(add_const) )    ALLOCATE(add_const(4,nat), source = 0.D0)
      IF ( .not. ALLOCATED(push_ids) )     ALLOCATE(push_ids(nat), source = 0)
      IF ( .not. ALLOCATED(push) )         ALLOCATE(push(3,nat), source = 0.D0)
      IF ( .not. ALLOCATED(eigenvec) )     ALLOCATE(eigenvec(3,nat), source = 0.D0)
      IF ( .not. ALLOCATED(eigen_saddle) ) ALLOCATE(eigen_saddle(3,nat), source = 0.D0)
      IF ( .not. ALLOCATED(tau_saddle) )   ALLOCATE(tau_saddle(3,nat), source = 0.D0)
      IF ( .not. ALLOCATED(tau_step) )   ALLOCATE(tau_step(3,nat), source = 0.D0)
      IF ( .not. ALLOCATED(force_step) )   ALLOCATE(force_step(3,nat), source = 0.D0)
      IF ( .not. ALLOCATED(force_old) ) ALLOCATE( force_old(3,nat), source = 0.D0 )
      IF ( .not. ALLOCATED(v_in) ) ALLOCATE( v_in(3,nat), source = 0.D0 )
      IF ( .not. ALLOCATED(elements) )     ALLOCATE(elements(300), source = "XXX")

      IF ( .not. ALLOCATED(delr) ) ALLOCATE( delr(3,nat), source = 0.D0 )

      ! ...Compute the size of ARTn lib
      mem = 0
      mem = mem + sizeof( add_const )
      mem = mem + sizeof( push_ids  )
      mem = mem + sizeof( push      )
      mem = mem + sizeof( eigenvec  )
      mem = mem + sizeof( eigen_saddle )
      mem = mem + sizeof( tau_saddle )
      mem = mem + sizeof( tau_step  )
      mem = mem + sizeof( force_step )
      mem = mem + sizeof( force_old )
      mem = mem + sizeof( v_in      )
      mem = mem + sizeof( elements  )
      mem = mem + sizeof( delr      )

      IF( verb )THEN
        print*, "* LIB-ARTn MEMORY: ", mem, "Bytes"
        print*, "* LIB-ARTn MEMORY: ", real(mem)/1.0e3, "KB"
        print*, "* LIB-ARTn MEMORY: ", real(mem)/1.0e6, "MB"
      ENDIF
      
      !
      ! read the ARTn input file
      OPEN( UNIT = iunartin, FILE = filnam, FORM = 'formatted', STATUS = 'unknown', IOSTAT = ios)
      READ( NML = artn_parameters, UNIT = iunartin)
      CLOSE ( UNIT = iunartin, STATUS = 'KEEP')
      lread_param = .true.
      !
      ! inital number of lanczos iterations
      !
      nlanc = lanc_mat_size
      !
      ! initialize lanczos matrices (user chooses wheter to change lanc_mat_size)
      !
      IF ( .not. ALLOCATED(H)) ALLOCATE( H(1:lanc_mat_size,1:lanc_mat_size), source = 0.D0 )
      IF ( .not. ALLOCATED(Vmat)) ALLOCATE( Vmat(3,nat,1:lanc_mat_size), source = 0.D0 )
      !
      ! initialize lanczos specific variables
    ENDIF


    !
    ! --- Read the counter file
    !
    !> min counter file
    ftmp = trim(prefix_min)//"counter"
    inquire( file=trim(ftmp), exist=file_exists )
    IF( file_exists )then
      open( newunit=ios, file=trim(ftmp), action="read" )
      read(ios,*) ctmp, ctmp, nmin 
      close( ios )
    endif

    !> saddle counter file
    ftmp = trim(prefix_sad)//"counter"
    inquire( file=trim(ftmp), exist=file_exists )
    IF( file_exists )then
      open( newunit=ios, file=trim(ftmp), action="read" )
      read(ios,*) ctmp, ctmp, nsaddle
      close( ios )
    endif


    !
    ! --- Define the Units conversion
    !
    call make_units( engine_units )


    if( verb )then
      write(*,2) repeat("*",50)
      write(*,2) "* Units:          ", trim(engine_units)
      write(*,1) "* dist_thr        = ", dist_thr
      write(*,1) "* init_forc_thr   = ", init_forc_thr
      write(*,1) "* forc_thr        = ", forc_thr
      write(*,1) "* fpara_thr       = ", fpara_thr
      write(*,1) "* eigval_thr      = ", eigval_thr
      write(*,1) "* frelax_ene_thr       = ", frelax_ene_thr
      !
      write(*,1) "* push_step_size  = ", push_step_size
      write(*,1) "* eigen_step_size = ", eigen_step_size
      write(*,1) "* dlanc           = ", dlanc
      write(*,2) repeat("*",50)
      1 format(x,a,x,g15.5)
      2 format(*(x,a))
    endif

    ! ...Convert the default values parameters from Engine_units
    !! For the moment the ARTn units is in a.u. (Ry, L, T)
    !! The default value are in ARTn units but the input values gives by the users
    !! are suppose in engine_units.
    !! We convert the USERS Values in ARTn units to be coherente:
    !! So we convert the value if it's differents from NAN initialized values

    if( dist_thr == NAN )then; dist_thr = def_dist_thr
    else;                      dist_thr = convert_length( dist_thr ); endif
    !
    if( init_forc_thr == NAN )then; init_forc_thr = def_init_forc_thr
    else;                           init_forc_thr = convert_force( init_forc_thr ); endif
    !convcrit_init = 1.0d-2
    if( forc_thr == NAN )     then;  forc_thr = def_forc_thr
    else;                            forc_thr = convert_force( forc_thr ); endif
    !convcrit_final = 1.0d-3
    if( fpara_thr == NAN )then; fpara_thr = def_fpara_thr
    else;                            fpara_thr = convert_force( fpara_thr ); endif
    !fpara_convcrit = 0.5d-2
    if( eigval_thr == NAN )then; eigval_thr = def_eigval_thr
    else;                        eigval_thr = convert_hessian( eigval_thr ); endif
    !eigval_thr = -0.01_DP ! in Ry/bohr^2 corresponds to 0.5 eV/Angs^2
    if( frelax_ene_thr == NAN )then; frelax_ene_thr = def_frelax_ene_thr
    else;                       frelax_ene_thr = convert_energy( frelax_ene_thr ); endif
    !relax_thr  = -0.01_DP ! in Ry; ( etot - etot_saddle ) < relax_thr
    !
    if( push_step_size == NAN )then; push_step_size = def_push_step_size
    else;                            push_step_size = convert_length( push_step_size ); endif
    !push_step_size = 0.3
    if( eigen_step_size == NAN )then; eigen_step_size = def_eigen_step_size
    else;                             eigen_step_size = convert_length( eigen_step_size ); endif
    !eigen_step_size = 0.2
    !
    if( dlanc == NAN )then; dlanc = def_dlanc
    else;                   dlanc = convert_length( dlanc ); endif
    !dlanc = 1.D-2

    if( verb )then
      write(*,2) repeat("*",50)
      write(*,2) "* Units:          ", trim(engine_units)
      write(*,1) "* dist_thr        = ", dist_thr
      write(*,1) "* init_forc_thr   = ", init_forc_thr
      write(*,1) "* forc_thr        = ", forc_thr
      write(*,1) "* fpara_thr       = ", fpara_thr
      write(*,1) "* eigval_thr      = ", eigval_thr
      write(*,1) "* frelax_ene_thr       = ", frelax_ene_thr
      !
      write(*,1) "* push_step_size  = ", push_step_size
      write(*,1) "* eigen_step_size = ", eigen_step_size
      write(*,1) "* dlanc           = ", dlanc
      write(*,2) repeat("*",50)
    endif


    ! ...Character verification
    converge_property = to_lower( converge_property )
    select case( converge_property )
      case( "norm", 'maxval' ); continue
      case default
        call warning( iunartout, "Initialize_artn",  &
          "converge_property has no good keyword (norm or maxval)" )
    end select

   CONTAINS

    !........................................................
    elemental Function to_lower( str )Result( string )
      !   ==============================
      !   Changes a string to lower case
      !   ==============================
      Implicit None
      Character(*), Intent(IN) :: str
      Character(LEN(str))      :: string
 
      Integer :: ic, i
 
      Character(26), Parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      Character(26), Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'
 
      !   Capitalize each letter if it is lowecase
      string = str
      do i = 1, LEN_TRIM(str)
          ic = INDEX(cap, str(i:i))
          if( ic > 0 )then
            string(i:i) = low(ic:ic)
          else
            string(i:i) = string(i:i)
          endif
      end do
    END FUNCTION to_lower


    !
  END SUBROUTINE initialize_artn


  !---------------------------------------------------------------------------
  SUBROUTINE Fill_param_step( nat, box, order, pos, etot, force )
    !
    !> @brief fill the *_step arrays on which ARTn works on.
    !! For parallel Engine each proc has list from 1 to natproc,
    !! pos_eng( i ) such as order( i ) = iat means pos( iat ) = pos_eng( i )
    !! So pos( order(i) ) = pos_eng( i )
    !
    !> @param[in]  nat    number of atoms
    !> @param[in]  box    box parameters
    !> @param[in]  order  index order of engine
    !> @param[in]  pos    atomic position
    !> @param[in]  etot   energy of the system
    !> @param[in]  force  atomic force
    !
    use units, only : convert_energy, convert_force
    
    INTEGER, INTENT(IN) :: nat, order(nat)!, types(nat)
    REAL(DP), INTENT(IN) :: box(3,3), etot, pos(3,nat), force(3,nat)

    lat = box
    tau_step(:,order(:)) = pos(:,:)
    etot_step = convert_energy( etot )
    force_step(:,order(:)) = convert_force( force(:,:) )

  END SUBROUTINE Fill_param_step



  !
  !---------------------------------------------------------------------------
  SUBROUTINE write_restart( filnres, nat )
    !
    !> @breif Subroutine that writes the minimum parameters required for restart of a calculation
    !! to a file
    !!
    ! LOGICAL FLAGS: linit, lperp, leigen, llanczos, lsaddle, lrelax
    ! COUNTERS : istep, iinit, ilanc, ieigen, ismooth, nlanc
    ! ARRAYS:
    !  - LANCZOS: eigenvec, H, Vmat, force_old, lowest_eigval
    !  - (INIT/SADDLE/STEP): etot, tau, force (, current_step_size, fpush_factor)
    !
    CHARACTER (LEN=255), INTENT(IN) :: filnres
    INTEGER, INTENT(IN) :: nat
    INTEGER :: ios
    INTEGER :: i,j
    OPEN( UNIT = iunartres, FILE = filnres, ACTION="WRITE", FORM = 'formatted', STATUS = 'unknown', IOSTAT = ios)
    !OPEN( UNIT = iunartres, FILE = filnres, FORM = 'unformatted', ACCESS="SEQUENTIAL",  &
    !      STATUS = 'unknown', ACTION="WRITE", IOSTAT = ios)

    !WRITE ( iunartres, * ) linit, lperp, leigen, llanczos, lsaddle, lrelax, &
    !     istep, iinit, ilanc, ieigen, ismooth, nlanc, nperp,  &
    !     etot_init, etot_step, lowest_eigval, etot_saddle, etot_final, de_back, &
    !     current_step_size, fpush_factor, &
    !     tau_step, force_step, push, eigenvec, H, Vmat, force_old, tau_saddle, eigen_saddle
    WRITE ( iunartres, * ) linit, lperp, leigen, llanczos, lsaddle, lrelax, &
         iartn, istep, iinit, ieigen, iperp, ilanc, irelax, ismooth,   &
         ninit, neigen, nlanc, lanc_mat_size, nperp, nmin, nsaddle, &
         etot_init, &
         etot_step, tau_step, force_step, current_step_size, fpush_factor, &    !> Actual step
         eigenvec, H, Vmat, force_old, lowest_eigval, &
         etot_saddle, tau_saddle 
    !IF( llanczos )  &
    !  WRITE( iunartres, * ) eigenvec, H, Vmat, force_old, lowest_eigval
    !IF( lsaddle )  &
    !  WRITE( iunartres, * ) etot_saddle, tau_saddle
    CLOSE ( UNIT = iunartres, STATUS = 'KEEP')

  END SUBROUTINE write_restart
  !
  !---------------------------------------------------------------------------
  SUBROUTINE read_restart( filnres, nat, order, ityp )
    !
    ! Subroutine that reads the restart file, if a restart is requested
    !
    ! LOGICAL FLAGS: lpush_init, lperp, leigen, llanczos, lsaddle, lrelax
    ! COUNTERS : istep, ipush, ilanc, ieigen, ismooth, nlanc
    ! ARRAYS: eigenvec, H, Vmat
    !
    CHARACTER (LEN=255), INTENT(IN) :: filnres
    LOGICAL :: file_exists
    INTEGER :: ios
    INTEGER :: i,nat, order(nat), ityp(nat)
    CHARACTER(LEN=255) :: fname

    INQUIRE ( file = filnres, exist = file_exists)
    IF ( file_exists ) THEN
       OPEN( UNIT = iunartres, FILE = filnres, ACTION="READ", FORM = 'formatted', STATUS = 'old', IOSTAT = ios)
       !OPEN( UNIT = iunartres, FILE = filnres, FORM = 'unformatted', ACTION="READ", STATUS = 'OLD', IOSTAT = ios)
       !READ ( iunartres, * ) linit, lperp, leigen, llanczos, lsaddle, lrelax, &
       !     istep, iinit, ilanc, ieigen, ismooth, nlanc, nperp, &
       !     etot_init, etot_step, lowest_eigval, etot_saddle, etot_final, de_back, &
       !     current_step_size, fpush_factor,  &
       !     tau_step, force_step, push, eigenvec, H, Vmat, force_old, tau_saddle, eigen_saddle
       READ( iunartres, * ) linit, lperp, leigen, llanczos, lsaddle, lrelax, &
         iartn, istep, iinit, ieigen, iperp, ilanc, irelax, ismooth,   &
         ninit, neigen, nlanc, lanc_mat_size, nperp, nmin, nsaddle, &
         etot_init, &
         etot_step, tau_step, force_step, current_step_size, fpush_factor, &   !> Actual step
         eigenvec, H, Vmat, force_old, lowest_eigval, &
         etot_saddle, tau_saddle
       !IF( llanczos )  &
       !  READ( iunartres, * ) eigenvec, H, Vmat, force_old, lowest_eigval
       !IF( lsaddle )  &
       !  READ( iunartres, * ) etot_saddle, tau_saddle
       CLOSE ( UNIT = iunartres, STATUS = 'KEEP')

       !> Maybe initialize de_back if needed

       ! ...Read the initial configuration => push, tau_init
       print*, "* RESTART:: init_structure file exist: ", trim(initpfname)
       if( .not.allocated(tau_init) )allocate( tau_init, source=tau_step)
       SELECT CASE( struc_format_out )
         CASE( 'xsf' )
           fname = TRIM(initpfname)//"."//TRIM(struc_format_out)
           CALL read_xsf( lat, nat, tau_init, order, elements, ityp, push, fname )

         CASE( 'xyz' )
           fname = TRIM(initpfname)//"."//TRIM(struc_format_out)
           CALL read_xyz( lat, nat, tau_init, order, elements, ityp, push, fname )
       END SELECT

    ELSE

       WRITE(iunartout,*) "ARTn: restart file does not exist, exiting ..."

    ENDIF


    !print*,"* RESTART::"
    !print*, nat 

  END SUBROUTINE read_restart


  !
  !---------------------------------------------------------------------------
  SUBROUTINE warning_nothing( u0, STEP, text )

    integer, intent( in ) :: u0
    character(*), intent( in ) :: STEP, text

    WRITE( u0,1 ) "* WARNING in ", STEP
    WRITE( u0,1 ) "* => ", text
    1 format(*(A))

  END SUBROUTINE warning_nothing

  SUBROUTINE warning_int( u0, STEP, text, intv )

    integer, intent( in ) :: u0, intv(:)
    character(*), intent( in ) :: STEP, text

    WRITE( u0,1 ) "* WARNING in ", STEP
    WRITE( u0,1 ) "* => ", text
    WRITE( u0,2 ) "* => ", intv
    1 format(*(A))
    2 format(A,*(x,i0))

  END SUBROUTINE warning_int


  SUBROUTINE warning_real( u0, STEP, text, realv )

    integer, intent( in ) :: u0
    REAL(DP), intent( in ) :: realv(:)
    character(*), intent( in ) :: STEP, text

    WRITE( u0,1 ) "* WARNING in ", STEP
    WRITE( u0,1 ) "* => ", text
    WRITE( u0,2 ) "* => ", realv
    1 format(*(A))
    2 format(A,*(x,f12.6))

  END SUBROUTINE warning_real


  !---------------------------------------------------------------------------
  subroutine flag_false()

    implicit none

    lrelax = .false.
    linit = .false.
    lbasin = .false.
    lperp = .false.
    llanczos = .false.
    leigen = .false.
    lsaddle = .false.

  end subroutine flag_false





  !---------------------------------------------------------------------------
  REAL(8) FUNCTION ran3( idum )
    !-------------------------------------------------------------------------
    !> @brief
    !!   Random number generator.
    !
    !> @param [in] idum   dummy integer: on first call to ran3, this is the seed,
    !                                    its value is put to 1 after the first
    !                                    call. If the calling program modifies it
    !                                    to a negative number, the generator is
    !                                    re-seeded.
    !> @return a real(8) ramdom number
    !
    !
    IMPLICIT NONE
    !
    SAVE
    !         implicit real*4(m)
    !         parameter (mbig=4000000.,mseed=1618033.,mz=0.,fac=2.5e-7)
    integer :: mbig, mseed, mz
    real(DP) :: fac
    parameter (mbig = 1000000000, mseed = 161803398, mz = 0, fac = 1.d-9)
   
    integer :: ma (55), iff, k, inext, inextp, ii, mj, idum, i, mk
    !inext = 0
    !inextp = 0
    !     common /ranz/ ma,inext,inextp
    data iff / 0 /
    if (idum.lt.0.or.iff.eq.0) then
       iff = 1
       mj = mseed-iabs (idum)
       mj = mod (mj, mbig)
       ma (55) = mj
       mk = 1
       do i = 1, 54
          ii = mod (21 * i, 55)
          ma (ii) = mk
          mk = mj - mk
          if (mk.lt.mz) mk = mk + mbig
          mj = ma (ii)
       enddo
       do k = 1, 4
          do i = 1, 55
           ma (i) = ma (i) - ma (1 + mod (i + 30, 55) )
           if (ma (i) .lt.mz) ma (i) = ma (i) + mbig
        enddo
     enddo
     inext = 0
     inextp = 31
     idum = 1
    endif
    inext = inext + 1
    if (inext.eq.56) inext = 1
    inextp = inextp + 1
    if (inextp.eq.56) inextp = 1
    mj = ma (inext) - ma (inextp)
    if (mj.lt.mz) mj = mj + mbig
    ma (inext) = mj
    ran3 = mj * fac
    return
  END FUNCTION ran3


  function dot_field( n, dx, dy )result( res )
    use units, only : DP
    implicit none
    integer, intent(in) :: n
    real(DP), intent(in) :: dx(*), dy(*)

    integer :: i
    real(DP) :: tmp, res

    res = 0.0_DP
    tmp = 0.0_DP
    do i = 1,n
       tmp = tmp + dx(i)*dy(i)
    enddo
    res = tmp
  end function dot_field
    
END MODULE artn_params
 
 

!......................................................................... FUNCTION
integer function get_iperp()
  !> @brief give the parameters IPERP
  !> @return the value of iperp
  USE artn_params, only : iperp
  get_iperp = iperp
end function get_iperp


integer function get_perp()
  !> @brief give the parameters PERP
  !> @return the value of perp
  USE artn_params, only : perp
  get_perp = perp
end function get_perp


integer function get_relx()
  !> @brief give the parameters RELX
  !> @return the value of RELX
  USE artn_params, only : relx
  get_relx = relx
end function get_relx



SUBROUTINE make_filename( f, prefix, n )

  character(*), intent(out) :: f
  character(*), intent(in) :: prefix
  integer,      intent(inout) :: n

  character(len=4) :: ctmp
  character(len=256) :: fcounter
  integer :: o0

  n = n + 1
  write( ctmp, '(I0.4)') n
  f = trim(prefix)//trim(ctmp)

  fcounter = trim(prefix)//"counter"
  open( newunit=o0, file=fcounter, action="write" )
    write( o0, '(x,"counter:",2(x,a))' ) trim(prefix), trim(ctmp)
  close( o0 )

END SUBROUTINE make_filename




!......................................................................... DEBUGING 

subroutine info_field( u0, n, v, txt )
  use units, only : DP
  implicit none

  integer, intent(in) :: u0, n
  real(DP), intent(in) :: v(3,n)
  character(*), intent(in) :: txt

  integer :: i
  real(DP) :: lmax, vtot

  write(u0,'("****** INFO FIELD:",x,a)') trim(txt)

  lmax = 0.0_DP
  do i = 1,n
     lmax = max( lmax, norm2(v(:,i)) )
  enddo
  call sum_force( v, n, vtot )

  write(u0,'("* total field: ",g12.4," | max norm: ",g12.4/)') vtot, lmax

end subroutine info_field



!subroutine max_component(n,)












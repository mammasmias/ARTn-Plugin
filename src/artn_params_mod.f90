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
  INTEGER, PARAMETER :: iunstruct = 556 ! fortran file unit for writing the structure
  INTEGER, PARAMETER :: iunrestart = 557 ! fortran file unit for writing the structure
  ! Constante move
  !INTEGER :: VOID, INIT, PERP, EIGN, LANC, RELX
  INTEGER :: VOID = 1, INIT = 2, PERP = 3, EIGN = 4, LANC = 5, RELX = 6
  !PARAMETER( VOID = 1, INIT = 2, PERP = 3, EIGN = 4, LANC = 5, RELX = 6 )
  CHARACTER(LEN=4) :: MOVE(6)
  PARAMETER( MOVE = [ 'void', 'init', 'perp', 'eign', 'lanc', 'relx' ])
  ! control flags
  LOGICAL :: lartn = .true. ! We use ARTn or not
  LOGICAL :: lpush_init ! initial push 
  LOGICAL :: lperp      ! perpendicular relax 
  LOGICAL :: leigen     ! push with lanczos eigenvector
  LOGICAL :: llanczos   ! lanczos algorithm
  LOGICAL :: lsaddle    ! saddle point obtained 
  ! counters
  integer :: istep
  INTEGER, target :: iperp      ! number of steps in perpendicular relaxation
  INTEGER :: ieigen     ! number of steps made with eigenvector
  INTEGER :: ipush      ! number of pushes made
  INTEGER :: ilanc      ! current lanczos iteration
  INTEGER :: nlanc      ! max number of lanczos iterations
  INTEGER :: ismooth    ! number of smoothing steps
  INTEGER :: if_pos_ct  ! counter used to determine the number of fixed coordinates 
  INTEGER :: zseed      ! random number generator seed
  ! lanczos variables
  REAL(DP) :: lowest_eigval
  !                                           ! 
  ! arrays that are needed by ARTn internally ! 
  !                                           ! 
  REAL(DP), ALLOCATABLE :: push(:,:)     ! initial push vector
  REAL(DP), ALLOCATABLE, target :: eigenvec(:,:) ! lanczos eigenvector
  REAL(DP), ALLOCATABLE :: tau_saddle(:,:) ! coordinates of saddle point  
  REAL(DP), ALLOCATABLE :: eigen_saddle(:,:) ! coordinates of saddle point  
  !
  ! stored total energies and energy differences 
  ! 
  REAL(DP) :: etot_init   !  the total energy of the initial state
  REAL(DP) :: etot_saddle !  the total energy of the saddle point
  REAL(DP) :: etot_final  !  the total energy of the next minimum along eigenvector 
  REAL(DP) :: de_saddle ! change in E from starting point  
  REAL(DP) :: de_back   ! backward barrier  
  REAL(DP) :: de_fwd    ! forward barrier 
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
  NAMELIST/artn_parameters/ lrelax, lpush_final, npush, neigen, nlanc_init, nsmooth, push_mode, dist_thr,  &
       convcrit_init,convcrit_final,fpara_convcrit, eigval_thr, relax_thr, &
       push_step_size, dlanc, eigen_step_size, &
       push_ids,add_const, engine_units, zseed
  ! 
  LOGICAL :: lrelax     ! do we want to relax to adjacent minima from the saddle point 
  LOGICAL :: lpush_final ! push to adjacent minimum along eigenvector 
  ! 
  INTEGER :: npush              ! number of initial pushes before lanczos start
  INTEGER :: neigen            ! number of steps made with eigenvector before perp relax 
  INTEGER :: nlanc_init        ! maximum number of lanczos iterations
  INTEGER :: nsmooth           ! number of smoothing steps from push to eigenvec 
  CHARACTER (LEN = 4) :: push_mode ! type of initial push (all , list or rad) 
  ! convergence criteria
  REAL(DP) :: dist_thr       ! distance threshold for push mode "rad" 
  REAL(DP) :: convcrit_init  ! initial perp force convergence criterion for perp relax
  REAL(DP) :: convcrit_final ! tightened force convergence criterion when near the saddle point
  REAL(DP) :: fpara_convcrit ! parallel force convergence criterion, used to determine when to tighten convcrit_final
  REAL(DP) :: eigval_thr     ! threshold for eigenvalue
  REAL(DP) :: relax_thr      ! threshold to start relaxation to adjacent minima
  ! step sizes
  REAL(DP) :: push_step_size ! step size of inital push in angstrom 
  REAL(DP) :: eigen_step_size      ! step size for a step with the lanczos eigenvector
  REAL(DP) :: current_step_size ! controls the current size of eigenvector step
  REAL(DP) :: fpush_factor     ! factor for the final push 
  REAL(DP), target :: dlanc         ! step size in the lanczos algorithm 
  ! arrays related to constraints 
  INTEGER,  ALLOCATABLE :: push_ids(:)    ! IDs of atoms to be pushed 
  REAL(DP), ALLOCATABLE :: add_const(:,:) ! constraints on initial push
  !
  CHARACTER(LEN=256) :: engine_units
  !CHARACTER(LEN=:), ALLOCATABLE :: engine_units
  
CONTAINS
  !
  SUBROUTINE initialize_artn( nat, iunartin, iunartout, filnam, filout )
    !
    ! sets defaults, reads input and creates ARTn output file
    ! 
    USE units
    IMPLICIT none
    INTEGER, INTENT(IN) :: nat,iunartin,iunartout
    CHARACTER (LEN=255), INTENT(IN) :: filnam,filout
    LOGICAL :: file_exists
    INTEGER :: ios


    lartn = .true.
    INQUIRE( file = filnam, exist = file_exists )
    
    IF( .not.file_exists )THEN

       WRITE(*,*) "ARTn: Input file does not exist!"
       lartn = .false.
       lpush_init = .false.
       lrelax = .true.
       RETURN 


    ELSE !%! FILE EXIST


    !
    ! set up defaults for flags and counters 
    !
    lrelax = .false. 
    lpush_init = .true.
    lperp = .false.
    llanczos = .false.
    leigen = .false.
    lsaddle = .false. 
    lpush_final = .false. 
    !
    istep = 0
    ipush = 0 
    iperp = 0
    ilanc = 0
    ieigen = 0
    ismooth = 1
    if_pos_ct = 0 
    !
    lowest_eigval = 0.D0
    !
    ! Defaults for input parameters
    ! 
    npush = 3    
    neigen = 1
    nsmooth = 1 
    !
    dist_thr = 0.0_DP 
    ! 
    convcrit_init = 1.0d-2
    convcrit_final = 1.0d-3
    fpara_convcrit = 0.5d-2
    eigval_thr = -0.01_DP ! in Ry/bohr^2 corresponds to 0.5 eV/Angs^2
    relax_thr  = -0.01_DP ! in Ry; ( etot - etot_saddle ) < relax_thr 
    ! 
    push_step_size = 0.3
    eigen_step_size = 0.2
    fpush_factor = 1.0
    !
    push_mode = 'all' 
    !
    dlanc = 1.D-2 
    nlanc_init = 16
    !
    ! Allocate the arrays
    !
    IF ( .not. ALLOCATED(add_const)) ALLOCATE(add_const(4,nat), source = 0.D0)
    IF ( .not. ALLOCATED(push_ids)) ALLOCATE(push_ids(nat), source = 0)
    IF ( .not. ALLOCATED(push)) ALLOCATE(push(3,nat), source = 0.D0)
    IF ( .not. ALLOCATED(eigenvec)) ALLOCATE(eigenvec(3,nat), source = 0.D0)
    IF ( .not. ALLOCATED(eigen_saddle)) ALLOCATE(eigen_saddle(3,nat), source = 0.D0)
    IF ( .not. ALLOCATED(tau_saddle)) ALLOCATE(tau_saddle(3,nat), source = 0.D0)
    ! 
    !IF ( file_exists ) THEN
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
       WRITE (iunartout,'(15X,"npush           = ", I6)') npush 
       WRITE (iunartout,'(15X,"convcrit_init   = ", F6.3)') convcrit_init 
       WRITE (iunartout,'(15X,"convcrit_final  = ", F6.3)') convcrit_final
       WRITE (iunartout,'(15X,"fpara_convcrit  = ", F6.3)') fpara_convcrit 
       WRITE (iunartout,'(15X,"eigval_thr      = ", F6.3)') eigval_thr  
       WRITE (iunartout,'(15X,"push_step_size  = ", F6.1)') push_step_size 
       WRITE (iunartout,'(15X,"eigen_step_size = ", F6.1)') eigen_step_size
       WRITE (iunartout,'(15X,"push_mode       = ", A6)') push_mode
       WRITE (iunartout,'(5X, "--------------------------------------------------")')
       WRITE (iunartout,'(5X, "Lanczos algorithm:")' )
       WRITE (iunartout,'(5X, "--------------------------------------------------")')
       WRITE (iunartout,'(15X, "nlanc_init     = ", I6)') nlanc_init
       WRITE (iunartout,'(15X, "dlanc          = ", F6.3)') dlanc
       WRITE (iunartout,'(5X, "--------------------------------------------------")')
       WRITE (iunartout,*) " "
       WRITE (iunartout,*) " "
       WRITE (iunartout,'(5X,"istep",4X,"ART_step",12X,"Etot",12X," Ftot ",9X," Fperp ",8X," Fpara ",8X,"eigval")') 
       WRITE (iunartout,'(34X, "[Ry]",15X,"-----------[Ry/a.u.]----------",10X,"Ry/a.u.^2")')
       CLOSE ( UNIT = iunartout, STATUS = 'KEEP')
    !ELSE
    !   WRITE(*,*) "ARTn: Input file does not exist!"
    !   lpush_init = .false.
    !   lrelax = .true.
    !   RETURN 
    ENDIF
    ! set initial number of lanczos iterations 
    nlanc = nlanc_init
    ! initialize lanczos specific variables
    IF ( .not. ALLOCATED(H)) ALLOCATE( H(1:nlanc,1:nlanc), source = 0.D0 )
    IF ( .not. ALLOCATED(Vmat)) ALLOCATE( Vmat(3,nat,1:nlanc), source = 0.D0 )
    IF ( .not. ALLOCATED(force_old) ) ALLOCATE( force_old(3,nat), source = 0.D0 )


    ! Define the UNits convertion
    call make_units( engine_units )

    ! ...Convert push/lanczos/eigenvec step size to bohr (because force units are in Ry/bohr) 
    !eigen_step_size = eigen_step_size/B2A
    !push_step_size = push_step_size/B2A
    !dlanc = dlanc/B2A
    !dist_thr = dist_thr/B2A
    eigen_step_size = convert_length( eigen_step_size )
    push_step_size = convert_length( push_step_size )
    dlanc = convert_length( dlanc )
    dist_thr = convert_length( dist_thr )
    ! 
  END SUBROUTINE initialize_artn
  !

!---------------------------------------------------------------------------
REAL(8) FUNCTION ran3( idum )
  !-------------------------------------------------------------------------
  !! Random number generator.
  !
  !USE kinds, ONLY : DP
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


END MODULE artn_params

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
  ! constants
  !INTEGER, PARAMETER ::  DP = selected_real_kind(14,200) ! double precision
  INTEGER, PARAMETER :: iunartin = 52    !> fortran file unit for ARTn input file
  INTEGER, PARAMETER :: iunartout = 53   !> fortran file unit for ARTn output file
  INTEGER, PARAMETER :: iunartres = 54   !> fortran file unit for ARTn restart file
  INTEGER, PARAMETER :: iunstruct = 556  !> fortran file unit for writing the structure
  INTEGER, PARAMETER :: iunrestart = 557 !> fortran file unit for writing the structure
  ! Constante move
  INTEGER :: VOID = 1, INIT = 2, PERP = 3, EIGN = 4, LANC = 5, RELX = 6
  CHARACTER(LEN=4) :: MOVE(6)
  PARAMETER( MOVE = [ 'void', 'init', 'perp', 'eign', 'lanc', 'relx' ])
  ! control flags
  LOGICAL :: lartn = .true. !> We use ARTn or not
  LOGICAL :: lpush_init !> initial push
  LOGICAL :: lperp      !> perpendicular relax
  LOGICAL :: leigen     !> push with lanczos eigenvector
  LOGICAL :: llanczos   !> lanczos algorithm
  LOGICAL :: lsaddle    !> saddle point obtained
  LOGICAL :: lbackward  !> backward saddle point obtained
  ! counters
  INTEGER :: istep
  INTEGER, target :: iperp      !> number of steps in perpendicular relaxation
  INTEGER :: ieigen     !> number of steps made with eigenvector
  INTEGER :: ipush      !> number of pushes made
  INTEGER :: ilanc      !> current lanczos iteration
  INTEGER :: nlanc      !> max number of lanczos iterations
  INTEGER :: ismooth    !> number of smoothing steps
  INTEGER :: if_pos_ct  !> counter used to determine the number of fixed coordinates
  INTEGER :: zseed      !> random number generator seed
  ! lanczos variables
  REAL(DP) :: lowest_eigval !> Lowest eigenvalues obtained by lanczos algorithm
  !                                           !
  ! arrays that are needed by ARTn internally !
  !                                           !
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
  REAL(DP) :: etot_step !>  the total energy in the current step
  REAL(DP) :: etot_saddle  !>  the total energy of the saddle point
  REAL(DP) :: etot_final   !>  the total energy of the next minimum along eigenvector
  REAL(DP) :: de_saddle !> change in E from starting point
  REAL(DP) :: de_back   !> backward barrier
  REAL(DP) :: de_fwd    !> forward barrier
  !                                               !
  ! arrays that are used by the Lanczos algorithm !
  !                                               !
  REAL(DP), ALLOCATABLE :: H(:,:)         !> tridiagonal matrix
  REAL(DP), ALLOCATABLE :: Vmat(:,:,:)    !> matrix containing the laczos vectors
  REAL(DP), ALLOCATABLE :: force_old(:,:) !> force in the previous step
  REAL(DP), ALLOCATABLE :: v_in(:,:)      !> first lanczos eigenvector
  !------------------------------------------------------------!
  ! variables that are read from the input  start here
  !------------------------------------------------------------!
  !
  LOGICAL :: lrestart    !> do we want to restart the calculation
  LOGICAL :: lrelax      !> do we want to relax to adjacent minima from the saddle point
  LOGICAL :: lpush_final !> push to adjacent minimum along eigenvector
  !
  INTEGER :: npush                !> number of initial pushes before lanczos start
  INTEGER :: neigen               !> number of steps made with eigenvector before perp relax
  INTEGER :: nlanc_init           !> maximum number of lanczos iterations
  INTEGER :: nsmooth              !> number of smoothing steps from push to eigenvec
  CHARACTER(LEN = 4) :: push_mode !> type of initial push (all , list or rad)
  ! convergence criteria
  REAL(DP) :: dist_thr       !> distance threshold for push mode "rad"
  REAL(DP) :: convcrit_init  !> initial perp force convergence criterion for perp relax
  REAL(DP) :: convcrit_final !> tightened force convergence criterion when near the saddle point
  REAL(DP) :: fpara_convcrit !> parallel force convergence criterion, used to determine when to tighten convcrit_final
  REAL(DP) :: eigval_thr     !> threshold for eigenvalue
  REAL(DP) :: relax_thr      !> threshold to start relaxation to adjacent minima
  ! step sizes
  REAL(DP) :: push_step_size        !> step size of inital push in angstrom
  REAL(DP) :: eigen_step_size       !> step size for a step with the lanczos eigenvector
  REAL(DP) :: current_step_size     !> controls the current size of eigenvector step
  REAL(DP) :: fpush_factor          !> factor for the final push
  REAL(DP), target :: dlanc         !> step size in the lanczos algorithm
  ! Default Values
  REAL(DP), PARAMETER :: NAN = HUGE( dlanc )  !! Biggest number in DP representation
  REAL(DP), PARAMETER :: def_dist_thr = 0.0_DP,       def_convcrit_init = 1.0d-2,   &
                         def_convcrit_final = 1.0d-3, def_fpara_convcrit = 0.5d-2,  &
                         def_eigval_thr = -0.01_DP,   def_relax_thr  = -0.01_DP,    &
                         def_push_step_size = 0.3,    def_eigen_step_size = 0.2,    &
                         def_dlanc = 1.D-2
  ! arrays related to constraints
  INTEGER,  ALLOCATABLE :: push_ids(:)    !> IDs of atoms to be pushed
  REAL(DP), ALLOCATABLE :: add_const(:,:) !> constraints on initial push
  !
  CHARACTER(LEN=256) :: engine_units
  CHARACTER(LEN=10) :: struc_format_out
  CHARACTER(LEN=3), ALLOCATABLE :: elements(:)
  !
  NAMELIST/artn_parameters/ lrestart, lrelax, lpush_final, npush, neigen, nlanc_init, nsmooth, push_mode, dist_thr,  &
       convcrit_init,convcrit_final,fpara_convcrit, eigval_thr, relax_thr, &
       push_step_size, dlanc, eigen_step_size, current_step_size, &
       push_ids,add_const, engine_units, zseed, struc_format_out, elements
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
    !> @param[in] iunartout Channel of Output
    !> @param[in] filnam Input file name
    !> @param[in] filout Ouput file name
    !
    USE units
    IMPLICIT none
    ! -- Arguments
    INTEGER,             INTENT(IN) :: nat,iunartin
    CHARACTER (LEN=255), INTENT(IN) :: filnam
    ! -- Local Variables
    LOGICAL :: file_exists, verbose
    INTEGER :: ios

    verbose = .true.
    verbose = .false.

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
      lbackward = .true.
      lrestart = .false.
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
      fpush_factor = 1.0
      !
      ! Defaults for input parameters
      !
      npush = 3
      neigen = 1
      nsmooth = 1
      !
      dist_thr = NAN
      !
      convcrit_init = NAN
      convcrit_final = NAN
      fpara_convcrit = NAN
      eigval_thr = NAN ! 0.1 Ry/bohr^2 corresponds to 0.5 eV/Angs^2
      relax_thr  = NAN ! in Ry; ( etot - etot_saddle ) < relax_thr
      !
      push_step_size = NAN
      eigen_step_size = NAN
      !
      push_mode = 'all'
      struc_format_out = 'xsf'
      !
      dlanc = NAN
      nlanc_init = 16
      !
      engine_units = 'qe'
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
      !
      ! read the ARTn input file
      OPEN( UNIT = iunartin, FILE = filnam, FORM = 'formatted', STATUS = 'unknown', IOSTAT = ios)
      READ( NML = artn_parameters, UNIT = iunartin)
      CLOSE ( UNIT = iunartin, STATUS = 'KEEP')
      !
      ! inital number of lanczos iterations
      !
      nlanc = nlanc_init
      !
      ! initialize lanczos matrices (user chooses wheter to change nlanc_init)
      !
      IF ( .not. ALLOCATED(H)) ALLOCATE( H(1:nlanc_init,1:nlanc_init), source = 0.D0 )
      IF ( .not. ALLOCATED(Vmat)) ALLOCATE( Vmat(3,nat,1:nlanc_init), source = 0.D0 )
      !
      ! initialize lanczos specific variables
    ENDIF


    ! Define the Units conversion
    call make_units( engine_units )


    !#! OLD VERSION
    ! ...Convert push/lanczos/eigenvec step size to bohr (because force units are in Ry/bohr)

    !eigen_step_size = convert_length( eigen_step_size )
    !push_step_size = convert_length( push_step_size )
    !dlanc = convert_length( dlanc )
    !dist_thr = convert_length( dist_thr )

    if( verbose )then
      write(*,2) repeat("*",50)
      write(*,2) "* Units:          ", trim(engine_units)
      write(*,1) "* dist_thr        = ", dist_thr
      write(*,1) "* convcrit_init   = ", convcrit_init
      write(*,1) "* convcrit_final  = ", convcrit_final
      write(*,1) "* fpara_convcrit  = ", fpara_convcrit
      write(*,1) "* eigval_thr      = ", eigval_thr
      write(*,1) "* relax_thr       = ", relax_thr
      !
      write(*,1) "* push_step_size  = ", push_step_size
      write(*,1) "* eigen_step_size = ", eigen_step_size
      write(*,1) "* dlanc           = ", dlanc
      write(*,2) repeat("*",50)
      !1 format(x,a,x,g15.5)
      !2 format(*(x,a))
    endif

    !#! NEW VERSION
    ! ...Convert the default values parameters from Engine_units
    !! For the moment the ARTn units is in a.u. (Ry, L, T)
    !! The default value are in ARTn units but the input values gives by the users
    !! are suppose in engine_units.
    !! We convert the USERS Values in ARTn units to be coherente:
    !! So we convert the value if it's differents from NAN initialized values

    if( dist_thr == NAN )then; dist_thr = def_dist_thr
    else;                      dist_thr = convert_length( dist_thr ); endif
    !
    if( convcrit_init == NAN )then; convcrit_init = def_convcrit_init
    else;                           convcrit_init = convert_force( convcrit_init ); endif
    !convcrit_init = 1.0d-2
    if( convcrit_final == NAN )then; convcrit_final = def_convcrit_final
    else;                            convcrit_final = convert_force( convcrit_final ); endif
    !convcrit_final = 1.0d-3
    if( fpara_convcrit == NAN )then; fpara_convcrit = def_fpara_convcrit
    else;                            fpara_convcrit = convert_force( fpara_convcrit ); endif
    !fpara_convcrit = 0.5d-2
    if( eigval_thr == NAN )then; eigval_thr = def_eigval_thr
    else;                        eigval_thr = convert_hessian( eigval_thr ); endif
    !eigval_thr = -0.01_DP ! in Ry/bohr^2 corresponds to 0.5 eV/Angs^2
    if( relax_thr == NAN )then; relax_thr = def_relax_thr
    else;                       relax_thr = convert_energy( relax_thr ); endif
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

    if( verbose )then
      write(*,2) repeat("*",50)
      write(*,2) "* Units:          ", trim(engine_units)
      write(*,1) "* dist_thr        = ", dist_thr
      write(*,1) "* convcrit_init   = ", convcrit_init
      write(*,1) "* convcrit_final  = ", convcrit_final
      write(*,1) "* fpara_convcrit  = ", fpara_convcrit
      write(*,1) "* eigval_thr      = ", eigval_thr
      write(*,1) "* relax_thr       = ", relax_thr
      !
      write(*,1) "* push_step_size  = ", push_step_size
      write(*,1) "* eigen_step_size = ", eigen_step_size
      write(*,1) "* dlanc           = ", dlanc
      write(*,2) repeat("*",50)
      1 format(x,a,x,g15.5)
      2 format(*(x,a))
    endif

    !
  END SUBROUTINE initialize_artn
  !
  SUBROUTINE write_initial_report(iunartout, filout)
    INTEGER,             INTENT(IN) :: iunartout
    CHARACTER (LEN=255), INTENT(IN) :: filout
    ! -- Local Variables
    INTEGER :: ios
    !
    ! Writes the header to the artn output file
    !
    OPEN ( UNIT = iunartout, FILE = filout, FORM = 'formatted', STATUS = 'unknown', IOSTAT = ios )
    WRITE (iunartout,'(5X, "--------------------------------------------------")')
    WRITE (iunartout,'(5X, "                ARTn plugin                       ")')
    WRITE (iunartout,'(5X, "--------------------------------------------------")')
    WRITE (iunartout,'(5X, " "                                                 )')
    WRITE (iunartout,'(5X, "               INPUT PARAMETERS                   ")')
    WRITE (iunartout,'(5X, "--------------------------------------------------")')
    WRITE (iunartout,'(5x, "engine_units:", *(x,A))') TRIM(engine_units)
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

  END SUBROUTINE write_initial_report
  !
  SUBROUTINE write_restart(filnres,nat)
    !
    ! Subroutine that writes the minimum parameters required for restart of a calculation
    ! to a file
    !
    ! LOGICAL FLAGS: lpush_init, lperp, leigen, llanczos, lsaddle, lrelax
    ! COUNTERS : istep, ipush, ilanc, ieigen, ismooth, nlanc
    ! ARRAYS: eigenvec, H, Vmat
    !
    CHARACTER (LEN=255), INTENT(IN) :: filnres
    INTEGER, INTENT(IN) :: nat
    INTEGER :: ios
    INTEGER :: i,j
    OPEN( UNIT = iunartres, FILE = filnres, FORM = 'formatted', STATUS = 'unknown', IOSTAT = ios)

    WRITE ( iunartres, * ) lpush_init, lperp, leigen, llanczos, lsaddle, lrelax, &
         istep, ipush, ilanc, ieigen, ismooth, nlanc,  &
         etot_init, etot_step, lowest_eigval, etot_saddle, etot_final, de_back, &
     current_step_size, fpush_factor, &
         tau_step, force_step, push, eigenvec, H, Vmat, force_old, tau_saddle, eigen_saddle
    CLOSE ( UNIT = iunartres, STATUS = 'KEEP')

  END SUBROUTINE write_restart
  !
  SUBROUTINE read_restart(filnres,nat)
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
    INTEGER :: i,nat
    INQUIRE ( file = filnres, exist = file_exists)
    IF ( file_exists ) THEN
       OPEN( UNIT = iunartres, FILE = filnres, FORM = 'formatted', STATUS = 'unknown', IOSTAT = ios)
       READ ( iunartres, * ) lpush_init, lperp, leigen, llanczos, lsaddle, lrelax, &
            istep, ipush, ilanc, ieigen, ismooth, nlanc , &
            etot_init, etot_step, lowest_eigval, etot_saddle, etot_final, de_back, &
            current_step_size, fpush_factor,  &
            tau_step, force_step, push, eigenvec, H, Vmat, force_old, tau_saddle, eigen_saddle
       CLOSE ( UNIT = iunartres, STATUS = 'KEEP')

    ELSE

       WRITE(iunartout,*) "ARTn: restart file does not exist, exiting ..."

    ENDIF

  END SUBROUTINE read_restart
  !
 !---------------------------------------------------------------------------
REAL(8) FUNCTION ran3( idum )
  !-------------------------------------------------------------------------
  !> @brief
  !!   Random number generator.
  !
  !> @param [in] idum   dummy integer
  !> @return a real(8) ramdom number
  !
  !#! WE HAVE TO INITIALIZE THE VARIABLE
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


END MODULE artn_params


!......................................................................... FUNCTION
!> @brief give the parameters IPERP
!> @return the value of iperp
integer function get_iperp()
  USE artn_params, only : iperp
  get_iperp = iperp
end function get_iperp


!> @brief give the parameters PERP
!> @return the value of perp
integer function get_perp()
  USE artn_params, only : perp
  get_perp = perp
end function get_perp


!> @brief give the parameters RELX
!> @return the value of RELX
integer function get_relx()
  USE artn_params, only : relx
  get_relx = relx
end function get_relx


















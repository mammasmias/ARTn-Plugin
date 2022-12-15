
!> @author Matic Poberznik,
!! @author Miha Gunde,
!! @author Nicolas Salles

!> @brief 
!!   UNITS module contains all the tool to reconize the Engine and its units
!!   to convert the energy/force/length/time in atomic units
!!   Atomic Units (au) in plugin-ARTn is the Rydberg-bohr-aut
!
!> @todo
!!   Change the unit philosophy: In principle ARTn could work whitout 
!!   to convert the quantities. 
!
Module units
  !
  PRIVATE

  PUBLIC :: DP, PI, Mass, B2A,  make_units,   &
            convert_length, unconvert_length,   &
            convert_force, unconvert_force,     &
            convert_hessian, unconvert_hessian, &
            convert_energy, unconvert_energy,   &
            convert_time, unconvert_time, strg_units, unit_char

  PUBLIC :: parser, lower, read_line
             
  

  INTEGER, PARAMETER ::  DP = selected_real_kind(14,200)           !< @brief double precision
  REAL(DP), PARAMETER :: PI     = 3.14159265358979323846_DP        !< @brief pi number 

  REAL(DP), PARAMETER :: H_PLANCK_SI      = 6.62607015E-34_DP      !< @brief J s
  REAL(DP), PARAMETER :: K_BOLTZMANN_SI   = 1.380649E-23_DP        !< @brief J K^-1 
  REAL(DP), PARAMETER :: ELECTRON_SI      = 1.602176634E-19_DP     !< @brief C
  REAL(DP), PARAMETER :: ELECTRONVOLT_SI  = 1.602176634E-19_DP     !< @brief J  
  REAL(DP), PARAMETER :: ELECTRONMASS_SI  = 9.1093837015E-31_DP    !< @brief Kg
  REAL(DP), PARAMETER :: HARTREE_SI       = 4.3597447222071E-18_DP !< @brief J
  REAL(DP), PARAMETER :: RYDBERG_SI       = HARTREE_SI/2.0_DP      !< @brief J
  REAL(DP), PARAMETER :: BOHR_RADIUS_SI   = 0.529177210903E-10_DP  !< @brief m
  REAL(DP), PARAMETER :: AMU_SI           = 1.66053906660E-27_DP   !< @brief Kg
  REAL(DP), PARAMETER :: C_SI             = 2.99792458E+8_DP       !< @brief m sec^-1
  REAL(DP), PARAMETER :: NA               = 6.022140857E+23_DP     !< @brief mol^-1

  REAL(DP), PARAMETER :: RY2EV            = 13.605691930242388_DP  !< @brief Ry to eV conversion 
  REAL(DP), PARAMETER :: RY2KCAL          = 5.2065348237317E-22_DP !< @brief Ry to kcal conversion 
  REAL(DP), PARAMETER :: RY2KJ            = 2.17987197E-21_DP      !< @brief Ry to kJoules conversion 
  REAL(DP), PARAMETER :: RY2KCALPMOL      = RY2KCAL*NA             !< @brief Ry to kcal/mole conversion 
  REAL(DP), PARAMETER :: RY2KJPMOL        = RY2KJ*NA               !< @brief Ry to kJoules per mole conversion 
  REAL(DP), PARAMETER :: B2A              = 0.529177210903_DP      !< @brief bohr to angstrom conversion (Used for QE engine)
  REAL(DP), PARAMETER :: AMU_RY2          = 911.44424310865645_DP  !< @brief calculated from QE using DP
  REAL(DP), PARAMETER :: ps2aut           = 41341.374575751 / 2.   !< @brief picosecond to atomic unit of time
  REAL(DP), PARAMETER :: aut2s            = 4.8278E-17_DP          !< @brief atomic units of times to second conversion (Ry atomic unit)

  REAL(DP), PARAMETER :: AMU_AU           = AMU_SI / ELECTRONMASS_SI  !< @brief Dimensionless Hartree
  REAL(DP), PARAMETER :: AMU_RY           = AMU_AU / 2.0_DP           !< @brief Dimensionless Rydberg

  !REAL(DP), PARAMETER :: AU_SEC           = H_PLANCK_SI/(2.*pi)/HARTREE_SI
  REAL(DP), PARAMETER :: AU_SEC           = H_PLANCK_SI/(2.*pi)/RYDBERG_SI    !< @brief Atomic unit of time to second
  REAL(DP), PARAMETER :: AU_PS            = AU_SEC * 1.0E+12_DP               !< @brief Atomic unit of time to picosecond
  REAL(DP), PARAMETER :: AU_FS            = AU_SEC * 1.0E+15_DP               !< @brief Atomic unit of time to femtosecond


  !! Units Character
  character(*), parameter :: AA = char(197)    !< @brief  Angstrom (ANSI code)
  !character(*), parameter :: to2 = char(178)  ! exponent 2
  character(*), parameter :: to2 = "**2"       !< @brief exponent 2

  character(:), allocatable :: cL, cE

  !! Units convertor
  CHARACTER(LEN=256) :: strg_units             !< @brief String containing the unit of the system with the output format 
  REAL(DP) :: Mass                             !< @brief Mass in Rydberg to buid the force - ARTn is in Rydberg (QE)
  REAL(DP) :: E2au, au2E, L2au, au2L, T2au, au2T, F2au, au2F, H2au, au2H

 contains

  !......................................................................................
  !> @brief  Convert an Array of Capital letter to lower case letter
  function lower( s1 )result( s2 )
    !
    !> @param[in]  s1   input string, contain some capital letter
    !! @return     s2   string with only lower case
    !
    character(*)       :: s1
    character(len(s1)) :: s2
    character          :: ch
    integer,parameter  :: duc = ichar('A') - ichar('a')
    integer            :: i

    do i = 1,len(s1)
       ch = s1(i:i)
       if (ch >= 'A'.and. ch <= 'Z') ch = char(ichar(ch)-duc)
       s2(i:i) = ch
    end do
  end function lower


  !................................................................................
  !> @brief 
  !!   Parse the instrg thank to the Field Separator FS and return 
  !!   the list of string and the number of element in the list
  integer function parser(instrg, FS, args )result( nargs )
    !
    !> @todo 
    !!   HAVE TO BE ADAPTED FOR MULTIPLE FS
    !
    !> @param[in]   instrg   input string 
    !> @param[in]   FS       Field Separator (one for the moment)
    !> @param[out]  args     arrays of string
    !> @return      nargs    number of string in output
    !
    implicit none
 
    ! -- ARGUMENT
    CHARACTER(len=*),              intent( in ) :: instrg
    character(len=1),              intent( in ) :: FS
    CHARACTER(len=:), allocatable, intent( inout ) :: args(:)
 
    ! -- LOCAL VAR
    character(len=:), allocatable :: str
    character(len=25) :: mot
    integer :: idx,leng
 
    ! +++ Copy in local variable the input_string
    str = adjustl(instrg)
    nargs = 0
    !
    ! +++ Repeat for each field
    do
    !   +++ Verification the length of sentence
        leng = len_TRIM( str )
        if( leng == 0 )exit
 
    !   +++ Find the Field Separator
        idx = SCAN( str, FS )
 
    !   +++ extract the word
        if( idx == 0 )then
          mot = trim(str)
        else
          mot = str( :idx-1 )
        endif
 
    !   +++ Add the word in args
        nargs = nargs + 1
        if( nargs == 1 )then
          args = [ mot ]
        else
          args = [ args(:), mot ]
        endif
 
    !   +++ cut the word
        if( idx == 0 )exit
        str = adjustl(str(idx+1:))
    !
    enddo
  end function parser
  
  ! .............................................................................
  !> @brief
  !!   Read a line, makes possible to use # for comment lines, skips empty lines, 
  !!   is pretty much a copy from QE.
  !!
  !> @details  Quantum ESPRESSO routine
  !
  !> @param[in]                 file descriptor
  !! @param[out]  line          what it reads
  !! @param[out]  end_of_file   logical to signal the EOF
  !
  subroutine read_line(fd, line, end_of_file)
    implicit none
    integer, intent(in) :: fd
    character(len=256), intent(out) :: line
    logical, optional, intent(out) :: end_of_file

    integer             :: ios
    logical :: tend

    !print*, "in read_line", fd

    tend = .false.
    101 read(unit=fd,fmt='(A256)',END=111, iostat=ios) line
    if (ios /= 0) then
       print*, " Reading Problem..."; stop; endif
    if(line == ' ' .or. line(1:1) == '#') go to 101
    go to 105
    111     tend = .true.
    !print*,"read_line", line
    go to 105
    105   continue

    if( present(end_of_file)) then
      end_of_file = tend
    endif
  end subroutine read_line


  !......................................................................................
  !> @brief 
  !!   Receive the keyword of Engine which contains the engine name and 
  !!   type of unit. Maybe we can also define the units for the output
  !
  !> @note
  !!   Important to know:
  !!   Hessian, Force, Position, Time are exchange with Engine
  !!   Energy is converted only for the ouput
  !!   Mass is needed for the fire integration. Defined in Ry can 
  !!   change depending the unit used.
  !
  !> WARNING: The mass in LJ is 1 but can be defined by the user so
  !!  we should take care about this
  !
  !> @param[in,out]  txt Name of the Engine 
  !
  subroutine make_units( txt )
    ! -- Arguments
    character(*), intent( inout ) :: txt
    ! -- Local variables
    character(:), allocatable :: engine, mode, words(:)
    integer :: n 

    logical :: verbose
    verbose = .true.
    !verbose = .false.


    ! ...Extract the Keyword from the engine_units
    engine = ""; mode = ""
    n = parser( trim(txt), "/",  words )
    if( n >= 1 )engine = lower( trim(words(1)) )
    if( n > 1 )mode = lower( trim(words(2)) )
    


    ! ...Initialization

    E2au = 1.0_DP
    au2E = 1.0_DP
    L2au = 1.0_DP
    au2L = 1.0_DP
    T2au = 1.0_DP
    au2T = 1.0_DP
    M2au = 1.0_DP
    au2M = 1.0_DP

    F2au = 1.0_DP
    au2F = 1.0_DP
    H2au = 1.0_DP
    au2H = 1.0_DP


    ! ...Select the units as function of engine and mode

    select case( engine )


      ! ---------------------------------------------- QE
      case( 'qe', 'quantum_espresso' )

        !! Energy: Rydberg   
        E2au = 1. !/ Ry2eV
        au2E = 1. !  Ry2eV

        !! Length: Bohr
        L2au = 1. ! / B2A
        au2L = 1. !  B2A

        !! Time: aut(Ry)
        T2au = 1.
        au2T = 1.

        !! Mass: au(Ry) AMU/2
        Mass = AMU_RY

        !! Force: Ry/au
        F2au = 1. !/ au2E / L2au
        au2F = 1. !/ F2au

        !! Hessian
        H2au = 1.0_DP 
        au2H = 1.0_DP 

        cE = "Ry"  ! "Ry"
        cL = "a.u." ! "bohr"
        !strg_units = '(27X, "[Ry]",17X,"-----------[Ry/a.u.]----------",3X,"Ry/a.u.^2")'

      ! ---------------------------------------------- LAMMPS
      case( 'lammps' )

        select case( mode )

          case( 'metal' )

            !! Energy: eV
            E2au = 1.0_DP / Ry2eV
            au2E = Ry2eV
 
            !! Length: Angstrom
            L2au = 1.0_DP / B2A
            au2L = B2A

            !! Time: picosecond
            T2au = 1.0_DP / AU_PS
            au2T = AU_PS

            !! Mass: gram/mol
            Mass = AMU_RY

            !! Force
            F2au = E2au / L2au
            au2F = 1.0_DP / F2au

            !! Hessian
            H2au = F2au / L2au
            au2H = 1.0_DP / H2au

            cE = "eV"
            !cL = AA
            cL = "Ang"

          case( 'lj' )
            !! Energy: 1
            E2au = 1.0_DP
            au2E = 1.0_DP
            !! Length: 1
            L2au = 1.0_DP
            au2L = 1.0_DP
            !! Mass: 1
            Mass = 1.0_DP
            !! Time: 1
            T2au = 1.0_DP
            au2T = 1.0_DP
            !! Force
            F2au = E2au / L2au
            au2F = 1.0_DP / F2au

            !! Hessian
            H2au = F2au / L2au
            au2H = 1.0_DP / H2au

            cE = "LJ"
            cL = "LJ"


          case( 'real' )
            !! Energy: Kcal/mol
            E2au = 1.0_DP / Ry2kcalPmol
            au2E = Ry2kcalPmol
 
            !! Length: Angstrom
            L2au = 1.0_DP / B2A
            au2L = B2A

            !! Time: femtosecond
            T2au = 1.0_DP / AU_FS
            au2T = AU_FS

            !! Mass: gram/mol
            Mass = AMU_RY

            !! Force
            F2au = E2au / L2au
            au2F = 1.0_DP / F2au

            !! Hessian
            H2au = F2au / L2au
            au2H = 1.0_DP / H2au

            cE = "Kcal/mol"
            !cL = AA
            cL = "Ang"

          !case( 'si' )
            !! Energy: J
            !! Length: metre
            !! Time: second
          !case( 'cgs' )
            !! Energy: ergs
            !! Length: cm
            !! Time: second
          !case( 'electron' )
            !! Energy: Hatree
            !! Length: Bohr
            !! Time: femtosecond
          !case( 'micro' )
            !! Energy: picogram-micrometer^2/microsecond^2
            !! Length: micrometer
            !! Time: microsecond
          !case( 'nano' )
            !! Energy: attogram-nanometer^2/nanosecond^2
            !! Length: nanometer
            !! Time: nanosecond

          case default
            print*, " * ARTn::WARNING::LAMMPS/mode not defined "

        end select


      ! ---------------------------------------------- OTHER
      case default
        print*, " * ARTn::WARNING::Engine not defined "

    end select   


    ! ...Define the output units string
    strg_units = '(27X, "['//cE//']",31X,"-----------['//cE//'/'//   &
                  cL//']-----------",2X,"['//cE//'/'//cL//to2//']   ['//cL//']")'

    if( verbose )then
      write(*,*) repeat("-",50)
      write(*,1) " * ARTn::UNITS::E2au::", E2au, "au2E", au2E
      write(*,1) " * ARTn::UNITS::L2au::", L2au, "au2L", au2L
      write(*,1) " * ARTn::UNITS::T2au::", T2au, "au2T", au2T
      write(*,1) " * ARTn::UNITS::F2au::", F2au, "au2F", au2F
      write(*,1) " * ARTn::UNITS::H2au::", H2au, "au2H", au2H
      write(*,1) " * ARTn::UNITS::Mass::", Mass
      write(*,*) repeat("-",50)
      1 format(*(x,a,x,g15.5))
    endif


  end subroutine make_units




  !......................................................................................
  ! FORCE

  !> @brief Convert the engine force to a.u.
  !> @param[in] f   force in engine unit
  !> @return a force in atomic units a.u.
  elemental pure function convert_force( f )result( fau )
    real(DP), intent( in ) :: f
    real(DP) :: fau
    fau = f * F2au
  end function

    !> @brief Convert the force in a.u. in engine units
    !> @param [in] fau   force in a.u.
    !> @return a force in engine units
  elemental pure function unconvert_force( fau )result( f )
    real(DP), intent( in ) :: fau
    real(DP) :: f
    f = fau * au2F
  end function



  !......................................................................................
  ! HESSIAN

    !> @brief Convert the engine hessian to a.u.
    !> @param [in] h   hessian in engine unit
    !> @return a hessain in atomic units a.u.
  elemental pure function convert_hessian( h )result( hau )
    real(DP), intent( in ) :: h
    real(DP) :: hau
    hau = h * H2au
  end function

    !> @brief Convert the force in a.u. in engine units
    !> @param [in] hau   force in a.u.
    !> @return a force in engine units
  elemental pure function unconvert_hessian( hau )result( h )
    real(DP), intent( in ) :: hau
    real(DP) :: h
    h = hau * au2H
  end function



  !......................................................................................
  ! LENGTH

    !> @brief Convert the engine length to a.u.
    !> @param [in] p   position in engine unit
    !> @return a position in a.u.
  elemental pure function convert_length( p )result( pau )
    real(DP), intent( in ) :: p
    real(DP) :: pau
    pau = p * L2au
  end function convert_length

    !> @brief Convert the a.u. length to engine unit
    !> @param [in] pau   position in a.u.
    !> @return  position in engine units
  elemental pure function unconvert_length( pau )result( p )
    real(DP), intent( in ) :: pau
    real(DP) :: p
    p = pau * au2L
  end function unconvert_length



  !......................................................................................
  ! ENERGY
  
    !> @brief Convert the engine energy to a.u.
    !> @param [in] e   enegy in engine unit
    !> @return an energy in a.u.
  elemental pure function convert_energy( e )result( eau )
    real(DP), intent( in ) :: e
    real(DP) :: eau
    eau = e * E2au
  end function convert_energy

    !> @brief Convert the a.u. energy to engine unit
    !> @param [in] eau   energy in a.u.
    !> @return an energy in engine units
  elemental pure function unconvert_energy( eau )result( e )
    real(DP), intent( in ) :: eau
    real(DP) :: e 
    e = eau * au2E
  end function unconvert_energy



  !......................................................................................
  ! TIME

    !> @brief Convert the engine time to a.u.
    !> @param [in] t   time in engine unit
    !> @return a time in a.u.
  elemental pure function convert_time( t )result( aut )
    real(DP), intent( in ) :: t
    real(DP) :: aut
    aut = t * T2au
  end function convert_time

    !> @brief Convert the a.u. TIME to engine unit
    !> @param [in] aut   time in a.u.
    !> @return a time in engine units
  elemental pure function unconvert_time( aut )result( t )
    real(DP), intent( in ) :: aut
    real(DP) :: t
    t = aut * au2T
  end function unconvert_time


  !......................................................................................
  ! Return UNIT
  !> @brief Return the unit in character of the quantity received 
  !> @param[in] quantity   (length, energy, force or hessian)
  !> @return correct unit in character
  function unit_char( quantity )result( uchar )
    character(*), intent(in) :: quantity
    character(:), allocatable :: uchar

    select case( quantity )
      case( 'length' );  uchar = cL
      case( 'energy' );  uchar = cE
      case( 'force' );   uchar = cE//'/'//cL
      case( 'hessian' ); uchar = cE//'/'//cL//to2
    end select 

  end function



end module units

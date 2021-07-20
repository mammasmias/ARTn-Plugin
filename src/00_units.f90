

module units
  !> @brief UNITS module contains all the tool to reconize the Engine and its units
  !! to convert the energy/force/length/time in atomic units
  !! Atomic Units (au) in plugin-ARTn is the Rydberg-borh-aut

  PRIVATE

  PUBLIC :: make_units,   &
            convert_length, unconvert_length,   &
            convert_force, unconvert_force,     &
            convert_energy, unconvert_energy,   &
            convert_time, unconvert_time
  
  !PUBLIC :: E2au, au2E, L2au, au2L, T2au, au2T

  INTEGER, PARAMETER ::  DP = selected_real_kind(14,200) ! double precision
  REAL(DP), PARAMETER :: PI     = 3.14159265358979323846_DP ! pi 

  REAL(DP), PARAMETER :: H_PLANCK_SI      = 6.62607015E-34_DP      ! J s
  REAL(DP), PARAMETER :: K_BOLTZMANN_SI   = 1.380649E-23_DP        ! J K^-1 
  REAL(DP), PARAMETER :: ELECTRON_SI      = 1.602176634E-19_DP     ! C
  REAL(DP), PARAMETER :: ELECTRONVOLT_SI  = 1.602176634E-19_DP     ! J  
  REAL(DP), PARAMETER :: ELECTRONMASS_SI  = 9.1093837015E-31_DP    ! Kg
  REAL(DP), PARAMETER :: HARTREE_SI       = 4.3597447222071E-18_DP ! J
  REAL(DP), PARAMETER :: RYDBERG_SI       = HARTREE_SI/2.0_DP      ! J
  REAL(DP), PARAMETER :: BOHR_RADIUS_SI   = 0.529177210903E-10_DP  ! m
  REAL(DP), PARAMETER :: AMU_SI           = 1.66053906660E-27_DP  ! Kg
  REAL(DP), PARAMETER :: C_SI             = 2.99792458E+8_DP    ! m sec^-1


  REAL(DP), PARAMETER :: RY2EV =  13.605691930242388_DP ! Ry to eV conversion 
  REAL(DP), PARAMETER :: B2A =  0.529177210903_DP ! bohr to angstrom conversion
  REAL(DP), PARAMETER :: AMU_RY2 = 911.44424310865645_DP ! calculated from QE using DP
  REAL(DP), PARAMETER :: ps2aut = 41341.374575751 / 2.
  REAL(DP), PARAMETER :: aut2s = 4.8278E-17_DP  !(Ry atomic unit)

  REAL(DP), PARAMETER :: AMU_AU           = AMU_SI / ELECTRONMASS_SI
  REAL(DP), PARAMETER :: AMU_RY           = AMU_AU / 2.0_DP

  !REAL(DP), PARAMETER :: AU_SEC           = H_PLANCK_SI/(2.*pi)/HARTREE_SI
  REAL(DP), PARAMETER :: AU_SEC           = H_PLANCK_SI/(2.*pi)/RYDBERG_SI
  REAL(DP), PARAMETER :: AU_PS            = AU_SEC * 1.0E+12_DP


  !! Units convertor
  REAL(DP) :: E2au, au2E, L2au, au2L, T2au, au2T, F2au, au2F

 contains

  !......................................................................................
  subroutine make_units( txt )
    !> @brief Receive the keyword of Engine which contains the engine name and 
    !! type of unit. Maybe we can also define the units for the output
    character(*), intent( in ) :: txt

    character(:), allocatable :: engine, mode
    integer :: i, n 


    ! ...Extract the Keyword from the engine_units

    n = LEN_TRIM(txt)
    i = SCAN(trim(txt), "/" )
    if( i /= 0.and. i < n )then
      engine = lower(trim(txt(1:i-1)))
      mode = lower(trim(txt(i+1:)))
    else if( i /= 0.and. i == n )then
      engine = lower(trim(txt(1:i-1)))
      mode = ""
    else
      engine = lower(trim(txt))
      mode = ""
    endif

   !print*, " * ARTn::UNITS::engine :", engine
   !print*, " * ARTn::UNITS::mode   :", mode
   !print*, " * ARTn::UNITS::AMU_RY  :", AMU_RY, AMU_RY2
   !print*, " * ARTn::UNITS::AU_SEC  :", AU_SEC, AU_PS, 1./ps2aut


    ! ...Initialization

    E2au = 1.0_DP
    au2E = 1.0_DP
    L2au = 1.0_DP
    au2L = 1.0_DP
    T2au = 1.0_DP
    au2T = 1.0_DP
    F2au = 1.0_DP
    au2F = 1.0_DP


    ! ...Select the units as function of engine and mode

    select case( engine )


      ! ---------------------------------------------- QE
      case( 'qe', 'quantum_espresso' )

        !! Energy: Rydberg
        E2au = 1. / Ry2eV
        au2E = Ry2eV

        !! Length: Bohr
        L2au = 1. / B2A
        au2L = B2A

        !! Time: aut(Ry)
        T2au = 1.
        au2T = 1.


      ! ---------------------------------------------- LAMMPS
      case( 'lammps' )

        select case( mode )

          case( 'metal' )

            !! Energy: eV
            E2au = 1. / Ry2eV
            au2E = Ry2eV
 
            !! Length: Angstrom
            L2au = 1. / B2A
            au2L = B2A

            !! Time: picosecond
            T2au = 1. / AU_PS
            au2T = AU_PS

          !case( 'real' )
            !! Energy: Kcal/mol
            !! Length: Angstrom
            !! Time: femtosecond
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
            print*, " * ARTn::WARNING::LAMMPS/mode not define "

        end select


      ! ---------------------------------------------- OTHER
      case default
        print*, " * ARTn::WARNING::Engine not define "

    end select   


    F2au = E2au / L2au
    au2F = 1. / F2au


    !! if( verbose )
    print*, " * ARTn::UNITS::E2au::", E2au, au2E
    print*, " * ARTn::UNITS::L2au::", L2au, au2L
    print*, " * ARTn::UNITS::T2au::", T2au, au2T
    print*, " * ARTn::UNITS::F2au::", F2au, au2F


   contains
    !................................................................................
    function lower( s1 )result( s2 )
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

  end subroutine make_units




  !......................................................................................
  subroutine units_convert_force( nat, tab )
    integer, intent( in ) :: nat
    real, intent( inout ) :: tab(3,nat)

    tab = tab * F2au 

  end subroutine units_convert_force

  elemental pure function convert_force( f )result( fau )
    real(DP), intent( in ) :: f
    real(DP) :: fau
    fau = f * F2au
  end function

  elemental pure function unconvert_force( fau )result( f )
    real(DP), intent( in ) :: fau
    real(DP) :: f
    f = fau * au2F
  end function



  !......................................................................................
  subroutine units_convert_length( nat, tab )
    integer, intent( in ) :: nat
    real, intent( inout ) :: tab(3,nat)
    tab = tab * L2au
  end subroutine units_convert_length

  elemental pure function convert_length( p )result( pau )
    real(DP), intent( in ) :: p
    real(DP) :: pau
    pau = p * L2au
  end function convert_length

  elemental pure function unconvert_length( pau )result( p )
    real(DP), intent( in ) :: pau
    real(DP) :: p
    p = pau * au2L
  end function unconvert_length



  !......................................................................................
  elemental pure function convert_energy( e )result( eau )
    real(DP), intent( in ) :: e
    real(DP) :: eau
    eau = e * E2au
  end function convert_energy

  elemental pure function unconvert_energy( eau )result( e )
    real(DP), intent( in ) :: eau
    real(DP) :: e 
    e = eau * au2E
  end function unconvert_energy



  !......................................................................................
  elemental pure function convert_time( t )result( aut )
    real(DP), intent( in ) :: t
    real(DP) :: aut
    aut = t * T2au
  end function convert_time

  elemental pure function unconvert_time( aut )result( t )
    real(DP), intent( in ) :: aut
    real(DP) :: t
    t = aut * au2T
  end function unconvert_time

end module units


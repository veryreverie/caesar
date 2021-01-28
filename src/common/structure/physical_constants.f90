! ======================================================================
! Global physical constants.
! ======================================================================
module caesar_physical_constants_module
  use caesar_utils_module
  implicit none
  
  private
  
  ! --------------------------------------------------
  ! Values taken from CODATA 2014, from NIST (physics.nist.gov) on 5/5/2017.
  ! --------------------------------------------------
  
  ! Twice Rydberg's constant in eV.
  real(dp), parameter, public :: EV_PER_HARTREE = 27.21138602_dp
  
  ! Boltzmann's constant in eV/K.
  real(dp), parameter, public :: KB_IN_EV_PER_K = 8.6173303e-5_dp
  
  ! Twice Rydberg's constant in J.
  real(dp), parameter, public :: JOULES_PER_HARTREE = 4.359744650e-18_dp
  
  ! Twice Rydberg's constant in inverse cm.
  real(dp), parameter, public :: INVERSE_CM_PER_HARTREE = 2.194746313702e4_dp
  
  ! Bohr radius in Angstrom.
  real(dp), parameter, public :: ANGSTROM_PER_BOHR = 0.52917721067_dp
  
  ! Electron mass in kg.
  real(dp), parameter, public :: KG_PER_ME = 9.10938356E-31_dp
  
  ! Atomic mass unit in kg.
  real(dp), parameter, public :: KG_PER_AMU = 1.660539040E-27_dp
  
  ! --------------------------------------------------
  ! Fundamental constants.
  ! --------------------------------------------------
  
  ! Rydberg's constant in Hartree.
  real(dp), parameter, public :: RYDBERG_PER_HARTREE = 2.0_dp
  
  real(dp), parameter, public :: RYDBERG_MASS_PER_ME = 0.5_dp
  
  ! --------------------------------------------------
  ! Derived values.
  ! --------------------------------------------------
  
  ! Boltzmann's constant in Hartree per Kelvin.
  real(dp), parameter, public :: KB_IN_AU = KB_IN_EV_PER_K / EV_PER_HARTREE
  
  ! AMU in ME.
  real(dp), parameter, public :: AMU_PER_ME = KG_PER_ME / KG_PER_AMU
  
  ! --------------------------------------------------
  ! Atomic symbols.
  ! --------------------------------------------------
  type(String), allocatable, private :: ATOMIC_SYMBOLS(:)
  public :: atomic_symbol_to_number
  public :: atomic_number_to_symbol
contains

subroutine initialise_atomic_symbols()
  implicit none
  
  if (.not. allocated(ATOMIC_SYMBOLS)) then
    ATOMIC_SYMBOLS = [ str('H'),  &
                     & str('He'), &
                     & str('Li'), &
                     & str('Be'), &
                     & str('B'),  &
                     & str('C'),  &
                     & str('N'),  &
                     & str('O'),  &
                     & str('F'),  &
                     & str('Ne'), &
                     & str('Na'), &
                     & str('Mg'), &
                     & str('Al'), &
                     & str('Si'), &
                     & str('P'),  &
                     & str('S'),  &
                     & str('Cl'), &
                     & str('Ar'), &
                     & str('K'),  &
                     & str('Ca'), &
                     & str('Sc'), &
                     & str('Ti'), &
                     & str('V'),  &
                     & str('Cr'), &
                     & str('Mn'), &
                     & str('Fe'), &
                     & str('Co'), &
                     & str('Ni'), &
                     & str('Cu'), &
                     & str('Zn'), &
                     & str('Ga'), &
                     & str('Ge'), &
                     & str('As'), &
                     & str('Se'), &
                     & str('Br'), &
                     & str('Kr'), &
                     & str('Rb'), &
                     & str('Sr'), &
                     & str('Y'),  &
                     & str('Zr'), &
                     & str('Nb'), &
                     & str('Mo'), &
                     & str('Tc'), &
                     & str('Ru'), &
                     & str('Rh'), &
                     & str('Pd'), &
                     & str('Ag'), &
                     & str('Cd'), &
                     & str('In'), &
                     & str('Sn'), &
                     & str('Sb'), &
                     & str('Te'), &
                     & str('I'),  &
                     & str('Xe'), &
                     & str('Cs'), &
                     & str('Ba'), &
                     & str('La'), &
                     & str('Ce'), &
                     & str('Pr'), &
                     & str('Nd'), &
                     & str('Pm'), &
                     & str('Sm'), &
                     & str('Eu'), &
                     & str('Gd'), &
                     & str('Tb'), &
                     & str('Dy'), &
                     & str('Ho'), &
                     & str('Er'), &
                     & str('Tm'), &
                     & str('Yb'), &
                     & str('Lu'), &
                     & str('Hf'), &
                     & str('Ta'), &
                     & str('W'),  &
                     & str('Re'), &
                     & str('Os'), &
                     & str('Ir'), &
                     & str('Pt'), &
                     & str('Au'), &
                     & str('Hg'), &
                     & str('Tl'), &
                     & str('Pb'), &
                     & str('Bi'), &
                     & str('Po'), &
                     & str('At'), &
                     & str('Rn'), &
                     & str('Fr'), &
                     & str('Ra'), &
                     & str('Ac'), &
                     & str('Th'), &
                     & str('Pa'), &
                     & str('U'),  &
                     & str('Np'), &
                     & str('Pu'), &
                     & str('Am'), &
                     & str('Cm'), &
                     & str('Bk'), &
                     & str('Cf'), &
                     & str('Es'), &
                     & str('Fm'), &
                     & str('Md'), &
                     & str('No'), &
                     & str('Lr'), &
                     & str('Rf'), &
                     & str('Db'), &
                     & str('Sg'), &
                     & str('Bh'), &
                     & str('Hs'), &
                     & str('Mt'), &
                     & str('Ds'), &
                     & str('Rg'), &
                     & str('Cn'), &
                     & str('Nh'), &
                     & str('Fl'), &
                     & str('Mc'), &
                     & str('Lv'), &
                     & str('Ts'), &
                     & str('Og')  ]
  endif
end subroutine

impure elemental function atomic_number_to_symbol(input) result(output)
  implicit none
  
  integer, intent(in) :: input
  type(String)        :: output
  
  call initialise_atomic_symbols()
  
  if (input<1 .or. input>size(ATOMIC_SYMBOLS)) then
    call print_line(ERROR//': '//input//' is not a valid atomic number.')
    call err()
  endif
  
  output = ATOMIC_SYMBOLS(input)
end function

impure elemental function atomic_symbol_to_number(input) result(output)
  implicit none
  
  type(String), intent(in) :: input
  integer                  :: output
  
  call initialise_atomic_symbols()
  
  output = first(ATOMIC_SYMBOLS==input, default=0)
end function
end module

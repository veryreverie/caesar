! ======================================================================
! A spectrum of energies.
! ======================================================================
module energy_spectrum_module
  use common_module
  implicit none
  
  private
  
  public :: EnergySpectrum
  public :: EnergySpectra
  public :: ThermodynamicData
  public :: converged
  
  ! A spectrum in a single space.
  type, extends(Stringable) :: EnergySpectrum
    real(dp),         allocatable          :: spectrum(:)
    integer,          allocatable          :: degeneracies(:)
    type(RealMatrix), allocatable, private :: stresses_(:)
  contains
    procedure, public :: min_energy => min_energy_EnergySpectrum
    procedure, public :: max_energy => max_energy_EnergySpectrum
    
    ! I/O.
    procedure, public :: read  => read_EnergySpectrum
    procedure, public :: write => write_EnergySpectrum
  end type
  
  interface EnergySpectrum
    module procedure new_EnergySpectrum
    module procedure new_EnergySpectrum_String
  end interface
  
  ! A product of spectra, spanning multiple spaces.
  type, extends(Stringsable) :: EnergySpectra
    type(EnergySpectrum), allocatable :: spectra(:)
  contains
    procedure, public :: min_energy => min_energy_EnergySpectra
    procedure, public :: max_energy => max_energy_EnergySpectra
    
    ! I/O.
    procedure, public :: read  => read_EnergySpectra
    procedure, public :: write => write_EnergySpectra
  end type
  
  interface EnergySpectra
    module procedure new_EnergySpectra
    module procedure new_EnergySpectra_Strings
    module procedure new_EnergySpectra_StringArray
  end interface
  
  interface size
    module procedure size_EnergySpectrum
    module procedure size_EnergySpectra
  end interface
  
  ! Construct ThermodynamicData from an EnergySpectra.
  interface ThermodynamicData
    module procedure new_ThermodynamicData_EnergySpectrum
    module procedure new_ThermodynamicData_EnergySpectra
  end interface
  
  ! Compare two spectra, checking if all energies are within tolerance.
  ! Returns whether or not the spectra are within tolerance.
  interface converged
    module procedure converged_EnergySpectrum_EnergySpectrum
    module procedure converged_EnergySpectra_EnergySpectra
  end interface
contains

! Constructors and size functions.
function new_EnergySpectrum(spectrum,degeneracies,stresses) result(this)
  implicit none
  
  real(dp),         intent(in)           :: spectrum(:)
  integer,          intent(in), optional :: degeneracies(:)
  type(RealMatrix), intent(in), optional :: stresses(:)
  type(EnergySpectrum)                   :: this
  
  integer :: i
  
  this%spectrum = spectrum
  
  if (present(degeneracies)) then
    if (size(degeneracies)/=size(spectrum)) then
      call print_line(CODE_ERROR//': spectrum and degeneracies do not match.')
      call err()
    endif
    this%degeneracies = degeneracies
  else
    this%degeneracies = [(1,i=1,size(spectrum))]
  endif
  
  if (present(stresses)) then
    if (size(stresses)/=size(spectrum)) then
      call print_line(CODE_ERROR//': spectrum and stresses do not match.')
      call err()
    endif
    this%stresses_ = stresses
  endif
end function

function size_EnergySpectrum(this) result(output)
  implicit none
  
  type(EnergySpectrum), intent(in) :: this
  integer                          :: output
  
  output = size(this%spectrum)
end function

function new_EnergySpectra(spectra) result(this)
  implicit none
  
  type(EnergySpectrum), intent(in) :: spectra(:)
  type(EnergySpectra)              :: this
  
  this%spectra = spectra
end function

function size_EnergySpectra(this) result(output)
  implicit none
  
  type(EnergySpectra), intent(in) :: this
  integer                         :: output
  
  output = size(this%spectra)
end function

! Minimum and maximum energy.
impure elemental function min_energy_EnergySpectrum(this) result(output)
  implicit none
  
  class(EnergySpectrum), intent(in) :: this
  real(dp)                          :: output
  
  output = minval(this%spectrum)
end function

impure elemental function max_energy_EnergySpectrum(this) result(output)
  implicit none
  
  class(EnergySpectrum), intent(in) :: this
  real(dp)                          :: output
  
  output = maxval(this%spectrum)
end function

impure elemental function min_energy_EnergySpectra(this) result(output)
  implicit none
  
  class(EnergySpectra), intent(in) :: this
  real(dp)                         :: output
  
  output = minval(this%spectra%min_energy())
end function

impure elemental function max_energy_EnergySpectra(this) result(output)
  implicit none
  
  class(EnergySpectra), intent(in) :: this
  real(dp)                         :: output
  
  output = maxval(this%spectra%max_energy())
end function

! Construct ThermodynamicData from an EnergySpectrum or EnergySpectra.
impure elemental function new_ThermodynamicData_EnergySpectrum( &
   & thermal_energy,spectrum) result(output)
  implicit none
  
  real(dp),             intent(in) :: thermal_energy
  type(EnergySpectrum), intent(in) :: spectrum
  type(ThermodynamicData)          :: output
  
  real(dp), allocatable :: energies(:)
  
  integer :: i,j
  
  energies = [real::]
  do i=1,size(spectrum)
    do j=1,spectrum%degeneracies(i)
      energies = [energies, spectrum%spectrum(i)]
    enddo
  enddo
  
  output = ThermodynamicData(thermal_energy, energies)
end function

impure elemental function new_ThermodynamicData_EnergySpectra(thermal_energy, &
   & spectra) result(output)
  implicit none
  
  real(dp),            intent(in) :: thermal_energy
  type(EnergySpectra), intent(in) :: spectra
  type(ThermodynamicData)         :: output
  
  output = sum(ThermodynamicData(thermal_energy, spectra%spectra))
end function

! Compare two spectra, up to a convergence threshold.
impure elemental function converged_EnergySpectrum_EnergySpectrum(this,that, &
   & convergence_threshold) result(output)
  implicit none
  
  type(EnergySpectrum), intent(in) :: this
  type(EnergySpectrum), intent(in) :: that
  real(dp),             intent(in) :: convergence_threshold
  logical                          :: output
  
  output = all(abs(this%spectrum-that%spectrum)<convergence_threshold)
end function

impure elemental function converged_EnergySpectra_EnergySpectra(this,that, &
   & convergence_threshold) result(output)
  implicit none
  
  type(EnergySpectra), intent(in) :: this
  type(EnergySpectra), intent(in) :: that
  real(dp),            intent(in) :: convergence_threshold
  logical                         :: output
  
  output = all(converged(this%spectra,that%spectra,convergence_threshold))
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_EnergySpectrum(this,input)
  implicit none
  
  class(EnergySpectrum), intent(out) :: this
  type(String),          intent(in)  :: input
  
  type(String), allocatable :: line(:)
  
  real(dp), allocatable :: spectrum(:)
  integer,  allocatable :: degeneracies(:)
  
  integer :: i,ialloc
  
  ! TODO
  
  select type(this); type is(EnergySpectrum)
    line = split_line(input)
    allocate( spectrum(size(line)),     &
            & degeneracies(size(line)), &
            & stat=ialloc); call err(ialloc)
    do i=1,size(line)
      spectrum(i) = dble(line(2))
      degeneracies(i) = int(line(1))
    enddo
    this = EnergySpectrum(spectrum, degeneracies)
  class default
    call err()
  end select
end subroutine

function write_EnergySpectrum(this) result(output)
  implicit none
  
  class(EnergySpectrum), intent(in) :: this
  type(String)                      :: output
  
  integer :: i
  
  select type(this); type is(EnergySpectrum)
    output = join([( this%degeneracies(i)//'x'//this%spectrum(i), &
                   & i=1,                                         &
                   & size(this)                                   )])
  class default
    call err()
  end select
end function

impure elemental function new_EnergySpectrum_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(EnergySpectrum)     :: this
  
  call this%read(input)
end function

subroutine read_EnergySpectra(this,input)
  implicit none
  
  class(EnergySpectra), intent(out) :: this
  type(String),         intent(in)  :: input(:)
  
  type(String), allocatable :: line(:)
  
  type(EnergySpectrum), allocatable :: spectra(:)
  
  integer :: i,ialloc
  
  select type(this); type is(EnergySpectra)
    allocate(spectra(size(input)), stat=ialloc); call err(ialloc)
    do i=1,size(input)
      line = split_line(input(i))
      spectra(i) = EnergySpectrum(dble(line(4:)))
    enddo
    this = EnergySpectra(spectra)
  class default
    call err()
  end select
end subroutine

function write_EnergySpectra(this) result(output)
  implicit none
  
  class(EnergySpectra), intent(in) :: this
  type(String), allocatable        :: output(:)
  
  integer :: i,ialloc
  
  select type(this); type is(EnergySpectra)
    allocate(output(size(this)), stat=ialloc); call err(ialloc)
    do i=1,size(this)
      output(i) = 'Spectrum '//left_pad(i,str(size(this)))//' : '// &
                & this%spectra(i)
    enddo
  class default
    call err()
  end select
end function

function new_EnergySpectra_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(EnergySpectra)      :: this
  
  call this%read(input)
end function

impure elemental function new_EnergySpectra_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(EnergySpectra)           :: this
  
  this = EnergySpectra(str(input))
end function
end module

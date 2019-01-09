! ======================================================================
! A spectrum of energies.
! ======================================================================
module energy_spectrum_module
  use common_module
  implicit none
  
  private
  
  public :: EnergySpectrum
  public :: EnergySpectra
  
  ! A spectrum in a single space.
  type, extends(Stringable) :: EnergySpectrum
    real(dp), allocatable :: spectrum(:)
  contains
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
contains

! Constructors and size functions.
function new_EnergySpectrum(spectrum) result(this)
  implicit none
  
  real(dp), intent(in) :: spectrum(:)
  type(EnergySpectrum) :: this
  
  this%spectrum = spectrum
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

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_EnergySpectrum(this,input)
  implicit none
  
  class(EnergySpectrum), intent(out) :: this
  type(String),          intent(in)  :: input
  
  select type(this); type is(EnergySpectrum)
    this = EnergySpectrum(dble(split_line(input)))
  class default
    call err()
  end select
end subroutine

function write_EnergySpectrum(this) result(output)
  implicit none
  
  class(EnergySpectrum), intent(in) :: this
  type(String)                      :: output
  
  select type(this); type is(EnergySpectrum)
    output = join(this%spectrum)
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

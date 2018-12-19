! ======================================================================
! A class for holding calculated energy and forces.
! ======================================================================
module electronic_structure_data_module
  use utils_module
  
  use structure_module
  use normal_mode_module
  
  use linear_response_module
  implicit none
  
  private
  
  public :: LinearResponse
  
  public :: ElectronicStructure
  
  type, extends(Stringsable) :: ElectronicStructure
    real(dp)                          :: energy
    type(CartesianForce)              :: forces
    type(RealMatrix),     allocatable :: virial
    type(LinearResponse), allocatable :: linear_response
  contains
    procedure, public :: read  => read_ElectronicStructure
    procedure, public :: write => write_ElectronicStructure
  end type
  
  interface ElectronicStructure
    module procedure new_ElectronicStructure
    module procedure new_ElectronicStructure_Strings
    module procedure new_ElectronicStructure_StringArray
  end interface
contains

! Constructor.
function new_ElectronicStructure(energy,forces,virial,linear_response) &
   & result(this)
  implicit none
  
  real(dp),             intent(in)           :: energy
  type(CartesianForce), intent(in)           :: forces
  type(RealMatrix),     intent(in), optional :: virial
  type(LinearResponse), intent(in), optional :: linear_response
  type(ElectronicStructure)                  :: this
  
  this%energy  = energy
  this%forces  = forces
  if (present(virial)) then
    this%virial = virial
  endif
  if (present(linear_response)) then
    this%linear_response = linear_response
  endif
end function

! I/O.
subroutine read_ElectronicStructure(this,input)
  implicit none
  
  class(ElectronicStructure), intent(out) :: this
  type(String),               intent(in)  :: input(:)
  
  real(dp)             :: energy
  type(CartesianForce) :: forces
  type(RealMatrix)     :: virial
  type(LinearResponse) :: linear_response
  
  integer :: virial_line
  integer :: linear_response_line
  integer :: i
  
  select type(this); type is(ElectronicStructure)
    ! Check if the virial is present.
    virial_line = 0
    do i=1,size(input)
      if (input(i)=='Virial') then
        virial_line = i
      endif
    enddo
    
    ! Check if linear response is present.
    linear_response_line = 0
    do i=1,size(input)
      if (input(i)=='Permittivity') then
        linear_response_line = i
      endif
    enddo
    
    if (virial_line==0 .and. linear_response_line==0) then
      energy = dble(input(2))
      forces = CartesianForce(input(4:))
      this = ElectronicStructure(energy,forces)
    elseif (linear_response_line==0) then
      energy = dble(input(2))
      forces = CartesianForce(input(4:virial_line-1))
      virial = RealMatrix(input(virial_line+1:))
      this = ElectronicStructure(energy,forces,virial)
    elseif (virial_line==0) then
      energy = dble(input(2))
      forces = CartesianForce(input(4:linear_response_line-1))
      linear_response = LinearResponse(input(linear_response_line:))
      this = ElectronicStructure(energy,forces,linear_response=linear_response)
    else
      energy = dble(input(2))
      forces = CartesianForce(input(4:virial_line-1))
      virial = RealMatrix(input(virial_line+1:linear_response_line-1))
      linear_response = LinearResponse(input(linear_response_line:))
      this = ElectronicStructure(energy,forces,virial,linear_response)
    endif
  class default
    call err()
  end select
end subroutine

function write_ElectronicStructure(this) result(output)
  implicit none
  
  class(ElectronicStructure), intent(in) :: this
  type(String), allocatable              :: output(:)
  
  select type(this); type is(ElectronicStructure)
    output = [ str('Energy (Hartree):'),      &
             & str(this%energy),              &
             & str('Forces (Hartree/Bohr):'), &
             & str(this%forces)               ]
    if (allocated(this%virial)) then
      output = [output, str('Virial'), str(this%virial)]
    endif
    if (allocated(this%linear_response)) then
      output = [output, str(this%linear_response)]
    endif
  class default
    call err()
  end select
end function

function new_ElectronicStructure_Strings(input) result(this)
  implicit none
  
  type(String), intent(in)  :: input(:)
  type(ElectronicStructure) :: this
  
  call this%read(input)
end function

impure elemental function new_ElectronicStructure_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(ElectronicStructure)     :: this
  
  this = ElectronicStructure(str(input))
end function
end module

! ======================================================================
! A class for holding calculated energy and forces.
! ======================================================================
module caesar_electronic_structure_data_module
  use caesar_utils_module
  
  use caesar_structure_module
  use caesar_normal_mode_module
  
  use caesar_linear_response_module
  implicit none
  
  private
  
  public :: LinearResponse
  
  public :: ElectronicStructure
  
  type, extends(Stringsable) :: ElectronicStructure
    real(dp),                            private :: energy_
    type(CartesianForce),   allocatable, private :: forces_
    type(CartesianHessian), allocatable, private :: hessian_
    type(RealMatrix),       allocatable, private :: stress_
    type(LinearResponse),   allocatable, private :: linear_response_
  contains
    procedure, public :: energy => energy_ElectronicStructure
    procedure, public :: has_forces => has_forces_ElectronicStructure
    procedure, public :: forces => forces_ElectronicStructure
    procedure, public :: has_hessian => has_hessian_ElectronicStructure
    procedure, public :: hessian => hessian_ElectronicStructure
    procedure, public :: has_stress => has_stress_ElectronicStructure
    procedure, public :: stress => stress_ElectronicStructure
    procedure, public :: has_linear_response => &
                       & has_linear_response_ElectronicStructure
    procedure, public :: linear_response => &
                       & linear_response_ElectronicStructure
    ! I/O.
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
function new_ElectronicStructure(energy,forces,hessian,stress, &
   & linear_response) result(this)
  implicit none
  
  real(dp),               intent(in)           :: energy
  type(CartesianForce),   intent(in), optional :: forces
  type(CartesianHessian), intent(in), optional :: hessian
  type(RealMatrix),       intent(in), optional :: stress
  type(LinearResponse),   intent(in), optional :: linear_response
  type(ElectronicStructure)                    :: this
  
  this%energy_ = energy
  
  if (present(forces)) then
    this%forces_ = forces
  endif
  
  if (present(hessian)) then
    this%hessian_ = hessian
  endif
  
  if (present(stress)) then
    this%stress_ = stress
  endif
  
  if (present(linear_response)) then
    this%linear_response_ = linear_response
  endif
end function

! Getters.
impure elemental function energy_ElectronicStructure(this) result(output)
  implicit none
  
  class(ElectronicStructure), intent(in) :: this
  real(dp)                               :: output
  
  output = this%energy_
end function

impure elemental function has_forces_ElectronicStructure(this) result(output)
  implicit none
  
  class(ElectronicStructure), intent(in) :: this
  logical                                :: output
  
  output = allocated(this%forces_)
end function

impure elemental function forces_ElectronicStructure(this) result(output)
  implicit none
  
  class(ElectronicStructure), intent(in) :: this
  type(CartesianForce)                   :: output
  
  if (this%has_forces()) then
    output = this%forces_
  else
    call print_line(ERROR//': Sample result does not contain forces.')
    call err()
  endif
end function

impure elemental function has_hessian_ElectronicStructure(this) result(output)
  implicit none
  
  class(ElectronicStructure), intent(in) :: this
  logical                                :: output
  
  output = allocated(this%hessian_)
end function

impure elemental function hessian_ElectronicStructure(this) result(output)
  implicit none
  
  class(ElectronicStructure), intent(in) :: this
  type(CartesianHessian)                 :: output
  
  if (this%has_hessian()) then
    output = this%hessian_
  else
    call print_line(ERROR//': Sample result does not contain the Hessian.')
    call err()
  endif
end function

impure elemental function has_stress_ElectronicStructure(this) result(output)
  implicit none
  
  class(ElectronicStructure), intent(in) :: this
  logical                                :: output
  
  output = allocated(this%stress_)
end function

impure elemental function stress_ElectronicStructure(this) result(output)
  implicit none
  
  class(ElectronicStructure), intent(in) :: this
  type(RealMatrix)                       :: output
  
  if (this%has_stress()) then
    output = this%stress_
  else
    call print_line(ERROR//': Sample result does not contain stress.')
    call err()
  endif
end function

impure elemental function has_linear_response_ElectronicStructure(this) &
   & result(output)
  implicit none
  
  class(ElectronicStructure), intent(in) :: this
  logical                                :: output
  
  output = allocated(this%linear_response_)
end function

impure elemental function linear_response_ElectronicStructure(this) &
   & result(output)
  implicit none
  
  class(ElectronicStructure), intent(in) :: this
  type(LinearResponse)                   :: output
  
  if (this%has_linear_response()) then
    output = this%linear_response_
  else
    call print_line(ERROR//': Sample result does not contain linear_response.')
    call err()
  endif
end function

! I/O.
subroutine read_ElectronicStructure(this,input)
  implicit none
  
  class(ElectronicStructure), intent(out) :: this
  type(String),               intent(in)  :: input(:)
  
  real(dp)                            :: energy
  type(CartesianForce),   allocatable :: forces
  type(CartesianHessian), allocatable :: hessian
  type(RealMatrix),       allocatable :: stress
  type(LinearResponse),   allocatable :: linear_response
  
  ! Line numbers.
  integer :: forces_start_line
  integer :: forces_end_line
  integer :: hessian_start_line
  integer :: hessian_end_line
  integer :: stress_start_line
  integer :: stress_end_line
  integer :: linear_response_start_line
  integer :: linear_response_end_line
  
  integer :: end_lines(5)
  
  integer :: i
  
  type(String), allocatable :: line(:)
  
  select type(this); type is(ElectronicStructure)
    ! Find the lines which start each section.
    forces_start_line = 0
    hessian_start_line = 0
    stress_start_line = 0
    linear_response_start_line = 0
    do i=1,size(input)
      line = split_line(lower_case(input(i)))
      if (size(line)>=1) then
        if (line(1)=='forces') then
          forces_start_line = i+1
        elseif (line(1)=='hessian') then
          hessian_start_line = i+1
        elseif (line(1)=='stress' .or. line(1)=='virial') then
          stress_start_line = i+1
        elseif (line(1)=='permittivity') then
          linear_response_start_line = i+1
        endif
      endif
    enddo
    
    ! Find the lines which end each section.
    end_lines = [ forces_start_line-2,          &
                & hessian_start_line-2,         &
                & stress_start_line-2,          &
                & linear_response_start_line-2, &
                & size(input)                   ]
    if (forces_start_line/=0) then
      forces_end_line = minval( end_lines,                          &
                              & 1,                                  &
                              & mask = end_lines>=forces_start_line )
    endif
    
    if (hessian_start_line/=0) then
      hessian_end_line = minval( end_lines,                           &
                               & 1,                                   &
                               & mask = end_lines>=hessian_start_line )
    endif
    
    if (stress_start_line/=0) then
      stress_end_line = minval( end_lines,                          &
                              & 1,                                  &
                              & mask = end_lines>=stress_start_line )
    endif
    
    if (linear_response_start_line/=0) then
      linear_response_end_line = minval(                &
         & end_lines,                                   &
         & 1,                                           &
         & mask = end_lines>=linear_response_start_line )
    endif
    
    ! Read file contents.
    energy = dble(input(2))
    if (forces_start_line/=0) then
      forces = CartesianForce(input(forces_start_line:forces_end_line))
    endif
    if (hessian_start_line/=0) then
      hessian = CartesianHessian(input(hessian_start_line:hessian_end_line))
    endif
    if (stress_start_line/=0) then
      stress = RealMatrix(input(stress_start_line:stress_end_line))
    endif
    if (linear_response_start_line/=0) then
      ! N.B. the linear response class expects its header line, hence the -1.
      linear_response = LinearResponse(input(                    &
         & linear_response_start_line-1:linear_response_end_line ))
    endif
    
    ! Construct output.
    this = ElectronicStructure( energy,         &
                              & forces,         &
                              & hessian,        &
                              & stress,         &
                              & linear_response )
  class default
    call err()
  end select
end subroutine

function write_ElectronicStructure(this) result(output)
  implicit none
  
  class(ElectronicStructure), intent(in) :: this
  type(String), allocatable              :: output(:)
  
  select type(this); type is(ElectronicStructure)
    output = [ str('Energy (Hartree):'), &
             & str(this%energy_)         ]
    
    if (allocated(this%forces_)) then
      output = [ output,                        &
               & str('Forces (Hartree/Bohr):'), &
               & str(this%forces_)              ]
    endif
    
    if (allocated(this%hessian_)) then
      output = [ output,                           &
               & str('Hessian (Hartree/Bohr^2):'), &
               & str(this%hessian_)                ]
    endif
    
    if (allocated(this%stress_)) then
      output = [ output,                          &
               & str('Stress (Hartree/Bohr^3):'), &
               & str(this%stress_)                ]
    endif
    
    if (allocated(this%linear_response_)) then
      output = [ output,                    &
               & str(this%linear_response_) ]
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

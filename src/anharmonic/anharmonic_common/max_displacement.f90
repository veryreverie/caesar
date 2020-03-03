! ======================================================================
! Generate the frequency of maximum displacement from
!    the maximum energy of displacement,
!    or vice versa.
! ======================================================================
module max_displacement_module
  use common_module
  implicit none
  
  private
  
  public :: MaxDisplacement
  
  type, extends(Stringsable) :: MaxDisplacement
    real(dp) :: maximum_weighted_displacement
    real(dp) :: frequency_of_max_displacement
    real(dp) :: max_energy_of_displacement
  contains
    ! I/O.
    procedure, public :: read  => read_MaxDisplacement
    procedure, public :: write => write_MaxDisplacement
  end type
  
  interface MaxDisplacement
    module procedure new_MaxDisplacement
    module procedure new_MaxDisplacement_displacement
    module procedure new_MaxDisplacement_Strings
    module procedure new_MaxDisplacement_StringArray
  end interface
contains

impure elemental function new_MaxDisplacement(maximum_weighted_displacement, &
   & frequency_of_max_displacement,max_energy_of_displacement) result(this)
  implicit none
  
  real(dp), intent(in), optional :: maximum_weighted_displacement
  real(dp), intent(in), optional :: frequency_of_max_displacement
  real(dp), intent(in), optional :: max_energy_of_displacement
  type(MaxDisplacement)          :: this
  
  if (count([ present(maximum_weighted_displacement), &
            & present(frequency_of_max_displacement), &
            & present(max_energy_of_displacement)     ])<2) then
    call print_line(ERROR//': At most one argument to MaxDisplacement can be &
       &reconstructed.')
    call err()
  endif
  
  ! Use E=0.5(wu)^2 to reconstruct a missing value if needed.
  if (present(maximum_weighted_displacement)) then
    this%maximum_weighted_displacement = maximum_weighted_displacement
  else
    this%maximum_weighted_displacement = sqrt(2*max_energy_of_displacement) &
                                     & / frequency_of_max_displacement
  endif
  
  if (present(frequency_of_max_displacement)) then
    this%frequency_of_max_displacement = frequency_of_max_displacement
  else
    this%frequency_of_max_displacement = sqrt(2*max_energy_of_displacement) &
                                     & / maximum_weighted_displacement
  endif
  
  if (present(max_energy_of_displacement)) then
    this%max_energy_of_displacement = max_energy_of_displacement
  else
    this%max_energy_of_displacement = 0.5_dp                          &
                                  & * ( frequency_of_max_displacement &
                                  &   * maximum_weighted_displacement )**2
  endif
end function

impure elemental function new_MaxDisplacement_displacement(        &
   & maximum_displacement,structure,frequency_of_max_displacement, &
   & max_energy_of_displacement) result(this)
  implicit none
  
  real(dp),            intent(in)           :: maximum_displacement
  type(StructureData), intent(in)           :: structure
  real(dp),            intent(in), optional :: frequency_of_max_displacement
  real(dp),            intent(in), optional :: max_energy_of_displacement
  type(MaxDisplacement)                     :: this
  
  real(dp) :: maximum_weighted_displacement
  
  ! Calculate the maximum mass-weighted displacement from the maximum
  !    displacement. This corresponds to a mode made entirely from the
  !    lightest element moving up to maximum_displacement.
  maximum_weighted_displacement = maximum_displacement &
                              & * sqrt(minval(structure%atoms%mass()))
  
  this = MaxDisplacement( maximum_weighted_displacement, &
                        & frequency_of_max_displacement, &
                        & max_energy_of_displacement     )
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_MaxDisplacement(this,input)
  implicit none
  
  class(MaxDisplacement), intent(out) :: this
  type(String),           intent(in)  :: input(:)
  
  real(dp) :: maximum_weighted_displacement
  real(dp) :: frequency_of_max_displacement
  real(dp) :: max_energy_of_displacement
  
  select type(this); type is(MaxDisplacement)
    maximum_weighted_displacement = dble(token(input(1),5))
    frequency_of_max_displacement = dble(token(input(2),6))
    max_energy_of_displacement = dble(token(input(3),6))
    
    this = MaxDisplacement( maximum_weighted_displacement, &
                          & frequency_of_max_displacement, &
                          & max_energy_of_displacement     )
  class default
    call err()
  end select
end subroutine

function write_MaxDisplacement(this) result(output)
  implicit none
  
  class(MaxDisplacement), intent(in) :: this
  type(String), allocatable          :: output(:)
  
  select type(this); type is(MaxDisplacement)
    output = [ 'Max. weighted displacement     : '// &
             &   this%maximum_weighted_displacement, &
             & 'Frequency of max. displacement : '// &
             &   this%frequency_of_max_displacement, &
             & 'Max. energy of displacement    : '// &
             &       this%max_energy_of_displacement ]
  class default
    call err()
  end select
end function

function new_MaxDisplacement_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(MaxDisplacement)    :: this
  
  call this%read(input)
end function

impure elemental function new_MaxDisplacement_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(MaxDisplacement)         :: this
  
  this = MaxDisplacement(str(input))
end function
end module

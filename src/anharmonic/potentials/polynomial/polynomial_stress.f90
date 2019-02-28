! ======================================================================
! A polynomial representation of the stress.
! ======================================================================
module polynomial_stress_module
  use common_module
  
  use states_module
  use anharmonic_common_module
  
  use coupling_stress_basis_functions_module
  implicit none
  
  private
  
  public :: startup_polynomial_stress
  
  public :: PolynomialStress
  
  type, extends(StressData) :: PolynomialStress
    type(RealMatrix), private :: reference_stress_
    type(CouplingStressBasisFunctions), allocatable, private :: &
       & basis_functions_(:)
  contains
    procedure, public, nopass :: representation => &
                               & representation_PolynomialStress
    
    procedure, public :: zero_stress => zero_stress_PolynomialStress
    
    procedure, public :: stress_RealModeDisplacement => &
                       & stress_RealModeDisplacement_PolynomialStress
    procedure, public :: stress_ComplexModeDisplacement => &
                       & stress_ComplexModeDisplacement_PolynomialStress
    
    procedure, public :: braket => braket_PolynomialStress
    
    procedure, public :: harmonic_expectation => &
                       & harmonic_expectation_PolynomialStress
    
    ! I/O.
    procedure, public :: read  => read_PolynomialStress
    procedure, public :: write => write_PolynomialStress
  end type
  
  interface PolynomialStress
    module procedure new_PolynomialStress
    module procedure new_PolynomialStress_Strings
    module procedure new_PolynomialStress_StringArray
  end interface
contains

! Startup procedure.
subroutine startup_polynomial_stress()
  implicit none
  
  type(PolynomialStress) :: stress
  
  call stress%startup()
end subroutine

! Constructor.
function new_PolynomialStress(reference_stress,basis_functions) result(this)
  implicit none
  
  type(RealMatrix),                   intent(in) :: reference_stress
  type(CouplingStressBasisFunctions), intent(in) :: basis_functions(:)
  type(PolynomialStress)                         :: this
  
  this%reference_stress_ = reference_stress
  this%basis_functions_  = basis_functions
end function

! Type representation.
impure elemental function representation_PolynomialStress() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'polynomial'
end function

! Set the undisplaced stress to zero.
impure elemental subroutine zero_stress_PolynomialStress(this)
  implicit none
  
  class(PolynomialStress), intent(inout) :: this
  
  this%reference_stress_ = this%reference_stress_ - this%undisplaced_stress()
end subroutine

! Calculate the stress at a given displacement.
impure elemental function stress_RealModeDisplacement_PolynomialStress(this, &
   & displacement) result(output)
  implicit none
  
  class(PolynomialStress),    intent(in) :: this
  type(RealModeDisplacement), intent(in) :: displacement
  type(RealMatrix)                       :: output
  
  if (size(this%basis_functions_)>0) then
    output = this%reference_stress_ &
         & + sum(this%basis_functions_%stress(displacement))
  else
    output = this%reference_stress_
  endif
end function

impure elemental function stress_ComplexModeDisplacement_PolynomialStress( &
   & this,displacement) result(output)
  implicit none
  
  class(PolynomialStress),       intent(in) :: this
  type(ComplexModeDisplacement), intent(in) :: displacement
  type(ComplexMatrix)                       :: output
  
  if (size(this%basis_functions_)>0) then
    output = cmplxmat(this%reference_stress_) &
         & + sum(this%basis_functions_%stress(displacement))
  else
    output = cmplxmat(this%reference_stress_)
  endif
end function

! Integrate the stress between two states.
function braket_PolynomialStress(this,bra,ket,subspace,subspace_basis, &
   & anharmonic_data) result(output)
  implicit none
  
  class(PolynomialStress), intent(in)           :: this
  class(SubspaceState),       intent(in)           :: bra
  class(SubspaceState),       intent(in), optional :: ket
  type(DegenerateSubspace),   intent(in)           :: subspace
  class(SubspaceBasis),       intent(in)           :: subspace_basis
  type(AnharmonicData),       intent(in)           :: anharmonic_data
  type(StressPointer)                           :: output
  
  type(PolynomialStress) :: stress
  
  logical, allocatable :: to_remove(:)
  
  integer :: i,j,k
  
  stress = this
  
  ! Integrate the reference stress (N.B. <i|e|j> = e<i|j> if e is a scalar.).
  stress%reference_stress_ = stress%reference_stress_ &
                        & * braket( bra,              &
                        &           ket,              &
                        &           subspace,         &
                        &           subspace_basis,   &
                        &           anharmonic_data )
  
  ! Integrate each basis function between the bra and the ket.
  to_remove = [(.false., i=1, size(stress%basis_functions_))]
  do i=1,size(stress%basis_functions_)
    j = first(stress%basis_functions_(i)%coupling%ids==subspace%id, default=0)
    if (j/=0) then
      do k=1,size(stress%basis_functions_(i))
        call stress%basis_functions_(i)%basis_functions(k)%braket( &
                                                 & bra,            &
                                                 & ket,            &
                                                 & subspace,       &
                                                 & subspace_basis, &
                                                 & anharmonic_data )
      enddo
      
      ! Simplify the stress.
      call stress%basis_functions_(i)%basis_functions%simplify()
      
      ! Update the coupling to remove the integrated subspace.
      stress%basis_functions_(i)%coupling%ids = [         &
         & stress%basis_functions_(i)%coupling%ids(:j-1), &
         & stress%basis_functions_(j)%coupling%ids(j+1:)  ]
      
      ! Check if the basis function is now a constant.
      ! If so, add the constant energy to the stress's reference energy,
      !    and flag the term for removal.
      ! Then check if a coupling is now the same as a previous coupling.
      ! If so, combine the two and flag the duplicate term for removal.
      if (size(stress%basis_functions_(i)%coupling)==0) then
        stress%reference_stress_ =           &
           &   stress%reference_stress_      &
           & + sum(stress%basis_functions_(i &
           &          )%basis_functions%undisplaced_stress())
        to_remove(i) = .true.
      else
        do k=1,i-1
          if (    size(stress%basis_functions_(k)%coupling%ids) &
             & == size(stress%basis_functions_(i)%coupling%ids) ) then
            if (all( stress%basis_functions_(k)%coupling%ids &
                & == stress%basis_functions_(i)%coupling%ids )) then
              stress%basis_functions_(k)%basis_functions = [   &
                 & stress%basis_functions_(k)%basis_functions, &
                 & stress%basis_functions_(i)%basis_functions  ]
              to_remove(i) = .true.
              exit
            endif
          endif
        enddo
      endif
    endif
  enddo
  
  ! Remove constant terms.
  stress%basis_functions_ = stress%basis_functions_(filter(.not.to_remove))
  
  output = StressPointer(stress)
end function

! Calculate the thermal expectation of the stress, <stress>, for a set of
!    harmonic states.
function harmonic_expectation_PolynomialStress(this,frequency, &
   & thermal_energy,no_states,subspace,anharmonic_data) result(output)
  implicit none
  
  class(PolynomialStress),  intent(in) :: this
  real(dp),                 intent(in) :: frequency
  real(dp),                 intent(in) :: thermal_energy
  integer,                  intent(in) :: no_states
  type(DegenerateSubspace), intent(in) :: subspace
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(RealMatrix)                     :: output
  
  output = this%reference_stress_                                          &
       & + sum(this%basis_functions_%harmonic_expectation( frequency,      &
       &                                                   thermal_energy, &
       &                                                   no_states,      &
       &                                                   subspace,       &
       &                                                   anharmonic_data ))
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_PolynomialStress(this,input)
  implicit none
  
  class(PolynomialStress), intent(out) :: this
  type(String),               intent(in)  :: input(:)
  
  type(RealMatrix)                                :: reference_stress
  type(CouplingStressBasisFunctions), allocatable :: basis_functions(:)
  
  select type(this); type is(PolynomialStress)
    reference_stress = RealMatrix(input(2:4))
    
    basis_functions = CouplingStressBasisFunctions(split_into_sections( &
                                       & input(8:),                     &
                                       & separating_line=repeat('=',50) ))
    
    this = PolynomialStress( reference_stress, &
                           & basis_functions   )
  class default
    call err()
  end select
end subroutine

function write_PolynomialStress(this) result(output)
  implicit none
  
  class(PolynomialStress), intent(in) :: this
  type(String), allocatable              :: output(:)
  
  select type(this); type is(PolynomialStress)
    output = [ str('Reference stress:'),                                  &
             & str(this%reference_stress_),                               &
             & str('Basis functions:'),                                   &
             & str(''),                                                   &
             & str(this%basis_functions_, separating_line=repeat('=',50)) ]
  class default
    call err()
  end select
end function

function new_PolynomialStress_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(PolynomialStress)   :: this
  
  call this%read(input)
end function

impure elemental function new_PolynomialStress_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(PolynomialStress)        :: this
  
  this = PolynomialStress(str(input))
end function
end module

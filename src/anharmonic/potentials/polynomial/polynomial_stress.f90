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
    
    procedure, public :: braket_SubspaceState  => &
                       & braket_SubspaceState_PolynomialStress
    procedure, public :: braket_BasisState  => &
                       & braket_BasisState_PolynomialStress
    procedure, public :: braket_BasisStates => &
                       & braket_BasisStates_PolynomialStress
    
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
subroutine braket_SubspaceState_PolynomialStress(this,bra,ket, &
   & anharmonic_data)
  implicit none
  
  class(PolynomialStress),  intent(inout)        :: this
  class(SubspaceState),     intent(in)           :: bra
  class(SubspaceState),     intent(in), optional :: ket
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  
  logical, allocatable :: to_remove(:)
  
  integer :: i,j,k
  
  ! Integrate the reference stress (N.B. <i|e|j> = e<i|j> if e is a scalar.).
  this%reference_stress_ = this%reference_stress_ &
                       & * bra%inner_product(ket,anharmonic_data)
  
  ! Integrate each basis function between the bra and the ket.
  do i=1,size(this%basis_functions_)
    call this%basis_functions_(i)%braket(bra,ket,anharmonic_data)
  enddo
  
  ! Check if the basis function is now a constant.
  ! If so, add the constant stress to the potential's reference stress,
  !    and flag the term for removal.
  ! Then check if a coupling is now the same as a previous coupling.
  ! If so, combine the two and flag the duplicate term for removal.
  to_remove = [(.false., i=1, size(this%basis_functions_))]
  do i=1,size(this%basis_functions_)
    if (size(this%basis_functions_(i)%coupling)==0) then
      this%reference_stress_ =      &
         &   this%reference_stress_ &
         & + this%basis_functions_(i)%undisplaced_stress()
      to_remove(i) = .true.
    else
      do j=1,i-1
        if (    this%basis_functions_(i)%coupling &
           & == this%basis_functions_(j)%coupling ) then
          if (to_remove(j)) then
            call err()
          endif
          call this%basis_functions_(j)%append( &
                     & this%basis_functions_(i) )
          to_remove(i) = .true.
        endif
      enddo
    endif
  enddo
  
  ! Remove constant and duplicate terms.
  this%basis_functions_ = this%basis_functions_(filter(.not.to_remove))
end subroutine

! Integrate the stress between two states.
subroutine braket_BasisState_PolynomialStress(this,bra,ket,subspace, &
   & subspace_basis,anharmonic_data)
  implicit none
  
  class(PolynomialStress),  intent(inout)        :: this
  class(BasisState),        intent(in)           :: bra
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  
  logical, allocatable :: to_remove(:)
  
  integer :: i,j,k
  
  ! Integrate the reference stress (N.B. <i|e|j> = e<i|j> if e is a scalar.).
  this%reference_stress_ = this%reference_stress_                        &
                       & * subspace_basis%inner_product( bra,            &
                       &                                 ket,            &
                       &                                 subspace,       &
                       &                                 anharmonic_data )
  
  ! Integrate each basis function between the bra and the ket.
  do i=1,size(this%basis_functions_)
    call this%basis_functions_(i)%braket( bra,            &
                                        & ket,            &
                                        & subspace,       &
                                        & subspace_basis, &
                                        & anharmonic_data )
  enddo
  
  ! Check if the basis function is now a constant.
  ! If so, add the constant stress to the potential's reference stress,
  !    and flag the term for removal.
  ! Then check if a coupling is now the same as a previous coupling.
  ! If so, combine the two and flag the duplicate term for removal.
  to_remove = [(.false., i=1, size(this%basis_functions_))]
  do i=1,size(this%basis_functions_)
    if (size(this%basis_functions_(i)%coupling)==0) then
      this%reference_stress_ =      &
         &   this%reference_stress_ &
         & + this%basis_functions_(i)%undisplaced_stress()
      to_remove(i) = .true.
    else
      do j=1,i-1
        if (    this%basis_functions_(i)%coupling &
           & == this%basis_functions_(j)%coupling ) then
          if (to_remove(j)) then
            call err()
          endif
          call this%basis_functions_(j)%append( &
                     & this%basis_functions_(i) )
          to_remove(i) = .true.
        endif
      enddo
    endif
  enddo
  
  ! Remove constant and duplicate terms.
  this%basis_functions_ = this%basis_functions_(filter(.not.to_remove))
end subroutine

subroutine braket_BasisStates_PolynomialStress(this,states,subspace, &
   & subspace_basis,anharmonic_data)
  implicit none
  
  class(PolynomialStress),  intent(inout) :: this
  class(BasisStates),       intent(in)    :: states
  type(DegenerateSubspace), intent(in)    :: subspace
  class(SubspaceBasis),     intent(in)    :: subspace_basis
  type(AnharmonicData),     intent(in)    :: anharmonic_data
  
  logical, allocatable :: to_remove(:)
  
  integer :: i,j,k
  
  ! Integrate each basis function between the bra and the ket.
  do i=1,size(this%basis_functions_)
    call this%basis_functions_(i)%braket( states,         &
                                        & subspace,       &
                                        & subspace_basis, &
                                        & anharmonic_data )
  enddo
  
  ! Check if the basis function is now a constant.
  ! If so, add the constant stress to the potential's reference stress,
  !    and flag the term for removal.
  ! Then check if a coupling is now the same as a previous coupling.
  ! If so, combine the two and flag the duplicate term for removal.
  to_remove = [(.false., i=1, size(this%basis_functions_))]
  do i=1,size(this%basis_functions_)
    if (size(this%basis_functions_(i)%coupling)==0) then
      this%reference_stress_ =      &
         &   this%reference_stress_ &
         & + this%basis_functions_(i)%undisplaced_stress()
      to_remove(i) = .true.
    else
      do j=1,i-1
        if (    this%basis_functions_(i)%coupling &
           & == this%basis_functions_(j)%coupling ) then
          if (to_remove(j)) then
            call err()
          endif
          call this%basis_functions_(j)%append( &
                     & this%basis_functions_(i) )
          to_remove(i) = .true.
        endif
      enddo
    endif
  enddo
  
  ! Remove constant and duplicate terms.
  this%basis_functions_ = this%basis_functions_(filter(.not.to_remove))
end subroutine

! Calculate the thermal expectation of the stress, <stress>, for a set of
!    harmonic states.
function harmonic_expectation_PolynomialStress(this,frequency, &
   & thermal_energy,subspace,anharmonic_data) result(output)
  implicit none
  
  class(PolynomialStress),  intent(in) :: this
  real(dp),                 intent(in) :: frequency
  real(dp),                 intent(in) :: thermal_energy
  type(DegenerateSubspace), intent(in) :: subspace
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(RealMatrix)                     :: output
  
  output = this%reference_stress_                                          &
       & + sum(this%basis_functions_%harmonic_expectation( frequency,      &
       &                                                   thermal_energy, &
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

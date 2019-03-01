! ======================================================================
! A set of stress basis functions spanning a given subspace coupling.
! ======================================================================
module coupling_stress_basis_functions_module
  use common_module
  
  use anharmonic_common_module
  
  use basis_function_module
  use stress_basis_function_module
  use sampling_points_module
  implicit none
  
  private
  
  public :: CouplingStressBasisFunctions
  public :: size
  public :: generate_stress_basis_functions
  
  type, extends(Stringsable) :: CouplingStressBasisFunctions
    type(SubspaceCoupling)                          :: coupling
    type(StressBasisFunction), allocatable, private :: basis_functions_(:)
  contains
    procedure, public :: basis_functions => &
                       & basis_functions_CouplingStressBasisFunctions
    
    generic,   public  ::                                          &
       & stress =>                                                 &
       & stress_RealModeDisplacement_CouplingStressBasisFunctions, &
       & stress_ComplexModeDisplacement_CouplingStressBasisFunctions
    procedure, private :: &
       & stress_RealModeDisplacement_CouplingStressBasisFunctions
    procedure, private :: &
       & stress_ComplexModeDisplacement_CouplingStressBasisFunctions
    
    procedure, public :: braket => braket_CouplingStressBasisFunctions
    
    procedure, public :: harmonic_expectation => &
                       & harmonic_expectation_CouplingStressBasisFunctions
    
    procedure, public :: undisplaced_stress => &
                       & undisplaced_stress_CouplingStressBasisFunctions
    
    procedure, public :: append => append_CouplingStressBasisFunctions
    
    ! I/O.
    procedure, public :: read  => read_CouplingStressBasisFunctions
    procedure, public :: write => write_CouplingStressBasisFunctions
  end type
  
  interface CouplingStressBasisFunctions
    module procedure new_CouplingStressBasisFunctions
    module procedure new_CouplingStressBasisFunctions_Strings
    module procedure new_CouplingStressBasisFunctions_StringArray
  end interface
  
  interface size
    module procedure size_CouplingStressBasisFunctions
  end interface
  
  interface generate_stress_basis_functions
    module procedure generate_stress_basis_functions_SubspaceCoupling
  end interface
contains

! Constructor and size function.
function new_CouplingStressBasisFunctions(coupling,basis_functions) &
   & result(this)
  implicit none
  
  type(SubspaceCoupling),    intent(in) :: coupling
  type(StressBasisFunction), intent(in) :: basis_functions(:)
  type(CouplingStressBasisFunctions)    :: this
  
  this%coupling = coupling
  this%basis_functions_ = basis_functions
end function

function size_CouplingStressBasisFunctions(this) result(output)
  implicit none
  
  class(CouplingStressBasisFunctions), intent(in) :: this
  integer                                         :: output
  
  output = size(this%basis_functions_)
end function

function basis_functions_CouplingStressBasisFunctions(this) result(output)
  implicit none
  
  class(CouplingStressBasisFunctions), intent(in) :: this
  type(StressBasisFunction), allocatable          :: output(:)
  
  output = this%basis_functions_
end function

impure elemental function                                           &
   & stress_RealModeDisplacement_CouplingStressBasisFunctions(this, &
   & displacement) result(output)
  implicit none
  
  class(CouplingStressBasisFunctions), intent(in) :: this
  class(RealModeDisplacement),         intent(in) :: displacement
  type(RealMatrix)                                :: output
  
  if (size(this%basis_functions_)>0) then
    output = sum(this%basis_functions_%stress(displacement))
  else
    output = dblemat(zeroes(3,3))
  endif
end function

impure elemental function                                              &
   & stress_ComplexModeDisplacement_CouplingStressBasisFunctions(this, &
   & displacement) result(output)
  implicit none
  
  class(CouplingStressBasisFunctions), intent(in) :: this
  class(ComplexModeDisplacement),      intent(in) :: displacement
  type(ComplexMatrix)                             :: output
  
  if (size(this%basis_functions_)>0) then
    output = sum(this%basis_functions_%stress(displacement))
  else
    output = cmplxmat(zeroes(3,3))
  endif
end function

impure elemental function braket_CouplingStressBasisFunctions(this,bra,ket, &
   & subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(CouplingStressBasisFunctions), intent(in)           :: this
  class(SubspaceState),                intent(in)           :: bra
  class(SubspaceState),                intent(in), optional :: ket
  type(DegenerateSubspace),            intent(in)           :: subspace
  class(SubspaceBasis),                intent(in)           :: subspace_basis
  type(AnharmonicData),                intent(in)           :: anharmonic_data
  type(CouplingStressBasisFunctions)                        :: output
  
  integer :: i,j
  
  output = this
  
  i = first(output%coupling%ids==subspace%id, default=0)
  if (i/=0) then
    do j=1,size(output)
      call output%basis_functions_(j)%braket( bra,            &
                                            & ket,            &
                                            & subspace,       &
                                            & subspace_basis, &
                                            & anharmonic_data )
    enddo
    
    ! Simplify the potential.
    call output%basis_functions_%simplify()
    
    ! Update the coupling to remove the integrated subspace.
    output%coupling%ids = [ output%coupling%ids(:i-1), &
                          & output%coupling%ids(i+1:)  ]
  endif
end function

impure elemental function harmonic_expectation_CouplingStressBasisFunctions( &
   & this,frequency,thermal_energy,no_states,subspace,anharmonic_data)       &
   & result(output)
  implicit none
  
  class(CouplingStressBasisFunctions), intent(in) :: this
  real(dp),                            intent(in) :: frequency
  real(dp),                            intent(in) :: thermal_energy
  integer,                             intent(in) :: no_states
  type(DegenerateSubspace),            intent(in) :: subspace
  type(AnharmonicData),                intent(in) :: anharmonic_data
  type(RealMatrix)                                :: output
  
  if (size(this%basis_functions_)==0) then
    output = dblemat(zeroes(3,3))
  else
    output = sum(this%basis_functions_%harmonic_expectation( frequency,      &
                                                           & thermal_energy, &
                                                           & no_states,      &
                                                           & subspace,       &
                                                           & anharmonic_data ))
  endif
end function

! Return the stress at zero displacement.
impure elemental function undisplaced_stress_CouplingStressBasisFunctions( &
   & this) result(output)
  implicit none
  
  class(CouplingStressBasisFunctions), intent(in) :: this
  type(RealMatrix)                                :: output
  
  type(RealModeDisplacement) :: zero_displacement
  
  zero_displacement = RealModeDisplacement([RealSingleDisplacement::])
  
  output = this%stress(zero_displacement)
end function

! Append another StressCouplingBasisFunctions to this.
subroutine append_CouplingStressBasisFunctions(this,that)
  implicit none
  
  class(CouplingStressBasisFunctions), intent(inout) :: this
  class(CouplingStressBasisFunctions), intent(in)    :: that
  
  if (this%coupling/=that%coupling) then
    call print_line(CODE_ERROR//': Appending incompatible basis functions.')
    call err()
  endif
  
  this%basis_functions_ = [this%basis_functions_, that%basis_functions_]
end subroutine

! Generate stress basis functions.
function generate_stress_basis_functions_SubspaceCoupling(coupling,     &
   & stress_expansion_order,structure,complex_modes,real_modes,qpoints, &
   & subspaces,degenerate_symmetries,vscf_basis_functions_only,logfile) &
   & result(output)
  implicit none
  
  type(SubspaceCoupling),   intent(in)    :: coupling
  integer,                  intent(in)    :: stress_expansion_order
  type(StructureData),      intent(in)    :: structure
  type(ComplexMode),        intent(in)    :: complex_modes(:)
  type(RealMode),           intent(in)    :: real_modes(:)
  type(QpointData),         intent(in)    :: qpoints(:)
  type(DegenerateSubspace), intent(in)    :: subspaces(:)
  type(DegenerateSymmetry), intent(in)    :: degenerate_symmetries(:)
  logical,                  intent(in)    :: vscf_basis_functions_only
  type(OFile),              intent(inout) :: logfile
  type(CouplingStressBasisFunctions)      :: output
  
  type(SubspaceMonomial), allocatable :: subspace_monomials(:)
  
  type(StressBasisFunction), allocatable :: coupling_basis_functions(:)
  type(StressBasisFunction), allocatable :: monomial_basis_functions(:)
  
  integer :: i
  
  ! Generate the set of subspace monomials corresponding to the subspace
  !    coupling.
  ! e.g. the coupling [1,2] might have monomials [1,2], [1,1,2] and [1,2,2].
  subspace_monomials = generate_subspace_monomials(     &
     & coupling,                                        &
     & subspaces,                                       &
     & minimum_expansion_order = 1,                     &
     & maximum_expansion_order = stress_expansion_order )
  
  ! Loop over the subspace monomials corresponding to the coupling.
  coupling_basis_functions = [StressBasisFunction::]
  do i=1,size(subspace_monomials)
    monomial_basis_functions = generate_stress_basis_functions( &
                                   & subspace_monomials(i),     &
                                   & structure,                 &
                                   & complex_modes,             &
                                   & real_modes,                &
                                   & qpoints,                   &
                                   & subspaces,                 &
                                   & degenerate_symmetries,     &
                                   & vscf_basis_functions_only, &
                                   & logfile                    )
    coupling_basis_functions = [ coupling_basis_functions, &
                               & monomial_basis_functions  ]
  enddo
  
  output = CouplingStressBasisFunctions( coupling,                &
                                       & coupling_basis_functions )
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_CouplingStressBasisFunctions(this,input)
  implicit none
  
  class(CouplingStressBasisFunctions), intent(out) :: this
  type(String),                        intent(in)  :: input(:)
  
  type(String), allocatable :: line(:)
  
  type(SubspaceCoupling)                 :: coupling
  type(StressBasisFunction), allocatable :: basis_functions(:)
  
  select type(this); type is(CouplingStressBasisFunctions)
    line = split_line(input(1))
    coupling = SubspaceCoupling(int(line(3:)))
    
    basis_functions = StressBasisFunction(split_into_sections( &
                              & input(3:),                     &
                              & separating_line=repeat('-',50) ))
    
    this = CouplingStressBasisFunctions( coupling        = coupling,       &
                                       & basis_functions = basis_functions )
  class default
    call err()
  end select
end subroutine

function write_CouplingStressBasisFunctions(this) result(output)
  implicit none
  
  class(CouplingStressBasisFunctions), intent(in) :: this
  type(String), allocatable                       :: output(:)
  
  select type(this); type is(CouplingStressBasisFunctions)
    output = [ 'Subspace Coupling: '//this%coupling,                      &
             & str(repeat('-',50)),                                       &
             & str(this%basis_functions_, separating_line=repeat('-',50)) ]
  class default
    call err()
  end select
end function

function new_CouplingStressBasisFunctions_Strings(input) result(this)
  implicit none
  
  type(String), intent(in)           :: input(:)
  type(CouplingStressBasisFunctions) :: this
  
  call this%read(input)
end function

impure elemental function new_CouplingStressBasisFunctions_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in)       :: input
  type(CouplingStressBasisFunctions)  :: this
  
  this = CouplingStressBasisFunctions(str(input))
end function
end module
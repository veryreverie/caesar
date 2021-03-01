! ======================================================================
! Reads harmonic forces, and generates the cartesian Hessian matrix.
! ======================================================================
module caesar_construct_supercell_hessian_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: construct_supercell_hessian
  
  interface
    ! ----------------------------------------------------------------------
    ! Read in the calculated forces, and use them with symmetry operations
    !    to construct the Hessian.
    ! ----------------------------------------------------------------------
    ! x is the collective vector of displacements, {xi}.
    ! f is the collective vector of forces, {fi}.
    ! Both are supercell%no_modes long.
    !
    ! Under the harmonic approximation, the hessian, F, is defined as:
    !    U = - x.F.x / 2
    !
    !    f = -dU/dx = F.x
    !
    ! Under the symmetry s:
    !    x -> x' = R.x
    !    f -> f' = R.f
    !
    ! F can be found by minimising L:
    !    L = sum(x,s)[ (f'-F.x')^2 ]
    !      = sum(x,s)[ f'.f' - 2*f'.F.x' + x'.F.F.x' ]
    !
    ! => 0 = dL/dF = -2 * sum(x,s)[ f'^x' - F.x'^x' ]
    ! 
    ! where ^ is the outer product.
    !
    ! => F = sum(x,s)[f'^x'] . inverse( sum(x,s)[x'^x'] )
    !
    ! sum(x,s)[x'^x'] is block diagonal, so can be inverted in 3x3 blocks.
    impure elemental module function construct_supercell_hessian(supercell, &
       & directory,calculation_reader,acoustic_sum_rule_forces,logfile)     &
       & result(output) 
      type(StructureData),     intent(in)    :: supercell
      type(String),            intent(in)    :: directory
      type(CalculationReader), intent(inout) :: calculation_reader
      logical,                 intent(in)    :: acoustic_sum_rule_forces
      type(OFile),             intent(inout) :: logfile
      type(CartesianHessian)                 :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Parse forces from electronic structure data,
    !    enforce acoustic sum rules if requires,
    !    and mass-reduce forces.
    ! ----------------------------------------------------------------------
    module function parse_forces(supercell,unique_directions, &
       & electronic_structure,acoustic_sum_rule_forces) result(output) 
      type(StructureData),       intent(in) :: supercell
      type(UniqueDirection),     intent(in) :: unique_directions(:)
      type(ElectronicStructure), intent(in) :: electronic_structure(:)
      logical,                   intent(in) :: acoustic_sum_rule_forces
      type(RealVector), allocatable         :: output(:,:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Construct and check sum[x'^x'].
    ! ----------------------------------------------------------------------
    module function construct_xx(unique_directions,supercell,logfile) &
       & result(output) 
      type(UniqueDirection), intent(in)    :: unique_directions(:)
      type(StructureData),   intent(in)    :: supercell
      type(OFile),           intent(inout) :: logfile
      type(RealMatrix), allocatable        :: output(:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Construct and check sum[f'^x'].
    ! ----------------------------------------------------------------------
    module function construct_fx(unique_directions,forces,supercell,logfile) &
       & result(output) 
      type(UniqueDirection), intent(in)    :: unique_directions(:)
      type(RealVector),      intent(in)    :: forces(:,:)
      type(StructureData),   intent(in)    :: supercell
      type(OFile),           intent(inout) :: logfile
      type(RealMatrix), allocatable        :: output(:,:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Construct and check F = s * sum[f'^x'] . inverse(sum[x'^x']).
    ! ----------------------------------------------------------------------
    ! sum[x'^x'] is block diagonal, so only the 3x3 diagonal blocks need inverting.
    ! N.B. the factor of s (sc_size) is to correct for the fact that only one
    !    atom in every s primitive cells is displaced, so the magnitude of x is
    !    effectively sqrt(s) times shorter than the equivalent displacement in
    !    the 1x1x1 supercell.
    module function construct_f(xx,fx,supercell,logfile) result(output) 
      type(RealMatrix),    intent(in)    :: xx(:)
      type(RealMatrix),    intent(in)    :: fx(:,:)
      type(StructureData), intent(in)    :: supercell
      type(OFile),         intent(inout) :: logfile
      type(CartesianHessian)             :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Checks the Hessian against the calculated forces and symmetries.
    ! ----------------------------------------------------------------------
    ! Checks that the elements corresponding to calculated forces have not been
    !    changed too much by symmetrisation.
    ! Checks that the Hessian has the correct symmetry properties.
    module subroutine check_hessian(hessian,forces,supercell, &
       & unique_directions,logfile) 
      type(CartesianHessian), intent(in)    :: hessian
      type(RealVector),       intent(in)    :: forces(:,:)
      type(StructureData),    intent(in)    :: supercell
      type(UniqueDirection),  intent(in)    :: unique_directions(:)
      type(OFile),            intent(inout) :: logfile
    end subroutine
  end interface
end module

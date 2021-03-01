! ======================================================================
! Calculate dynamical matrices at each q-point from the Hessians of the
!    non-diagonal supercells.
! ======================================================================
module caesar_calculate_dynamical_matrices_module
  use caesar_utils_module
  
  use caesar_structure_module
  use caesar_normal_mode_module
  
  use caesar_dynamical_matrix_module
  implicit none
  
  private
  
  public :: MatrixAndModes
  public :: calculate_dynamical_matrices
  
  ! Result type for calculate_dynamical_matrices.
  type :: MatrixAndModes
    type(DynamicalMatrix)          :: matrix
    type(ComplexMode), allocatable :: modes(:)
  end type
  
  interface calculate_dynamical_matrices
    ! Calculate dynamical matrices at each q-point from the Hessians of each
    !    non-diagonal supercell.
    module function calculate_dynamical_matrices_hessians(structure, &
       & supercells,supercell_hessians,qpoints,logfile) result(output) 
      type(StructureData),    intent(in)    :: structure
      type(StructureData),    intent(in)    :: supercells(:)
      type(CartesianHessian), intent(in)    :: supercell_hessians(:)
      type(QpointData),       intent(in)    :: qpoints(:)
      type(OFile),            intent(inout) :: logfile
      type(MatrixAndModes), allocatable     :: output(:)
    end function
  end interface
  
  interface
    ! --------------------------------------------------
    ! Construct the dynamical matrix at the specified q-point from the given
    !    Hessian.
    ! Considers all supercell, S, and symmetry, R, combinations such that
    !    S.R.q is a vector of integers, i.e. the q-point is a G-vector of the
    !    supercell.
    ! --------------------------------------------------
    module function construct_dynamical_matrix(qpoint,supercells,hessian, &
       & structure,subspace_id,logfile) result(output) 
      type(QpointData),       intent(in)    :: qpoint
      type(StructureData),    intent(in)    :: supercells(:)
      type(CartesianHessian), intent(in)    :: hessian(:)
      type(StructureData),    intent(in)    :: structure
      integer,                intent(in)    :: subspace_id
      type(OFile),            intent(inout) :: logfile
      type(MatrixAndModes)                  :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Compares two dynamical matrices.
    ! ----------------------------------------------------------------------
    module subroutine compare_dynamical_matrices(a,b,logfile) 
      type(DynamicalMatrix), intent(in)    :: a
      type(DynamicalMatrix), intent(in)    :: b
      type(OFile),           intent(inout) :: logfile
    end subroutine
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Transform a dynamical matrix and set of normal modes
    !    from one q-point to another.
    ! ----------------------------------------------------------------------
    ! Construct data at q_new from data at q_old, where
    !    R . q_old = q_new.
    ! N.B. the provided q-point should be q_new not q_old.
    ! The symmetry, S, maps equilibrium position ri to rj+R, and q-point q to q'.
    ! The +R needs to be corrected for, by multiplying the phase by exp(-iq'.R).
    module function transform_modes(input,symmetry,qpoint_from,qpoint_to) &
       & result(output) 
      type(DynamicalMatrix),  intent(in)    :: input
      type(SymmetryOperator), intent(in)    :: symmetry
      type(QpointData),       intent(in)    :: qpoint_from
      type(QpointData),       intent(in)    :: qpoint_to
      type(DynamicalMatrix)                 :: output
    end function
  end interface
  
  interface
    module function transform_dynamical_matrix(input,symmetry,qpoint_from, &
       & qpoint_to) result(output) 
      type(ComplexMatrix),    intent(in) :: input(:,:)
      type(SymmetryOperator), intent(in) :: symmetry
      type(QpointData),       intent(in) :: qpoint_from
      type(QpointData),       intent(in) :: qpoint_to
      type(ComplexMatrix), allocatable   :: output(:,:)
    end function
  end interface
end module

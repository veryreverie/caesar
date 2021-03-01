! ======================================================================
! Methods for generating supercells.
! ======================================================================
module caesar_supercell_module
  use caesar_utils_module
  
  use caesar_atom_module
  
  use caesar_qpoint_module
  use caesar_structure_data_module
  implicit none
  
  private
  
  public :: construct_supercell
  public :: construct_supercell_matrix
  
  interface
    ! ----------------------------------------------------------------------
    ! Takes three linearly independent (LI) vectors, a, b and c.
    ! The longest vector is replaced with the shortest linear combination
    !    of the three s.t. the three output vectors are still LI.
    ! The linear combinations considered are a+b+c, a+b-c, a-b+c and -a+b+c.
    ! ----------------------------------------------------------------------
    module function reduce_once(input,metric) result(output) 
      type(IntVector),  intent(in) :: input(3)
      type(RealMatrix), intent(in) :: metric
      type(IntVector)              :: output(3)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Takes a lattice matrix, L, whose rows are three vectors, a, b and c.
    ! Returns a matrix whose rows are linear combinations of a, b and c s.t.:
    !    - |A|/=0, i.e. the rows or the output are linearly independent.
    !    - Each vector is shorter than any linear combination of those vectors.
    ! ----------------------------------------------------------------------
    module function reduce(input,structure) result(output) 
      type(IntMatrix),     intent(in) :: input
      type(StructureData), intent(in) :: structure
      type(IntMatrix)                 :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Find a supercell matrix, S, s.t. S.q is a vector of integers for all q.
    ! ----------------------------------------------------------------------
    ! S is found s.t. |S| is as small as possible (whilst being >0).
    ! Returns answer in Hermite Normal Form.
    module function find_hnf_supercell_matrix(qpoints) result(output) 
      type(QpointData), intent(in) :: qpoints(:)
      type(IntMatrix)              :: output
    end function
  end interface
  
  interface construct_supercell_matrix
    ! ----------------------------------------------------------------------
    ! Construct the smallest matrix S s.t. S.q is a G-vector for all q.
    ! Smallest means:
    !    - |S| is as small as possible.
    !    - The vectors defined in fractional co-ordinates by the rows of S
    !         are as short as possible in cartesian co-ordinates.
    ! ----------------------------------------------------------------------
    module function construct_supercell_matrix_qpoints(qpoints,structure) &
       & result(output) 
      type(QpointData),    intent(in) :: qpoints(:)
      type(StructureData), intent(in) :: structure
      type(IntMatrix)                 :: output
    end function
  
    module function construct_supercell_matrix_qpoint(qpoint,structure) &
       & result(output) 
      type(QpointData),    intent(in) :: qpoint
      type(StructureData), intent(in) :: structure
      type(IntMatrix)                 :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Checks supercells.
    ! ----------------------------------------------------------------------
    ! Throws an error if there is a problem.
    impure elemental module subroutine check_supercell(supercell,structure) 
      type(StructureData), intent(in) :: supercell
      type(StructureData), intent(in) :: structure
    end subroutine
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Calculates the set of vectors which are not related to one another by
    !    lattice vectors.
    ! ----------------------------------------------------------------------
    ! Either calculates the R-vectors of the primitive cell which are unique
    !    in the supercell, or the G-vectors of the reciprocal supercell which are
    !    unique in the reciprocal primitive cell.
    module function calculate_unique_vectors(supercell,centre_on_origin) &
       & result(output) 
      ! Supercell vectors are the rows of supercell.
      type(IntMatrix), intent(in)  :: supercell
      logical,         intent(in)  :: centre_on_origin
      type(IntVector), allocatable :: output(:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Makes a supercell from a given primitive cell and supercell matrix.
    ! ----------------------------------------------------------------------
    module function construct_supercell(structure,supercell_matrix) &
       & result(output) 
      type(StructureData), intent(in) :: structure
      type(IntMatrix),     intent(in) :: supercell_matrix
      type(StructureData)             :: output
    end function
  end interface
end module

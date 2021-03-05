! ======================================================================
! Reads in normal modes from a separate calculation.
! ======================================================================
module caesar_read_normal_modes_module
  use caesar_common_module
  
  use caesar_harmonic_data_module
  implicit none
  
  private
  
  public :: read_normal_modes_mode
  
  interface
    ! ----------------------------------------------------------------------
    ! Generate keywords and helptext.
    ! ----------------------------------------------------------------------
    module function read_normal_modes_mode() result(output)
      type(ProgramMode) :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Main program.
    ! ----------------------------------------------------------------------
    module subroutine read_normal_modes_subroutine(arguments) 
      type(Dictionary), intent(in) :: arguments
    end subroutine
  end interface
  
  interface
    module function read_dynamical_matrices_castep(seedname,structure, &
       & qpoints) result(output) 
      type(String),        intent(in)    :: seedname
      type(StructureData), intent(in)    :: structure
      type(QpointData),    intent(in)    :: qpoints(:)
      type(DynamicalMatrix), allocatable :: output(:)
    end function
  end interface
  
  interface
    module function read_dynamical_matrices_qe(seedname,structure, &
       & large_supercell,qpoints) result(output) 
      type(String),        intent(in)    :: seedname
      type(StructureData), intent(in)    :: structure
      type(StructureData), intent(in)    :: large_supercell
      type(QpointData),    intent(in)    :: qpoints(:)
      type(DynamicalMatrix), allocatable :: output(:)
    end function
  end interface
  
  interface
    module function calculate_normal_modes(structure,qpoints, &
       & dynamical_matrices) result(output) 
      type(StructureData),   intent(in) :: structure
      type(QpointData),      intent(in) :: qpoints(:)
      type(DynamicalMatrix), intent(in) :: dynamical_matrices(:)
      type(ComplexMode), allocatable    :: output(:,:)
    end function
  end interface
end module

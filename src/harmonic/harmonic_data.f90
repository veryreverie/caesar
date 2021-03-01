! ======================================================================
! Data common to all harmonic calculations.
! ======================================================================
module caesar_harmonic_data_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: HarmonicData
  
  type, extends(NoDefaultConstructor) :: HarmonicData
    type(StructureData)           :: structure
    type(StructureData)           :: large_supercell
    type(QpointData), allocatable :: qpoints(:)
  end type
  
  interface HarmonicData
    module function new_HarmonicData(structure,large_supercell,qpoints) &
       & result(this) 
      type(StructureData), intent(in) :: structure
      type(StructureData), intent(in) :: large_supercell
      type(QpointData),    intent(in) :: qpoints(:)
      type(HarmonicData)              :: this
    end function
  
    ! Construct the harmonic data from the input arguments to setup_harmonic or
    !    read_normal_modes.
    ! This module function reads the input file.
    module function new_HarmonicData_arguments(file_type,seedname,grid, &
       & symmetry_precision,snap_to_symmetry,loto_direction) result(this) 
      type(String),         intent(in)           :: file_type
      type(String),         intent(in)           :: seedname
      integer,              intent(in)           :: grid(3)
      real(dp),             intent(in)           :: symmetry_precision
      logical,              intent(in)           :: snap_to_symmetry
      type(Fractionvector), intent(in), optional :: loto_direction
      type(HarmonicData)                         :: this
    end function
  end interface
end module

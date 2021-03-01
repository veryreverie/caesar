! ======================================================================
! Converges the cutoff energy.
! ======================================================================
! Runs multiple single-point calculations with different cutoff energies.
module caesar_converge_cutoff_and_kpoints_module
  use caesar_common_module
  
  use caesar_castep_output_file_module
  implicit none
  
  private
  
  public :: converge_cutoff_and_kpoints
  
  interface
    ! ----------------------------------------------------------------------
    ! Generate keywords and helptext.
    ! ----------------------------------------------------------------------
    module function converge_cutoff_and_kpoints() result(output) 
      type(CaesarMode) :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Main program.
    ! ----------------------------------------------------------------------
    module subroutine converge_cutoff_and_kpoints_subroutine(arguments) 
      type(Dictionary), intent(in) :: arguments
    end subroutine
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Runs CASTEP with a particular k-point spacing and cutoff.
    ! ----------------------------------------------------------------------
    module function run_castep(cutoff,kpoint_spacing,wd,dir,seedname, &
       & run_script,no_cores,cell_file,param_file,structure,          &
       & symmetry_precision) result(output) 
      integer,             intent(in) :: cutoff
      real(dp),            intent(in) :: kpoint_spacing
      type(String),        intent(in) :: wd
      type(String),        intent(in) :: dir
      type(String),        intent(in) :: seedname
      type(String),        intent(in) :: run_script
      integer,             intent(in) :: no_cores
      type(IFile),         intent(in) :: cell_file
      type(IFile),         intent(in) :: param_file
      type(StructureData), intent(in) :: structure
      real(dp),            intent(in) :: symmetry_precision
      type(CastepOutputFile)          :: output
    end function
  end interface
end module

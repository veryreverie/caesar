! ======================================================================
! Converges harmonic free energies w/r/t q-point spacing.
! Runs multiple harmonic caesar calculations with different q-point grids.
! ======================================================================
module caesar_converge_harmonic_qpoints_module
  use caesar_common_module
  
  implicit none
  
  private
  
  public :: converge_harmonic_qpoints_mode
  
  interface
    ! ----------------------------------------------------------------------
    ! Generates keywords and helptext.
    ! ----------------------------------------------------------------------
    module function converge_harmonic_qpoints_mode() result(output)
      type(ProgramMode) :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Main program.
    ! ----------------------------------------------------------------------
    module subroutine converge_harmonic_qpoints_subroutine(arguments) 
      type(Dictionary), intent(in) :: arguments
    end subroutine
  end interface
  
  interface
    module function calculate_free_energies(directory,qpoint_grid,  &
       & input_file_name,repeat_calculations,random_seed,arguments) &
       & result(output) 
      type(String),     intent(in) :: directory
      type(KpointGrid), intent(in) :: qpoint_grid
      type(String),     intent(in) :: input_file_name
      logical,          intent(in) :: repeat_calculations
      integer,          intent(in) :: random_seed
      type(Dictionary), intent(in) :: arguments
      type(RealVector)             :: output
    end function
  end interface
  
  interface
    module subroutine copy_file(input,output) 
      type(String), intent(in) :: input
      type(String), intent(in) :: output
    end subroutine
  end interface
  
  interface
    module subroutine write_output_file(structure,energy_tolerance, &
       & min_temperature,max_temperature,qpoint_spacings,free_energies) 
      type(StructureData), intent(in) :: structure
      real(dp),            intent(in) :: energy_tolerance
      real(dp),            intent(in) :: min_temperature
      real(dp),            intent(in) :: max_temperature
      real(dp),            intent(in) :: qpoint_spacings(:)
      type(RealVector),    intent(in) :: free_energies(:)
    end subroutine
  end interface
  
  interface
    ! Find the maximum element-wise difference between two vectors.
    impure elemental module function maximum_difference(this,that) &
       & result(output) 
      type(RealVector), intent(in) :: this
      type(RealVector), intent(in) :: that
      real(dp)                     :: output
    end function
  end interface
end module

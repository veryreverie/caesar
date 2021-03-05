! ======================================================================
! Converges harmonic frequencies and free energies w/r/t
!    DFT cutoff energy, k-point spacing and electronic smearing.
! Runs multiple harmonic caesar calculations with different
!    cutoff energies and k-point spacings.
! ======================================================================
module caesar_converge_harmonic_frequencies_module
  use caesar_common_module

  implicit none
  
  private
  
  public :: converge_harmonic_frequencies_mode
  
  interface
    ! ----------------------------------------------------------------------
    ! Generate keywords and helptext
    ! ----------------------------------------------------------------------
    module function converge_harmonic_frequencies_mode() result(output)
      type(ProgramMode) :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Main program.
    ! ----------------------------------------------------------------------
    module subroutine converge_harmonic_frequencies_subroutine(arguments) 
      type(Dictionary), intent(in) :: arguments
    end subroutine
  end interface
  
  interface
    module subroutine write_output_file(structure,cutoffs,        &
       & cutoff_frequencies,cutoff_free_energies,kpoint_spacings, &
       & kpoint_frequencies,kpoint_free_energies,smearings,       &
       & smearing_frequencies,smearing_free_energies) 
      type(StructureData), intent(in) :: structure
      real(dp),            intent(in) :: cutoffs(:)
      type(RealVector),    intent(in) :: cutoff_frequencies(:)
      type(RealVector),    intent(in) :: cutoff_free_energies(:)
      real(dp),            intent(in) :: kpoint_spacings(:)
      type(RealVector),    intent(in) :: kpoint_frequencies(:)
      type(RealVector),    intent(in) :: kpoint_free_energies(:)
      real(dp),            intent(in) :: smearings(:)
      type(RealVector),    intent(in) :: smearing_frequencies(:)
      type(RealVector),    intent(in) :: smearing_free_energies(:)
    end subroutine
  end interface
  
  interface
    module function calculate_frequencies(directory,cutoff,kpoint_grid, &
       & smearing,seedname,no_qpoints,file_type,arguments) result(output) 
      type(String),     intent(in)           :: directory
      real(dp),         intent(in), optional :: cutoff
      type(KpointGrid), intent(in), optional :: kpoint_grid
      real(dp),         intent(in), optional :: smearing
      type(String),     intent(in)           :: seedname
      integer,          intent(in)           :: no_qpoints
      type(String),     intent(in)           :: file_type
      type(Dictionary), intent(in)           :: arguments
      type(RealVector)                       :: output
    end function
  end interface
  
  interface
    module function calculate_free_energies(directory,repeat_calculations, &
       & random_seed,arguments) result(output) 
      type(String),     intent(in) :: directory
      logical,          intent(in) :: repeat_calculations
      integer,          intent(in) :: random_seed
      type(Dictionary), intent(in) :: arguments
      type(RealVector)             :: output
    end function
  end interface
  
  interface
    module subroutine write_castep_files(directory,seedname,cutoff, &
       & kpoint_grid,smearing) 
      type(String),     intent(in)           :: directory
      type(String),     intent(in)           :: seedname
      real(dp),         intent(in), optional :: cutoff
      type(KpointGrid), intent(in), optional :: kpoint_grid
      real(dp),         intent(in), optional :: smearing
    end subroutine
  end interface
  
  interface
    module subroutine write_qe_file(directory,seedname,cutoff,kpoint_grid, &
       & smearing) 
      type(String),     intent(in)           :: directory
      type(String),     intent(in)           :: seedname
      real(dp),         intent(in), optional :: cutoff
      type(KpointGrid), intent(in), optional :: kpoint_grid
      real(dp),         intent(in), optional :: smearing
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

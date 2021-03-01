! ======================================================================
! Allows user to monitor convergence of harmonic frequencies and
!    free energies wrt qpoint grid
! Runs multiple harmonic caesar calculations with different qpoint grids
! ======================================================================
module caesar_converge_qpoint_grid_module
  use caesar_common_module

  implicit none
  
  private
  
  public :: converge_qpoint_grid
  
  interface
    ! ----------------------------------------------------------------------
    ! Generate keywords and helptext
    ! ----------------------------------------------------------------------
    module function converge_qpoint_grid() result(output) 
      type(CaesarMode) :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Main program
    ! ----------------------------------------------------------------------
    module subroutine converge_qpoint_grid_subroutine(arguments) 
      type(Dictionary), intent(in) :: arguments
    end subroutine
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! FUNCTIONS
    ! ----------------------------------------------------------------------
    
    ! ----------------------------------------------------------------------
    ! Write setup_harmonic input file in specified directory
    ! ----------------------------------------------------------------------
    module function write_setupgrid(file_type,seedname,grid, &
       & symmetry_precision,harmonic_displacement,dir)       &
       & result(setup_harmonic_inputs) 
      type(String),        intent(in) :: dir, seedname, file_type, grid
      real(dp),            intent(in) :: symmetry_precision, harmonic_displacement
      
      type(OFile)  :: setup_harmonic_inputs
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Write run_harmonic input file in specified directory
    ! ----------------------------------------------------------------------
    module function write_runinput( supercell_file,run_script,no_cores, &
       & no_nodes,dir) result(run_harmonic_inputs) 
      type(IFile),         intent(in) :: supercell_file
      type(String),        intent(in) :: dir
      type(String),        intent(in) :: run_script
      integer,             intent(in) :: no_cores, no_nodes
    
      type(OFile)  :: run_harmonic_inputs
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Write calculate_normal_modes input file in specified directory
    ! ----------------------------------------------------------------------
    module function write_normmode( acoustic_sum_rule,dir) &
       & result(normal_mode_inputs) 
      type(String),        intent(in) :: dir, acoustic_sum_rule
    
      type(OFile)  :: normal_mode_inputs
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Write calculate_harmonic_observables input file in specified directory
    ! ----------------------------------------------------------------------
    module function write_obs_input( min_temperature,max_temperature, &
       & no_temperature_steps,min_frequency,path,no_dos_samples,dir)  &
       & result(harm_obs_input) 
      real(dp),            intent(in) :: min_temperature
      real(dp),            intent(in) :: max_temperature
      integer,             intent(in) :: no_temperature_steps
      real(dp),            intent(in) :: min_frequency
      type(String),        intent(in) :: path
      integer,             intent(in) :: no_dos_samples
      type(String),        intent(in) :: dir
    
      type(OFile)  :: harm_obs_input
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Index module function for two strings
    ! ----------------------------------------------------------------------
    module function findindex_character(string1,string2) result(index_found) 
      character(*), intent(in)     :: string1, string2
      integer      :: index_found
    end function
  end interface
  
  interface
    module function findindex_string(string1,string2) result(index_found) 
      type(String),        intent(in) :: string1, string2
      integer      :: index_found
    end function
  end interface
end module

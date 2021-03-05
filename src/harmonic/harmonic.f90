! ======================================================================
! The operations of Caesar within the harmonic approximation.
! ======================================================================
! This module is simply an interface for the various harmonic modules.
module caesar_harmonic_module
  use caesar_common_module
  
  use caesar_hartree_to_ev_module
  use caesar_snap_to_symmetry_module
  use caesar_setup_harmonic_module
  use caesar_run_harmonic_module
  use caesar_calculate_normal_modes_module
  use caesar_read_normal_modes_module
  use caesar_plot_normal_modes_module
  use caesar_calculate_harmonic_observables_module
  use caesar_plot_dos_and_dispersion_module
  use caesar_plot_thermodynamic_variables_module
  use caesar_converge_harmonic_qpoints_module
  use caesar_plot_harmonic_qpoint_convergence_module
  implicit none
  
  interface
    module function harmonic_modes() result(output)
      type(ProgramMode), allocatable :: output(:)
    end function
  end interface
end module

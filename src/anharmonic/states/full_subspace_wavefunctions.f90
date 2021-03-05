! ======================================================================
! Wavefunctions spanning the full subspace.
! ======================================================================
module caesar_full_subspace_wavefunctions_module
  use caesar_common_module
  
  use caesar_subspace_wavefunctions_module
  implicit none
  
  private
  
  public :: FullSubspaceWavefunctions
  
  type, extends(SubspaceWavefunctions) :: FullSubspaceWavefunctions
    integer                   :: subspace_id
    integer,      allocatable :: mode_ids(:)
    integer,      allocatable :: paired_mode_ids(:)
    type(String)              :: harmonic_ground_state
    real(dp),     allocatable :: energies(:)
    type(String), allocatable :: wavefunctions(:)
  contains
    ! Type representation.
    procedure, public, nopass :: representation => &
                               & representation_FullSubspaceWavefunctions
    
    ! I/O.
    procedure, public :: read  => read_FullSubspaceWavefunctions
    procedure, public :: write => write_FullSubspaceWavefunctions
  end type
  
  interface FullSubspaceWavefunctions
    module function new_FullSubspaceWavefunctions(subspace_id,mode_ids, &
        paired_mode_ids,harmonic_ground_state,energies,wavefunctions) &
          & result(this) 
      integer,      intent(in)        :: subspace_id
      integer,      intent(in)        :: mode_ids(:)
      integer,      intent(in)        :: paired_mode_ids(:)
      type(String), intent(in)        :: harmonic_ground_state
      real(dp),     intent(in)        :: energies(:)
      type(String) ,intent(in) :: wavefunctions(:) 
      type(FullSubspaceWavefunctions) :: this
    end function
  end interface
  
  interface
    impure elemental module function representation_FullSubspaceWavefunctions() result(output) 
      type(String) :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_FullSubspaceWavefunctions(this,input) 
      class(FullSubspaceWavefunctions) ,intent(out) :: this
      type(String),                     intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_FullSubspaceWavefunctions(this) result(output) 
      class(FullSubspaceWavefunctions) ,intent(in) :: this
      type(String), allocatable                    :: output(:)
    end function
  end interface
  
  interface FullSubspaceWavefunctions
    module function new_FullSubspaceWavefunctions_Strings(input) result(this) 
      type(String), intent(in)        :: input(:)
      type(FullSubspaceWavefunctions) :: this
    end function
  
    impure elemental module function new_FullSubspaceWavefunctions_StringArray(input) result(this) 
      type(StringArray), intent(in)   :: input
      type(FullSubspaceWavefunctions) :: this
    end function
  end interface
end module

! ======================================================================
! Wavefunctions which treat q-points separately.
! ======================================================================
module caesar_split_qpoints_wavefunctions_module
  use caesar_common_module
  
  use caesar_subspace_wavefunctions_module
  implicit none
  
  private
  
  public :: SplitQpointsWavefunctions
  
  type, extends(SubspaceWavefunctions) :: SplitQpointsWavefunctions
    integer                   :: subspace_id
    integer,      allocatable :: mode_ids(:)
    integer,      allocatable :: paired_mode_ids(:)
    type(String)              :: harmonic_ground_state
    real(dp),     allocatable :: energies(:)
    type(String), allocatable :: wavefunctions(:)
  contains
    ! Type representation.
    procedure, public, nopass :: representation => &
                               & representation_SplitQpointsWavefunctions
    
    ! I/O.
    procedure, public :: read  => read_SplitQpointsWavefunctions
    procedure, public :: write => write_SplitQpointsWavefunctions
  end type
  
  interface SplitQpointsWavefunctions
    module function new_SplitQpointsWavefunctions(subspace_id,mode_ids, &
        paired_mode_ids,harmonic_ground_state,energies,wavefunctions) &
          & result(this) 
      integer,      intent(in)        :: subspace_id
      integer,      intent(in)        :: mode_ids(:)
      integer,      intent(in)        :: paired_mode_ids(:)
      type(String), intent(in)        :: harmonic_ground_state
      real(dp),     intent(in)        :: energies(:)
      type(String) ,intent(in) :: wavefunctions(:) 
      type(SplitQpointsWavefunctions) :: this
    end function
  end interface
  
  interface
    impure elemental module function representation_SplitQpointsWavefunctions() result(output) 
      type(String) :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_SplitQpointsWavefunctions(this,input) 
      class(SplitQpointsWavefunctions) ,intent(out) :: this
      type(String),                     intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_SplitQpointsWavefunctions(this) result(output) 
      class(SplitQpointsWavefunctions) ,intent(in) :: this
      type(String), allocatable                    :: output(:)
    end function
  end interface
  
  interface SplitQpointsWavefunctions
    module function new_SplitQpointsWavefunctions_Strings(input) result(this) 
      type(String), intent(in)        :: input(:)
      type(SplitQpointsWavefunctions) :: this
    end function
  
    impure elemental module function new_SplitQpointsWavefunctions_StringArray(input) result(this) 
      type(StringArray), intent(in)   :: input
      type(SplitQpointsWavefunctions) :: this
    end function
  end interface
end module

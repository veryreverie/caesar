! ======================================================================
! Takes a potential or stress, and a list of subspaces,
!    and generates the single-subspace potentials or stresses.
! ======================================================================
module caesar_generate_subspace_potentials_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  implicit none
  
  private
  
  public :: generate_subspace_potentials
  public :: generate_subspace_stresses
  
  interface
    ! ----------------------------------------------------------------------
    ! Takes a potential V and an array of subspace states {|i>}, and generates the
    !    set of single-subspace potentials {V_i}, defined by
    !    V_i = (prod_{j/=i}<j|)V(prod_{j/=i}|j>)
    !        - (prod_i<i|)V(prod_i|i>) * (n-1)/n
    ! ----------------------------------------------------------------------
    ! The naive method of calculating {V_i} for n subspaces takes
    !    n(n-1) operations.
    ! This can be accelerated using a bisection method, outlined below.
    ! 
    ! V0(1) = V, the input potential.
    ! 
    ! The first iteration splits the states into two intervals,
    !    [1,s-1] and [s,n], where s=n/2, and two potentials are calculated:
    !    - V1(1) = (<s|<s+1|...<n|)V0(1)(|s>|s+1>...|n>)
    !    - V1(s) = (<1|<2|...<s-1|)V0(1)(|1>|2>...|s-1>)
    ! These intervals are recorded in terms of their min and max values:
    !   mins = [1  , s]
    !   maxs = [s-1, n]
    !
    ! The next iteration splits each of the intervals into two intervals,
    !    copies the potential to both intervals, and integrates the potential
    !    corresponding to each interval over the states in the other interval.
    ! This method takes O(n.log(n)) operations.
    module function generate_subspace_potentials(potential,subspaces, &
       & subspace_bases,subspace_states,old_subspace_potentials,      &
       & anharmonic_data) result(output) 
      class(PotentialData),     intent(in)           :: potential
      type(DegenerateSubspace), intent(in)           :: subspaces(:)
      class(SubspaceBasis),     intent(in)           :: subspace_bases(:)
      class(BasisStates),       intent(inout)        :: subspace_states(:)
      type(PotentialPointer),   intent(in), optional :: old_subspace_potentials(:)
      type(AnharmonicData),     intent(in)           :: anharmonic_data
      type(PotentialPointer), allocatable            :: output(:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Breaks a general stress into single-subspace stresses,
    !    in a manner analogous to generate_subspace_potentials.
    ! ----------------------------------------------------------------------
    module function generate_subspace_stresses(stress,subspaces, &
       & subspace_bases,subspace_states,anharmonic_data) result(output) 
      class(StressData),        intent(in)    :: stress
      type(DegenerateSubspace), intent(in)    :: subspaces(:)
      class(SubspaceBasis),     intent(in)    :: subspace_bases(:)
      class(BasisStates),       intent(inout) :: subspace_states(:)
      type(AnharmonicData),     intent(in)    :: anharmonic_data
      type(StressPointer), allocatable        :: output(:)
    end function
  end interface
end module

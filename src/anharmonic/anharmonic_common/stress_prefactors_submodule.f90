submodule (caesar_stress_prefactors_module) caesar_stress_prefactors_submodule
  use caesar_anharmonic_common_module
contains

module procedure new_StressPrefactors
  this%subspace_id_ = subspace_id
  this%mode_ids_    = mode_ids
  this%pair_locs_   = pair_locs
  this%prefactors_  = prefactors
end procedure

module procedure new_StressPrefactors_subspace
  type(ComplexMode), allocatable :: subspace_modes(:)
  
  integer, allocatable :: pair_locs(:)
  
  type(RealMatrix), allocatable :: prefactors(:,:)
  
  integer :: i,j,ialloc
  
  subspace_modes = subspace%modes(modes)
  
  allocate( prefactors(size(subspace_modes),size(subspace_modes)), &
          & stat=ialloc); call err(ialloc)
  prefactors = dblemat(zeroes(3,3))
  do i=1,size(subspace_modes)
    do j=1,size(subspace_modes)
      if (subspace_modes(i)%qpoint_id==subspace_modes(j)%paired_qpoint_id) then
        prefactors(i,j) = subspace_modes(i)%stress_prefactor(subspace_modes(j))
      endif
    enddo
  enddo
  
  pair_locs = [( first(subspace_modes%id==subspace_modes(i)%paired_id), &
               & i=1,                                                   &
               & size(subspace_modes)                                   )]
  
  this = StressPrefactors( subspace%id,       &
                         & subspace_modes%id, &
                         & pair_locs,         &
                         & prefactors         )
end procedure

module procedure prefactor_StressPrefactors_modes
  output = this%prefactors_( first(this%mode_ids_==mode1%id), &
                           & first(this%mode_ids_==mode2%id)  )
end procedure

module procedure prefactor_StressPrefactors_ids
  output = this%prefactors_( first(this%mode_ids_==mode1), &
                           & first(this%mode_ids_==mode2)  )
end procedure

module procedure average_prefactor
  integer :: no_modes
  
  integer :: i
  
  no_modes = size(this%mode_ids_)
  
  output = sum([(this%prefactors_(i,this%pair_locs_(i)), i=1, no_modes)]) &
       & / no_modes
end procedure
end submodule

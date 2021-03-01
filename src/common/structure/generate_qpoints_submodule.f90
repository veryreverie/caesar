submodule (caesar_generate_qpoints_module) caesar_generate_qpoints_submodule
  use caesar_structure_module
contains

module procedure generate_qpoints
  ! Working variables
  integer, allocatable :: paired_qpoints(:)
  
  ! Temporary variables
  integer :: i,j,ialloc
  
  ! --------------------------------------------------
  ! Construct q-points from G-vectors of large supercell.
  ! --------------------------------------------------
  allocate(output(large_supercell%sc_size), stat=ialloc); call err(ialloc)
  do i=1,large_supercell%sc_size
    output(i)%id     = i
    output(i)%qpoint = transpose(large_supercell%recip_supercell) &
                   & * large_supercell%gvectors(i)
  enddo
  
  ! --------------------------------------------------
  ! Find paired q-points.
  ! --------------------------------------------------
  ! qpoint + paired_qpoint = G, for a primitive-cell G-vector.
  allocate( paired_qpoints(large_supercell%sc_size), &
          & stat=ialloc); call err(ialloc)
  paired_qpoints = 0
  do i=1,size(output)
    do j=1,size(output)
      if (is_int(output(i)%qpoint+output(j)%qpoint)) then
        if (paired_qpoints(i)==0 .and. paired_qpoints(j)==0) then
          paired_qpoints(i) = j
          paired_qpoints(j) = i
        else
          if (paired_qpoints(i)/=j .or. paired_qpoints(j)/=i) then
            call print_line(CODE_ERROR//': error pairing q-points.')
          endif
        endif
      endif
    enddo
  enddo
  
  if (any(paired_qpoints==0)) then
    call print_line(CODE_ERROR//': q-points were not succesfully paired up.')
    call err()
  endif
  
  output%paired_qpoint_id = paired_qpoints
end procedure
end submodule

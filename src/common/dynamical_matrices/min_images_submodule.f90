submodule (caesar_min_images_module) caesar_min_images_submodule
  use caesar_dynamical_matrices_module
contains

module procedure size_MinImages
  output = size(this%image_rvectors)
end procedure

module procedure calculate_min_images
  type(AtomData) :: atom_i
  type(AtomData) :: atom_j
  integer        :: rvec
  
  type(MinImages), allocatable :: min_images(:,:,:)
  
  integer :: i,j,ialloc

  ! Calculate minimum images for atom 2 in the primitive cell.
  allocate( min_images( supercell%sc_size,        &
          &             supercell%no_atoms_prim,  &
          &             supercell%no_atoms_prim), &
          & stat=ialloc); call err(ialloc)
  do i=1,supercell%no_atoms
    atom_i = supercell%atoms(i)
    do j=1,supercell%no_atoms_prim
      atom_j = supercell%atoms(j)
      min_images(atom_i%rvec_id(),atom_j%id(),atom_i%prim_id()) = &
         & min_images_brute_force(supercell,atom_j%id(),atom_i%id())
    enddo
  enddo
  
  ! Copy across minimum images to pairs of atoms related by R-vectors.
  allocate( output(supercell%no_atoms_prim, supercell%no_atoms), &
          & stat=ialloc); call err(ialloc)
  do i=1,supercell%no_atoms
    atom_i = supercell%atoms(i)
    rvec   = atom_i%rvec_id()
    do j=1,supercell%no_atoms_prim
      atom_j = supercell%atoms(j)
      output(atom_j%id(),atom_i%id()) = min_images( rvec,             &
                                                  & atom_j%prim_id(), &
                                                  & atom_i%prim_id()  )
    enddo
  enddo
end procedure

module procedure min_images_brute_force
  ! Number of "shells" of supercell R-vectors to check.
  integer,parameter :: check_shell=3
  
  ! Maximum number of images.
  ! There is no real penalty for having this be larger than strictly necessary.
  integer, parameter :: max_images = 27
  
  ! Positions in fractional supercell co-ordinates.
  type(RealVector) :: a_frac_pos
  type(RealVector) :: b_frac_pos
  
  ! The R-vector from a to b ignoring supercell R-vectors.
  type(IntVector)  :: rvector
  
  ! Tolerance parameters.
  real(dp), parameter :: tol=1.0e-8_dp
  real(dp)            :: tol_l2
  
  ! Distances in cartesian co-ordinates.
  type(RealVector) :: displacement
  real(dp)         :: l2_distance
  real(dp)         :: min_l2_distance
  
  ! Output data.
  integer                      :: no_images
  type(IntVector), allocatable :: supercell_rvectors(:)
  
  ! Temporary variables.
  integer :: i,j,k,ialloc
  
  rvector = supercell%rvectors(supercell%atoms(b)%rvec_id()) &
        & - supercell%rvectors(supercell%atoms(a)%rvec_id())
  
  tol_l2 = tol*supercell%volume**(2/3.0_dp)
  
  ! Get the fractional position of both atoms.
  a_frac_pos = supercell%atoms(a)%fractional_position()
  b_frac_pos = supercell%atoms(b)%fractional_position()
  
  ! Identify R-vectors.
  allocate(supercell_rvectors(max_images), stat=ialloc); call err(ialloc)
  min_l2_distance = huge(1.0_dp)
  no_images = 0
  do i=-check_shell,check_shell
    do j=-check_shell,check_shell
      do k=-check_shell,check_shell
        displacement = transpose(supercell%lattice) &
                   & * (b_frac_pos-a_frac_pos+vec([i,j,k]))
        l2_distance = displacement * displacement
        
        if (l2_distance < min_l2_distance-tol_l2) then
          ! The new R-vector makes a+R and b closer than any previous R-vector.
          min_l2_distance = l2_distance
          no_images = 1
          supercell_rvectors(no_images) = vec([i,j,k])
        elseif (l2_distance <= min_l2_distance+tol_l2) then
          ! The new R-vector makes a+R and b similarly close to
          !    previous R-vectors.
          no_images = no_images+1
          if (no_images>max_images) then
            call print_line(ERROR//': max_images too small.')
            call err()
          endif
          supercell_rvectors(no_images) = vec([i,j,k])
        endif
      enddo
    enddo
  enddo
  
  ! Check that at least one minimum image has been found.
  if (no_images==0) then
    call print_line(ERROR//': bug in min_images_brute_force.')
    call err()
  endif
  
  ! Construct output.
  allocate(output%image_rvectors(no_images), stat=ialloc); call err(ialloc)
  do i=1,no_images
    ! Convert the supercell R-vector from fractional supercell co-ords
    !    to fractional primitive co-ords,
    !    and add on the primitive R-vector.
    output%image_rvectors(i) = transpose(supercell%supercell) &
                           & * supercell_rvectors(i)          &
                           & + rvector
  enddo
  
  ! Check that all output R-vectors are equal to
  !    rvector + a supercell R-vector.
  do i=1,size(output%image_rvectors)
    if (.not. is_int( supercell%recip_supercell &
                  & * (output%image_rvectors(i)-rvector))) then
      call print_line(CODE_ERROR//': A minimum image R-vector is not &
         &equivalent to the input R-vector modulo a supercell R-vector.')
      call print_line('Image R-vector: '//output%image_rvectors(i))
      call print_line('Input R-vector: '//rvector)
      call print_line('Supercell:')
      call print_lines(supercell%supercell)
      call print_line('Recip Supercell:')
      call print_lines(supercell%recip_supercell)
      call err()
    endif
  enddo
end procedure
end submodule

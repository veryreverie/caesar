module coupled_atomic_units_module
  use string_module
  use io_module
  use constants_module, only : dp
contains

subroutine coupled_atomic_units()
  use constants_module, only : ev_per_hartree
  implicit none
  
  integer :: no_modes,no_indep_data,no_data,no_cells,no_couplings,no_params,order
  real(dp) :: energy_ref
  real(dp), allocatable :: harmonic(:,:),positions_indep(:,:),positions(:,:,:),energy(:,:),energy_indep(:,:)
  real(dp), allocatable :: params(:,:)
  integer, allocatable :: symmetry_mode(:),symmetry_ref(:),working_modes(:,:)
  type(String) :: functional
  
  integer  :: modes_read
  real(dp) :: output
  
  ! Files.
  type(String), allocatable :: input_file(:)
  type(String), allocatable :: symmetry_file(:)
  type(String), allocatable :: symmetry_coupling_file(:)
  type(String), allocatable :: energy_file(:)
  type(String), allocatable :: coupled_energy_file(:)
  integer                   :: output_file
  
  ! Temporary variables.
  type(String), allocatable :: line(:)
  integer :: i,j,ialloc
  
  ! Read in inputs.
  input_file = read_lines('post_input_coupled.dat')
  line = split(input_file(1))
  no_modes = int(line(1)) - 3
  no_indep_data = int(line(2))
  no_data = int(line(3))
  no_cells = int(line(4))
  no_couplings = int(line(5))
  energy_ref = dble(line(5))
  functional = line(6)
  order = int(line(7))
  
  if(functional == 'polynomial')then
    no_params=order-1
  else if(functional == 'sine')then
    no_params=order+1
  else
    write(*,*)'incorrect functional form.'
    call err()
  endif

  allocate( harmonic(no_couplings,2),                &
          & positions_indep(no_modes,no_indep_data), &
          & energy_indep(no_modes,no_indep_data),    &
          & positions(no_couplings,no_data,2),       &
          & energy(no_couplings,no_data),            &
          & working_modes(no_couplings,2),           &
          & symmetry_mode(no_modes),                 &
          & symmetry_ref(no_modes),                  &
          & params(no_modes,no_params),              &
          & stat=ialloc); call err(ialloc)

  ! Read in symmetry.
  symmetry_file = read_lines('symmetry.dat')
  do i=1,no_modes
    line = split(symmetry_file(i+3))
    symmetry_mode(i) = int(line(1))
    symmetry_ref(i) = int(line(2))
  enddo

  ! Read in symmetry couplings.
  symmetry_coupling_file = read_lines('symmetry_coupling.dat')
  do i=1,no_couplings
    line = split(symmetry_coupling_file(i))
    working_modes(i,:) = int(line)-3
  enddo
  
  ! Read in independent energy.
  modes_read = 0
  energy_file = read_lines('energy.dat')
  do i=1,no_modes
    if (symmetry_mode(i)==symmetry_ref(i)) then
      do j=1,no_indep_data
        line = split(energy_file(modes_read*(no_indep_data+3)+j+1))
        positions_indep(i,j) = dble(line(1))
        energy_indep(i,j) = dble(line(2))
      enddo
    else
      positions_indep(i,:) = positions_indep(symmetry_ref(i)-3,:)
      energy_indep(i,:) = energy_indep(symmetry_ref(i)-3,:)
    endif
  enddo
  
  ! Read in coupled energy.
  coupled_energy_file = read_lines('coupled_energy.dat')
  do i=1,no_couplings
    line = split(coupled_energy_file((i-1)*(no_data+3)+1))
    harmonic(i,:) = dble(line(2:3))
    do j=1,no_data
      line = split(coupled_energy_file((i-1)*(no_data+3)+1+j))
      positions(i,j,:) = dble(line(1:2))
      energy(i,j) = dble(line(3))
    enddo
  enddo
  
  ! Write out result.
  output_file = open_write_file('coupled_energy_au.dat')
  do i=1,no_couplings
    call print_line(output_file, '# '//harmonic(i,:))
    do j=1,no_data
      output = (energy(i,j)-energy_ref)/(no_cells*ev_per_hartree)
      call print_line(output_file, positions(i,j,:)//' '//output)
    enddo
  enddo
  close(output_file)
end subroutine

! ----------------------------------------------------------------------
! Single mode potential
! ----------------------------------------------------------------------
function v_indep(i,no_modes,no_params,params,q,functional,order) result(output)
  implicit none
  
  integer,      intent(in) :: i
  integer,      intent(in) :: no_modes
  integer,      intent(in) :: no_params
  integer,      intent(in) :: order
  real(dp),     intent(in) :: params(no_modes,no_params)
  real(dp),     intent(in) :: q
  type(String), intent(in) :: functional
  real(dp) :: output

  if(functional=='polynomial')then
    if(order==6)then
      output = params(i,1)*q**2 + params(i,2)*q**3 + &
       &params(i,3)*q**4 + params(i,4)*q**5 + &
       &params(i,5)*q**6
    else if(order==8)then
      output = params(i,1)*q**2 + params(i,2)*q**3 + &
       &params(i,3)*q**4 + params(i,4)*q**5 + &
       &params(i,5)*q**6 + params(i,6)*q**7 + &
       &params(i,7)*q**8
    else if(order==10)then
      output = params(i,1)*q**2 + params(i,2)*q**3 + &
       &params(i,3)*q**4 + params(i,4)*q**5 + &
       &params(i,5)*q**6 + params(i,6)*q**7 + &
       &params(i,7)*q**8 + params(i,8)*q**9 + &
       &params(i,9)*q**10
    else if(order==12)then
      output = params(i,1)*q**2 + params(i,2)*q**3 + &
       &params(i,3)*q**4 + params(i,4)*q**5 + &
       &params(i,5)*q**6 + params(i,6)*q**7 + &
       &params(i,7)*q**8 + params(i,8)*q**9 + &
       &params(i,9)*q**10 + params(i,10)*q**11 + &
       &params(i,11)*q**12
    endif
  else if(functional=='sine')then
    if(order==6)then
      output = params(i,1)*q**2 + params(i,2)*q**3 + &
       &params(i,3)*q**4 + params(i,4)*q**5 + &
       &params(i,5)*q**6 + &
       &params(i,6)*sin(params(i,7)*q)**2
    else if(order==8)then
      output = params(i,1)*q**2 + params(i,2)*q**3 + &
       &params(i,3)*q**4 + params(i,4)*q**5 + &
       &params(i,5)*q**6 + params(i,6)*q**7 + &
       &params(i,7)*q**8 + &
       &params(i,8)*sin(params(i,9)*q)**2
    else if(order==10)then
      output = params(i,1)*q**2 + params(i,2)*q**3 + &
       &params(i,3)*q**4 + params(i,4)*q**5 + &
       &params(i,5)*q**6 + params(i,6)*q**7 + &
       &params(i,7)*q**8 + params(i,8)*q**9 + &
       &params(i,9)*q**10 + &
       &params(i,10)*sin(params(i,11)*q)**2
    endif
  else
    call print_line('Error: the functional form of the independent potential &
       &is incorrect!')
    call err()
  endif
end function
end module
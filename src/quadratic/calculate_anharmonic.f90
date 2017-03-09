module calculate_anharmonic_module
  use constants, only : dp
  implicit none
  
contains

subroutine calculate_anharmonic(multiplicity, no_modes, Nbasis, harmonic,&
    &eigenvals, result_file)
  use constants, only : kB
  use string_module
  use file_module
  implicit none
  
  integer,  intent(in) :: multiplicity(:)
  integer,  intent(in) :: no_modes
  integer,  intent(in) :: Nbasis
  real(dp), intent(in) :: harmonic(:,:,:)
  real(dp), intent(in) :: eigenvals(:,:,:)
  integer,  intent(in) :: result_file
  
  ! Temperature parameters
  real(dp), parameter :: dtemperature    = 50.0_dp ! delta T (K)
  integer,  parameter :: no_temperatures = 21
  
  ! Temperatures
  real(dp) :: betas(no_temperatures)    ! {1/kB*T}
  
  ! number of kpoints, = size(multiplicity)
  integer :: no_kpoints
  
  ! Working variables
  integer :: i,j,k,l
  integer :: total_kpoints
  real(dp) :: renormalised_eigenvals,renormalised_harmonic
  real(dp),allocatable :: part_fn(:,:,:),har_part_fn(:,:,:)
  
  ! Get kpoint data
  no_kpoints = size(multiplicity)
  total_kpoints = sum(multiplicity)
  
  ! Calculate thermal energies
  do i=1,no_temperatures
    betas(i) = 1.0_dp/((i-1)*dtemperature*kB)
  enddo

  ! Calculate partition function
  allocate(part_fn(no_modes,no_kpoints,no_temperatures))
  allocate(har_part_fn(no_modes,no_kpoints,no_temperatures))
  part_fn=0.0_dp
  har_part_fn=0.0_dp
  do i=1,no_temperatures
    do j=1,no_kpoints
      do k=1,no_modes
        do l=1,Nbasis
          part_fn(k,j,i)=part_fn(k,j,i)+exp(-eigenvals(l,k,j)*betas(i))
          har_part_fn(k,j,i)=har_part_fn(k,j,i)+exp(-harmonic(l,k,j)*betas(i))
        enddo
      enddo
    enddo
  enddo

  ! Calculate anharmonic vibrational correction
  do i=1,no_temperatures
    renormalised_eigenvals=0.0_dp
    renormalised_harmonic=0.0_dp
    if (i==1) then ! T=0
      do j=1,no_kpoints
        do k=1,no_modes
          renormalised_eigenvals=renormalised_eigenvals+&
           &eigenvals(1,k,j)*multiplicity(j)/total_kpoints
          renormalised_harmonic=renormalised_harmonic+&
           &harmonic(1,k,j)*multiplicity(j)/total_kpoints
        enddo ! k
      enddo ! j
    else ! T /= 0
       do j=1,no_kpoints
         do k=1,no_modes
           renormalised_eigenvals = renormalised_eigenvals &
            & -log(part_fn(k,j,i))*multiplicity(j)/(total_kpoints*betas(i))
           renormalised_harmonic = renormalised_harmonic &
            & -log(har_part_fn(k,j,i))*multiplicity(j)/(total_kpoints*betas(i))
         enddo
       enddo
    endif
    call print_line(result_file, 1.0_dp/(kB*betas(i))  //' '// &
                               & renormalised_harmonic //' '// &
                               & renormalised_eigenvals)
  enddo
end subroutine
end module

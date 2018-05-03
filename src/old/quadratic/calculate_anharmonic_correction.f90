module calculate_anharmonic_correction_module
  use common_module
  implicit none
contains

subroutine calculate_anharmonic_correction(structure,structure_grid,qpoints, &
   & Nbasis,harmonic,eigenvals,result_file)
  implicit none
  
  type(StructureData), intent(in)    :: structure
  type(StructureData), intent(in)    :: structure_grid
  type(QpointData),    intent(in)    :: qpoints(:)
  integer,             intent(in)    :: Nbasis
  real(dp),            intent(in)    :: harmonic(:,:,:)
  real(dp),            intent(in)    :: eigenvals(:,:,:)
  type(OFile),         intent(inout) :: result_file
  
  ! Temperature parameters
  real(dp), parameter :: dtemperature    = 50.0_dp ! delta T (K)
  integer,  parameter :: no_temperatures = 21
  
  ! Temperatures
  real(dp) :: betas(no_temperatures)    ! {1/kB*T}
  
  integer :: no_modes
  integer :: no_qpoints
  
  ! Working variables
  integer :: i,j,k,l
  integer :: total_qpoints
  real(dp) :: renormalised_eigenvals,renormalised_harmonic
  real(dp),allocatable :: part_fn(:,:,:),har_part_fn(:,:,:)
  
  ! Get qpoint data
  no_modes = structure%no_modes
  no_qpoints = size(qpoints)
  total_qpoints = structure_grid%sc_size
  
  ! Calculate thermal energies
  do i=1,no_temperatures
    betas(i) = 1.0_dp/((i-1)*dtemperature*KB_IN_AU)
  enddo

  ! Calculate partition function
  allocate(part_fn(no_modes,no_qpoints,no_temperatures))
  allocate(har_part_fn(no_modes,no_qpoints,no_temperatures))
  part_fn=0.0_dp
  har_part_fn=0.0_dp
  do i=1,no_temperatures
    do j=1,no_qpoints
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
      do j=1,no_qpoints
        do k=1,no_modes
          renormalised_eigenvals=renormalised_eigenvals+&
           &eigenvals(1,k,j)*size(qpoints(j)%gvectors)/total_qpoints
          renormalised_harmonic=renormalised_harmonic+&
           &harmonic(1,k,j)*size(qpoints(j)%gvectors)/total_qpoints
        enddo ! k
      enddo ! j
    else ! T /= 0
       do j=1,no_qpoints
         do k=1,no_modes
           renormalised_eigenvals = renormalised_eigenvals &
            & -log(part_fn(k,j,i))*size(qpoints(j)%gvectors)/(total_qpoints*betas(i))
           renormalised_harmonic = renormalised_harmonic &
            & -log(har_part_fn(k,j,i))*size(qpoints(j)%gvectors)/(total_qpoints*betas(i))
         enddo
       enddo
    endif
    call result_file%print_line( 1.0_dp/(KB_IN_AU*betas(i))  //' '// &
                               & renormalised_harmonic       //' '// &
                               & renormalised_eigenvals)
  enddo
end subroutine
end module

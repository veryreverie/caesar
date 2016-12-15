module calculate_anharmonic_module
  use constants, only : dp
  implicit none
  
contains

subroutine calculate_anharmonic(multiplicity, no_modes, Nbasis, harmonic,&
    &eigenvals, result_file_unit)
  use utils,     only : i2s
  use constants, only : kB
  implicit none
  
  integer,  intent(in) :: multiplicity(:)
  integer,  intent(in) :: no_modes
  integer,  intent(in) :: Nbasis
  real(dp), intent(in) :: harmonic(:,:,:)
  real(dp), intent(in) :: eigenvals(:,:,:)
  integer,  intent(in) :: result_file_unit
  
  integer :: no_kpoints ! number of kpoints, = size(multiplicity)
  
  ! temperature variables
  real(dp), parameter :: dtemperature=50.0d0       ! delta T (K)
  integer,  parameter :: no_temperatures=21
  real(dp)            :: thermals(no_temperatures) ! {kB*T}
  real(dp)            :: thermal

  ! Working variables
  integer :: i,j,k,l
  integer :: total_kpoints
  real(dp) :: renormalised_eigenvals,renormalised_harmonic
  real(dp),allocatable :: part_fn(:,:,:),har_part_fn(:,:,:)
  
  ! get the number of kpoints
  no_kpoints = size(multiplicity)
  
  ! calculate total kpoints
  total_kpoints = sum(multiplicity)

  ! Allocate various arrays
  allocate(part_fn(no_kpoints,no_modes,21))
  allocate(har_part_fn(no_kpoints,no_modes,21))
  
  ! calculate thermal energies
  do i=1,size(thermals)
    thermals(i) = (i-1)*dtemperature*kB
  enddo

  ! Calculate partition function
  part_fn=0.d0
  har_part_fn=0.d0
  do l=1,size(thermals)
    thermal = thermals(l)
    do i=1,no_kpoints
      do j=1,no_modes
        do k=1,Nbasis
          part_fn(i,j,l)=part_fn(i,j,l)+exp(-eigenvals(i,j,k)/thermal)
          har_part_fn(i,j,l)=har_part_fn(i,j,l)+exp(-harmonic(i,j,k)/thermal)
        enddo ! k
      enddo ! j
    enddo ! i
  enddo ! l

  ! Calculate anharmonic vibrational correction
  do l=1,size(thermals)
    thermal = thermals(l)
    renormalised_eigenvals=0.d0
    renormalised_harmonic=0.d0
    if (l==1) then ! T=0
      do i=1,no_kpoints
        do j=1,no_modes
          renormalised_eigenvals=renormalised_eigenvals+&
           &eigenvals(i,j,1)*multiplicity(i)/total_kpoints
          renormalised_harmonic=renormalised_harmonic+&
           &harmonic(i,j,1)*multiplicity(i)/total_kpoints
        enddo ! j
      enddo ! i
    else ! T /= 0
       do i=1,no_kpoints
         do j=1,no_modes
           renormalised_eigenvals=renormalised_eigenvals+&
            &(-thermal)*log(part_fn(i,j,l))*multiplicity(i)/total_kpoints
           renormalised_harmonic=renormalised_harmonic+&
            &(-thermal)*log(har_part_fn(i,j,l))*multiplicity(i)/total_kpoints
         enddo ! j
       enddo ! i
    endif ! T>0
    write(result_file_unit,*) thermal/kB,&
                            & renormalised_harmonic,&
                            & renormalised_eigenvals
  enddo ! l
end subroutine
end module

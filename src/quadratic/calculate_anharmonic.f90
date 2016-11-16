module calculate_anharmonic_module
  implicit none
contains

subroutine calculate_anharmonic()
  implicit none
  integer,parameter :: dp=kind(1.d0)
  real :: temperature=0.0,dtemperature=50.0,thermal
  ! Hard coded numerical parameters
  ! If this value is changed, it also needs to be changed in 'vscf_1d.f90'
  integer :: Nbasis=20 

  ! Input variables
  integer :: no_kpoints,no_modes
  integer,allocatable :: multiplicity(:)
  real(dp),allocatable :: eigenvals(:,:,:),harmonic(:,:,:)

  ! Working variables
  integer :: i,j,k,l
  integer :: total_kpoints,dump_int
  real(dp) :: dump_real,renormalised_eigenvals,renormalised_harmonic
  real(dp),allocatable :: part_fn(:,:,:),har_part_fn(:,:,:)
  character(80) :: filename,c_kp,c_at
  logical :: file_exists
  

  ! Read in general input
  open(1,file='input.dat')
  read(1,*)no_kpoints,no_modes
  close(1)

  ! Allocate various arrays
  allocate(multiplicity(no_kpoints))
  allocate(eigenvals(no_kpoints,no_modes,Nbasis))
  allocate(harmonic(no_kpoints,no_modes,Nbasis))
  allocate(part_fn(no_kpoints,no_modes,21))
  allocate(har_part_fn(no_kpoints,no_modes,21))

  ! Read in multiplicity
  open(1,file='ibz.dat')
  do i=1,no_kpoints
    read(1,*)dump_real,dump_real,dump_real,multiplicity(i)
  enddo ! i
  close(1)
  total_kpoints=0
  do i=1,no_kpoints
    total_kpoints=total_kpoints+multiplicity(i)
  enddo ! i

  ! Read in anharmonic eigenvals
  harmonic=0.d0; eigenvals=0.d0
  do i=1,no_kpoints
    if(i==1)then
      do j=1,no_modes
        write(c_kp,*)i; c_kp=adjustl(c_kp)
        write(c_at,*)j; c_at=adjustl(c_at)
        filename=trim('eigenvals.')//trim(c_kp)//trim('.')//trim(c_at)//trim('.dat')
        inquire(file=filename,exist=file_exists)
        if(file_exists)then
          open(1,file=filename)
          do k=1,Nbasis
            read(1,*)dump_int,harmonic(i,j,k),eigenvals(i,j,k)
          enddo ! k
          close(1)
        endif
      enddo ! j
    else
      do j=1,no_modes
        write(c_kp,*)i; c_kp=adjustl(c_kp)
        write(c_at,*)j; c_at=adjustl(c_at)
        filename=trim('eigenvals.')//trim(c_kp)//trim('.')//trim(c_at)//trim('.dat')
        open(1,file=filename)
        do k=1,Nbasis
          read(1,*)dump_int,harmonic(i,j,k),eigenvals(i,j,k)
        enddo ! k
        close(1)
      enddo ! j
    endif ! i==1 (acoustic modes)
  enddo ! i

  ! Calculate partition function
  temperature=-dtemperature
  part_fn=0.d0
  har_part_fn=0.d0
  do l=1,21
    temperature=temperature+dtemperature
    thermal=temperature*8.6173324E-5
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
  open(1,file='anharmonic_correction.dat')
  temperature=-dtemperature
  do l=1,21 ! loop over temperature
    renormalised_eigenvals=0.d0; renormalised_harmonic=0.d0
    temperature=temperature+dtemperature
    if(temperature<1.d-5)then
      do i=1,no_kpoints
        do j=1,no_modes
          renormalised_eigenvals=renormalised_eigenvals+&
           &eigenvals(i,j,1)*multiplicity(i)/total_kpoints
          renormalised_harmonic=renormalised_harmonic+&
           &harmonic(i,j,1)*multiplicity(i)/total_kpoints
        enddo ! j
      enddo ! i
    else
       thermal=temperature*8.6173324E-5
       do i=1,no_kpoints
         do j=1,no_modes
           renormalised_eigenvals=renormalised_eigenvals+&
            &(-thermal)*log(part_fn(i,j,l))*multiplicity(i)/total_kpoints
           renormalised_harmonic=renormalised_harmonic+&
            &(-thermal)*log(har_part_fn(i,j,l))*multiplicity(i)/total_kpoints
         enddo ! j
       enddo ! i
    endif ! T>0
    write(1,*)temperature,renormalised_harmonic,renormalised_eigenvals
  enddo ! l
  close(1)
end subroutine
end module

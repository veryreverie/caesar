module calculate_bs_module
  implicit none
contains

subroutine calculate_bs()
  implicit none
  integer,parameter :: dp=kind(1.d0)
  ! Hard-coded variables
  real :: temperature=0.0,dtemperature=50.0,thermal
  ! Input variables
  integer :: no_kpoints,no_modes,degeneracy
  real,allocatable :: bands(:,:),multiplicity(:),frequency(:,:),deformation(:,:)
  real :: mapping_amplitude
  ! Working variables
  integer :: i,j,k,total_kpoints
  real :: renormalised_band,renormalised_band_kpoint
  real,allocatable :: db1(:),db2(:)
  real :: b1,b2,b0,dump_r
  character(80) :: filename1,filename2,c_kp,c_at

  ! Read in general input
  open(1,file='input.dat')
  read(1,*)no_kpoints,no_modes,degeneracy
  close(1)
  allocate(bands(no_kpoints,no_modes))
  allocate(multiplicity(no_kpoints))
  allocate(frequency(no_kpoints,no_modes))
  allocate(deformation(no_kpoints,no_modes))
  bands=0.0
  frequency=1.0
  deformation=0.0
  allocate(db1(degeneracy),db2(degeneracy))
  
  open(1,file='mapping.dat')
  read(1,*)mapping_amplitude
  close(1)
  !write(*,*)'The mapping amplitude is',mapping_amplitude
  
  open(1,file='ibz.dat')
  do i=1,no_kpoints 
    read(1,*)dump_r,dump_r,dump_r,multiplicity(i)
  enddo ! i
  close(1)
  total_kpoints=0
  do i=1,no_kpoints
    total_kpoints=total_kpoints+multiplicity(i)
  enddo ! i
  !write(*,*)'The total number of k-points is:',total_kpoints

  ! Read in bands
  do i=1,no_kpoints
    if(i==1)then
      do j=4,no_modes
        write(c_kp,*)i; c_kp=adjustl(c_kp)
        write(c_at,*)j; c_at=adjustl(c_at)
        filename1=trim('bs.')//trim(c_kp)//trim('.')//trim(c_at)//trim('.dat')
        filename2=trim('frequency.')//trim(c_kp)//trim('.')//trim(c_at)//trim('.dat')
        open(1,file=filename2)
        read(1,*)frequency(i,j)
        close(1)
        open(1,file=filename1)
        b1=0.0
        do k=1,degeneracy
          read(1,*)db1(k) 
          b1=b1+(db1(k)/degeneracy)
        enddo ! k
        b2=0.0
        do k=1,degeneracy
          read(1,*)db2(k) 
          b2=b2+(db2(k)/degeneracy)
        enddo ! k
        read(1,*)b0 
        close(1)
        bands(i,j)=(b1+b2)/2.0-b0      
        !write(*,*)i,j,bands(i,j),frequency(i,j),multiplicity(i)
      enddo ! j
    else
      do j=1,no_modes
        write(c_kp,*)i; c_kp=adjustl(c_kp)
        write(c_at,*)j; c_at=adjustl(c_at)
        filename1=trim('bs.')//trim(c_kp)//trim('.')//trim(c_at)//trim('.dat')
        filename2=trim('frequency.')//trim(c_kp)//trim('.')//trim(c_at)//trim('.dat')
        open(1,file=filename2)
        read(1,*)frequency(i,j)
        close(1)
        open(1,file=filename1)
        b1=0.0
        do k=1,degeneracy
          read(1,*)db1(k) 
          b1=b1+(db1(k)/degeneracy)
        enddo ! k
        b2=0.0
        do k=1,degeneracy
          read(1,*)db2(k) 
          b2=b2+(db1(k)/degeneracy)
        enddo ! k
        read(1,*)b0 
        close(1)
        bands(i,j)=(b1+b2)/2.0-b0      
      enddo ! j
    endif ! skip acoustic modes
  enddo ! i

  ! Calculate deformation potential
  do i=1,no_kpoints
    do j=1,no_modes
      deformation(i,j)=bands(i,j)/((mapping_amplitude/sqrt(2.0*abs(frequency(i,j))))**2)/(2.0*abs(frequency(i,j))) !/27.21138602
      !write(*,*)i,j,deformation(i,j)
    enddo ! j
  enddo ! i

  ! Calculate quadratic vibrational correction
  open(1,file='band_gap_correction.dat')
  open(2,file='bg_correction_kp.dat')
  temperature=-dtemperature
  do k=1,21  ! loop over temperature
    renormalised_band=0.0
    temperature=temperature+dtemperature
    if(temperature<1.d-5)then
      do i=1,no_kpoints
        write(2,*)'k-point',i
        renormalised_band_kpoint=0.0
        do j=1,no_modes
          renormalised_band=renormalised_band+deformation(i,j)*multiplicity(i)/total_kpoints
          renormalised_band_kpoint=renormalised_band_kpoint+deformation(i,j)
          write(2,*)i,j,deformation(i,j)
        enddo ! j
        write(2,*)i,renormalised_band_kpoint
      enddo ! i
    else
      thermal=temperature*8.6173324E-5!/3.1577464E5
      do i=1,no_kpoints
        do j=1,no_modes
          renormalised_band=renormalised_band+&
          &deformation(i,j)*(1.0+2.0/(exp(frequency(i,j)/thermal)-1))*multiplicity(i)/total_kpoints
        enddo ! j
      enddo ! i
    endif ! temperature
    write(1,*)temperature,renormalised_band 
  enddo ! k
  close(1)
  close(2)
end subroutine
end module

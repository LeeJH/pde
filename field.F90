! Copyright Li Jiahui@ifts <jiahuili@zju.edu.cn>
!
! This code is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details:
! <http://www.gnu.org/licenses/>.

subroutine init_field

  use parameters
  use input_parameters
  use topology
  use grid
  use udata
  use solver 
  use field
  
  implicit none
  
  integer :: ijk,ilen
  character(len=128) :: file_name
  character(len=12) :: varname
  character(len=4) :: cn
  real(rk), dimension(:), allocatable :: tmp
  
  ilen=mxl*myl
  
  allocate(rho(mxl*myl*mzl),stat=ierr)
  if(ierr .ne. 0) then
  write(*,*) "Memory allocating error."
  stop
  endif
  
  allocate(phi(mxsub,mysub,mzsub),stat=ierr)
  if(ierr .ne. 0) then
  write(*,*) "Memory allocating error."
  stop
  endif

#ifdef _test_poisson
  call solve_phi
  call fucomm(phi,1)
  if(ibndryx==0 .and. npx==1) call periodx(phi, 1)
  if(ibndryy==0 .and. npy==1) call periody(phi, 1)
  if(ibndryz==0 .and. npz==1) call periodz(phi, 1)


  allocate(tmp(mxl*myl*mzl),stat=ierr)
  if(ierr .ne. 0) then
  write(*,*) "Memory allocating error."
  stop
  endif

  ijk=1
  do k=mzs,mze
  do j=mys,mye
  do i=mxs,mxe
  tmp(ijk)=phi(i,j,k)
  ijk=ijk+1
  enddo
  enddo
  enddo

  file_name="pois"//cn(irank)//cn(0)//".nc"
  varname="ux"
  call dump1d(file_name,tmp,varname,mxl*myl*mzl)
  deallocate(tmp,stat=ierr)
#endif

  return
    
end subroutine init_field


subroutine solve_phi

  use parameters
  use input_parameters
  use topology
  use grid
  use udata
  use solver 
  use field
  
  implicit none
  
  integer :: ijk,ilen
  
  ilen=mxl*myl

  do k=mzs,mze
  ijk=1

  if(ipx==-1) then !Set zeros on the boundary
  do j=mys,mye
  rho(ijk)=0.d0
  ijk=ijk+1
  do i=mxs+1,mxe
  rho(ijk)=x(i,j,k,1)
  ijk=ijk+1
  enddo
  enddo
 
  else if(ipx==-1) then
  do j=mys,mye
  do i=mxs,mxe-1
  rho(ijk)=x(i,j,k,1)
  ijk=ijk+1
  enddo
  rho(ijk)=0.d0
  ijk=ijk+1
  enddo

  else
  do j=mys,mye
  do i=mxs,mxe
  rho(ijk)=x(i,j,k,1)
  ijk=ijk+1
  enddo
  enddo

  endif
  
  call solve_poisson(rho,indglxy,ilen)

  ijk=1
  do j=mys,mye
  do i=mxs,mxe
  phi(i,j,k)=rho(ijk)
  ijk=ijk+1
  enddo
  enddo

  enddo

  return
    
end subroutine solve_phi


subroutine solve3d

  use parameters
  use input_parameters
  use topology
  use grid
  use udata
  use solver 
  use field
  
  implicit none
  
  integer :: ijk
  
  ijk=1
  do i=mxs,mxe
  do j=mys,mye
  do k=mzs,mze
  rho(ijk)=x(i,j,k,1)
  ijk=ijk+1
  enddo
  enddo

  enddo

  call fftpoisson(rho,mxl,myl,mzl)  

  ijk=1
  do i=mxs,mxe
  do j=mys,mye
  do k=mzs,mze
  phi(i,j,k)=rho(ijk)
  ijk=ijk+1
  enddo
  enddo
  enddo
  
  return
    
end subroutine solve3d

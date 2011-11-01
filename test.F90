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

!Initialize initial profile and solvers
subroutine test_poisson

  use parameters
  use input_parameters
  use topology
  use grid
  use udata
  use solver 
  use field
  
  implicit none
  
  integer :: isol,ilen
  character(len=128) :: file_name
  character(len=12) :: varname
  character(len=4) :: cn
    
  ilen=mxl*myl
 
  do k=mzs,mze
  isol=1
  do j=mys,mye
  do i=mxs,mxe
  rho(isol)=x(i,j,k,1)
  isol=isol+1
  enddo
  enddo
  
  call solve_poisson(rho,indglxy,ilen)
  enddo

  file_name="pois"//cn(irank)//".nc"
  varname="phi"
  call dump1d(file_name,rho,varname,ilen)
  
  return
    
end subroutine test_poisson


subroutine test_field

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
  do j=mys,mye
  do i=mxs,mxe
  rho(ijk)=x(i,j,k,1)
  ijk=ijk+1
  enddo
  enddo
  
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
    
end subroutine test_field

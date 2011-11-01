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

subroutine dumpgrid

  use grid
  use input_parameters
  use topology
  
  implicit none
  include 'netcdf.inc'
  
  character(len=128) :: file_name
  character(len=4) :: cn
  
  integer, parameter :: NDIMS=1, topolen=7
  
  integer :: ncid, varid, dimids(NDIMS), dimids0(1)
  integer :: varidx, varidy, varidz, varid0
  integer :: x_dimid,x_dimid0
  
  integer :: topo(topolen)
  real(rk), dimension(:), allocatable :: tmp
  integer :: II
  file_name='grid'//cn(irank)//'.nc'
  
  allocate(tmp(mxl*myl*mzl),stat=ierr)
  if(ierr .ne. 0) then
  write(*,*) "Memory allocating error."
  stop
  endif
  
  topo(2)=mxl
  topo(3)=myl
  topo(4)=mzl
  topo(1)=isize
  topo(5)=mx
  topo(6)=my
  topo(7)=mz
  
  call check(nf_create(FILE_NAME, NF_CLOBBER, ncid))

  call check(nf_def_dim(ncid, "x", mxl*myl*mzl, x_dimid))
  call check(nf_def_dim(ncid, "x0", topolen, x_dimid0))
  
  dimids=(/x_dimid/)
  dimids0(1)=x_dimid0
  
  call check(nf_def_var(ncid, "topo", NF_INT, NDIMS, dimids0, varid0))
  call check(nf_def_var(ncid, "index", NF_INT, NDIMS, dimids, varid))
  call check(nf_def_var(ncid, "xgrid", NF_DOUBLE, NDIMS, dimids, varidx))
  call check(nf_def_var(ncid, "ygrid", NF_DOUBLE, NDIMS, dimids, varidy))
  call check(nf_def_var(ncid, "zgrid", NF_DOUBLE, NDIMS, dimids, varidz))
  call check(nf_enddef(ncid))
  
  call check(nf_put_var_int(ncid, varid0, topo))
  call check(nf_put_var_int(ncid, varid, indgl3d))
  
  II=1
  do k=mzs,mze
  do j=mys,mye
  do i=mxs,mxe
  tmp(II)=xx(i)
  enddo
  enddo
  enddo
  call check(nf_put_var_double(ncid, varidx, tmp))
  
  II=1
  do k=mzs,mze
  do j=mys,mye
  do i=mxs,mxe
  tmp(II)=yy(j)
  enddo
  enddo
  enddo
  call check(nf_put_var_double(ncid, varidy, tmp))
  
  II=1
  do k=mzs,mze
  do j=mys,mye
  do i=mxs,mxe
  tmp(II)=zz(k)
  enddo
  enddo
  enddo
  call check(nf_put_var_double(ncid, varidz, tmp))
  
  call check(nf_close(ncid))
  
  deallocate(tmp)
  
  return
  
end subroutine dumpgrid
  
subroutine dump1d(file_name,var,varname,ixdim)
  
  implicit none
  include 'netcdf.inc'
  
  character(len=128) :: file_name
  character(len=12) :: varname
  
  integer, parameter :: NDIMS=1
  
  integer :: ncid, varid, dimids(NDIMS)
  integer :: x_dimid
  integer :: ixdim
  real*8 var(ixdim)
  
  
  call check(nf_create(FILE_NAME, NF_CLOBBER, ncid))

  call check(nf_def_dim(ncid, "x", ixdim, x_dimid))
  
  dimids=(/x_dimid/)
  
  call check(nf_def_var(ncid, varname, NF_DOUBLE, NDIMS, dimids, varid))
  call check(nf_enddef(ncid))
  
  call check(nf_put_var_double(ncid, varid, var))
  
  
  call check(nf_close(ncid))
  
end subroutine dump1d

subroutine dump2d(file_name,var,varname,ixdim,iydim)
  
  implicit none
  include 'netcdf.inc'
  
  character(len=128) :: file_name
  character(len=12) :: varname
  
  integer, parameter :: NDIMS=2
  
  integer :: ncid, varid, dimids(NDIMS)
  integer :: x_dimid,y_dimid
  integer :: ixdim,iydim
  real*8  :: var(ixdim,iydim)
  
  
  call check(nf_create(FILE_NAME, NF_CLOBBER, ncid))
  
  call check(nf_def_dim(ncid, "y", iydim, y_dimid))
  call check(nf_def_dim(ncid, "x", ixdim, x_dimid))
  
  dimids=(/x_dimid, y_dimid/)
  
  call check(nf_def_var(ncid, varname, NF_DOUBLE, NDIMS, dimids, varid))
  call check(nf_enddef(ncid))
  
  call check(nf_put_var_double(ncid, varid, var))
  
  call check(nf_close(ncid))
  
end subroutine dump2d

    
subroutine dump4d(file_name)

 ! use netcdf
  use grid
  use topology
  use udata
  
  implicit none
  include 'netcdf.inc'

  character(len=128):: FILE_NAME
  character(len=4) :: cn

  integer, parameter :: NDIMS=4
  
  integer :: ncid, varid, dimids(NDIMS)
  integer :: x_dimid, y_dimid, z_dimid, var_dimid

  call check(nf_create(FILE_NAME, NF_CLOBBER, ncid))
  call check(nf_def_dim(ncid, "var", nvar, var_dimid))
  call check(nf_def_dim(ncid, "z", mzsub, z_dimid))
  call check(nf_def_dim(ncid, "y", mysub, y_dimid))
  call check(nf_def_dim(ncid, "x", mxsub, x_dimid))
  
  dimids=(/x_dimid, y_dimid, z_dimid, var_dimid/)
  
  call check(nf_def_var(ncid, "data", NF_DOUBLE, NDIMS, dimids, varid))
  call check(nf_enddef(ncid))
 
  call check(nf_put_var_double(ncid, varid, x))
  
  call check(nf_close(ncid))
  
end subroutine dump4d

subroutine check(status)
include 'netcdf.inc'

  integer, intent(in) :: status
  
  if(status /= nf_noerr) then
    print *, 'Error', nf_strerror(status)
    stop 2
  endif
 
end subroutine check

character*4 function cn(n)
!
!-----assume that n is no greater than 999
!
!
!-----separate the digits
!
  implicit none
  integer :: n, n1, n2, n3, n4
      n1=n/1000
      n2=(n-1000*n1)/100
      n3=(n-1000*n1-100*n2)/10
      n4=n-1000*n1-100*n2-10*n3
!
!-----stick together cn using char function
!
      n1=n1+48
      n2=n2+48
      n3=n3+48
      n4=n4+48
      cn(1:1)=char(n1)
      cn(2:2)=char(n2)
      cn(3:3)=char(n3)
      cn(4:4)=char(n4)
!
      return
      end


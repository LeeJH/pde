
      program main
      implicit none
      
      include 'netcdf.inc'
      
      character(len=128) :: file_name
      character(len=4) :: cn, filein
      character(len=128) :: output
      
      integer, parameter :: topolen=7
      integer :: ncid, varid,varidx,varidy,varidz
      integer :: i,j,k,l,II,JJ,maxlen
      integer :: topo(topolen), isize, varnum, fileid

      integer*8 ilen
      integer, dimension(:,:), allocatable :: topotol
      integer*4, dimension(:,:), allocatable :: indgl
      integer*4, dimension(:), allocatable :: ind
      real*8, dimension(:), allocatable :: xxlocal,xlocal
      real*8, dimension(:,:), allocatable:: grid,x
      
      write(*,*) "Collecting grid ..."
      file_name='grid'//cn(0)//'.nc'
      
      call check(nf_open(file_name,NF_NOWRITE,ncid))
      
      call check(nf_inq_varid(ncid,"topo",varid))
      
      call check(nf_get_var_int(ncid,varid,topo))
      
      call check(nf_close(ncid))
      
      write(*,*) topo
      
      isize=topo(1)

      ilen=topo(5)*topo(6)*topo(7)
      
      allocate(topotol(isize,topolen))
      maxlen=0
      !*******
      do i=1,isize
      file_name='grid'//cn(i-1)//'.nc'
      
      call check(nf_open(file_name,NF_NOWRITE,ncid))
      
      call check(nf_inq_varid(ncid,"topo",varid))
      call check(nf_get_var_int(ncid,varid,topo))
      
      call check(nf_close(ncid))
      
      do II=1,topolen
      topotol(i,II)=topo(II)
      enddo
      
      maxlen=max(maxlen,topo(2)*topo(3)*topo(4))
      enddo
      !********
      allocate(indgl(isize,maxlen),grid(ilen,3))
      
      do i=1,isize
      allocate(ind(topotol(i,2)*topotol(i,3)*topotol(i,4)))
      file_name='grid'//cn(i-1)//'.nc'
      
      call check(nf_open(file_name,NF_NOWRITE,ncid))
      
      call check(nf_inq_varid(ncid,"index",varid))
      
      call check(nf_get_var_int(ncid,varid,ind))
      do II=1,topotol(i,2)*topotol(i,3)*topotol(i,4)
      indgl(i,II)=ind(II)
      enddo
      
      call check(nf_close(ncid))
      deallocate(ind)
      enddo
      
      do i=1,isize
      allocate(xxlocal(topotol(i,2)*topotol(i,3)*topotol(i,4)))
      file_name='grid'//cn(i-1)//'.nc'
      
      call check(nf_open(file_name,NF_NOWRITE,ncid))
      
      call check(nf_inq_varid(ncid,"xgrid",varidx))
      call check(nf_inq_varid(ncid,"ygrid",varidy))
      call check(nf_inq_varid(ncid,"zgrid",varidz))
      
!      call check(nf_inq_varid(ncid,"xgrid",varid))
      call check(nf_get_var_double(ncid,varidx,xxlocal))
      do II=1,topotol(i,2)*topotol(i,3)*topotol(i,4)
      grid(indgl(i,II),1)=xxlocal(II)
      enddo
      
!      call check(nf_inq_varid(ncid,"ygrid",varid))
      call check(nf_get_var_double(ncid,varidy,xxlocal))
      do II=1,topotol(i,2)*topotol(i,3)*topotol(i,4)
      grid(indgl(i,II),2)=xxlocal(II)
      enddo
      
!      call check(nf_inq_varid(ncid,"zgrid",varid))
      call check(nf_get_var_double(ncid,varidz,xxlocal))
      do II=1,topotol(i,2)*topotol(i,3)*topotol(i,4)
      grid(indgl(i,II),3)=xxlocal(II)
      enddo
      
      call check(nf_close(ncid))
      deallocate(xxlocal)
      
      enddo
      !************
      write(*,*) "Input file id and number of variables:"
      read(*,*) fileid, varnum

      allocate(x(ilen,varnum))
      
      do i=1,isize
      allocate(xlocal(topotol(i,2)*topotol(i,3)*topotol(i,4)*varnum))
      file_name=cn(fileid)//'/data'//cn(i-1)//'.nc'
      
      call check(nf_open(file_name,NF_NOWRITE,ncid))
      
      call check(nf_inq_varid(ncid,"ux",varid))
      
      call check(nf_get_var_double(ncid,varid,xlocal))
      
      do II=1,topotol(i,2)*topotol(i,3)*topotol(i,4)
      do JJ=1,varnum
!      x(indgl(i,II),1)=xlocal(II)
      x(indgl(i,II),JJ)=xlocal(II+topotol(i,2)*topotol(i,3)*topotol(i,4)*(JJ-1))
      enddo
      enddo
      
      
      call check(nf_close(ncid))
      deallocate(xlocal)
      enddo
      
      call dump2d('init.nc',x,'x',ilen,varnum)
      return
      
      end

subroutine dump1d(file_name,var,varname,ilen)
  
  use digit
  
  implicit none
  include 'netcdf.inc'
  
  character(len=128) :: file_name
  character(len=12) :: varname
  
  integer, parameter :: NDIMS=1
  
  integer :: ncid, varid, dimids(NDIMS)
  integer :: x_dimid
  real(rk) var(*)
  integer :: ilen
  
  call check(nf_create(FILE_NAME, NF_CLOBBER, ncid))

  call check(nf_def_dim(ncid, "x", ilen, x_dimid))
  
  dimids=(/x_dimid/)
  
  call check(nf_def_var(ncid, varname, NF_DOUBLE, NDIMS, dimids, varid))
  call check(nf_enddef(ncid))
  
  call check(nf_put_var_double(ncid, varid, var))
  
  
  call check(nf_close(ncid))
  
end subroutine dump1d

subroutine dump2d(file_name,var,varname,ilen,varnum)
  
  
  implicit none
  include 'netcdf.inc'
  
  character(len=128) :: file_name
  character(len=12) :: varname
  
  integer :: varnum
  integer :: ncid, varid, dimids(varnum),ilen
  integer :: x_dimid,y_dimid
  real*8 :: var(ilen,varnum)
  
  call check(nf_create("data.nc", NF_CLOBBER, ncid))
  call check(nf_def_dim(ncid, "y", varnum, y_dimid))
  call check(nf_def_dim(ncid, "x", ilen, x_dimid))
  
  dimids=(/x_dimid, y_dimid/)
  
  call check(nf_def_var(ncid, "ux", NF_DOUBLE, varnum, dimids, varid))
  call check(nf_enddef(ncid))
  
  call check(nf_put_var_double(ncid, varid, var))
  
  call check(nf_close(ncid))
  
end subroutine dump2d

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

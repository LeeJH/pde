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
subroutine init

  use parameters
  use input_parameters
  use topology
  use grid
  use udata
  use solver 
  
  implicit none
  
  integer :: isol,ilen
  character(len=128) :: file_name
  character(len=12) :: varname
  character(len=4) :: cn
  
  ilen=mxl*myl
  if(irank.eq.0) write(*,*) "Initializing grid ..."
  call init_grid
  
  if(irank.eq.0) write(*,*) "Initializing profile ..."
  call init_profile
  
  if(irank.eq.0) write(*,*) "Initializing poisson solver ..."
  call init_poisson(indglxy,ilen,mx,my,dx(mxs),dy(mys),ptol,commz)
  
  if(irank.eq.0) write(*,*) "Initializing time integrator ..."
  call init_solver

  if(irank.eq.0) write(*,*) "Initializing fields ..."
  call init_field

  if(irank.eq.0) write(*,*) "Initializing fftw ..."
  call init_fft
!poisson3d  call init_fft2d
  
  call MPI_Barrier(comm,ierr)
  
  call dumpgrid
  
  if(irank.eq.0) write(*,*) "Initialization clear."
  return
    
end subroutine init


subroutine init_grid

  use parameters
  use input_parameters
  use topology
  use grid
  
  implicit none
  
  real(rk) :: dxtmp, dytmp, dztmp
  !For fft2d.c
  integer :: commxc,commyzc
  integer :: npxc, npyc, npzc
  integer :: mxc, myc, mzc
  integer :: mxlc, mylc, mzlc
  integer :: ipxc, ipyc, ipzc
  integer :: ixglobalc, iyglobalc, izglobalc
  real(rk) :: dx0,dy0,dz0
  
  common /topo/ commxc,commyzc,npxc,npyc,npzc,mxc,myc,mzc,&
                 mxlc,mylc,mzlc,ipxc,ipyc,ipzc,&
                 ixglobalc,iyglobalc,izglobalc,&
                 dx0,dy0,dz0
                 !Pay attention to the order of the variables
  !End
  
  allocate(xx(mxl+2*mgx),dx(mxl+2*mgx),yy(myl+2*mgy),&
  	   dy(myl+2*mgy),zz(mzl+2*mgz),dz(mzl+2*mgz), stat=ierr)
  if(ierr .ne. 0) then
  write(*,*) "Memory allocating error."
  stop
  endif
  
  if(ibndryx>0) then
    dxtmp=(xmax-xmin)/(mx-1)
  else
    dxtmp=(xmax-xmin)/mx
  endif
    
  if(ibndryy>0) then
    dytmp=(ymax-ymin)/(my-1)
  else
    dytmp=(ymax-ymin)/my
  endif
  
  if(ibndryz>0) then
    dztmp=(zmax-zmin)/(mz-1)
  else
    dztmp=(zmax-zmin)/mz
  endif
  
  do i=1, mxl+2*mgx
    dx(i)=dxtmp
    xx(i)=xmin+dxtmp*(ixglobal-mgx+i-1-1)
  enddo
  
  do i=1, myl+2*mgy
    dy(i)=dytmp
    yy(i)=ymin+dytmp*(iyglobal-mgy+i-1-1)
  enddo
  
  do i=1, mzl+2*mgz
    dz(i)=dztmp
    zz(i)=zmin+dztmp*(izglobal-mgz+i-1-1)
  enddo
  
  !For fft2d.c
  commxc=commx
  commyzc=commyz
  npxc=npx
  npyc=npy
  npzc=npz
  mxc=mx
  myc=my
  mzc=mz
  mxlc=mxl
  mylc=myl
  mzlc=mzl
  ipxc=ipx
  ipyc=ipy
  ipzc=ipz
  ixglobalc=ixglobal
  iyglobalc=iyglobal
  izglobalc=izglobal
  dx0=dx(1)
  dy0=dy(1)
  dz0=dz(1)
  !End
  
  call MPI_Barrier(comm, ierr)

#ifdef debug
if(irank.eq.0) write(*,*) "Grid initialized."
#endif

  return
    
end subroutine init_grid

subroutine init_profile

  use parameters
  use input_parameters
  use topology
  use grid
  use udata 
  use field
  
  implicit none
  
  integer :: lxy
  lxy=mxl*myl
  
  allocate(x(mxsub,mysub,mzsub,nvar),&
           xdot(mxsub,mysub,mzsub,nvar), stat=ierr)
  if(ierr .ne. 0) then
  write(*,*) "Memory allocating error."
  stop
  endif
  
  !Block frequently modified
  allocate(Ni0(mxsub,mysub),B(mxsub,mysub,mzsub,4),&
           Te0(mxsub,mysub), uepar(mxsub,mysub,mzsub),stat=ierr)
  if(ierr .ne. 0) then
  write(*,*) "Memory allocating error."
  stop
  endif
  
  do k=1, mze+mgz
    do j=1, mye+mgy
      do i=1, mxe+mgx
!      x(i,j,k,1)=dsin(2*pi*xx(i)/(xmax-xmin))&
                !*dsin(2*pi*yy(j)/(ymax-ymin))!&
!                *dsin(2*pi*zz(k)/(zmax-zmin))

!      x(i,j,k,1)=dsin(2*pi*xx(i)/(xmax-xmin))*(1./(dcosh(zz(k)/zmax/0.1)**2))
      x(i,j,k,1)=(1./(dcosh(xx(i)/xmax/0.1)**2))*(1./(dcosh(zz(k)/zmax/0.1)**2))
!      x(i,j,k,1)=(1./(dcosh(zz(k)/zmax/0.1)**2))
!      x(i,j,k,1)=(1./(dcosh(xx(i)/xmax/0.1)**2))*(dsin(2*pi*zz(k)/(zmax-zmin)))
      x(i,j,k,2)=0.0*(1./(dcosh(xx(i)/xmax/0.2)**2))*(1./(dcosh(zz(k)/zmax/0.2)**2))

      
      xdot(i,j,k,1)=0.d0 
      xdot(i,j,k,2)=0.d0
      ! Initialize it, because the boundary is never used.
      enddo
    enddo
  enddo
  
  do j=mys, mye
    do i=mxs, mxe
      Ni0(i,j)=1.0
      Te0(i,j)=1.0
    enddo
  enddo
  
  call MPI_Barrier(comm, ierr)
  
  call fucomm(x, nvar)

  call MPI_Barrier(comm, ierr)

#ifdef debug
if(irank.eq.0) write(*,*) "Profile initialized."
#endif
  
  return

end subroutine init_profile


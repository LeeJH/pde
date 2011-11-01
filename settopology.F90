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

!Imformation of the division of the simulation 
!domain over all processes
subroutine settopology

  use input_parameters
  use topology
  
  implicit none
  
  integer*4 :: II
  integer*4 :: ixtmp, iytmp, iztmp
  
  ipz=irank/(npx*npy)
  ipy=(irank-ipz*npx*npy)/npx
  ipx=irank-ipz*npx*npy-ipy*npx
  
  if(ipz>=npz) then
    write(*,*) "ipz error!"
    stop
  endif
  
  if(ipy>=npy) then 
    write(*,*) "ipy error!"
    stop
  endif
  
  if(ipx>=npx) then
    write(*,*) "ipx error!"
    stop
  endif
  
  if(ipx<(npx-1)) then
    ixup=irank+1
  else
    ixup=irank-(npx-1)
  endif
  if((ipx .eq. npx-1) .and. ibndryx>0) ixup=MPI_PROC_NULL
  
  if(ipx>0) then
    ixdown=irank-1
  else
    ixdown=irank+(npx-1)
  endif
  if((ipx .eq. 0) .and. ibndryx>0) ixdown=MPI_PROC_NULL
  
  if(ipy<(npy-1)) then
    iyup=irank+npx
  else
    iyup=irank-(npy-1)*npx
  endif
  if((ipy .eq. npy-1) .and. ibndryy>0) iyup=MPI_PROC_NULL
  
  if(ipy>0) then
    iydown=irank-npx
  else
    iydown=irank+(npy-1)*npx
  endif
  if((ipy .eq. 0) .and. ibndryy>0) iydown=MPI_PROC_NULL
  
  if(ipz<(npz-1)) then
    izup=irank+npx*npy
  else
    izup=irank-(npz-1)*npx*npy
  endif
  if((ipz .eq. npz-1) .and. ibndryz>0) izup=MPI_PROC_NULL
  
  if(ipz>0) then
    izdown=irank-npx*npy
  else
    izdown=irank+(npz-1)*npx*npy
  endif
  if((ipy .eq. 0) .and. ibndryz>0) izdown=MPI_PROC_NULL
  
  if(mod(mx,npx)==0) then
    mxl=mx/npx
    mxs=mgx+1
    mxe=mxl+mgx
  else
    mxl=mx/npx+1
    if(ipx==(npx-1)) mxl=mx-mxl*(npx-1)
    mxs=mgx+1
    mxe=mgx+mxl
  endif
  mxsub=mxl+2*mgx
  
  if(mod(my,npy)==0) then
    myl=my/npy
    mys=mgy+1
    mye=myl+mgy
  else
    myl=my/npy+1
    if(ipy==(npy-1)) myl=my-myl*(npy-1)
    mys=mgy+1
    mye=mgy+myl
  endif
  mysub=myl+2*mgy
  
  if(mod(mz,npz)==0) then
    mzl=mz/npz
    mzs=mgz+1
    mze=mzl+mgz
  else
    mzl=mz/npz+1
    if(ipz==(npz-1)) mzl=mz-mzl*(npz-1)
    mzs=mgz+1
    mze=mgz+mzl
  endif
  mzsub=mzl+2*mgz
  
  if(ipx==0) then
    ixglobal=1
  else
    ixglobal=(ipx-1)*mxl+(mx-mxl*(npx-1))+1
  endif
  
  if(ipy==0) then
    iyglobal=1
  else
    iyglobal=(ipy-1)*myl+(my-myl*(npy-1))+1
  endif
  
  if(ipz==0) then
    izglobal=1
  else
    izglobal=(ipz-1)*mzl+(mz-mzl*(npz-1))+1
  endif
  
  call MPI_Barrier(comm, ierr)
  
  !Set index for 2d poisson solver
  allocate(indglxy(mxl*myl), stat=ierr)
  if(ierr .ne. 0) then
  write(*,*) "Memory allocating error."
  stop
  endif
  
  II=1
  do j=1,myl
  do i=1,mxl
  ixtmp=ixglobal+i-1
  iytmp=iyglobal+j-1
  indglxy(II)=ixtmp+mx*(iytmp-1)-1
  II=II+1
  enddo
  enddo
  
  !Set index for output
  allocate(indgl3d(mxl*myl*mzl),stat=ierr)
  if(ierr .ne. 0) then
  write(*,*) "Memory allocating error."
  stop
  endif
  
  II=1
  do k=1,mzl
  do j=1,myl
  do i=1,mxl
  ixtmp=ixglobal+i-1
  iytmp=iyglobal+j-1
  iztmp=izglobal+k-1
  indgl3d(II)=ixtmp+mx*(iytmp-1)+mx*my*(iztmp-1)
  II=II+1
  enddo
  enddo
  enddo
  
  !Split the MPI_COMM_WORLD according to 'ipz'
  call MPI_COMM_SPLIT(comm,ipz,irank,commz,ierr)
  call MPI_COMM_RANK(commz, irankz, ierr)
  call MPI_COMM_SIZE(commz, isizez, ierr)
  
  !Split the MPI_COMM_WORLD according to 'ipx'
  call MPI_COMM_SPLIT(comm,ipx,irank,commx,ierr)
  call MPI_COMM_RANK(commx, irankx, ierr)
  call MPI_COMM_SIZE(commx, isizex, ierr)
  
  call MPI_COMM_SPLIT(comm,ipx+ipy*npx,ipz,commxy,ierr)
  call MPI_COMM_RANK(commxy, irankxy, ierr)
  call MPI_COMM_SIZE(commxy, isizexy, ierr)

  call MPI_COMM_SPLIT(comm,ipy+ipz*npy,ipx,commyz,ierr)
  call MPI_COMM_RANK(commyz, irankyz, ierr)
  call MPI_COMM_SIZE(commyz, isizeyz, ierr)
  
#ifdef debug
  write(*,*) irank,ipz,irankz  
#endif

  return
end subroutine settopology

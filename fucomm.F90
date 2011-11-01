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

subroutine fucomm(tmp, ndim) ! Corner ghost points are not setted

  use parameters
  use input_parameters
  use topology
  
  implicit none
  
  integer :: ndim
  integer :: iscg ! send corner ghost points
  !integer(kind=8) :: icount, ibufsize, ibufx, ibufy, ibufz
  integer :: icount, ibufsize, ibufx, ibufy, ibufz
  integer :: istatus(MPI_STATUS_SIZE)
  real(rk) :: tmp(mxl+2*mgx,myl+2*mgy,mzl+2*mgz,ndim)
  real(rk), dimension(:), allocatable :: sendbuf, recvbuf
   
              
  iscg=1
  
  ibufx=0
  ibufy=0
  ibufz=0
  
  
  if(npx>1) ibufx=(myl+iscg*2*mgy)*(mzl+iscg*2*mgz)*mgx*ndim
  if(npy>1) ibufy=(mxl+iscg*2*mgx)*(mzl+iscg*2*mgz)*mgy*ndim
  if(npz>1) ibufz=(mxl+iscg*2*mgx)*(myl+iscg*2*mgy)*mgz*ndim
  
  ibufsize=max(ibufx,ibufy,ibufz)
  allocate(sendbuf(ibufsize),recvbuf(ibufsize),stat=ierr)
  if(ierr .ne. 0) then
  write(*,*) "Memory allocating error."
  stop
  endif
  
! Communication in x direction
  if(npx>1) then
  icount=1
  do l=1,ndim
    do k=mzs-iscg*mgz,mze+iscg*mgz  
      do j=mys-iscg*mgy,mye+iscg*mgy
        do i=mxs,mxs+mgx-1        
        sendbuf(icount)=tmp(i,j,k,l)
        icount=icount+1
        enddo
      enddo
    enddo
  enddo

  call MPI_SENDRECV(sendbuf,ibufx,mpiR,ixdown,irank,recvbuf,&
               ibufx,mpiR,ixup,ixup,comm,istatus,ierr)
               
  icount=1
  do l=1,ndim
    do k=mzs-iscg*mgz,mze+iscg*mgz
      do j=mys-iscg*mgy,mye+iscg*mgy
        do i=mxe+1,mxe+mgx        
        tmp(i,j,k,l)=recvbuf(icount)
        icount=icount+1
        enddo
      enddo
    enddo
  enddo
  
  !
  icount=1
  do l=1,ndim
    do k=mzs-iscg*mgz,mze+iscg*mgz
      do j=mys-iscg*mgy,mye+iscg*mgy
        do i=mxe-mgx+1,mxe
        sendbuf(icount)=tmp(i,j,k,l)
        icount=icount+1
        enddo
      enddo
    enddo
  enddo

  call MPI_SENDRECV(sendbuf,ibufx,mpiR,ixup,irank,recvbuf,&
               ibufx,mpiR,ixdown,ixdown,comm,istatus,ierr)
               
  icount=1
  do l=1,ndim
    do k=mzs-iscg*mgz,mze+iscg*mgz
      do j=mys-iscg*mgy,mye+iscg*mgy  
        do i=mxs-mgx,mxs-1
        tmp(i,j,k,l)=recvbuf(icount)
        icount=icount+1
        enddo
      enddo
    enddo
  enddo
    
  endif ! npx>1
  
! Communication in y direction
  if(npy>1) then
  
  icount=1  
  do l=1,ndim
    do k=mzs-iscg*mgz,mze+iscg*mgz
      do j=mys,mys+mgy-1
        do i=mxs-iscg*mgx,mxe+iscg*mgx
        sendbuf(icount)=tmp(i,j,k,l)
        icount=icount+1
        enddo
      enddo
    enddo
  enddo

  call MPI_SENDRECV(sendbuf,ibufy,mpiR,iydown,irank,recvbuf,&
               ibufy,mpiR,iyup,iyup,comm,istatus,ierr)
               
  icount=1
  do l=1,ndim
    do k=mzs-iscg*mgz,mze+iscg*mgz
      do j=mye+1,mye+mgy
        do i=mxs-iscg*mgx,mxe+iscg*mgx
        tmp(i,j,k,l)=recvbuf(icount)
        icount=icount+1
        enddo
      enddo
    enddo
  enddo
  
  icount=1
  do l=1,ndim
    do k=mzs-iscg*mgz,mze+iscg*mgz
      do j=mye-mgy+1,mye
        do i=mxs-iscg*mgx,mxe+iscg*mgx      
        sendbuf(icount)=tmp(i,j,k,l)
        icount=icount+1
        enddo
      enddo
    enddo
  enddo

  call MPI_SENDRECV(sendbuf,ibufy,mpiR,iyup,irank,recvbuf,&
               ibufy,mpiR,iydown,iydown,comm,istatus,ierr)
               
  icount=1
  do l=1,ndim
    do k=mzs-iscg*mgz,mze+iscg*mgz
      do j=mys-mgy,mys-1
        do i=mxs-iscg*mgx,mxe+iscg*mgx
        tmp(i,j,k,l)=recvbuf(icount)
        icount=icount+1
        enddo
      enddo
    enddo
  enddo
    
  endif ! npy>1
  
! Communication in z direction
  if(npz>1) then
  
  icount=1  
  do l=1,ndim
    do k=mzs,mzs+mgz-1
      do j=mys-iscg*mgy,mye+iscg*mgy
        do i=mxs-iscg*mgx,mxe+iscg*mgx    
        sendbuf(icount)=tmp(i,j,k,l)
        icount=icount+1
        enddo
      enddo
    enddo
  enddo

  call MPI_SENDRECV(sendbuf,ibufz,mpiR,izdown,irank,recvbuf,&
               ibufz,mpiR,izup,izup,comm,istatus,ierr)
               
  icount=1
  do l=1,ndim
    do k=mze+1,mze+mgz
      do j=mys-iscg*mgy,mye+iscg*mgy
        do i=mxs-iscg*mgx,mxe+iscg*mgx        
        tmp(i,j,k,l)=recvbuf(icount)
        icount=icount+1
        enddo
      enddo
    enddo
  enddo
  
  icount=1
  do l=1,ndim
    do k=mze-mgz+1,mze
      do j=mys-iscg*mgy,mye+iscg*mgy
        do i=mxs-iscg*mgx,mxe+iscg*mgx
        sendbuf(icount)=tmp(i,j,k,l)
        icount=icount+1
        enddo
      enddo
    enddo
  enddo

  call MPI_SENDRECV(sendbuf,ibufz,mpiR,izup,irank,recvbuf,&
               ibufz,mpiR,izdown,izdown,comm,istatus,ierr)
               
  icount=1
  do l=1,ndim
    do k=mzs-mgz,mzs-1
      do j=mys-iscg*mgy,mye+iscg*mgy
        do i=mxs-iscg*mgx,mxe+iscg*mgx
        tmp(i,j,k,l)=recvbuf(icount)
        icount=icount+1
        enddo
      enddo
    enddo
  enddo
    
  endif ! npz>1
  
  deallocate(sendbuf)
  deallocate(recvbuf)
  
end subroutine fucomm

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

subroutine periodx(vb,ndim)

  use input_parameters
  use topology
  
  implicit none

#define d2fc2(fm,f0,fp,h0) (fm+fp-2.d0*f0)/h0/h0
#define d2fc4(fm2,fm,f0,fp,fp2,h0) (-fm2+16*fm-30*f0+16*fp-fp2)/(12.d0*h0*h0)
#define d2fp4(fm,f0,fp,fp2,fp3,h0) (11.d0*fm-20.d0*f0+6.d0*fp+4.d0*fp2-fp3)/(12.d0*h0*h0)    
#define d2fm4(fm3,fm2,fm,f0,fp,h0) (-fm3+4.d0*fm2+6.d0*fm-20.d0*f0+11.d0*fp)/(12.d0*h0*h0)   
  
  integer :: ndim
  real(rk) :: vb(mxsub, mysub, mzsub, ndim)
  
  if(npx .ne. 1) then
  write(*,*) "boundary error."
  stop 3
  endif
  
  do l=1,ndim
    do k=1,mzsub
      do j=1,mysub        
        do i=1,mgx
        vb(i,j,k,l)=vb(mxe-mgx+i,j,k,l)
        vb(mxe+i,j,k,l)=vb(mxs+i-1,j,k,l)
        enddo
      enddo
    enddo
  enddo
  
  return
  
end subroutine periodx


subroutine periody(vb,ndim)

  use input_parameters
  use topology
  
  implicit none

#define d2fc2(fm,f0,fp,h0) (fm+fp-2.d0*f0)/h0/h0
#define d2fc4(fm2,fm,f0,fp,fp2,h0) (-fm2+16*fm-30*f0+16*fp-fp2)/(12.d0*h0*h0)
#define d2fp4(fm,f0,fp,fp2,fp3,h0) (11.d0*fm-20.d0*f0+6.d0*fp+4.d0*fp2-fp3)/(12.d0*h0*h0)    
#define d2fm4(fm3,fm2,fm,f0,fp,h0) (-fm3+4.d0*fm2+6.d0*fm-20.d0*f0+11.d0*fp)/(12.d0*h0*h0)   
  
  integer :: ndim
  real(rk) :: vb(mxsub, mysub, mzsub, ndim)
  
  if(npy .ne. 1) then
  write(*,*) "boundary error."
  stop 3
  endif
  
  do l=1,ndim
    do k=1,mzsub
      do j=1,mgy
        do i=1,mxsub
        vb(i,j,k,l)=vb(i,mye-mgy+j,k,l)
        vb(i,mye+j,k,l)=vb(i,mys+j-1,k,l)
        enddo
      enddo
    enddo
  enddo
  
  return
  
end subroutine periody

subroutine periodz(vb,ndim)

  use input_parameters
  use topology
  
  implicit none

#define d2fc2(fm,f0,fp,h0) (fm+fp-2.d0*f0)/h0/h0
#define d2fc4(fm2,fm,f0,fp,fp2,h0) (-fm2+16*fm-30*f0+16*fp-fp2)/(12.d0*h0*h0)
#define d2fp4(fm,f0,fp,fp2,fp3,h0) (11.d0*fm-20.d0*f0+6.d0*fp+4.d0*fp2-fp3)/(12.d0*h0*h0)    
#define d2fm4(fm3,fm2,fm,f0,fp,h0) (-fm3+4.d0*fm2+6.d0*fm-20.d0*f0+11.d0*fp)/(12.d0*h0*h0)   
  
  integer :: ndim
  real(rk) :: vb(mxsub, mysub, mzsub, ndim)
  
  if(npz .ne. 1) then
  write(*,*) "boundary error."
  stop 3
  endif
  
  do l=1,ndim
    do k=1,mgz
      do j=1,mysub
        do i=1,mxsub     
        vb(i,j,k,l)=vb(i,j,mze-mgz+k,l)
        vb(i,j,mze+k,l)=vb(i,j,mzs+k-1,l)
        enddo
      enddo
    enddo
  enddo
  
  return
  
end subroutine periodz
!*******************************************************************************************
  
subroutine bndryx(vb, vbdot,ndim)

  use input_parameters
  use topology
  
  implicit none

#define d2fc2(fm,f0,fp,h0) (fm+fp-2.d0*f0)/h0/h0
#define d2fc4(fm2,fm,f0,fp,fp2,h0) (-fm2+16*fm-30*f0+16*fp-fp2)/(12.d0*h0*h0)
#define d2fp4(fm,f0,fp,fp2,fp3,h0) (11.d0*fm-20.d0*f0+6.d0*fp+4.d0*fp2-fp3)/(12.d0*h0*h0)    
#define d2fm4(fm3,fm2,fm,f0,fp,h0) (-fm3+4.d0*fm2+6.d0*fm-20.d0*f0+11.d0*fp)/(12.d0*h0*h0)   
  
  integer :: ndim
  real(rk) :: vb(mxsub,mysub,mzsub,ndim), vbdot(mxsub,mysub,mzsub,ndim)
  
  select case(ibndryx)
  
  case(1)
  if(ipx==0) then
  do l=1,ndim
    do k=1,mzsub
      do j=1,mysub
        do i=1,mgx
        vb(i,j,k,l)=vb(mxs,j,k,l)
        enddo
      enddo
    enddo
  enddo
  endif
  
  if(ipx==npx-1) then
  do l=1,ndim
    do k=1,mzsub        
      do j=1,mysub
        do i=1,mgx
        vb(mxe+i,j,k,l)=vb(mxe,j,k,l)
        enddo
      enddo
    enddo
  enddo
  endif
   
  case(2)
  if(ipx==0) then
  do l=1,ndim
    do k=1,mzsub
      do j=1,mysub
        vbdot(mxs,j,k,l)=0
!        vbdot(mxs+1,j,k,l)=(vbdot(mxs,j,k,l)+vbdot(mxs+2,j,k,l))/2
      enddo
    enddo
  enddo
  endif
  
  if(ipx==npx-1) then
  do l=1,ndim
    do k=1,mzsub        
      do j=1,mysub
        vbdot(mxe,j,k,l)=0
!        vbdot(mxe-1,j,k,l)=(vbdot(mxe,j,k,l)+vbdot(mxe-2,j,k,l))/2.d0
      enddo
    enddo
  enddo
  endif
  
  case default
  write(*,*) "Boundary not set."
  stop
  
  end select 
  
end subroutine bndryx

subroutine bndryy(vb, vbdot, ndim)

  use input_parameters
  use topology
  
  implicit none
  
  integer :: ndim
  real(rk) :: vb(mxsub,mysub,mzsub,ndim),vbdot(mxsub,mysub,mzsub,ndim)
  
  select case(ibndryy)
  
  case(1)
  if(ipy==0) then
  do l=1,ndim
    do k=1,mzsub
      do j=1,mgy
        do i=1,mxsub
        vb(i,j,k,l)=vb(i,mys,k,l)
        enddo
      enddo
    enddo
  enddo
  endif
  
  if(ipy==npy-1) then
  do l=1,ndim
    do k=1,mzsub
      do j=1,mgy
        do i=1,mxsub       
        vb(i,mye+j,k,l)=vb(i,mye,k,l)
        enddo
      enddo
    enddo
  enddo
  endif
  
  case default
  write(*,*) "Boundary not set."
  stop
  
  end select
  
end subroutine bndryy

subroutine bndryz(vb, ndim)

  use input_parameters
  use topology
  
  implicit none
  
  integer :: ndim
  real(rk) :: vb(mxsub,mysub,mzsub,ndim),vbdot(mxsub,mysub,mzsub,ndim)
  
  select case(ibndryz)

  case(1)
  if(ipz==0) then
  do l=1,ndim
    do k=1,mgz
      do j=1,mysub    
        do i=1,mxsub    
        vb(i,j,k,l)=vb(i,j,mzs,l)
        enddo
      enddo
    enddo
  enddo
  endif
  
  if(ipz==npz-1) then
  do l=1,ndim
    do k=1,mgz
      do j=1,mysub
        do i=1,mxsub 
        vb(i,j,mze+k,l)=vb(i,j,mze,l)
        enddo
      enddo
    enddo
  enddo
  endif
  
  case default
  write(*,*) "Boundary not set."
  stop
  
  end select
  
end subroutine bndryz

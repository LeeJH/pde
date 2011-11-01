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

subroutine fcvfun(t, ux, uxdot, ipar, rpar, ierror)

  use input_parameters
  use topology
  use grid
  use udata
  use field
  
  implicit none
  
#define d1fc2(fm,f0,fp,h0) 0.5d0*(fp-fm)/h0
#define d1fc4(fm2,fm,f0,fp,fp2,h0) (fm2-8.d0*fm+8.d0*fp-fp2)/12.d0/h0
#define d2fc2(fm,f0,fp,h0) (fm+fp-2.d0*f0)/h0/h0
#define d2fc4(fm2,fm,f0,fp,fp2,h0) (-fm2+16*fm-30*f0+16*fp-fp2)/(12.d0*h0*h0)
#define d2fp4(fm,f0,fp,fp2,fp3,h0) (11.d0*fm-20.d0*f0+6.d0*fp+4.d0*fp2-fp3)/(12.d0*h0*h0)    
#define d2fm4(fm3,fm2,fm,f0,fp,h0) (-fm3+4.d0*fm2+6.d0*fm-20.d0*f0+11.d0*fp)/(12.d0*h0*h0)   
#define indtran(i,j,k,l) ((i-mgx)+mxl*(j-mgy-1)+mxl*myl*(k-mgz-1)+mxl*myl*mzl*(nvar-1))
  !integer(kind=8) :: ijk
  integer :: ijk
  integer :: ierror
  integer(kind=8) :: ipar(*)
  real(rk) :: t
  real(rk) :: ux(*), uxdot(*), rpar(*)
  
  !Copy data from time integrator to grid
  ijk=1
  do l=1, nvar
  do k=mzs, mze
  do j=mys, mye
  do i=mxs, mxe
    x(i,j,k,l)=ux(ijk)
    ijk=ijk+1
  enddo
  enddo
  enddo
  enddo
  
  
  
  call fucomm(x,nvar)
  if(ibndryx==0 .and. npx==1) call periodx(x, nvar)
  if(ibndryy==0 .and. npy==1) call periody(x, nvar)
  if(ibndryz==0 .and. npz==1) call periodz(x, nvar)
  !*************************************************
  
!poisson3d  call solve3d
  call solve_phi
  call fftfilter
  call fucomm(phi,1)
  if(ibndryx==0 .and. npx==1) call periodx(phi, 1)
  if(ibndryy==0 .and. npy==1) call periody(phi, 1)
  if(ibndryz==0 .and. npz==1) call periodz(phi, 1)
  
  do k=mzs, mze
  do j=mys, mye
  do i=mxs, mxe
    uepar(i,j,k)=2.0*(d2fc4(x(i-2,j,k,2),x(i-1,j,k,2),x(i,j,k,2),x(i+1,j,k,2),x(i+2,j,k,2),dx(i))&
      +d2fc4(x(i,j-2,k,2),x(i,j-1,k,2),x(i,j,k,2),x(i,j+1,k,2),x(i,j+2,k,2),dy(j)))/beta
!    uepar(i,j,k)=2.0*(d2fc2(x(i-1,j,k,2),x(i,j,k,2),x(i+1,j,k,2),dx(i))&
!      +d2fc2(x(i,j-1,k,2),x(i,j,k,2),x(i,j+1,k,2),dy(j)))/beta
  enddo
  enddo
  enddo
!Inner boundary treatment
  if(ipx==-1) then
  i=mxs+1
  do k=mzs, mze
  do j=mys, mye
!    uepar(i,j,k)=2.0*(d2fp4(x(i-1,j,k,2),x(i,j,k,2),x(i+1,j,k,2),x(i+2,j,k,2),x(i+3,j,k,2),dx(i))+d2fc4(x(i,j-2,k,2),x(i,j-1,k,2),x(i,j,k,2),x(i,j+1,k,2),x(i,j+2,k,2),dy(j)))/beta
    uepar(mxs,j,k)=0.d0
  enddo
  enddo
  endif

  if(ipx==-1) then
  i=mxe-1
  do k=mzs, mze
  do j=mys, mye
!    uepar(i,j,k)=2.0*(d2fm4(x(i-3,j,k,2),x(i-2,j,k,2),x(i-1,j,k,2),x(i,j,k,2),x(i+1,j,k,2),dx(i))+d2fc4(x(i,j-2,k,2),x(i,j-1,k,2),x(i,j,k,2),x(i,j+1,k,2),x(i,j+2,k,2),dy(j)))/beta
    uepar(mxe,j,k)=0.d0
  enddo
  enddo
  endif
!End of inner boundary treatment

  call fucomm(uepar,1)
  if(ibndryx==0 .and. npx==1) call periodx(uepar, 1)
  if(ibndryy==0 .and. npy==1) call periody(uepar, 1)
  if(ibndryz==0 .and. npz==1) call periodz(uepar, 1)

#ifdef debug
if(irank==0) write(*,*) t,x(mxe/4,mye/2,mze/4,1)
#endif

  do k=mzs, mze
  do j=mys, mye
  do i=mxs, mxe
    xdot(i,j,k,1)=-d1fc4(uepar(i,j,k-2),uepar(i,j,k-1),uepar(i,j,k),uepar(i,j,k+1),uepar(i,j,k+2),dz(k))
    xdot(i,j,k,2)=tau*d1fc4(x(i,j,k-2,1),x(i,j,k-1,1),x(i,j,k,1),x(i,j,k+1,1),x(i,j,k+2,1),dz(k))&
      -d1fc4(phi(i,j,k-2),phi(i,j,k-1),phi(i,j,k),phi(i,j,k+1),phi(i,j,k+2),dz(k))
!    xdot(i,j,k,1)=-d1fc2(uepar(i,j,k-1),uepar(i,j,k),uepar(i,j,k+1),dz(k))
!    xdot(i,j,k,2)=tau*d1fc2(x(i,j,k-1,1),x(i,j,k,1),x(i,j,k+1,1),dz(k))&
!      -d1fc2(phi(i,j,k-1),phi(i,j,k),phi(i,j,k+1),dz(k))
    
  enddo
  enddo
  enddo

  if(ibndryx>0) call bndryx(x, xdot, nvar)
  if(ibndryy>0) call bndryy(x, xdot, nvar)
  if(ibndryz>0) call bndryz(x, xdot, nvar)
  
  !Copy data to time integrator
  ijk=1
  do l=1, nvar
  do k=mzs, mze
  do j=mys, mye   
  do i=mxs, mxe
    uxdot(ijk)=xdot(i,j,k,l)
    ijk=ijk+1
  enddo
  enddo
  enddo
  enddo
  
  ierror=0  !important
  return
  
end subroutine fcvfun

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

subroutine init_rk

  use input_parameters
  use topology
  use udata
  use solver
  use runge_kutta
 
  implicit none
  
  !integer(kind=8) :: ijk  
  integer :: ijk
  
  allocate(f1(neq), f2(neq),&
     f3(neq),f4(neq),f5(neq),ftmp(neq),stat=ierr)  
  if(ierr .ne. 0) then
  write(*,*) "Memory allocating error."
  stop
  endif
  
  ijk=1
  do l=1, nvar
  do k=mzs, mze
  do j=mys, mye   
  do i=mxs, mxe

    ux(ijk)=x(i,j,k,l)
    uxdot(ijk)=0.d0

    ijk=ijk+1
  
  enddo
  enddo
  enddo
  enddo
  
  call MPI_Barrier(comm, ierr)  
#ifdef debug
if(irank==0) write(*,*) "Solver initialized."
#endif


end subroutine init_rk


subroutine rk45f(fcvfun, t, h, ux, ipar, rpar, est)

  use input_parameters
  use topology
  use runge_kutta
  
  implicit none
  
  external fcvfun
  
  integer :: ierror
  integer(kind=8) :: neq, ik
  integer(kind=8) :: ipar(*)
  real(rk) :: h, t
  real(rk) :: ux(*), est(*), rpar(*)
  
  neq = ipar(1)
  
  call fcvfun(t, ux, f1, ipar, rpar, ierror)
  do ik=1,neq
    f1(ik)=h*f1(ik)
  enddo
  
  do ik=1,neq
    ftmp(ik)=ux(ik)+a21*f1(ik)
  enddo
  
  call fcvfun(t+c2*h, ftmp, f2, ipar, rpar, ierror)
  do ik=1,neq
    f2(ik)=h*f2(ik)
  enddo
  
  do ik=1,neq
    ftmp(ik)=ux(ik)+a31*f1(ik)+a32*f2(ik)
  enddo
  
  call fcvfun(t+c3*h, ftmp, f3, ipar, rpar, ierror)
  do ik=1,neq
    f3(ik)=h*f3(ik)
  enddo
  
  do ik=1,neq
    ftmp(ik)=ux(ik)+a41*f1(ik)+a42*f2(ik)+a43*f3(ik)
  enddo
  
  call fcvfun(t+c4*h, ftmp, f4, ipar, rpar, ierror)
  do ik=1,neq
    f4(ik)=h*f4(ik)
  enddo
  
  do ik=1,neq
    ftmp(ik)=ux(ik)+a51*f1(ik)+a52*f2(ik)+a53*f3(ik)+a54*f4(ik)
  enddo
  
  call fcvfun(t+c5*h, ftmp, f5, ipar, rpar, ierror)
  do ik=1,neq
    f5(ik)=h*f5(ik)
  enddo
  
  do ik=1,neq
    ux(ik)=ux(ik)+b1*f1(ik)+b2*f2(ik)+b3*f3(ik)+b4*f4(ik)+b5*f5(ik)
    est(ik)=h*((b1-d1)*f1(ik)+(b2-d2)*f2(ik)+(b3-d3)*f3(ik)&
           +(b4-d4)*f4(ik)+(b5-d5)*f5(ik))
  enddo
  
  
end subroutine rk45f

subroutine rk2(fcvfun, t, h, ux, ipar, rpar, est)

  use input_parameters
  use topology
  use runge_kutta
  
  implicit none
  
  external fcvfun
  
  integer :: ierror
  integer(kind=8) :: neq, ik
  integer(kind=8) :: ipar(*)
  real(rk) :: h, t
  real(rk) :: ux(*), rpar(*), est(*)
  
  neq=ipar(1)
  
  call fcvfun(t, ux, f1, ipar, rpar, ierror)
  do ik=1,neq
    f1(ik)=h*f1(ik)
  enddo
  
  do ik=1,neq
    ftmp(ik)=ux(ik)+2.d0/3.d0*f1(ik)
  enddo
  
  call fcvfun(t+2./3.*h, ftmp, f2, ipar, rpar, ierror)
  do ik=1,neq
    f2(ik)=h*f2(ik)
  enddo
  
  do ik=1,neq
    ux(ik)=ux(ik)+f1(ik)/4.d0+f2(ik)*3.d0/4.d0
  enddo
  
end subroutine rk2


subroutine rk4(fcvfun, t, h, ux, ipar, rpar, est)

  use input_parameters
  use topology
  use runge_kutta
  
  implicit none
  
  external fcvfun
  
  integer :: ierror
  integer(kind=8) :: neq, ik
  integer(kind=8) :: ipar(*)
  real(rk) :: h, t
  real(rk) :: ux(*), rpar(*), est(*)
  
  neq=ipar(1)
  
  call fcvfun(t, ux, f1, ipar, rpar, ierror)
  do ik=1,neq
    f1(ik)=h*f1(ik)
  enddo
  
  do ik=1,neq
    ftmp(ik)=ux(ik)+0.5d0*f1(ik)
  enddo
  
  call fcvfun(t+0.5d0*h, ftmp, f2, ipar, rpar, ierror)
  do ik=1,neq
    f2(ik)=h*f2(ik)
  enddo
  
  do ik=1,neq
    ftmp(ik)=ux(ik)+0.5d0*f2(ik)
  enddo
  
  call fcvfun(t+0.5d0*h, ftmp, f3, ipar, rpar, ierror)
  do ik=1,neq
    f3(ik)=h*f3(ik)
  enddo
  
  do ik=1,neq
    ftmp(ik)=ux(ik)+f3(ik)
  enddo
  
  call fcvfun(t+h, ftmp, f4, ipar, rpar, ierror)
  do ik=1,neq
    f4(ik)=h*f4(ik)
  enddo
  
  do ik=1,neq
    ux(ik)=ux(ik)+(f1(ik)+2.d0*f2(ik)+2.d0*f3(ik)+f4(ik))/6.d0
  enddo
  
end subroutine rk4

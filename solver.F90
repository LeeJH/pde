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

subroutine init_solver
  
  use input_parameters
  use topology
  use solver
  use udata, only:nvar
  
  implicit none

  neq = mxl*myl*mzl*nvar
  
  allocate(ux(neq),uxdot(neq),est(neq),stat=ierr)
  if(ierr .ne. 0) then
  write(*,*) "Memory allocating error."
  stop
  endif
  
!  allocate(ipar(2),stat=ierr)
!  if(ierr .ne. 0) then
!  write(*,*) "Memory allocating error."
!  stop
!  endif
  
  allocate(rpar(2),stat=ierr)
  if(ierr .ne. 0) then
  write(*,*) "Memory allocating error."
  stop
  endif
  
  
  dt=hstart 
  tnow=tstart
  
  ipar(1)=neq
  ipar(2)=irank

#ifdef _use_cvode
  call init_cvode
#else  
  call init_rk
#endif


#ifdef debug
if(irank.eq.0) write(*,*) "Solver initialized."
#endif
  
  return
  
end subroutine init_solver


subroutine solve(tout)

  use input_parameters
  use topology
  use solver
  
  implicit none

#ifndef _use_cvode  
  external fcvfun
#endif

  integer :: it
  real(rk) :: tout

#ifdef _use_cvode
  
  call fcvode(tout,tnow,ux,itask,ierr)

#ifdef debug
if(irank==0) write(*,*) "Time=",tnow, "Flag=",ierr
#endif
  
  if(ierr .ne. 0) then
  write(*,*) "CVode error=",ierr
  call MPI_ABORT(comm, 1, ierr)
  endif
  
#else  
  do while(tnow<tout)

#ifdef _use_rk2
    call rk2(fcvfun, tnow, dt, ux, ipar, rpar, est)
#else 
    call rk4(fcvfun, tnow, dt, ux, ipar, rpar, est)
#endif

    tnow=tnow+dt
    call nextdt
  enddo
#endif

  tout=tnow
  
  return
  
end subroutine solve


subroutine nextdt

  use input_parameters
  use topology
  use solver
  
  dt=hstart
  
  call MPI_Barrier(comm, ierr)
  
  return

end subroutine nextdt
    


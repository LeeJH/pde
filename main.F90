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

program main

  use input_parameters
  use topology
  use udata
  use solver
  use fft_filter

  implicit none
  
  integer :: iout
  real(rk) :: tout
  character(len=128) :: file_name
  character(len=12) :: varname
  character(len=4) :: cn
  character(len=5) :: dirname
  
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
  
  comm=PETSC_COMM_WORLD
  call MPI_COMM_RANK(comm, irank, ierr)
  call MPI_COMM_SIZE(comm, isize, ierr)
  
  call input
  call settopology
  call init
  
  call MPI_Barrier(comm, ierr)
   
  do iout=1, nout

  dirname=cn(iout-1)//char(0)
  if(irank==0) then
  call createdir(dirname)
  endif
  call MPI_Barrier(comm, ierr)

  file_name=cn(iout-1)//"/data"//cn(irank)//".nc"
  varname="ux"
  call dump1d(file_name,ux,varname,mxl*myl*mzl*nvar)

  call MPI_Barrier(comm, ierr)
    
    tout=tstart+iout*(tend-tstart)/nout
    if(irank.eq.0) write(*,*) "Step ",iout,", Next time=",tout
    !call UT(f,t,tnow,ux,uxdot,uxmax,work,uflag)
    call solve(tout)

!  dirname=cn(iout)//char(0)
!  if(irank==0) then
!  call createdir(dirname)
!  endif
!  call MPI_Barrier(comm, ierr)

!  file_name=cn(iout)//"/data"//cn(irank)//".nc"
!  varname="ux"
!  call dump1d(file_name,ux,varname,mxl*myl*mzl*nvar)
!   
  enddo
  call MPI_Barrier(comm, ierr)
  


  call final
  call PetscFinalize(ierr)
  
end

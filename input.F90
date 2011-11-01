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

!Reads parameters form the file "input.in"
!The Rank 0 process reads the parameters and then 
!broadcast them to other process

subroutine input

  use input_parameters
  use topology
  
  implicit none
  
  logical file_exist
  
  namelist /domainsize/ mx, my, mz, mgx, mgy, mgz,&
  			xmin, xmax, ymin, ymax, zmin, zmax
  namelist /bndry/ ibndryx, ibndryy, ibndryz
  namelist /procs/ npx, npy, npz
  namelist /solver/ method,nout,tstart, tend, tol, hstart,&
                    meth, itmeth, iatol, ipre, igs, itask,&
                    maxl, mudq, mldq, mu, ml, delt, dqrely,&
                    rtol, atol, ptol
                      
  if(irank==0) then
  
  inquire(file='input.in', exist=file_exist)
  if(file_exist) then
    write(*,*) "Reading input parameters from input.in ..."
    open(21, file='input.in', status='old')
    read(21, nml=procs)
    read(21, nml=domainsize)
    read(21, nml=bndry)
    read(21, nml=solver)
    close(21)
  else
    write(*,*) " input.in does not exist!"
    stop
  endif
  
  endif
  
  call MPI_Barrier(comm, ierr)
  call broadcast_input_parameters
  
end subroutine input

!Broadcast parameters from process 0 to other processes
subroutine broadcast_input_parameters

  use input_parameters
  use topology
  
  implicit none
    
  integer, parameter :: n_integer=25, n_real=15
  integer :: intpack(n_integer)
  real(rk) :: realpack(n_real)
  
  if(irank==0) then
  intpack(1)=mx
  intpack(2)=my
  intpack(3)=mz
  intpack(4)=mgx
  intpack(5)=mgy
  intpack(6)=mgz
  intpack(7)=npx
  intpack(8)=npy
  intpack(9)=npz
  intpack(10)=ibndryx
  intpack(11)=ibndryy
  intpack(12)=ibndryz
  
  !solver
  intpack(13)=method
  intpack(14)=nout
  
  intpack(15)=meth
  intpack(16)=itmeth
  intpack(17)=iatol
  intpack(18)=ipre
  intpack(19)=igs
  intpack(20)=maxl
  intpack(21)=mudq
  intpack(22)=mldq
  intpack(23)=mu
  intpack(24)=ml
  intpack(25)=itask
  
  realpack(1)=xmin
  realpack(2)=xmax
  realpack(3)=ymin
  realpack(4)=ymax
  realpack(5)=zmin
  realpack(6)=zmax
  
  !solver
  realpack(7)=tstart
  realpack(8)=tend
  realpack(9)=tol
  realpack(10)=hstart
  
  realpack(11)=rtol
  realpack(12)=atol
  realpack(13)=delt
  realpack(14)=dqrely
  realpack(15)=ptol
  
  endif
  
  call MPI_BCAST(intpack,n_integer,MPI_INTEGER,0,comm,ierr)
  call MPI_BCAST(realpack,n_real,mpiR,0,comm,ierr)
  
  if(irank>0) then
  mx=intpack(1)
  my=intpack(2)
  mz=intpack(3)
  mgx=intpack(4)
  mgy=intpack(5)
  mgz=intpack(6)
  npx=intpack(7)
  npy=intpack(8)
  npz=intpack(9)
  ibndryx=intpack(10)
  ibndryy=intpack(11)
  ibndryz=intpack(12)
  
  !solver
  method=intpack(13)
  nout=intpack(14)
  
  meth=intpack(15)
  itmeth=intpack(16)
  iatol=intpack(17)
  ipre=intpack(18)
  igs=intpack(19)
  maxl=intpack(20)
  mudq=intpack(21)
  mldq=intpack(22)
  mu=intpack(23)
  ml=intpack(24)
  itask=intpack(25)
  
  xmin=realpack(1)
  xmax=realpack(2)
  ymin=realpack(3)
  ymax=realpack(4)
  zmin=realpack(5)
  zmax=realpack(6)
  
  !solver
  tstart=realpack(7)
  tend=realpack(8)
  tol=realpack(9)
  hstart=realpack(10)
  
  rtol=realpack(11)
  atol=realpack(12)
  delt=realpack(13)
  dqrely=realpack(14)
  ptol=realpack(15)
  
  endif
  
  call MPI_Barrier(comm, ierr)
  
end subroutine broadcast_input_parameters


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

subroutine init_cvode

  use input_parameters
  use topology
  use udata
  use solver
  use cvodef
  
  integer :: ijk
  
  nlocal=neq
  nglobal=mx*my*mz*nvar
  
  CALL FNVINITP(comm, key, NLOCAL, NGLOBAL, IERR)
  if(ierr .ne. 0) then
  write(*,*) "FNVINITP err."
  call MPI_ABORT(comm, 1, ierr)
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
  
!  allocate(iout(21),stat=ierr)
!  if(ierr .ne. 0) then
!  write(*,*) "Memory allocating error."
!  stop
!  endif
  
!  allocate(rout(6),stat=ierr)
!  if(ierr .ne. 0) then
!  write(*,*) "Memory allocating error."
!  stop
!  endif
  
  allocate(atolv(nlocal),stat=ierr)
  if(ierr .ne. 0) then
  write(*,*) "Memory allocating error."
  stop
  endif
  
  CALL FCVMALLOC(tstart, ux, METH, ITMETH, IATOL, RTOL, ATOL,&
                 IOUT, ROUT, IPAR, RPAR, IERR)
  if(ierr .ne. 0) then
  write(*,*) "FCVMALLOC err."
  call MPI_ABORT(comm, 1, ierr)
  endif

  CALL FCVSETIIN('MAX_NSTEPS', 10000, IERR)
  if(ierr .ne. 0) then
  write(*,*) "FCVSETIIN err."
  call MPI_ABORT(comm, 1, ierr)
  endif
  
  CALL FCVSPGMR(IPRE, IGS, MAXL, delt, IERR)
  if(ierr .ne. 0) then
  write(*,*) "FCVSPGMR err."
  call MPI_ABORT(comm, 1, ierr)
  endif
  
  CALL FCVBBDINIT(NLOCAL, MUDQ, MLDQ, MU, ML, dqrely, IERR)
  if(ierr .ne. 0) then
  write(*,*) "FCVBBDINIT err."
  call MPI_ABORT(comm, 1, ierr)
  endif
  
!  CALL FCVSPILSSETJAC(0, IERR)
  if(ierr .ne. 0) then
  write(*,*) "FCVSPILSSETJAC err."
  call MPI_ABORT(comm, 1, ierr)
  endif
  
  if(irank.eq.0) write(*,*) "CVODE initialized."
  
  return

end subroutine init_cvode

SUBROUTINE FCVGLOCFN(NLOC, T, YLOC, GLOC, IPAR, RPAR, IERROR)
!     Routine to define local approximate function g, here the same as f. 
  
  use digit
  IMPLICIT NONE
!
  INTEGER(kind=8) :: NLOC, IPAR(*)
  integer :: IERror
  real(rk) :: T, YLOC(*), GLOC(*), RPAR(*)
!     
  CALL FCVFUN(T, YLOC, GLOC, IPAR, RPAR, IERRor)
!
  RETURN
END
!
!     ------------------------------------------------------------------------
!      
SUBROUTINE FCVCOMMFN(NLOC, T, YLOC, IPAR, RPAR, IERROR)
!  Routine to perform communication required for evaluation of g.
  
  use digit
  implicit none
  
  integer :: ierror
  integer(kind=8) :: NLOC, IPAR(*)
  real(rk) :: t, yloc(*), rpar(*)
  
  IERROR = 0
  
  RETURN
END

SUBROUTINE FCVJTIMES (V, FJV, T, Y, FY, H, IPAR, RPAR, WORK, IERROR)

  use digit
  implicit none
  
  integer :: ierror
  integer(kind=8) :: ipar(*)
  real(rk) :: t, h
  real(rk) ::  V(*), FJV(*), Y(*), FY(*), RPAR(*), WORK(*)
  
  
  ierror=1
  return
  
end subroutine fcvjtimes

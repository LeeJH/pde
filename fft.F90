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

subroutine init_fft

  use parameters
  use input_parameters
  use topology
  use fft_filter
  
  implicit none
  include 'fftw3.f'

  if(mod(mz,npz) .ne. 0) then
  write(*,*) "mod(mz,npz) should equal to zero."
  write(*,*) "Processor number error. Aborting ... "
  call MPI_Abort(comm,1,ierr)
  endif

  allocate(gathsub(mxl*myl*mzl), stat=ierr)
  if(ierr .ne. 0) then
  write(*,*) "Memory allocating error."
  stop
  endif

  if(ipz==0) then
  allocate(fftin(mz),fftout(mz),gathtol(mz*mxl*myl),stat=ierr)
  if(ierr .ne. 0) then
  write(*,*) "Memory allocating error."
  stop
  endif

  do k=1,mz
    fftin(k)=(0.0,0.0)
    fftout(k)=(0.0,0.0)
  enddo

  call dfftw_plan_dft_1d(planF,mz,fftin,fftout,FFTW_FORWARD,FFTW_MEASURE)
  call dfftw_plan_dft_1d(planR,mz,fftout,fftin,FFTW_BACKWARD,FFTW_MEASURE)
  
  endif
  
  return
  
end subroutine init_fft

subroutine fft_final

  use topology
  use fft_filter
  
  implicit none
  include 'fftw3.f'
  
  deallocate(gathsub,stat=ierr)
  
  if(ipz==0) then
  call dfftw_destroy_plan(planF)
  call dfftw_destroy_plan(planR)
  
  deallocate(fftin, fftout, gathtol,stat=ierr)
  endif
  
  return 
  
end subroutine fft_final

subroutine fftfilter

  use parameters
  use input_parameters
  use topology
  use udata
  use field
  use fft_filter
  
  implicit none
  include 'fftw3.f'

  integer*8 :: ijk, ij
  
  ijk=1
  do j=mys,mye
  do i=mxs,mxe
     do k=mzs,mze
     gathsub(ijk)=phi(i,j,k)
     ijk=ijk+1
     enddo
  enddo
  enddo
   
  call MPI_Gather(gathsub,mxl*myl*mzl,mpiR,&
         gathtol,mxl*myl*mzl,mpiR,0,commxy,ierr)
  
  if(irankxy==0) then       

  do j=1,myl
  do i=1,mxl
    do ij=0,isizexy-1
    ijk=ij*mxl*myl*mzl
    do k=1,mzl
    fftin(ij*mzl+k)=cmplx(gathtol(ijk+((i-1)+(j-1)*mxl)*mzl+k),0.d0)
    enddo
    enddo
    call dfftw_execute_dft(planF,fftin,fftout)
   
!    do k=4,mz-2
    do k=mz/4,mz-mz/4+2
    fftout(k)=(0.d0,0.d0)
    enddo
 
    call dfftw_execute_dft(planR,fftout,fftin)
    
    do ij=0,isizexy-1
    ijk=ij*mxl*myl*mzl
    do k=1,mzl
    gathtol(ijk+((i-1)+(j-1)*mxl)*mzl+k)=real(fftin(ij*mzl+k))/mz
    enddo
    enddo
    
  enddo
  enddo
  
  endif
  call MPI_Scatter(gathtol,mxl*myl*mzl,mpiR,gathsub,mxl*myl*mzl,mpiR,0, commxy,ierr)
  ijk=1
  do j=mys,mye
  do i=mxs,mxe
     do k=mzs,mze
     phi(i,j,k)=gathsub(ijk)
     ijk=ijk+1
     enddo
  enddo
  enddo
         
  return

end subroutine fftfilter
    

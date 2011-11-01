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

subroutine final

  use input_parameters
  use topology
  use grid
  use udata
  use solver 
  use field
  use runge_kutta
  use cvodef
  
  implicit none
  
  call final_poisson
  call fft_final
!  call final_fft2d

  deallocate(xx,dx,yy,dy,zz,dz,stat=ierr)
  deallocate(x,xdot,stat=ierr)
  deallocate(Ni0,Te0,uepar,B,stat=ierr)
  deallocate(rho,phi,stat=ierr)
#ifdef _use_cvode
  deallocate(atolv,stat=ierr)
#else
  deallocate(f1,f2,f3,f4,f5,ftmp,stat=ierr)
#endif
  deallocate(indglxy,indgl3d,stat=ierr)
  deallocate(ux,uxdot,est,stat=ierr)
!  deallocate(ipar,rpar,stat=ierr)
  
  return
end subroutine final

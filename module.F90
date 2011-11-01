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

module digit
#include "finclude/petscsys.h"

  integer, parameter :: doublekind=selected_real_kind(12),&
    singlekind=selected_real_kind(6),&
    defaultkind=kind(0.0)
    
  integer, parameter :: rk=doublekind,mpiR=MPI_DOUBLE_PRECISION,&
                        mpiC=MPI_DOUBLE_COMPLEX
  
  integer :: i, j, k, l, ierr
  
save
end module digit

! Module for parameters, mainly constants
module parameters
  use digit
  real(rk), parameter :: Pi = 3.1415926535898
  !Problem part
  
save
end module parameters

! Module for problem
module input_parameters
  use digit
  use parameters

  integer :: mx, my, mz
  integer :: mgx, mgy, mgz
  integer :: ibndryx, ibndryy, ibndryz
  real(rk) :: xmin, xmax, ymin, ymax, zmin, zmax
  
  !Solver part
  !Rk
  integer :: method
  integer :: nout
  real(rk) :: tstart, tend, tol, hstart
  
  !COVDE
  integer :: METH, ITMETH, IATOL, IPRE, IGS,&
             MAXL, MUDQ, MLDQ, MU, ML, itask
  real(rk) :: rtol, atol, delt, dqrely
  
  !2D Poisson solver
  real(rk) :: ptol
  
save
end module input_parameters

! Module for topology of processes
module topology
  use digit
  
  integer :: comm, commx, commz, commxy, commyz
  integer :: irank, isize, irankz, isizez, irankxy, isizexy,isizeyz,irankyz
  integer :: irankx, isizex
  integer :: npx, npy, npz
  integer :: mxs, mxe, mxl, mys, mye, myl, mzs, mze, mzl
  integer :: mxsub, mysub, mzsub
  integer :: ipx, ipy, ipz
  integer :: ixup, ixdown, iyup, iydown, izup, izdown
  integer :: ixglobal, iyglobal, izglobal
  integer*4 ,dimension(:), allocatable :: indglxy !Canot use integer*8
  integer*4 ,dimension(:), allocatable :: indgl3d
  
  
  
save
end module topology

! Module for grid
module grid
  use digit

  real(rk),dimension(:), allocatable :: xx, dx, yy, dy, zz, dz
  real(rk),parameter :: beta=0.01, tau=0.1
save
end module grid

! Module for main variables
module udata
  use digit

  integer, parameter :: nvar=2
  real(rk), dimension(:,:,:,:), allocatable :: x, xdot
  
  !Problem part
  real(rk), dimension(:,:), allocatable :: Ni0,Te0
  real(rk), dimension(:,:,:), allocatable :: uepar
  real(rk), dimension(:,:,:,:), allocatable :: B
  
save
end module udata

! Module for fields
module field
  use digit

  real(rk), dimension(:), allocatable :: rho
  real(rk), dimension(:,:,:), allocatable :: phi
  
save
end module field

! Module for solver
module solver
  use digit
  
  !rk solver
!*****************************************************
!  integer :: OUTCH
!  real(rk) :: MCHPES,DWARF
!  
!  integer :: neq, lenwrk
!  logical :: errass, mesage
!  character(len=2) :: task ! CT or UT
!  
!  real(rk), allocatable, dimension(:) :: work, thres, ystart
!  
!  integer :: uflag
!  real(rk) :: dt
!  real(rk), dimension(:), allocatable :: uxmax
!  real(rk), allocatable, dimension(:) :: ux, uxdot
!*****************************************************  
  real(rk) :: dt, tnow
  integer(kind=8) :: neq
  integer(kind=8) :: ipar(2)
  real(rk), allocatable, dimension(:) :: ux, uxdot, rpar, est
  
save
end module solver

! Module for runge-kutta
module runge_kutta
  use digit
  
  !rk solver
  
  real(rk), parameter :: c2=0.25d0, c3=3.d0/8.d0, c4=12.d0/13.d0,&
           c5=1.d0, c6=0.5d0, b1=16.d0/135.d0, b2=0.d0,&
           b3=6656.d0/12825.d0, b4=28561.d0/56430.d0, b5=-9.d0/50.d0,&
           b6=2.d0/55.d0, d1=25.d0/216.d0, d2=0.d0, d3=1408.d0/2565.d0,&
           d4=2197.d0/4104.d0, d5=-0.2d0, d6=0.d0
  real(rk), parameter :: a21=0.25d0, a31=3.d0/32.d0, a32=9.d0/32.d0,&
           a41=1932.d0/2197.d0, a42=-7200.d0/2197.d0, a43=7296.d0/2197.d0,&
           a51=439.d0/216.d0, a52=8.d0, a53=3680.d0/513.d0, &
           a54=-845.d0/4104.d0, a61=-8.d0/27.d0, a62=2.d0,&
           a63=-3544.d0/2565.d0, a64=1859.d0/4104.d0, a65=-11.d0/40.d0
  
!  real(rk), allocatable, dimension(:) :: fres
  real(rk), allocatable, dimension(:) :: f1, f2, f3, f4, f5, ftmp
  
  
save
end module runge_kutta

module cvodef
  use digit
  
  integer, parameter :: key=1
  integer(kind=8) :: nlocal, nglobal

  integer :: iout(21)
  real(rk) :: rout(6)
  real(rk), dimension(:), allocatable :: atolv
  
save
end module cvodef

module fft_filter
  use digit
  
  integer(kind=8) :: planF, planR
  complex(rk), dimension(:), allocatable :: fftin, fftout
  real(rk), dimension(:), allocatable :: gathsub, gathtol
  
save
end module fft_filter
 

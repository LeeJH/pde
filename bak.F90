
subroutine init_rk

  use input_parameters
  use topology
  use udata
  use solver
 
  implicit none
  
  external envirn, setup
  
  integer :: isol
  
  call envirn(OUTCH,MCHPES,DWARF)
  
  neq = mxl*myl*mzl*nvar

  task  ="UT"
  errass=.false.
  mesage=.true.

  lenwrk=32*neq
  
  dt=(tend-tstart)/nout
  
  allocate(ystart(neq), thres(neq), work(lenwrk))
  
  allocate(ux(neq),uxdot(neq),uxmax(neq))
  
  isol=1
  do i=mxs, mxe
  do j=mys, mye
  do k=mzs, mze
  do l=1, nvar
    
    ystart(isol)=x(i,j,k,l)
    
    ux(isol)=x(i,j,k,l)
    uxdot(isol)=0.d0
    
    thres(isol) =1.0e-8
    isol=isol+1
  
  enddo
  enddo
  enddo
  enddo
  
  call setup(neq,tstart,ystart,tend,tol,thres,method,&
             task,errass,hstart,work,lenwrk,mesage)
  
#ifdef debug
write(*,*) "Solver initialized."
#endif


end subroutine init_rk


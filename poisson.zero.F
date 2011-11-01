
      subroutine init_poisson(indgl,ilen,m,n,dx,dy,tol,comm) 
      
      implicit none
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscpc.h"
#include "finclude/petscksp.h"
#include "finclude/petscis.h"

      integer i,j,ione,ierr
      integer m,n,comm
      integer*4 indgl(*),ilen
      
      double precision dx, dy, tol
      integer*4 Istart, Iend, ilow, ihigh
      double precision dxs, dys, dydxs, dxdys
      Vec x,b,xlocal
      IS sfrom, sto
      VecScatter scatter
      Mat A
      PC pc
      KSP ksp
      
      integer*4 II,JJ
      double precision vm4,vm3,vm2,vm1,v,vp1,vp2,vp3,vp4,vb,
     1     one,neg_one,tmp
      
      
      common /coeff/ dxs, dys, dxdys, dydxs
      common /range/ Istart,Iend,ilow,ihigh
      common /pvec/ x,b,xlocal,sfrom,sto,scatter
      common /pmat/ A
      common /ppc/ pc
      common /pksp/ ksp
      
      !problem related part
      dxs=dx*dx
      dys=dy*dy
      
      dxdys=dxs/dys
      dydxs=dys/dxs
      
      ione=1
      one=1.0d0
      neg_one=-1.0
      vm4=-1.d0/12.d0
      vm3=16.d0/12.d0
      vm2=-1.d0/12.d0*dydxs
      vm1=16.d0/12.d0*dydxs
      v=-(1.d0+dydxs)*30.d0/12.d0
      vp1=16.d0/12.d0*dydxs
      vp2=-1.d0/12.d0*dydxs
      vp3=16.d0/12.d0
      vp4=-1.d0/12.d0
      
      call MatCreate(comm,A,ierr)
      call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m*n,m*n,ierr)
      call MatSetType(A, MATAIJ,ierr)
      call MatSetFromOptions(A,ierr)
      call MatMPIAIJSetPreallocation(A,9,PETSC_NULL_INTEGER,9,  
     &                     PETSC_NULL_INTEGER,ierr)
      call MatSeqAIJSetPreallocation(A,9,PETSC_NULL_INTEGER,ierr)
      
      call MatGetOwnershipRange(A,Istart,Iend,ierr)
      
      do II=Istart,Iend-1
      
        j=II/m
        i=II-j*m
        
        if ((j.gt.1).and.(j.lt.n-2).and.(i.gt.0) .and. (i.lt.m-1)) then
          JJ = II - 2*m
          call MatSetValues(A,ione,II,ione,JJ,vm4,INSERT_VALUES,ierr)
          JJ = II - m
          call MatSetValues(A,ione,II,ione,JJ,vm3,INSERT_VALUES,ierr)
          JJ = II + m
          call MatSetValues(A,ione,II,ione,JJ,vp3,INSERT_VALUES,ierr)
          JJ = II + 2*m
          call MatSetValues(A,ione,II,ione,JJ,vp4,INSERT_VALUES,ierr)
        endif
        
        if ((i.gt.1) .and. (i.lt.m-2)) then
          JJ = II - 2
          call MatSetValues(A,ione,II,ione,JJ,vm2,INSERT_VALUES,ierr)
          JJ = II - 1
          call MatSetValues(A,ione,II,ione,JJ,vm1,INSERT_VALUES,ierr)  
          JJ = II + 1
          call MatSetValues(A,ione,II,ione,JJ,vp1,INSERT_VALUES,ierr)
          JJ = II + 2
          call MatSetValues(A,ione,II,ione,JJ,vp2,INSERT_VALUES,ierr)
        endif
        
        if ((j.eq.0).and.(i.gt.0) .and. (i.lt.m-1)) then
          JJ = II + m
          call MatSetValues(A,ione,II,ione,JJ,vp3,INSERT_VALUES,ierr)
          JJ = II + 2*m
          call MatSetValues(A,ione,II,ione,JJ,vp4,INSERT_VALUES,ierr)
          JJ = II + m*(n-1)
          call MatSetValues(A,ione,II,ione,JJ,vm3,INSERT_VALUES,ierr)
          JJ = II + m*(n-2)
          call MatSetValues(A,ione,II,ione,JJ,vm4,INSERT_VALUES,ierr)
        endif
        
        if ((j.eq.1).and.(i.gt.0) .and. (i.lt.m-1)) then
          JJ = II + m
          call MatSetValues(A,ione,II,ione,JJ,vp3,INSERT_VALUES,ierr)
          JJ = II + 2*m
          call MatSetValues(A,ione,II,ione,JJ,vp4,INSERT_VALUES,ierr)
          JJ = II - m
          call MatSetValues(A,ione,II,ione,JJ,vm3,INSERT_VALUES,ierr)
          JJ = II + m*(n-2)
          call MatSetValues(A,ione,II,ione,JJ,vm4,INSERT_VALUES,ierr)
        endif
        
        if ((j.eq.(n-2)).and.(i.gt.0) .and. (i.lt.m-1)) then
          JJ = II - 2*m
          call MatSetValues(A,ione,II,ione,JJ,vm4,INSERT_VALUES,ierr)
          JJ = II - m
          call MatSetValues(A,ione,II,ione,JJ,vm3,INSERT_VALUES,ierr)
          JJ = II + m
          call MatSetvalues(A,ione,II,ione,JJ,vp3,INSERT_VALUES,ierr)
          JJ = II - m*(n-2)
          call MatSetValues(A,ione,II,ione,JJ,vp4,INSERT_VALUES,ierr)
        endif
        
        if ((j.eq.(n-1)).and.(i.gt.0) .and. (i.lt.m-1)) then
          JJ = II - 2*m
          call MatSetValues(A,ione,II,ione,JJ,vm4,INSERT_VALUES,ierr)
          JJ = II - m
          call MatSetValues(A,ione,II,ione,JJ,vm3,INSERT_VALUES,ierr)
          JJ = II - m*(n-1)
          call MatSetvalues(A,ione,II,ione,JJ,vp3,INSERT_VALUES,ierr)
          JJ = II - m*(n-2)
          call MatSetValues(A,ione,II,ione,JJ,vp4,INSERT_VALUES,ierr)
        endif
        
        call  MatSetValues(A,ione,II,ione,II,v,INSERT_VALUES,ierr)
        
        if (i.eq.1) then
          JJ = II - 1
          vb = 11.d0*dydxs/12.d0
          call MatSetValues(A,ione,II,ione,JJ,vb,INSERT_VALUES,ierr)
          JJ = II
          vb = -30.d0/12.d0-5.d0/3.d0*dydxs
          call MatSetValues(A,ione,II,ione,JJ,vb,INSERT_VALUES,ierr)
          JJ = II + 1
          vb = 0.5d0*dydxs
          call MatSetValues(A,ione,II,ione,JJ,vb,INSERT_VALUES,ierr)
          JJ = II + 2
          vb = 1.d0/3.d0*dydxs
          call MatSetValues(A,ione,II,ione,JJ,vb,INSERT_VALUES,ierr)
          JJ = II + 3
          vb = -1.d0/12.d0*dydxs
          call MatSetValues(A,ione,II,ione,JJ,vb,INSERT_VALUES,ierr)
        endif
        
        if (i.eq.(m-2)) then
          JJ = II + 1
          vb = 11.d0*dydxs/12.d0
          call MatSetValues(A,ione,II,ione,JJ,vb,INSERT_VALUES,ierr)
          JJ = II
          vb = -30.d0/12.d0-5.d0/3.d0*dydxs
          call MatSetValues(A,ione,II,ione,JJ,vb,INSERT_VALUES,ierr)
          JJ = II - 1
          vb = 0.5d0*dydxs
          call MatSetValues(A,ione,II,ione,JJ,vb,INSERT_VALUES,ierr)
          JJ = II - 2
          vb = 1.d0/3.d0*dydxs
          call MatSetValues(A,ione,II,ione,JJ,vb,INSERT_VALUES,ierr)
          JJ = II - 3
          vb = -1.d0/12.d0*dydxs
          call MatSetValues(A,ione,II,ione,JJ,vb,INSERT_VALUES,ierr)
        endif
  
      enddo
      
      call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
      
      call VecCreateMPI(comm,PETSC_DECIDE,m*n,x,ierr)
      call VecDuplicate(x,b,ierr)
      call VecGetOwnershipRange(x,ilow,ihigh,ierr)
      call VecCreateSeq(PETSC_COMM_SELF,ilen,xlocal,ierr)

#ifdef debug
      if(Istart.eq.0) then
      write(*,*) "Vectors and Matrix for poisson ready."
      write(*,*) "Poisson Vector layout:"
      endif
      call MPI_Barrier(comm,ierr)
      write(*,*) Istart-ilow,Iend-ihigh
      call MPI_Barrier(comm,ierr)
#endif
      
      call ISCreateStride(comm,ilen,0,1,sto,ierr)
      call ISCreateGeneral(comm,ilen,indgl,sfrom,ierr)
      call VecScatterCreate(x,sfrom,xlocal,sto,scatter,ierr)
      
      call KSPCreate(comm,ksp,ierr)
      
      call KSPSetOperators(ksp,A,A,SAME_NONZERO_PATTERN,ierr)
      
      call KSPSetTolerances(ksp,tol,PETSC_DEFAULT_DOUBLE_PRECISION,    
     &     PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_INTEGER,ierr)
      
      call KSPGetPC(ksp,pc,ierr)
      call PCSetType(pc,PCJACOBI,ierr)
      call KSPSetFromOptions(ksp,ierr)

#ifdef debug
      if(Istart.eq.0) write(*,*) "Poisson initialized."
#endif
        
      return

      end subroutine init_poisson
      
      subroutine solve_poisson(rho,indgl,ilen)
      
      implicit none
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscpc.h"
#include "finclude/petscksp.h"
#include "finclude/petscis.h"

      
      integer m,n
      integer*4 IJ, mn
      integer*4 ilen, indgl(ilen)
      
      double precision rho(ilen)
      double precision dxs, dys, dydxs, dxdys
      
      integer*4 Istart, Iend, ilow, ihigh
      PetscScalar xtmp(1)
      PetscOffset ioff
      
      Vec x,b,xlocal
      IS sfrom, sto
      VecScatter scatter
      Mat A
      PC pc
      KSP ksp
      
      integer i,j,ierr
      
      common /coeff/ dxs, dys, dxdys, dydxs
      common /range/ Istart,Iend,ilow,ihigh
      common /pvec/ x,b,xlocal,sfrom,sto,scatter
      common /pmat/ A
      common /ppc/ pc
      common /pksp/ ksp
      
      do i=1,ilen
      rho(i)=rho(i)*dys
      enddo
      call VecSetValues(b,ilen,indgl,rho,INSERT_VALUES,ierr)
      call VecAssemblyBegin(b, ierr)
      call VecAssemblyEnd(b,ierr)
      
      call KSPSolve(ksp,b,x,ierr)
      
      call VecScatterBegin(scatter,x,xlocal,INSERT_VALUES
     &     ,SCATTER_FORWARD,ierr)
      call VecScatterEnd(scatter,x,xlocal,INSERT_VALUES
     &     ,SCATTER_FORWARD,ierr)
      
      call VecGetArray(xlocal,xtmp,ioff,ierr)
      
      do i=1,ilen
      rho(i)=xtmp(ioff+i)
      enddo
      
      call VecRestoreArray(xlocal,xtmp,ioff,ierr)
      
      return
      
      end subroutine solve_poisson
      
      
      subroutine final_poisson
      implicit none
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscpc.h"
#include "finclude/petscksp.h"
#include "finclude/petscis.h"
      
      Vec x,b,xlocal
      IS sfrom, sto
      VecScatter scatter
      Mat A
      PC pc
      KSP ksp
      
      integer ierr
      
      common /pvec/ x,b,xlocal,sfrom,sto,scatter
      common /pmat/ A
      common /ppc/ pc
      common /pksp/ ksp
      
      call VecDestroy(x,ierr)
      call VecDestroy(b,ierr)
      call VecDestroy(xlocal,ierr)
      call MatDestroy(a,ierr)
      call ISDestroy(sfrom,ierr)
      call ISDestroy(sto,ierr)
      call KSPDestroy(ksp,ierr)
      
      return
      
      end subroutine final_poisson

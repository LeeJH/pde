/********************************************
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
*********************************************/

#include <stdio.h>
#include "petscksp.h"
#include "fftw3.h"
#include "mpi.h"
#include "fftw3-mpi.h"

#define CHKERRQ(ierr)
#define PI 3.1415926535898
typedef struct
{
  //Pay attention to the order of the variables
  int commxc,commyzc;
  int npxc, npyc, npzc;
  int mxc, myc, mzc;
  int mxlc, mylc, mzlc;
  int ipxc, ipyc, ipzc;
  int ixglobalc, iyglobalc, izglobalc;
  double dx0, dy0, dz0;
} topology;

  extern topology topo_;
  MPI_Comm commx, commyz;
  int irankx, isizex, irankyz, isizeyz;
  
  const int rnk=2;
  const ptrdiff_t myz[]={64,128};
  ptrdiff_t howmany;
  
  const ptrdiff_t mx=128, my=64, mz=128;
  fftw_plan mplanF, mplanR;
  fftw_complex *minp, *mout;

  ptrdiff_t ly, lys, alloc_ly;
  
  double ***freq, ***data;

//**** Solver part ********

  Vec xr, xi, br, bi;
  Mat A;
  KSP ksp;
  PC  pc;
  PetscErrorCode ierr;
  PetscInt Ii,J,Istart,Iend,its;
  PetscTruth flg = PETSC_FALSE;
  double dxs,dys,dzs;

void init_fft2d_(void );
void final_fft2d_(void);
void fftpoisson_(double *rho, int *mxl, int *myl, int *mzl);

// Memory functions
double *rvector(int size);
int *ivector(int size);
void free_rvector(double *rv);
void free_ivector(int *iv);
double **rmatrix(int xsize, int ysize);
int **imatrix(int xsize, int ysize);
void free_rmatrix(double **m);
void free_imatrix(int **m);
double ***r3tensor(int nrow, int ncol, int ndep);
void free_r3tensor(double ***t);


void init_fft2d_(void )
{
  int i,j,k;
  double vm2,vm1,v,vp1,vp2,vb;
  
  commx = MPI_Comm_f2c(topo_.commxc);
  MPI_Comm_rank(commx, &irankx);
  MPI_Comm_size(commx, &isizex);

  commyz = MPI_Comm_f2c(topo_.commyzc);
  MPI_Comm_rank(commyz, &irankyz);
  MPI_Comm_size(commyz, &isizeyz);

  fftw_mpi_init();

  howmany = topo_.mxlc;
//***********  
  alloc_ly = fftw_mpi_local_size_2d(my, mz, commx, &ly, &lys);
/*  alloc_ly=fftw_mpi_local_size_many(rnk, myz, howmany,
           FFTW_MPI_DEFAULT_BLOCK, commx, &ly, &lys);
*/
//***********
  
  if(((ly-topo_.mylc)!=0) || topo_.npzc>1) {
  printf("Error,npz should equal to 1, or %d\t%d\n",irankx,ly-topo_.mylc);
  MPI_Abort(commx,1);
  }

  minp = fftw_alloc_complex(alloc_ly);
  mout = fftw_alloc_complex(alloc_ly);
  
  if( !(freq = r3tensor(topo_.mxlc, topo_.mylc*topo_.mzlc, 2)) ) 
  	printf("Malloc error!\n");  
  if( !(data = r3tensor(topo_.mxlc, topo_.mylc*topo_.mzlc, 2)) ) 
  	printf("Malloc error!\n");
/*  if( !(dar = r3tensor(topo_.mxlc, topo_.mylc, topo_.mzlc)) ) 
  	printf("Malloc error!\n");
  if( !(dai = r3tensor(topo_.mxlc, topo_.mylc, topo_.mzlc)) ) 
  	printf("Malloc error!\n");
*/
//***********
  mplanF = fftw_mpi_plan_dft_2d(my, mz, minp, mout, commx, FFTW_FORWARD, FFTW_MEASURE);
  mplanR = fftw_mpi_plan_dft_2d(my, mz, minp, mout, commx, FFTW_BACKWARD, FFTW_MEASURE);
/*  mplanF = fftw_mpi_plan_many_dft(rnk, myz, howmany,
        FFTW_MPI_DEFAULT_BLOCK,FFTW_MPI_DEFAULT_BLOCK,
        minp, mout, commx, FFTW_FORWARD, FFTW_MEASURE);
  mplanR = fftw_mpi_plan_many_dft(rnk, myz, howmany,
        FFTW_MPI_DEFAULT_BLOCK,FFTW_MPI_DEFAULT_BLOCK,
        minp, mout, commx, FFTW_BACKWARD, FFTW_MEASURE);
*/
//***********
//***** Solver part ******

  dxs = topo_.dx0*topo_.dx0;
  dys = topo_.dy0*topo_.dy0;
  dzs = topo_.dz0*topo_.dz0;

  vm2=-1.0/12.0;
  vm1=16.0/12.0;
  v  =-30.0/12.0;
  vp1=16.0/12.0;
  vp2=-1.0/12.0;

  MatCreateMPIAIJ(commyz, PETSC_DECIDE, PETSC_DECIDE, mx, mx, 5, PETSC_NULL, 5, PETSC_NULL, &A);

  ierr = MatGetOwnershipRange(A,&Istart,&Iend);

  
  for (Ii=Istart; Ii<Iend; Ii++) {
    i = Ii; j = Ii;
    if ((i>1)&&(i<mx-2))   {

	J = Ii - 2; MatSetValues(A,1,&Ii,1,&J,&vm2,INSERT_VALUES);
	J = Ii - 1; ierr = MatSetValues(A,1,&Ii,1,&J,&vm1,INSERT_VALUES);CHKERRQ(ierr);
	J = Ii; ierr = MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);CHKERRQ(ierr);
	J = Ii + 1; ierr = MatSetValues(A,1,&Ii,1,&J,&vp1,INSERT_VALUES);CHKERRQ(ierr);
	J = Ii + 2; ierr = MatSetValues(A,1,&Ii,1,&J,&vp2,INSERT_VALUES);CHKERRQ(ierr);
    }

    if (i==0) {
	J = Ii; ierr = MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);CHKERRQ(ierr);}
    if (i==1)   {
	J = Ii - 1; vb = 11.0/12.0; ierr = MatSetValues(A,1,&Ii,1,&J,&vb,INSERT_VALUES);CHKERRQ(ierr);
	J = Ii ; vb = -5.0/3.0; ierr = MatSetValues(A,1,&Ii,1,&J,&vb,INSERT_VALUES);
	J = Ii + 1; vb = 0.5; ierr = MatSetValues(A,1,&Ii,1,&J,&vb,INSERT_VALUES);
	J = Ii + 2; vb = 1.0/3.0; ierr = MatSetValues(A,1,&Ii,1,&J,&vb,INSERT_VALUES);
	J = Ii + 3; vb = -1.0/12.0; ierr = MatSetValues(A,1,&Ii,1,&J,&vb,INSERT_VALUES);
}

    if (i==mx-2) {
        J = Ii + 1; vb = 11.0/12.0; ierr = MatSetValues(A,1,&Ii,1,&J,&vb,INSERT_VALUES);CHKERRQ(ierr);
        J = Ii ; vb = -5.0/3.0; ierr = MatSetValues(A,1,&Ii,1,&J,&vb,INSERT_VALUES);
        J = Ii - 1; vb = 0.5; ierr = MatSetValues(A,1,&Ii,1,&J,&vb,INSERT_VALUES);
        J = Ii - 2; vb = 1.0/3.0; ierr = MatSetValues(A,1,&Ii,1,&J,&vb,INSERT_VALUES);
        J = Ii - 3; vb = -1.0/12.0; ierr = MatSetValues(A,1,&Ii,1,&J,&vb,INSERT_VALUES);
}
    if (i==mx-1) {J = Ii; ierr = MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);CHKERRQ(ierr);} 

    }

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  ierr = VecCreate(commyz,&br);CHKERRQ(ierr);
  ierr = VecSetSizes(br,PETSC_DECIDE,mx);CHKERRQ(ierr);
  ierr = VecSetFromOptions(br);CHKERRQ(ierr);
  ierr = VecDuplicate(br,&xr);CHKERRQ(ierr);
  ierr = VecDuplicate(br,&bi);CHKERRQ(ierr);
  ierr = VecDuplicate(br,&xi);CHKERRQ(ierr);


  ierr = KSPCreate(commyz,&ksp);CHKERRQ(ierr);

  ierr = KSPSetOperators(ksp,A,A,SAME_PRECONDITIONER);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);
  PCSetType(pc,PCJACOBI);
  ierr = KSPSetTolerances(ksp,1.e-7,1.e-50,PETSC_DEFAULT,
                          PETSC_DEFAULT);CHKERRQ(ierr);

//******* End  *******
}

void final_fft2d_(void)
{

  fftw_destroy_plan(mplanF);
  fftw_destroy_plan(mplanR);
  
  free_r3tensor(freq);
  free_r3tensor(data);

  ierr = KSPDestroy(ksp);CHKERRQ(ierr);
  ierr = VecDestroy(xr);CHKERRQ(ierr);
  ierr = VecDestroy(br);CHKERRQ(ierr);
  ierr = VecDestroy(xi);CHKERRQ(ierr);
  ierr = VecDestroy(bi);CHKERRQ(ierr);
  ierr = MatDestroy(A);CHKERRQ(ierr);

}

void fftpoisson_(double *rho, int *mxl, int *myl, int *mzl)
{

  int i,j,k,l;
  int lyz,igx,igy,igz,igyh,igzh;
  double vr, vi, v0, v,vy,vz;
  PetscScalar *brp, *bip;
  double view[64][2048][2],tmp[64];

  lyz = *myl*(*mzl);

  for(i=0;i<*mxl;++i) {
  
  for(l=0;l<lyz;++l) {
   minp[l][0]=rho[i*lyz + l];
   minp[l][1]=0.0;
  } 
  fftw_execute(mplanF);
  
  for(l=0;l<lyz;++l) {
    freq[i][l][0] = mout[l][0];
    freq[i][l][1] = mout[l][1];
    view[i][l][0] = freq[i][l][0];
    view[i][l][1] = freq[i][l][1];
    }
  
  }

// ********* Solver part *********
  vy = 4*PI*PI*dxs/(dys*myz[0]*myz[0]);
  vz = 4*PI*PI*dxs/(dzs*myz[1]*myz[1]);
  igyh = my/2;
  igzh = mz/2;

  for(j=0;j<*myl;++j) {
  for(k=0;k<*mzl;++k) {
// Calculate 'ky' and 'kz'
  igy = j + topo_.iyglobalc - 1;
  igz = k + topo_.izglobalc - 1;

    if((igy>igyh)&&(igz>igzh)) {
      v0=-(my-igy)*(my-igy)*vy-(mz-igz)*(mz-igz)*vz;}
    else if((igy>igyh)&&(igz<=igzh)) {
      v0=-(my-igy)*(my-igy)*vy-igz*igz*vz;}
    else if((igy<=igyh)&&(igz>igzh)) {
      v0=-igy*igy*vy-(mz-igz)*(mz-igz)*vz;}
    else {
      v0=-igy*igy*vy-igz*igz*vz;}

    VecGetArray(br, &brp);
    VecGetArray(bi, &bip);

    for(i=0;i<*mxl;++i) {
      igx = i + topo_.ixglobalc - 1;
      if((igx==1) || (igx==topo_.mxc-2)) 
	v = v0 - 5.0/3.0;
      else 
        v = v0 - 5.0/2.0;
      MatSetValues(A,1,&igx,1,&igx,&v,INSERT_VALUES); // BUG Overlap with previous ones, use INSERT
      vr = freq[i][j*(*mzl)+k][0]*dxs/(my*mz);
      vi = freq[i][j*(*mzl)+k][1]*dxs/(my*mz);
      brp[i] = vr;
      bip[i] = vi;
      tmp[i] = freq[i][j*(*mzl)+k][0]/(my*mz);
      if(igx==0 || igx==mx-1) {
        brp[i]=0;
        bip[i]=0;}
    }
    VecRestoreArray(br, &brp);
    VecRestoreArray(bi, &bip);

    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

    ierr = KSPSolve(ksp,br,xr);CHKERRQ(ierr);
    ierr = KSPSolve(ksp,bi,xi);CHKERRQ(ierr);

MatMult(A,xr,xi);
//MatMult(A,bi,xi);
VecView(br,PETSC_VIEWER_STDOUT_WORLD);
VecView(xi,PETSC_VIEWER_STDOUT_WORLD);
   MPI_Abort(commyz,1);
// Copy result to minp

    VecGetArray(xr, &brp);
    VecGetArray(xi, &bip);
    for(i=0;i<*mxl;++i) {
      data[i][j*(*mzl)+k][0]=brp[i];
      data[i][j*(*mzl)+k][1]=bip[i];
      view[i][j*(*mzl)+k][0]=brp[i];
      view[i][j*(*mzl)+k][1]=bip[i];
      tmp[i]=brp[i];
    }
    VecRestoreArray(xr, &brp);
    VecRestoreArray(xi, &bip);
  }
  }
//********************************

  for(i=0;i<*mxl;++i){
    for(l=0;l<lyz;++l) {
    minp[l][0]=data[i][l][0];
    minp[l][1]=data[i][l][1];
    }
  
  fftw_execute(mplanR);
  
  for(l=0;l<lyz;++l) {
    rho[i*lyz+l]=mout[l][0];
  }
  
  }
}


/**************************************************************************

 * Memory allocation routines
 * From BOUT++
 *************************************************************************/

double *rvector(int size)
{
  return (double*) malloc(sizeof(double)*size);
}

int *ivector(int size)
{
  return (int*) malloc(sizeof(int)*size);
}

void free_rvector(double *rv)
{
  free(rv);
}

void free_ivector(int *iv)
{
  free(iv);
}

double **rmatrix(int xsize, int ysize)
{
  long i;
  double **m;

  if((m = (double**) malloc(xsize*sizeof(double*))) == (double**) NULL) {

    printf("Error: could not allocate memory:%d\n", xsize);

    exit(1);

  }

  if((m[0] = (double*) malloc(xsize*ysize*sizeof(double))) == (double*) NULL) {

    printf("Error: could not allocate memory\n");

    exit(1);

  }

  for(i=1;i!=xsize;i++) {

    m[i] = m[i-1] + ysize;

  }

  return(m);

}



int **imatrix(int xsize, int ysize)

{

  long i;

  int **m;

  

  if((m = (int**) malloc(xsize*sizeof(int*))) == (int**) NULL) {

    printf("Error: could not allocate memory:%d\n", xsize);

    exit(1);

  }



  if((m[0] = (int*) malloc(xsize*ysize*sizeof(int))) == (int*) NULL) {

    printf("Error: could not allocate memory\n");

    exit(1);

  }

  for(i=1;i!=xsize;i++) {

    m[i] = m[i-1] + ysize;

  }



  return(m);

}



void free_rmatrix(double **m)

{

  free(m[0]);

  free(m);

}



void free_imatrix(int **m)

{

  free(m[0]);

  free(m);

}



double ***r3tensor(int nrow, int ncol, int ndep)

{

  int i,j;

  double ***t;



  /* allocate pointers to pointers to rows */

  t=(double ***) malloc((size_t)(nrow*sizeof(double**)));



  /* allocate pointers to rows and set pointers to them */

  t[0]=(double **) malloc((size_t)(nrow*ncol*sizeof(double*)));



  /* allocate rows and set pointers to them */

  t[0][0]=(double *) malloc((size_t)(nrow*ncol*ndep*sizeof(double)));



  for(j=1;j!=ncol;j++) t[0][j]=t[0][j-1]+ndep;

  for(i=1;i!=nrow;i++) {

    t[i]=t[i-1]+ncol;

    t[i][0]=t[i-1][0]+ncol*ndep;

    for(j=1;j!=ncol;j++) t[i][j]=t[i][j-1]+ndep;

  }



  /* return pointer to array of pointers to rows */

  return t;

}



void free_r3tensor(double ***t)

{

  free(t[0][0]);

  free(t[0]);

  free(t);

}





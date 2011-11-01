#
CC=mpicc
FC=mpif90

CFLAG= -O3
CCFLAG= -c -O3

FFLAG=-g
FCFLAG= -g -Ddebug -c -D_use_cvode #-D_test_poisson
#FCFLAG=-g -c -Ddebug -D_use_rk2 #-D_test_poisson

LIB=-Wl,-rpath,/home/jhli/mylib/petsc_mpich2_run/lib -L/home/jhli/mylib/petsc_mpich2_run/lib -lpetsc    -lX11 -Wl,-rpath,/home/jhli/mylib/petsc_mpich2_run/lib -L/home/jhli/mylib/petsc_mpich2_run/lib -lHYPRE -Wl,-rpath,/home/chess_usr/pgi/linux86-64/9.0-3/lib -lmpichcxx -lstd -lC -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lblacs -lsuperlu_dist_2.4 -lparmetis -lmetis -Wl,-rpath,/home/jhli/mylib/lib -L/home/jhli/mylib/cvode/lib -lsundials_fcvode -lsundials_fnvecparallel -lsundials_fnvecserial -lsundials_cvode -lsundials_nvecserial -lsundials_nvecparallel -Wl,-rpath,/home/chess_usr/lib/lapack-3.2.1/pgi -L/home/chess_usr/lib/lapack-3.2.1/pgi -llapack -lblas -Wl,-rpath,/home/chess_usr/lib/netcdf/lib -L/home/chess_usr/lib/netcdf/lib -lnetcdf -lnsl -lrt -L/home/chess_usr/mpi2/pgi/mpich2/lib -L/home/chess_usr/pgi/linux86-64/9.0-3/lib -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2 -ldl -lmpich -lopa -lmpl -lrt -lpthread -Wl,-rpath,/home/chess_usr/pgi/linux86-64/9.0-3/lib -lnspgc -lpgc -lmpichf90 -lpgf90 -lpgf90_rpm1 -lpgf902 -lpgftnrtl -lpgf90rtl -lm -lmpichcxx -lstd -lC -ldl -lmpich -lopa -lmpl -lrt -lpthread -lnspgc -lpgc -ldl -L/home/jhli/mylib/fftw3/lib -lfftw3_mpi -lfftw3
INC=-I/home/jhli/mylib/petsc_mpich2_run/include -I/home/jhli/mylib/petsc_mpich2_run/include/finclude -I/home/chess_usr/lib/netcdf/include -I/home/jhli/mylib/cvode/include/cvode -I/home/jhli/mylib/fftw3/include


FSOURCE=main.F90 poisson.F input.F90 fucomm.F90 solver.F90 init.F90 settopology.F90 dump.F90 rk.F90 boundary.F90 rhs.F90 cvodef.F90 final.F90 field.F90 fft.F90


OBJ=main.o poisson.o input.o fucomm.o solver.o init.o settopology.o dump.o rk.o boundary.o rhs.o cvodef.o final.o field.o fft.o

MODULE=module.o
EXEC=pde

$(EXEC): $(MODULE) $(OBJ) fortrandir.o fft2d.o
	$(FC) $(FLFLAG) -o $(EXEC) $(OBJ) fortrandir.o fft2d.o $(LIB)

$(MODULE): module.F90
	$(FC) $(FCFLAG)  $(INC) module.F90

${OBJ}: ${FSOURCE} 
	$(FC) $(FCFLAG) $(FSOURCE) $(INC)

fortrandir.o: fortrandir.c
	mpicc -c fortrandir.c
fft2d.o: fft2d.c
	mpicc -c -g fft2d.c $(INC)

read: read.F90
	$(FC) read.F90 $(INC) $(LIB) -o read

clean: 
	rm $(EXEC) *.mod *.o


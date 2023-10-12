
echo "Compiling MPI test"

lapack_path=/usr/lib/x86_64-linux-gnu/lapack
fishpack_path=$HOME/projects/def-lpnadeau/celiz2/fishpack-master/lib

mpif90 -O3 -o  exec_testmpi mpi.f90 -L$fishpack_path -lfishpack 
#mpifort -O3 -o  exec_testmpi mpi.f90 
#mpiifort -O3 -o  exec_testmpi mpi.f90 

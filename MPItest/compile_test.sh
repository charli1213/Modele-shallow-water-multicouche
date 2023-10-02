
echo "Compiling MPI test"

lapack_path=/usr/lib/x86_64-linux-gnu/lapack
fishpack_path=/home/charlesedouard/Desktop/Travail/fishpack/lib

mpif90 -O3 -o  exec_testmpi mpi.f90 
#mpifort -O3 -o  exec_testmpi mpi.f90 
#mpiifort -O3 -o  exec_testmpi mpi.f90 

To compile on the Hercules cluster, run

./configure --prefix=${HOME} CC=mpiicc CXX=mpiicpc F77=mpiifort FC=mpiifort   LDFLAGS="-L${MKLROOT}/lib/intel64  " LIBS="-lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -lpthread -lm  -lmkl_avx2  -lmkl_def -qopenmp -liomp5 -DMKL_ILP64  -lstdc++ -lifcore -limf -lifport"

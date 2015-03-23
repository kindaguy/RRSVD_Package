RRSVD_DIR = ../RRSVD
DATA_DIR = ../Data

all:	timewrap.x

dyn_array.o: dyn_array.f90 
	ifort  -c -I$(MKL_HOME)/include/intel64/lp64 -I$(MKL_HOME)/include dyn_array.f90


timewrap.x: dyn_array.o timewrap.f90
	ifort -L$(MKL_HOME)/lib/intel64 -L./ -I$(MKL_HOME)/include/intel64/lp64 -I$(MKL_HOME)/include -I./ \
	timewrap.f90 dyn_array.o $(RRSVD_DIR)/libdecompositionsmkl.so -o timewrap.x \
	$(MKL_HOME)/lib/intel64/libmkl_lapack95_lp64.a -m64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -liomp5 -lpthread -lm


exec: timewrap.x 
	./timewrap.x $(DATA_DIR)/theta_2500_0.dat

.PHONY: clean

clean:
	rm -f *.x *.o *.mod

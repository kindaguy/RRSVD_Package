BASEDIR  = $(realpath ..)

INCDIR = $(BASEDIR)/Incs

LIBDIR = $(BASEDIR)/Libs

LIB =  $(LIBDIR)/liblapacke.a $(LIBDIR)/libcblas.a $(LIBDIR)/librefblas.a $(LIBDIR)/liblapack.a -L$(LIBDIR) -lgfortran  -lrt  -lstdc++


RRSVD_LIB =  $(BASEDIR)/RRSVD/libdecompositions.so
RRSVD_INC = $(BASEDIR)/RRSVD/Include

DATA_DIR = ../Data

$(TARGET): $(OBJS)
	    $(CXX) $(CXXFLAGS) $(OBJS) $(RRSVD_LIB) $(LIB) -L.  -o $(TARGET)

TEBD_Trials.o: TEBD_Trials.cpp
	    $(CXX) $(CXXFLAGS) -I$(RRSVD_INC) -I$(INCDIR)  -c  TEBD_Trials.cpp


all:	timewrap.x

dyn_array.o: dyn_array.f90 
	$(FC)  -c -I$(RRSVD_INC) -I$(INCDIR) dyn_array.f90


timewrap.x: dyn_array.o timewrap.f90
	$(FC)  $(RRSVD_LIB)  $(LIB) -L./ -I$(RRSVD_INC) -I$(INCDIR)  -I./ \
	timewrap.f90 dyn_array.o -o timewrap.x


exec: timewrap.x 
	./timewrap.x $(DATA_DIR)/theta_900_0.dat

.PHONY: clean

clean:
	rm -f *.x *.o *.mod *.dat

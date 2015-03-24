CXXFLAGS =	-O2 -std=c++11  

TARGET = TEBD_Trials.x	

OBJS = TEBD_Trials.o 

BASEDIR  = $(realpath ..)

INCDIR = $(BASEDIR)/Incs

LIBDIR = $(BASEDIR)/Libs

LIB =  $(LIBDIR)/liblapacke.a $(LIBDIR)/libcblas.a $(LIBDIR)/librefblas.a $(LIBDIR)/liblapack.a -L$(LIBDIR) -lgfortran  -lrt  -lstdc++


RRSVD_LIB =  $(BASEDIR)/RRSVD/libdecompositions.so
RRSVD_INC = $(BASEDIR)/RRSVD/Include

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(RRSVD_LIB) $(LIB) -L.  -o $(TARGET)

TEBD_Trials.o: TEBD_Trials.cpp
	$(CXX) $(CXXFLAGS) -I$(RRSVD_INC) -I$(INCDIR)  -c  TEBD_Trials.cpp

exec: $(TARGET)
	LD_LIBRARY_PATH=./:$(LD_LIBRARY_PATH) ./TEBD_Trials.x ../Data/theta_900_0.dat 100 2 1e-4 1 1 1 
	
buildSample.x: buildSample.cpp
	$(CXX) $(CXXFLAGS) -I$(RRSVD_INC) -I$(INCDIR) -c buildSample.cpp 
	$(CXX) $(CXXFLAGS) buildSample.o $(RRSVD_LIB) $(LIB) -o buildSample.x 
generateSample: buildSample.x
	./buildSample.x 4900 4900 d3.dat ../../RRSVD_MKL_Release_20150317/decSA_0.dat

clean:
	rm -f $(OBJS) $(TARGET) *.dat



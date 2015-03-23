CXX = icc
CXXFLAGS =  -O2 -std=c++11

TARGET = TEBD_Trials.x

OBJS = TEBD_Trials.o

#Change the path here if needed
RRSVD_LIB = ../RRSVD

#Change the MKL home directory environment variable, or set the MKL library path if needed
MKL_FLAGS = -L${MKL_HOME}lib/intel64 -lmkl_intel_thread -liomp5

#MKL openMP num threads:
NUM_THREADS=1

all:
	make $(TARGET)
	make buildSample.x

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(RRSVD_LIB)/libdecompositionsmkl.so $(MKL_FLAGS) -o $(TARGET)
	rm $(OBJS)

TEBD_Trials.o: TEBD_Trials.cpp
	$(CXX) $(CXXFLAGS) -I$(RRSVD_LIB)/Include -I${MKL_INC} -c  -DMKL_ENABLED TEBD_Trials.cpp

exec: TEBD_Trials.x
	export MKL_NUM_THREADS=$(NUM_THREADS)
	./TEBD_Trials.x ../Data/theta_900_0.dat 100 2 1e-4 1 1 1

buildSample.x: buildSample.cpp
	$(CXX) $(CXXFLAGS) -c buildSample.cpp -DMKL_ENABLED -I../RRSVD/Include -I${MKL_INC}
	$(CXX) $(CXXFLAGS)  buildSample.o $(RRSVD_LIB)/libdecompositionsmkl.so  $(MKL_FLAGS) -L${MKL_HOME}lib/intel64 -lmkl_intel_lp64  -o buildSample.x
	rm buildSample.o

generateSample: buildSample.x
		./buildSample.x 1000 1000 1 criticalSample.dat

.PHONY: clean
clean:
	rm  $(TARGET) buildSample.x 


#Garbage

	#icc buildSample.o ../RRSVD/libdecompositionsmkl.so   /soft/local/modules/compiler/intel/compxe.2013.sp1.4.211/mkl/lib/intel64/libmkl_intel_lp64.so /soft/local/modules/compiler/intel/compxe.2013.sp1.4.211/lib/intel64/libintlc.so.5 -liomp5 -o buildSample.x
	

CXX = icc
CXXFLAGS =  -O2 -std=c++11

TARGET = libdecompositionsmkl.so

OBJS = decompositions.o

INCDIR = ./Include

SRCDIR = ./Src

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS)  -funroll-loops -m64  -Wl,--no-as-needed -L${MKL_HOME}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread  -shared -lgfortran -lstdc++ -o  libdecompositionsmkl.so decompositions.o
	rm -f decompositions.o


decompositions.o: $(SRCDIR)/decompositions.cpp $(INCDIR)/decompositions.h
	$(CXX) $(CXXFLAGS)  -DMKL_ENABLED  -export -funroll-loops -fPIC -I$(INCDIR) -I${MKL_INC} -c $(SRCDIR)/decompositions.cpp

#libdecompositionsmkl.so: decompositions.o
#	$(CXX) $(CXXFLAGS)  -funroll-loops -m64  -Wl,--no-as-needed -L${MKL_HOME}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread  -shared -lgfortran -lstdc++ -o  libdecompositionsmkl.so decompositions.o

.PHONY: clean
clean:
	rm -f $(OBJS) $(TARGET)



#decompositions.o: decompositions.cpp decompositions.h
#	g++ -std=c++11 -O3 -DMKL_ENABLED -fPIC  -c decompositions.cpp

#libdecompositionsmkl.so: decompositions.o
#	g++ -std=c++11 -v -O3 -m64  -Wl,--no-as-needed -L$(MKL_HOME)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread  -shared -lgfortran -lstdc++ -o libdecompositionsmkl.so decompositions.o

exec: TEBD_Trials.x
	export MKL_NUM_THREADS=$(NUM_THREADS)
	LD_LIBRARY_PATH=./:${LD_LIBRARY_PATH} ./TEBD_Trials.x 1500 750 100 2 1 1 theta_4900_0.dat


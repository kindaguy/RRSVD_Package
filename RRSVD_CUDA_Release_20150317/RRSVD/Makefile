CXXFLAGS =	-O2 -std=c++11 -DCUDA_ENABLED 

TARGET = libdecompositionscuda.so	

OBJS = decompositions.o		

ROOTDIR = $(realpath ..)

INCDIR = $(ROOTDIR)/Incs

BASEDIR = $(realpath .)

INCLUDE = $(BASEDIR)/Include

SRC = $(BASEDIR)/Src



CULA_LIB_PATH_64 = ${HOME}/cula/lib64
CULA_INC_PATH = ${HOME}/cula/include
CUDA_INC_PATH = ${CUDA_INC}
CUDA_LIB_PATH_64 = ${CUDA_LIB}
#./liblapacke.a ./libcblas.a ./librefblas.a 


LIBS = -lgfortran -lstdc++  -L$(CUDA_LIB_PATH_64) -lcuda -lcudart -lcurand  -lcublas -L$(CULA_LIB_PATH_64) -lcula_lapack 


$(TARGET): $(OBJS) 
	$(CXX) -o $(TARGET) -shared -Wl,--no-as-needed $(OBJS) $(LIBS)
	rm $(OBJS)

decompositions.o: $(SRC)/decompositions.cpp $(INCLUDE)/decompositions.h
	$(CXX) $(CXXFLAGS) -c -export -fPIC -I./ -I$(CUDA_INC_PATH) -I$(CULA_INC_PATH) -I$(INCDIR) -I$(INCLUDE) $(SRC)/decompositions.cpp

#libdecompositionscuda.so: decompositions.o
#	$(CXX) $(CXXFLAGS) -shared -Wl,--no-as-needed decompositions.o $(LIBS) -o libdecompositionscuda.so
#all:	$(TARGET) 

clean:
	rm -f $(OBJS) $(TARGET)

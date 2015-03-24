CXXFLAGS =	-O2 -std=c++11 -DCUDA_ENABLED 


SOURCE = TEBD_Trials.cpp
TARGET = TEBD_Trials.x	
OBJS = TEBD_Trials.o		

ROOTDIR = $(realpath ..)
INCDIR = $(ROOTDIR)/Incs

BASEDIR = $(realpath .)

RRSVD_LIB =  $(ROOTDIR)/RRSVD/libdecompositionscuda.so
RRSVD_INC = $(ROOTDIR)/RRSVD/Include

CULA_LIB_PATH_64 = ${HOME}/cula/lib64
CULA_INC_PATH = ${HOME}/cula/include
CUDA_INC_PATH = ${CUDA_INC}
CUDA_LIB_PATH_64 = ${CUDA_LIB}
#./liblapacke.a ./libcblas.a ./librefblas.a 
LIBS = -lgfortran -lstdc++  -L$(CUDA_LIB_PATH_64) -lcuda -lcudart -lcurand  -lcublas -L$(CULA_LIB_PATH_64) -lcula_lapack 
$(TARGET):  $(OBJS) 
	$(CXX) -o $(TARGET) $(OBJS) $(PATHS) -Wl,--no-as-needed $(RRSVD_LIB) $(LIBS)

$(OBJS): $(SOURCE) 
	$(CXX) $(CXXFLAGS) -c -I$(CUDA_INC_PATH) -I$(CULA_INC_PATH) -I$(INCDIR) -I$(RRSVD_INC) $(SOURCE)
all:	$(TARGET) 

exec: $(TARGET)
	LD_LIBRARY_PATH=./:$(CULA_LIB_PATH_64):$(LD_LIBRARY_PATH) $(TARGET)  $(ROOTDIR)/Data/theta_900_0.dat 100 2 1e-3 1 1 5

clean:
	rm -f $(OBJS) $(TARGET)

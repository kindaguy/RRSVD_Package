CXXFLAGS =	-O2 -std=c++11  

TARGET = libdecompositions.so

OBJS =  decompositions.o


ROOTDIR = $(realpath ..)

INCDIR = $(ROOTDIR)/Incs

LIBDIR = $(ROOTDIR)/Libs

BASEDIR = $(realpath .)

INCLUDE  = $(BASEDIR)/Include

SRC = $(BASEDIR)/Src

LIB =  $(LIBDIR)/liblapacke.a $(LIBDIR)/libcblas.a $(LIBDIR)/librefblas.a $(LIBDIR)/liblapack.a -lgfortran  -lrt 


$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -shared  $(OBJS) $(LIB) -o $(TARGET)
	rm -f $(OBJS)

$(OBJS): $(SRC)/decompositions.cpp $(INCLUDE)/decompositions.h
	$(CXX) $(CXXFLAGS) -c -fPIC -Wl,--no-as-needed  -I$(INCDIR) -I$(INCLUDE) $(LIB) $(SRC)/decompositions.cpp 	


clean:
	rm -f $(TARGET)
	echo $(ROOTDIR)

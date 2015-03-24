#Reduced-Rank Randomized SVD package (RRSVD-Package) 
This repository contains the implementation of the method described in the paper 'Fast Time-Evolving Block-Decimation algorithm through Reduced-Rank Radomized Singular Value Decomposition'

[![Linux Build Status](https://travis-ci.org/qubit-ulm/ebs.svg?branch=master)](https://travis-ci.org/qubit-ulm/ebs)

# Introduction 

The computational complexity of the Singular Value Decomposition of an m x n matrix **A** is *O*(m n^2). For large matrices, the SVD can therefore require a significative amount of time. In many situations the full SVD is not needed: only the largest singular values (and corresponding left- and right-singular vectors) are actually needed.  It is possible to avoid the *full* SVD of **A** and compute only the first **k** singular values and corresponding singular vectors by using *Truncated SVD* methods; such methods are standard tools in data-classification algorithms, signal-processing and other research fields. The  Implicitly Restarted Arnoldi Method  and the Lanczos-Iteration algorithms, both belonging to the Krylov-subspace iterative methods are two examples.  The Reduced-Rank Singular Value Decomposition (**RRSVD**),  originally presented in by [N. Halko *et al.*](http://arxiv.org/abs/0909.4061), is a *randomized* truncated SVD. It is particularly suited to decompose structured matrices, such as the one appearing in some simulation of non-critical quantum systems (e.g. the Time-Evolving Block-Decimation algorithm, based on the Matrix Product States formalism). Most interestingly, the algorithm is insensitive to the quality of the random number generator used, delivers highly accurate results and is, despite its random nature, very stable: the probability of failure can be made arbitrarily small with a minor impact on the computational resource required.

# Usage
The RRSVD routine comes in three versions: 

* CBLAS-LAPACKE-based (BL-RRSVD)
* Intel MKL-based (MKL-RRSVD)
* Nvidia CUDA-based (CUDA-RRSVD). 

The interface to the library is uniform for all the different implementations. For every implementation, we provide single/double and real/complex RRSVD. 

The RRSVD-Package contains the three implementations of RRSVD. For the sake of modularity we provide them in three different folders. In each folder we put the library source files (RRSVD folder), some examplifying C++ code showing how to use the library (Examples), a Fortran wrapper that allows to call the RRSVD routine directly from Fortran code (FortranWrapper) and some sample matrices to use for testing (Data). All the folders and subfolders have README files illustrating their content and use. Where necessary, a Makefile is provided as well.

## General remarks
The RRSVD routine, for the time being, works only on square (n x n) and tall (m x n, m>n) matrices. In order to handle wide matrices (m x n, n>m) matrices, you must transpose (or adjoin) them in your own code and then manipulate the output provided by RRSVD as needed.
The RRSVD code, and most of the example programs, are written in C++11 and use some of the new libraries that come with this version of the language. 

# Building

## Requirements
RRSVD-Package has the following requirements on Unix systems

* BL-RRSVD
  * gcc >= 4.7
  * CBLAS and LAPAcke libraries [here](http://www.netlib.org/lapack/lapacke.html)

* MKL-RRSVD
  * icc >= 14.0
  * Intel MKL library (11.1.4)

* CUDA-RRSVD
  * CUDA >= 5.0
  * CUDA capable device with compute capability >= 2.0
  * CULA (Release 17)
  * gcc >=4.7

#Warning
While the building of the MKL- and CUDA-RRSVD is rather standard, the creation of the BL-RRSVD library requires a little bit of attention. In the RRSVD_BL_Release folder we include a precompiled version of LAPacke and CBLAS. The libraries have been compiled with the GNU gcc 4.7 compiler and should work on most Linux distributions with gcc 4.7 installed. Before compiling this libraries we have set complex types to be represented as a struct, with the real and imaginary part accessible via the field names .r and .i. This is one of the options provided by the LAPacke installation. If you recompile the LAPacke and CBLAS libraries pay attention to set the right complex numbers representation in the appropriate configuration files.

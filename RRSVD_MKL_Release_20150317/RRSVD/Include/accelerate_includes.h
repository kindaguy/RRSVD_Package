/*
 * accelerate_includes.h
 *
 *  Created on: Sep 29, 2014
 *      Author: Tama
 */

#ifndef ACCELERATE_INCLUDES_H_
#define ACCELERATE_INCLUDES_H_

#ifdef __linux__
/*LINUX IMPLEMENTATION*/
/*the libraries:

*liblapacke.a
*liblapack.a
*libcblas.a
*librefblas.a

are supposed to be in some standard place*/
//!!! The constant MKL_ENABLED must be defined at compile time
#ifdef CUDA_ENABLED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas.h>
#include <cula_lapack.h>
#include <cula_lapack_device.h>
#include <curand.h>
#endif

#ifdef MKL_ENABLED
//#include <mkl_lapacke.h>
#include <mkl_cblas.h>
#include <mkl.h>

//#include <lapacke_mangling.h>
#else
#include <lapacke.h>
#include <cblas.h>
//#include <lapacke_mangling.h>
#endif

typedef lapack_int integer;

typedef float __lpkreal;
typedef double __lpkdoublereal;
typedef lapack_complex_float __lpkcomplex;
typedef lapack_complex_double __lpkdoublecomplex;

template <typename T>
struct constants{
//	const __lpkcomplex complex_one ={1.f,0.f};
//		const __lpkcomplex complex_zero = {0.f,0.f};
//		const __lpkreal real_one = 1.f;
//		const __lpkreal real_zero = 0.f;
//		const __lpkreal typed_one = 1.f;
//		const __lpkreal typed_zero = 0.f;
//		const __lpkreal typed_mone = -1.f;
};

template <>
struct constants<__lpkreal>{
	const __lpkcomplex complex_one ={1.f,0.f};
	const __lpkcomplex complex_zero = {0.f,0.f};
	const __lpkreal real_one = 1.f;
	const __lpkreal real_zero = 0.f;
	const __lpkreal typed_one = 1.f;
	const __lpkreal typed_zero = 0.f;
	const __lpkreal typed_mone = -1.f;
	const char typeCheck = 's';
};

template <>
struct constants<__lpkdoublereal>{
	const __lpkdoublecomplex complex_one ={1.,0.};
	const __lpkdoublecomplex complex_zero = {0.,0.};
	const __lpkdoublereal real_one = 1.;
	const __lpkdoublereal real_zero = 0.;
	const __lpkdoublereal typed_one = 1.;
	const __lpkdoublereal typed_zero = 0.;
	const __lpkdoublereal typed_mone = -1.;
	const char typeCheck = 'd';
};

template <>
struct constants<__lpkcomplex>{
	const __lpkcomplex complex_one ={1.f,0.f};
	const __lpkcomplex complex_zero = {0.f,0.f};
	const __lpkreal real_one = 1.f;
	const __lpkreal real_zero = 0.f;
	const __lpkcomplex typed_one = {1.f,0.f};
	const __lpkcomplex typed_zero = {0.f,0.f};
	const __lpkcomplex typed_mone = {-1.f,0.f};
	const char typeCheck = 'c';
};


template <>
struct constants<__lpkdoublecomplex>{
	const __lpkdoublecomplex complex_one ={1.,0.};
	const __lpkdoublecomplex complex_zero = {0.,0.};
	const __lpkdoublereal real_one = 1.;
	const __lpkdoublereal real_zero = 0.;
	const __lpkdoublecomplex typed_one = {1.,0.};
	const __lpkdoublecomplex typed_zero = {0.,0.};
	const __lpkdoublecomplex typed_mone = {-1.,0.};
const char typeCheck = 'z';
};


//#ifdef SINGLEPREC
//typedef float lpkreal;
//typedef lapack_complex_float lpkcomplex;
//const lpkcomplex complex_one = {1.f,0.f};
//const lpkcomplex complex_zero = {0.f,0.f};
//const lpkreal real_one = 1.0f;
//const lpkreal real_zero = 0.0f;
////Real Norm 2 single prec
//#define NORM2REAL cblas_snrm2
//
//#ifdef REAL
////SINGLE PREC REAL
//const lpkreal typed_one = 1.0f;
//const lpkreal typed_zero = 0.0f;
//const lpkreal typed_mone = -1.0f;
//#define GEQRF_ sgeqrf_
//#define UNGQR_ sungqr_
//#define GEMM cblas_sgemm
//#define AXPY cblas_saxpy
//#define DOT cblas_sdotc
//#define COPY cblas_scopy
//#define GESVD_ sgesvd_
//#define DAG CblasTrans
//
//#else
//
//const lpkcomplex typed_one = {1.0f,0.f};
//const lpkcomplex typed_mone = {-1.0f,0.f};
//const lpkcomplex typed_zero = {0.0f,0.0f};
////SINGLE PREC COMPLEX
//#define GEQRF_ cgeqrf_
//#define UNGQR_ cungqr_
//#define GEMM cblas_cgemm
//#define AXPY cblas_caxpy
//#define DOT cblas_cdotc_sub
//#define COPY cblas_ccopy
//#define GESVD_ cgesvd_
//#define DAG CblasConjTrans
//#endif //REAL
//
//
//#else //DOUBLEPREC
//typedef double  lpkreal;
//typedef lapack_complex_double lpkcomplex;
//const lpkcomplex complex_one = {1.,0.};
//const lpkcomplex complex_zero = {0.,0.};
//const lpkreal real_one = 1.0;
//const lpkreal real_zero = 0.0;
//
////Real Norm 2 double prec
//#define NORM2REAL cblas_dnrm2
//#ifdef REAL
////DOULBE PREC REAL
//const lpkreal typed_one = 1.0;
//const lpkreal typed_mone = -1.0;
//const lpkreal typed_zero = 0.0;
//#define GEQRF_ dgeqrf_
//#define UNGQR_ dungqr_
//#define GEMM cblas_dgemm
//#define AXPY cblas_daxpy
//#define DOT cblas_ddotc
//#define COPY cblas_dcopy
//#define GESVD_ dgesvd_
//#define DAG CblasTrans
//#else
////DOUBLE PREC COMPLEX
//const lpkcomplex typed_one = {1.0,0.};
//const lpkcomplex typed_mone = {-1.0,0.};
//const lpkcomplex typed_zero = {0.0,0.0};
//
//
//#define GEQRF_ zgeqrf_
//#define UNGQR_ zungqr_
//#define GEMM cblas_zgemm
//#define AXPY cblas_zaxpy
//#define DOT cblas_zdotc_sub
//#define COPY cblas_zcopy
//#define GESVD_ zgesvd_
//#define DAG CblasConjTrans
//#endif
//
//#endif //END DOUBLE PREC


#endif //LINUX




#ifdef __APPLE__
#include "/System/Library/Frameworks/Accelerate.framework/Headers/Accelerate.h"
//#include <Accelerate.framework/Headers/Accelerate.h>
#include "/System/Library/Frameworks/Accelerate.framework/Frameworks/vecLib.framework/Headers/clapack.h"
#include "/System/Library/Frameworks/Accelerate.framework/Frameworks/vecLib.framework/Headers/cblas.h"



/******APPLE IMPLEMENTATION******/
typedef __CLPK_integer integer;

typedef __CLPK_real __lpkreal;
typedef __CLPK_doublereal __lpkdoublereal;
typedef __CLPK_complex __lpkcomplex;
typedef __CLPK_doublecomplex __lpkdoublecomplex;

template <typename T>
struct constants{
//	const __lpkcomplex complex_one ={1.f,0.f};
//		const __lpkcomplex complex_zero = {0.f,0.f};
//		const __lpkreal real_one = 1.f;
//		const __lpkreal real_zero = 0.f;
//		const __lpkreal typed_one = 1.f;
//		const __lpkreal typed_zero = 0.f;
//		const __lpkreal typed_mone = -1.f;
};

template <>
struct constants<__lpkreal>{
	const __lpkcomplex complex_one ={1.f,0.f};
	const __lpkcomplex complex_zero = {0.f,0.f};
	const __lpkreal real_one = 1.f;
	const __lpkreal real_zero = 0.f;
	const __lpkreal typed_one = 1.f;
	const __lpkreal typed_zero = 0.f;
	const __lpkreal typed_mone = -1.f;
	const char typeCheck = 's';
};

template <>
struct constants<__lpkdoublereal>{
	const __lpkdoublecomplex complex_one ={1.,0.};
	const __lpkdoublecomplex complex_zero = {0.,0.};
	const __lpkdoublereal real_one = 1.;
	const __lpkdoublereal real_zero = 0.;
	const __lpkdoublereal typed_one = 1.;
	const __lpkdoublereal typed_zero = 0.;
	const __lpkdoublereal typed_mone = -1.;
	const char typeCheck = 'd';
};

template <>
struct constants<__lpkcomplex>{
	const __lpkcomplex complex_one ={1.f,0.f};
	const __lpkcomplex complex_zero = {0.f,0.f};
	const __lpkreal real_one = 1.f;
	const __lpkreal real_zero = 0.f;
	const __lpkcomplex typed_one = {1.f,0.f};
	const __lpkcomplex typed_zero = {0.f,0.f};
	const __lpkcomplex typed_mone = {-1.f,0.f};
	const char typeCheck = 'c';
};


template <>
struct constants<__lpkdoublecomplex>{
	const __lpkdoublecomplex complex_one ={1.,0.};
	const __lpkdoublecomplex complex_zero = {0.,0.};
	const __lpkdoublereal real_one = 1.;
	const __lpkdoublereal real_zero = 0.;
	const __lpkdoublecomplex typed_one = {1.,0.};
	const __lpkdoublecomplex typed_zero = {0.,0.};
	const __lpkdoublecomplex typed_mone = {-1.,0.};
const char typeCheck = 'z';
};

//#ifdef SINGLEPREC
//typedef __CLPK_real lpkreal;
//typedef __CLPK_complex lpkcomplex;
//const lpkcomplex complex_one = {1.f,0.f};
//const lpkcomplex complex_zero = {0.f,0.f};
//const lpkreal real_one = 1.0f;
//const lpkreal real_zero = 0.0f;
////Real Norm 2 single prec
//#define NORM2REAL cblas_snrm2
//
//#ifdef REAL
////SINGLE PREC REAL
//const lpkreal typed_one = 1.0f;
//const lpkreal typed_zero = 0.0f;
//const lpkreal typed_mone = -1.0f;
//
//#define GEQRF_ sgeqrf_
//#define UNGQR_ sungqr_
//#define GEMM cblas_sgemm
//#define AXPY cblas_saxpy
//#define DOT cblas_sdotc
//#define COPY cblas_scopy
//#define GESVD_ sgesvd_
//
//#else
//
//const lpkcomplex typed_one = {1.0f,0.f};
//const lpkcomplex typed_mone = {-1.0f,0.f};
//const lpkcomplex typed_zero = {0.0f,0.0f};
////SINGLE PREC COMPLEX
//#define GEQRF_ cgeqrf_
//#define UNGQR_ cungqr_
//#define GEMM cblas_cgemm
//#define AXPY cblas_caxpy
//#define DOT cblas_cdotc_sub
//#define COPY cblas_ccopy
//#define GESVD_ cgesvd_
//#endif //REAL


//#else //DOUBLEPREC
//typedef __CLPK_doublereal lpkreal;
//typedef __CLPK_doublecomplex lpkcomplex;
//const lpkcomplex complex_one = {1.,0.};
//const lpkcomplex complex_zero = {0.,0.};
//const lpkreal real_one = 1.0;
//const lpkreal real_zero = 0.0;
//
////Real Norm 2 double prec
////#define NORM2REAL cblas_dnrm2
////#ifdef REAL
////DOULBE PREC REAL
//template <>
//struct constants<lpkreal>{
//	const lpkcomplex complex_one ={1.,0.};
//	const lpkcomplex complex_zero = {0.,0.};
//	const lpkreal real_one = 1.;
//	const lpkreal real_zero = 0.;
//	const lpkreal typed_one = 1.;
//	const lpkreal typed_zero = 0.;
//	const lpkreal typed_mone = -1.;
//	const char typeCheck = 'd';
//};
//const lpkreal typed_one = 1.0;
//const lpkreal typed_mone = -1.0;
//const lpkreal typed_zero = 0.0;
//#define GEQRF_ dgeqrf_
//#define UNGQR_ dorgqr_
//#define GEMM cblas_dgemm
//#define AXPY cblas_daxpy
//#define DOT cblas_ddotc
//#define COPY cblas_dcopy
//#define GESVD_ dgesvd_

//#else
//DOUBLE PREC COMPLEX
//Questa struttura dovra` essere modificata alla fine
//template<>
//struct constants<lpkcomplex>{
//	const lpkcomplex complex_one ={1.,0.};
//	const lpkcomplex complex_zero = {0.,0.};
//	const lpkreal real_one = 1.;
//	const lpkreal real_zero = 0.;
//	const lpkcomplex typed_one = {1.,0.};
//	const lpkcomplex typed_zero = {0.,0.};
//	const lpkcomplex typed_mone = {-1.,0.};
//	const char typeCheck = 'z';
//};



//const lpkcomplex typed_one = {1.0,0.};
//const lpkcomplex typed_mone = {-1.0,0.};
//const lpkcomplex typed_zero = {0.0,0.0};


//#define GEQRF_ zgeqrf_
//#define UNGQR_ zungqr_
//#define GEMM cblas_zgemm
//#define AXPY cblas_zaxpy
//#define DOT cblas_zdotc_sub
//#define COPY cblas_zcopy
//#define GESVD_ zgesvd_

//#endif

//#endif //END DOUBLE PREC
#endif //APPLE
/******Tama's constants******/

#endif /* ACCELERATE_INCLUDES_H_ */

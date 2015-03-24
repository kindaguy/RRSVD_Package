/*
 * decompositions.cpp
 *
 *  Created on: Oct 28, 2014
 *      Author: Tama
 */

/***IMPLEMENTATIONS***/

#include "decompositions.h"

template <typename T>
QRfactLaPack<T>::QRfactLaPack(int m,int n){
	nrows = m;
	ncols = n;
	mQ = NULL;
	mR = NULL;
}

template <typename T>
QRfactLaPack<T>::~QRfactLaPack(){

	//if (mQ!=NULL)
	//	delete  mQ;
	if (mR != NULL)
		delete mR;


}

//cublasHandle_t *handle;

/************************************************************************/
/*NEW VERSION: "nasty functions templated"
 *
 */

/*GEMM: the problem comes with the parameters \alpha and \beta: if the product involves real matrices, they must be passed to the function
 * otherwise they ust be referenced by the function.*/
template <typename T>
void TamaGEMM(char opA, char opB, tamaClss::genMatrix<T> &A, tamaClss::genMatrix<T> &B, tamaClss::genMatrix<T> &C, const T c1, const T c2){


}
template<>
void TamaGEMM<__lpkdoublecomplex>(char opA, char opB, tamaClss::genMatrix<__lpkdoublecomplex> &A, tamaClss::genMatrix<__lpkdoublecomplex> &B, tamaClss::genMatrix<__lpkdoublecomplex> &C, const __lpkdoublecomplex c1, const __lpkdoublecomplex c2){


#ifdef MKL_ENABLED
CBLAS_TRANSPOSE tA,tB;
#else
enum CBLAS_TRANSPOSE tA,tB;
#endif

	if(opA == 'N'){
		tA = CblasNoTrans;

	}
	else{//
		tA = CblasConjTrans;

	}
	if(opB == 'N'){
		tB = CblasNoTrans;

	}
	else{//
		tB = CblasConjTrans;

	}

#ifdef CUDA_ENABLED
    if(opA=='T')
        opA='C';
    if(opB == 'T')
        opB = 'C';

    cuDoubleComplex c1dev,c2dev;
    c1dev.x=lapack_complex_double_real(c1);
    c1dev.y=lapack_complex_double_imag(c1);
    c2dev.x=lapack_complex_double_real(c2);
    c2dev.y=lapack_complex_double_imag(c2);

    cublasZgemm(opA,opB,C.nrows,C.ncols,(opA=='N'?A.ncols:A.nrows), c1dev,(cuDoubleComplex*)(A.raw),A.nrows,(cuDoubleComplex*) (B.raw),B.nrows,c2dev,(cuDoubleComplex*)(C.raw),C.nrows);

#else

    cblas_zgemm(CblasColMajor,tA,tB,C.nrows,C.ncols,(opA=='N'?A.ncols:A.nrows),&(c1),A.raw,A.nrows,B.raw,B.nrows,&(c2),C.raw,C.nrows);

#endif
}


template<>
void TamaGEMM<__lpkdoublereal>(char opA, char opB, tamaClss::genMatrix<__lpkdoublereal> &A, tamaClss::genMatrix<__lpkdoublereal> &B, tamaClss::genMatrix<__lpkdoublereal> &C, const __lpkdoublereal c1, const __lpkdoublereal c2){

#ifdef MKL_ENABLED
CBLAS_TRANSPOSE tA,tB;
#else
enum CBLAS_TRANSPOSE tA,tB;
#endif

	if(opA == 'N'){
		tA = CblasNoTrans;
	}
	else{//
		tA = CblasTrans;
	}
	if(opB == 'N'){
			tB = CblasNoTrans;
		}
		else{//
			tB = CblasTrans;
		}
#ifdef CUDA_ENABLED
        //cublasInit();
        cout << "\nCall to cublasDgemm: \n";
 cublasDgemm(opA,opB,C.nrows,C.ncols,(opA=='N'?A.ncols:A.nrows), c1, A.raw,A.nrows, B.raw,B.nrows,c2,C.raw,C.nrows);

cout <<"\ncublasDgemm completed!!!\n";

#else

	cblas_dgemm(CblasColMajor,tA,tB,C.nrows,C.ncols,(opA=='N'?A.ncols:A.nrows),(c1),A.raw,A.nrows,B.raw,B.nrows,(c2),C.raw,C.nrows);
#endif
}



template<>
void TamaGEMM<__lpkcomplex>(char opA, char opB, tamaClss::genMatrix<__lpkcomplex> &A, tamaClss::genMatrix<__lpkcomplex> &B, tamaClss::genMatrix<__lpkcomplex> &C, const __lpkcomplex c1, const __lpkcomplex c2){

	#ifdef MKL_ENABLED
	CBLAS_TRANSPOSE tA,tB;
#else
	enum CBLAS_TRANSPOSE tA,tB;
#endif

	if(opA == 'N'){
		tA = CblasNoTrans;

	}
	else{//
		tA = CblasConjTrans;

	}
	if(opB == 'N'){
		tB = CblasNoTrans;

	}
	else{//
		tB = CblasConjTrans;

	}
#ifdef CUDA_ENABLED
    if(opA=='T')
        opA='C';
    if(opB == 'T')
        opB = 'C';

    cuComplex c1dev,c2dev;
    c1dev.x=lapack_complex_float_real(c1);
    c1dev.y=lapack_complex_float_imag(c1);
    c2dev.x=lapack_complex_float_real(c2);
    c2dev.y=lapack_complex_float_imag(c2);

    cublasCgemm(opA,opB,C.nrows,C.ncols,(opA=='N'?A.ncols:A.nrows), c1dev,(cuComplex*)A.raw,A.nrows,(cuComplex*) B.raw,B.nrows,c2dev,(cuComplex*)C.raw,C.nrows);

#else

	cblas_cgemm(CblasColMajor,tA,tB,C.nrows,C.ncols,(opA=='N'?A.ncols:A.nrows),&(c1),A.raw,A.nrows,B.raw,B.nrows,&(c2),C.raw,C.nrows);

#endif

}

template<>
void TamaGEMM<__lpkreal>(char opA, char opB, tamaClss::genMatrix<__lpkreal> &A, tamaClss::genMatrix<__lpkreal> &B, tamaClss::genMatrix<__lpkreal> &C, const __lpkreal c1, const __lpkreal c2){
#ifdef MKL_ENABLED
	CBLAS_TRANSPOSE tA,tB;
#else
	enum CBLAS_TRANSPOSE tA,tB;
#endif
	if(opA == 'N'){
		tA = CblasNoTrans;
	}
	else{//
		tA = CblasTrans;
	}
	if(opB == 'N'){
			tB = CblasNoTrans;
		}
		else{//
			tB = CblasTrans;
		}
#ifdef CUDA_ENABLED
        //cublasInit();
 cublasSgemm(opA,opB,C.nrows,C.ncols,(opA=='N'?A.ncols:A.nrows), c1, A.raw,A.nrows, B.raw,B.nrows,c2,C.raw,C.nrows);


#else


	cblas_sgemm(CblasColMajor,tA,tB,C.nrows,C.ncols,(opA=='N'?A.ncols:A.nrows),(c1),A.raw,A.nrows,B.raw,B.nrows,(c2),C.raw,C.nrows);
#endif
}
//Are missing


/************************************************************************/
/*Here the problem is the name of the function that provides the matrix Q:
 * for complexes it is ?UNGQR
 * for reals it is ?ORGQR
 * Moreover we must pay attention to the way we extract the lwork size*/

template <typename T>
//GEQRF_(&m,&n,mat,&m,tau,lwork,&lworkSize,&info);
void TamaGEQRF(int *m, int *n, T *mat, int *lda, T * tau, T * lwork, int *lworkSize, int *info){

}

template <>
void TamaGEQRF<__lpkdoublecomplex>(int *m, int *n, __lpkdoublecomplex *mat, int *lda, __lpkdoublecomplex * tau, __lpkdoublecomplex * lwork, int *lworkSize, int *info){
#ifdef CUDA_ENABLED
    culaDeviceZgeqrf(*m,*n,(culaDeviceDoubleComplex *) mat,*lda,(culaDeviceDoubleComplex *) tau);
    *info = culaGetErrorInfo();
#elif defined MKL_ENABLED
	*info = LAPACKE_zgeqrf(LAPACK_COL_MAJOR,*m,*n,mat,*lda,tau);
#else
	zgeqrf_(m,n,mat,lda,tau,lwork,lworkSize,info);
#endif
}


template <>
void TamaGEQRF<__lpkcomplex>(int *m, int *n, __lpkcomplex *mat, int *lda, __lpkcomplex * tau, __lpkcomplex * lwork, int *lworkSize, int *info){
#ifdef CUDA_ENABLED
    culaDeviceCgeqrf(*m,*n,(culaDeviceFloatComplex *) mat,*lda,(culaDeviceFloatComplex *) tau);
    *info = culaGetErrorInfo();

#elif defined MKL_ENABLED
	*info = LAPACKE_cgeqrf(LAPACK_COL_MAJOR,*m,*n,mat,*lda,tau);
#else
	cgeqrf_(m,n,mat,lda,tau,lwork,lworkSize,info);
#endif
}

template <>
void TamaGEQRF<__lpkdoublereal>(int *m, int *n, __lpkdoublereal *mat, int *lda, __lpkdoublereal * tau, __lpkdoublereal * lwork, int *lworkSize, int *info){

#ifdef CUDA_ENABLED
    //_DEBUG
   // cout << endl <<"Debug: " <<endl;
   // cout << *m << endl;
   // cout << *n << endl;
   //// cout << mat[0] << " " << mat[1] << endl;
   // cout  << *lda << endl;
   // cout << tau << endl;
   // {
   //     int pippo;
   //     cin >> pippo;
   // }

    culaDeviceDgeqrf(*m,*n, mat,*lda, tau);
    *info = culaGetErrorInfo();
    cout << "\nInfo = " << *info << endl;

#elif defined MKL_ENABLED
	//printf("\nm= %d, n=%d, lda = %d\n",*m,*n,*lda);
	*info = LAPACKE_dgeqrf(LAPACK_COL_MAJOR,*m,*n,mat,*lda,tau);
	if(*info !=0)
			printf("\nProblem in dgeqrf\n");
#else
	dgeqrf_(m,n,mat,lda,tau,lwork,lworkSize,info);
#endif
}

template <>
void TamaGEQRF<__lpkreal>(int *m, int *n, __lpkreal *mat, int *lda, __lpkreal * tau, __lpkreal * lwork, int *lworkSize, int *info){
#ifdef CUDA_ENABLED
    
    culaDeviceSgeqrf(*m,*n, mat,*lda, tau);
    *info = culaGetErrorInfo();

#elif defined MKL_ENABLED
	*info = LAPACKE_sgeqrf(LAPACK_COL_MAJOR,*m,*n,mat,*lda,tau);
#else
	sgeqrf_(m,n,mat,lda,tau,lwork,lworkSize,info);
#endif
}
template <typename T>
void TamaUNGQR(int *m, int *n, int * kk,  T *mat, int *lda, T * tau, T * lwork, int *lworkSize, int *info){

}

template <>
void TamaUNGQR<__lpkdoublecomplex>(int *m, int *n, int * kk, __lpkdoublecomplex *mat, int *lda, __lpkdoublecomplex * tau, __lpkdoublecomplex * lwork, int *lworkSize, int *info){
#ifdef CUDA_ENABLED

    culaDeviceZungqr(*m,*n,*kk,(culaDeviceDoubleComplex *) mat, *lda, (culaDeviceDoubleComplex *) tau);
    *info = culaGetErrorInfo();

#elif defined MKL_ENABLED

	*info = LAPACKE_zungqr(LAPACK_COL_MAJOR,*m,*n,*kk,mat,*lda,tau);
#else

	zungqr_(m,n,kk,mat,lda,tau,lwork,lworkSize,info);

#endif
}

template <>
void TamaUNGQR<__lpkcomplex>(int *m, int *n, int * kk, __lpkcomplex *mat, int *lda, __lpkcomplex * tau, __lpkcomplex * lwork, int *lworkSize, int *info){
#ifdef CUDA_ENABLED

    culaDeviceCungqr(*m,*n,*kk,(culaDeviceFloatComplex *) mat, *lda, (culaDeviceFloatComplex *) tau);
    *info = culaGetErrorInfo();


#elif defined MKL_ENABLED
	*info = LAPACKE_cungqr(LAPACK_COL_MAJOR,*m,*n,*kk,mat,*lda,tau);
#else
	cungqr_(m,n,kk,mat,lda,tau,lwork,lworkSize,info);
#endif
}


template <>
void TamaUNGQR<__lpkdoublereal>(int *m, int *n, int * kk, __lpkdoublereal *mat, int *lda, __lpkdoublereal * tau, __lpkdoublereal * lwork, int *lworkSize, int *info){
#ifdef CUDA_ENABLED

    culaDeviceDorgqr(*m,*n,*kk, mat, *lda, tau);
    *info = culaGetErrorInfo();

#elif defined MKL_ENABLED
	*info = LAPACKE_dorgqr(LAPACK_COL_MAJOR,*m,*n,*kk,mat,*lda,tau);
	if(*info !=0)
		printf("\nProblem in dorgqr\n");
#else
	dorgqr_(m,n,kk,mat,lda,tau,lwork,lworkSize,info);
#endif
}

template <>
void TamaUNGQR<__lpkreal>(int *m, int *n, int * kk, __lpkreal *mat, int *lda, __lpkreal * tau, __lpkreal * lwork, int *lworkSize, int *info){
#ifdef CUDA_ENABLED

    culaDeviceSorgqr(*m,*n,*kk, mat, *lda, tau);
    *info = culaGetErrorInfo();

#elif defined MKL_ENABLED

	*info = LAPACKE_sorgqr(LAPACK_COL_MAJOR,*m,*n,*kk,mat,*lda,tau);

#else

	sorgqr_(m,n,kk,mat,lda,tau,lwork,lworkSize,info);

#endif
}
/************************************************************************/

/************************************************************************/
/*GESVD is another beast! The real version requires only the lwork array. The complex one requires a rwork array as well. */
template <typename T>
void *TamaGESVD(int * nrA, int * ncA, T * rawA,int * lda, void * matrS, T * matrU, int * ldu, T* matrVT,int * ldvt, int *info){
	//Empty body
}

template <>
void *  TamaGESVD<__lpkdoublecomplex>(int * nrA, int * ncA, __lpkdoublecomplex * rawA,int * lda, void * matrS, __lpkdoublecomplex * matrU, int * ldu, __lpkdoublecomplex* matrVT,int * ldvt, int *info){

	char jobu = 'S';
	char jobvt = 'S';

#ifdef CUDA_ENABLED
    __lpkdoublereal *arrayRWork=NULL;
    culaDeviceZgesvd(jobu,jobvt,*nrA, *ncA,(culaDeviceDoubleComplex *) rawA,*lda,(__lpkdoublereal *) matrS, (culaDeviceDoubleComplex *) matrU, *ldu,(culaDeviceDoubleComplex *) matrVT,*ldvt);
    *info = culaGetErrorInfo();
    //If we have an error on exit we set the arrayRWork value returned by the function to 1;
    if(*info !=0){
        arrayRWork = (__lpkdoublereal *)1;
    }
    //else it remains NULL;
    
#elif defined  MKL_ENABLED

	__lpkdoublereal *arrayRWork = __lpkdoublereal[*ncA];
	*info = LAPACKE_zgesvd(LAPACK_COL_MAJOR,jobu,jobvt,*nrA,*ncA,rawA,*lda, (__lpkdoublereal *) matrS,matrU,*ldu,matrVT,*ldvt,arrayRWork);
	if((*info)!=0){
		printf("\nProblem in MKL SVD setup!\n");
		return NULL;
	}
#else


	 integer rworkDim = 5* (*ldvt); //Set following the documentation
	 __lpkdoublereal *arrayRWork =  new __lpkdoublereal[rworkDim];


	integer lworkDim = -1; //
	__lpkdoublecomplex *arrayLWork = new __lpkdoublecomplex[10];

	//For complexes we allocate the arrayRwork with dimensions fixed following the documentation

	 //Workspace query


	zgesvd_(&jobu,&jobvt,nrA,ncA,rawA,lda, (__lpkdoublereal *) matrS,matrU,ldu,matrVT,ldvt,arrayLWork,&lworkDim, arrayRWork,info);

	//arrayLWork contains the right size

	if((*info)!=0){
		printf("\nProblem in SVD setup!\n");
		return NULL;
	}

	//Set the WORK array size
	lworkDim = (integer) *((__lpkdoublereal *)arrayLWork);

	//Remove the old array
	delete [] arrayLWork;

	//Allocate the new array
	arrayLWork = new __lpkdoublecomplex[lworkDim];

	zgesvd_(&jobu,&jobvt,nrA,ncA,rawA,lda, (__lpkdoublereal *) matrS,matrU,ldu,matrVT,ldvt,arrayLWork,&lworkDim, arrayRWork,info);

	delete []arrayLWork;

#endif

	return (void *) arrayRWork;

}

template <>
void *  TamaGESVD<__lpkcomplex>(int * nrA, int * ncA, __lpkcomplex * rawA,int * lda, void * matrS, __lpkcomplex * matrU, int * ldu, __lpkcomplex* matrVT,int * ldvt, int *info){


	char jobu = 'S';
	char jobvt = 'S';
#ifdef CUDA_ENABLED
    __lpkreal *arrayRWork=NULL;
    culaDeviceCgesvd(jobu,jobvt,*nrA, *ncA,(culaDeviceFloatComplex *) rawA,*lda,(__lpkreal *) matrS, (culaDeviceFloatComplex *) matrU, *ldu,(culaDeviceFloatComplex *) matrVT,*ldvt);
    *info = culaGetErrorInfo();
    //If we have an error on exit we set the arrayRWork value returned by the function to 1;
    if(*info !=0){
        arrayRWork = (__lpkreal *)1;
    }
    //else it remains NULL;
    
#elif defined  MKL_ENABLED

	__lpkreal *arrayRWork = new __lpkreal[*ncA];
	*info = LAPACKE_cgesvd(LAPACK_COL_MAJOR,jobu,jobvt,*nrA,*ncA,rawA,*lda, (__lpkreal *) matrS,matrU,*ldu,matrVT,*ldvt,arrayRWork);
	if((*info)!=0){
		printf("\nProblem in MKL SVD setup!\n");
		return NULL;
	}


#else


	integer lworkDim = -1; //
	__lpkcomplex *arrayLWork = new __lpkcomplex[10];

	//For complexes we allocate the arrayRwork with dimensions fixed following the documentation
	integer rworkDim = 5* (*ldvt); //Set following the documentation
	 __lpkreal *arrayRWork =  new __lpkreal[rworkDim];

//	 char jobu = 'S';
//	 char jobvt = 'S';
//	 //Workspace query

	cgesvd_(&jobu,&jobvt,nrA,ncA,rawA,lda, (__lpkreal *) matrS,matrU,ldu,matrVT,ldvt,arrayLWork,&lworkDim, arrayRWork,info);
	//arrayLWork contains the right size

	if((*info)!=0){
		printf("\nProblem in SVD setup!\n");
		return NULL;
	}

	//Set the WORK array size
	lworkDim = (integer) *((__lpkreal *)arrayLWork);

	//Remove the old array
	delete [] arrayLWork;

	//Allocate the new array
	arrayLWork = new __lpkcomplex[lworkDim];

	cgesvd_(&jobu,&jobvt,nrA,ncA,rawA,lda, (__lpkreal *) matrS,matrU,ldu,matrVT,ldvt,arrayLWork,&lworkDim, arrayRWork,info);

	delete []arrayLWork;

#endif

	return (void *) arrayRWork;

}

template <>
void *  TamaGESVD<__lpkdoublereal>(int * nrA, int * ncA, __lpkdoublereal * rawA,int * lda, void * matrS, __lpkdoublereal * matrU, int * ldu, __lpkdoublereal* matrVT,int * ldvt, int *info){


	char jobu = 'S';
	char jobvt = 'S';

 #ifdef CUDA_ENABLED
    __lpkdoublereal *arrayLWork=NULL;
    culaDeviceDgesvd(jobu,jobvt,*nrA, *ncA, rawA,*lda,(__lpkdoublereal *) matrS, matrU, *ldu, matrVT,*ldvt);
    *info = culaGetErrorInfo();
    //If we have an error on exit we set the arrayRWork value returned by the function to 1;
    if(*info !=0){
        printf("\nProblem in culaDeviceDgesvd!!!\n");
        arrayLWork = (__lpkdoublereal *)1;
    }
    //else it remains NULL;
    
#elif defined  MKL_ENABLED

	__lpkdoublereal *arrayLWork = new __lpkdoublereal[*ncA];
	*info = LAPACKE_dgesvd(LAPACK_COL_MAJOR,jobu,jobvt,*nrA,*ncA,rawA,*lda, (__lpkdoublereal *) matrS,matrU,*ldu,matrVT,*ldvt,arrayLWork);
	if((*info)!=0){
		printf("\nProblem in MKL SVD setup!\n");
		return NULL;
	}
#else

	integer lworkDim = -1; //
	__lpkdoublereal *arrayLWork = new __lpkdoublereal[10];


//	 char jobu = 'S';
//	 char jobvt = 'S';
	 //Workspace query

	dgesvd_(&jobu,&jobvt,nrA,ncA,rawA,lda, (__lpkdoublereal *) matrS,matrU,ldu,matrVT,ldvt,arrayLWork,&lworkDim,info);
	//arrayLWork contains the right size

	if((*info)!=0){
		printf("\nProblem in SVD setup!\n");
		return NULL;
	}

	//Set the WORK array size
	lworkDim = (integer) *((__lpkdoublereal *)arrayLWork);

	//Remove the old array
	delete [] arrayLWork;

	//Allocate the new array
	arrayLWork = new __lpkdoublereal[lworkDim];

	dgesvd_(&jobu,&jobvt,nrA,ncA,rawA,lda, (__lpkdoublereal *) matrS,matrU,ldu,matrVT,ldvt,arrayLWork,&lworkDim,info);


#endif

	return (void *) arrayLWork;

}

template <>
void *  TamaGESVD<__lpkreal>(int * nrA, int * ncA, __lpkreal * rawA,int * lda, void * matrS, __lpkreal * matrU, int * ldu, __lpkreal* matrVT,int * ldvt, int *info){

	char jobu = 'S';
	char jobvt = 'S';
 #ifdef CUDA_ENABLED
    __lpkreal *arrayLWork=NULL;
    culaDeviceSgesvd(jobu,jobvt,*nrA, *ncA, rawA,*lda,(__lpkreal *) matrS, matrU, *ldu, matrVT,*ldvt);
    *info = culaGetErrorInfo();
    //If we have an error on exit we set the arrayRWork value returned by the function to 1;
    if(*info !=0){
        arrayLWork = (__lpkreal *)1;
    }
    //else it remains NULL;
    
#elif defined  MKL_ENABLED

	__lpkreal *arrayLWork = new __lpkreal[*ncA];
	*info = LAPACKE_sgesvd(LAPACK_COL_MAJOR,jobu,jobvt,*nrA,*ncA,rawA,*lda, (__lpkreal *) matrS,matrU,*ldu,matrVT,*ldvt,arrayLWork);
	if((*info)!=0){
		printf("\nProblem in MKL SVD setup!\n");
		return NULL;
	}
#else

	integer lworkDim = -1; //
	__lpkreal *arrayLWork = new __lpkreal[10];


//	 char jobu = 'S';
//	 char jobvt = 'S';
	 //Workspace query

	sgesvd_(&jobu,&jobvt,nrA,ncA,rawA,lda, (__lpkreal *) matrS,matrU,ldu,matrVT,ldvt,arrayLWork,&lworkDim,info);
	//arrayLWork contains the right size

	if((*info)!=0){
		printf("\nProblem in SVD setup!\n");
		return NULL;
	}

	//Set the WORK array size
	lworkDim = (integer) *((__lpkreal *)arrayLWork);

	//Remove the old array
	delete [] arrayLWork;

	//Allocate the new array
	arrayLWork = new __lpkreal[lworkDim];

	sgesvd_(&jobu,&jobvt,nrA,ncA,rawA,lda, (__lpkreal *) matrS,matrU,ldu,matrVT,ldvt,arrayLWork,&lworkDim,info);

#endif

	return (void *) arrayLWork;

}

template <typename T>
__lpkdoublereal TamaNORM2REAL(integer dim, T *data){
	//General function not implemented
}

template <>
__lpkdoublereal TamaNORM2REAL<__lpkdoublecomplex>(integer dim, __lpkdoublecomplex * data){
#ifdef CUDA_ENABLED
    
   return cublasDnrm2(2*dim,(__lpkdoublereal *)data,1);

#else
	return cblas_dnrm2(2*dim, (__lpkdoublereal *) data,1);
#endif
    //appo2=NORM2REAL(dimRatio*(AO.nrows),((lpkreal *)(AO.raw)) + dimRatio*i*(AO.nrows),1);
}

template <>
__lpkdoublereal TamaNORM2REAL<__lpkcomplex>(integer dim, __lpkcomplex * data){
#ifdef CUDA_ENABLED
    
   return cublasSnrm2(2*dim,(__lpkreal *)data,1);

#else

	return cblas_snrm2(2*dim, (__lpkreal *) data,1);

#endif
	//appo2=NORM2REAL(dimRatio*(AO.nrows),((lpkreal *)(AO.raw)) + dimRatio*i*(AO.nrows),1);
}

template <>
__lpkdoublereal TamaNORM2REAL<__lpkdoublereal>(integer dim, __lpkdoublereal * data){

#ifdef CUDA_ENABLED
    
   return cublasDnrm2(dim,data,1);

#else

	return cblas_dnrm2(dim, data,1);
#endif

	//appo2=NORM2REAL(dimRatio*(AO.nrows),((lpkreal *)(AO.raw)) + dimRatio*i*(AO.nrows),1);
}

template <>
__lpkdoublereal TamaNORM2REAL<__lpkreal>(integer dim, __lpkreal * data){

#ifdef CUDA_ENABLED
    return cublasSnrm2(dim,data,1);
#else
	return cblas_snrm2(dim,data,1);
#endif
	//appo2=NORM2REAL(dimRatio*(AO.nrows),((lpkreal *)(AO.raw)) + dimRatio*i*(AO.nrows),1);
}
//appo2=NORM2REAL(dimRatio*(AO.nrows),((lpkreal *)(AO.raw)) + dimRatio*i*(AO.nrows),1);





template <typename T>
QRfactLaPack<T> * reorthogonalizeLaPack( tamaClss::genMatrix<T>  &matObj,char r){

	integer m = matObj.nrows;
	integer n = matObj.ncols;

	constants<T> costanti;

	integer kk = min(m,n);
#ifdef CUDA_ENABLED
    T *tau;
    cudaMalloc((void **)&tau,kk*sizeof(T));
#else
	T *tau = new T[kk];
#endif
	T * mat = matObj.raw;

	//doublecomplex * mat = new doublecomplex[m*n];
	//
	//	//Copy the original matrix
	//	cblas_zcopy((const int) n*m, (void *) (matObj.raw),1,(void *)mat,1);

	integer info;

	integer lworkSize = -1;

	T *lwork=new T[1];
#if defined(CUDA_ENABLED) ||defined(MKL_ENABLED)

//We need to call the GEQRF and GEUNGQR functions only once
	TamaGEQRF(&m,&n,mat,&m,tau,lwork,&lworkSize,&info);
	if(info!=0){
		printf("\nproblem with ??geqrf!\n");
		return NULL;
	}

	TamaUNGQR(&m,&n,&kk,mat,&m,tau,lwork,&lworkSize,&info);

#else

	//tamaClss::genMatrix<T> *matR=NULL;

	//Query for workspace

	//GEQRF_(&m,&n,mat,&m,tau,lwork,&lworkSize,&info);
	TamaGEQRF(&m,&n,mat,&m,tau,lwork,&lworkSize,&info);
	if(info!=0){
		printf("\nproblem with zgeqrf!\n");
		return NULL;
	}

	if (costanti.typeCheck == 'z'){
		lworkSize =(integer )  *((__lpkdoublereal *) lwork);
	}
	else if (costanti.typeCheck == 'c'){
		lworkSize =(integer )  *((__lpkreal *) lwork);
	}
	else if (costanti.typeCheck == 'd'){
		lworkSize =(integer )  *((__lpkdoublereal *) lwork);

	}
	else //costanti.typeCheck == 's'
		lworkSize =(integer )  *((__lpkreal *) lwork);
	delete [] lwork;


	lwork = new T[lworkSize];
	//Compute the first step of the QR decomposition
	TamaGEQRF(&m,&n,mat,&m,tau,lwork,&lworkSize,&info);

	//__OBSERVATION
	//This part of the function is still specialized on complexes. We NEVER use it, so I do not touch it.
	//For the time being, I silence this part of the code. In case it becomes necessary to have the matrix R
	//remind to modify it
	if(r=='r'){//Create the object matR and assign its adress to *rmatrix;

		printf("Option 'r' not enabled. I do not do anything. Take care of unallocated memory");
		//Since we will call this very rarely (and mainly fo r debugging purposes) I do not devote to much attention to optimization
//		matR = new tamaClss::genMatrix<T>(min(m,n),n);
//		for(int j=0;j<matR->nrows;j++){//Here we copy the upper triangular part
//			for(int i=0; i<=j;i++){
//				matR->raw[j*(matR->nrows)+i].r= mat[j*m+i].r;
//				matR->raw[j*(matR->nrows)+i].i= mat[j*m+i].i;
//			}
//			for(int i=j+1; i < matR->nrows;i++){
//				matR->raw[j*(matR->nrows)+i].r= 0.;
//				matR->raw[j*(matR->nrows)+i].i= 0.;
//			}
//		}
//		for(int j=matR->nrows;j<matR->ncols;j++){//Here we copy the rectangular part
//			for(int i=0; i<=j;i++){
//				matR->raw[j*(matR->nrows)+i].r= mat[j*m+i].r;
//				matR->raw[j*(matR->nrows)+i].i= mat[j*m+i].i;
//			}
//		}

	}

	//Cleaning up
	delete [] lwork;
	lworkSize = -1;
	lwork = new T[10];

	//Construct the actual orthogonal matrix

	//1st: query for optimal lwork dim

	TamaUNGQR(&m,&n,&kk,mat,&m,tau,lwork,&lworkSize,&info);

	if (costanti.typeCheck == 'z'){
			lworkSize =(integer )  *((__lpkdoublereal *) lwork);
		}
		else if (costanti.typeCheck == 'c'){
			lworkSize =(integer )  *((__lpkreal *) lwork);
		}
		else if (costanti.typeCheck == 'd'){
			lworkSize =(integer )  *((__lpkdoublereal *) lwork);

		}
		else //costanti.typeCheck == 's'
			lworkSize =(integer )  *((__lpkreal *) lwork);

	delete [] lwork;

	lwork = new T[lworkSize];

	//2nd: construct the matrix
	TamaUNGQR(&m,&n,&kk,mat,&m,tau,lwork,&lworkSize,&info);



//	QRfactLaPack<T> *result = new QRfactLaPack<T>(m,n);
//	result -> mQ = &matObj;
	//result -> tau = tau;
//	if (r=='r'){
//		result -> mR = matR;
//	}
#endif

	QRfactLaPack<T> *result = new QRfactLaPack<T>(m,n);
	result -> mQ = &matObj;

	delete [] lwork;
#ifdef CUDA_ENABLED
    cudaFree(tau);
#else
	delete [] tau;
#endif
	return result;

}

//SUCCESS: NULL returned;
//FAILURE (no convergence): address of the RWORK array returned




/************************************************************************/
template <typename T>
QRfactLaPack<T> * onlineCheck(tamaClss::genMatrix<T> & matA, tamaClss::genMatrix<T> &Q, __lpkdoublereal tolerance, integer incr){

	//__DEBUG
//	printf("\nSono in onlineCheck...\n");

	QRfactLaPack<T> * newQ = new QRfactLaPack<T>(Q.nrows,Q.ncols);
	newQ->mQ = &Q;

	constants<T> costanti;

   // tamaClss::genMatrix<T> Qappo(Q.nrows,Q.ncols);
   // cudaMemcpy(Qappo.raw,Q.raw,Q.nrows * Q.ncols *sizeof(T), cudaMemcpyDeviceToHost);
   // Qappo.toBinFile("Qappo.dat");

	//__DEBUG
	// Q.toBinFile("provaOrthoPrima");

	//T norms[incr];
#ifdef CUDA_ENABLED
    curandGenerator_t cuGen;
    curandCreateGenerator(&cuGen,CURAND_RNG_PSEUDO_MTGP32)/*MT19937 if cc>3.2*/;
    curandSetPseudoRandomGeneratorSeed(cuGen,time(NULL));

#else
	mt19937_64 generator;
	//Set rndGen seed
	generator.seed(time(NULL));
	//Normal distribution
	//__ATTENTION: minimal waste of space and time
	normal_distribution<__lpkdoublereal> norm_distribution_double(0.,1.);
	normal_distribution<__lpkreal> norm_distribution_single(0.,1.);
#endif

	__lpkdoublereal max;

	//__OBSERVATION
	//The norm of a n-component complex vector can be computed as the norm of a real vector with 2*n components


	do{
#ifdef CUDA_ENABLED
		//Prepare the additional columns Omega matrix
		tamaClss::genMatrix<T> addOmega;
        addOmega.nrows = matA.ncols;
        addOmega.ncols = incr;
        cudaMalloc((void **)&(addOmega.raw),addOmega.nrows *addOmega.ncols *sizeof(T));

		//A.addOmega
		tamaClss::genMatrix<T> AO;
        AO.nrows = matA.nrows;
        AO.ncols = addOmega.ncols;
        cudaMalloc((void **)&(AO.raw),AO.nrows*AO.ncols*sizeof(T));
		//QdagaO
		tamaClss::genMatrix<T> QdagaAO;
        QdagaAO.nrows = Q.ncols;
        QdagaAO.ncols = incr;
        cudaMalloc((void **) &(QdagaAO.raw),QdagaAO.nrows*QdagaAO.ncols*sizeof(T));
#else
		//Prepare the additional columns Omega matrix
		tamaClss::genMatrix<T> addOmega(matA.ncols,incr);
		//A.addOmega
		tamaClss::genMatrix<T> AO(matA.nrows,addOmega.ncols);
		//QdagaO
		tamaClss::genMatrix<T> QdagaAO(Q.ncols,incr);
#endif
		//Generate the additional columns
		//__OBSERVATION
		//We treat always the raw vector as a vector of reals.
		//The number of components depends on the real type: n if real, 2n if complex

		if (costanti.typeCheck == 'z'){
			//printf("\nTipo complesso doppia precisione\n");
#ifdef CUDA_ENABLED
            if(CURAND_STATUS_SUCCESS != curandGenerateNormalDouble(cuGen,(__lpkdoublereal *) addOmega.raw, 2* addOmega.nrows * addOmega.ncols,0.,1.)){
                    
            printf("\nCuRand addOmega problem\n.");        
            }
            //__DEBUG
           // tamaClss::genMatrix<__lpkdoublecomplex> addOmegaProva(addOmega.nrows,addOmega.ncols);
           // cudaMemcpy(addOmegaProva.raw,addOmega.raw, addOmega.nrows *addOmega.ncols * sizeof(__lpkdoublecomplex),cudaMemcpyDeviceToHost);
           // addOmegaProva.toBinFile("addOmegaProva.dat");

#else
			for(int i=0; i<2 * addOmega.nrows * addOmega.ncols ;i++){
				*(((__lpkdoublereal *) (addOmega.raw)) +i)=norm_distribution_double(generator);
			}
#endif
		}
		else if (costanti.typeCheck == 'c'){
#ifdef CUDA_ENABLED
            curandGenerateNormal(cuGen,(__lpkreal *) addOmega.raw, 2* addOmega.nrows * addOmega.ncols,0.,1.);
#else

			for(int i=0; i<2 * addOmega.nrows * addOmega.ncols ;i++){
				*(((__lpkreal *) (addOmega.raw)) +i)=norm_distribution_single(generator);
			}
#endif
		}
		else if(costanti.typeCheck == 'd'){

#ifdef CUDA_ENABLED
            curandGenerateNormalDouble(cuGen,(__lpkdoublereal *) addOmega.raw, addOmega.nrows * addOmega.ncols,0.,1.);
#else

			//	printf("\ndouble real type detected!\n");
			for(int i=0; i< addOmega.nrows * addOmega.ncols ;i++){
				*(( (__lpkdoublereal *)(addOmega.raw)) +i)=norm_distribution_double(generator);
			}
#endif
		}
		else{
#ifdef CUDA_ENABLED
            curandGenerateNormal(cuGen,(__lpkreal *) addOmega.raw, addOmega.nrows * addOmega.ncols,0.,1.);
#else

			for(int i=0; i< addOmega.nrows * addOmega.ncols ;i++){
				*(((__lpkreal *) (addOmega.raw)) +i)=norm_distribution_single(generator);
			}
#endif
		}




		//addOmega.toBinFile("addOmega.dat");

		//Compute A.addOmega [m x incr]
		//		printf("\nPrima moltiplicazione...\n");

		TamaGEMM('N','N',matA,addOmega,AO,costanti.typed_one,costanti.typed_zero);
        //tamaClss::genMatrix<T> AOappo(AO.nrows,AO.ncols);
        //cudaMemcpy(AOappo.raw,AO.raw, AO.nrows *AO.ncols *sizeof(T),cudaMemcpyDeviceToHost);
        //AOappo.toBinFile("AOappo.dat");
		//We save the AO matrix for future use
		//AO.toBinFile("AO1.dat");
#ifdef CUDA_ENABLED
        tamaClss::genMatrix<T>AO2;
        AO2.nrows = AO.nrows;
        AO2.ncols = AO.ncols;
        cudaMalloc((void **) &(AO2.raw),AO2.nrows * AO2.ncols *sizeof(T));
        cudaMemcpy(AO2.raw,AO.raw,AO.nrows*AO.ncols*sizeof(T),cudaMemcpyDeviceToDevice);

#else
		tamaClss::genMatrix<T> AO2 = AO;
#endif
        ////		{int pippo;
		////		cin >> pippo;
		////		}
		//		//Compute QdagaAO
		////		printf("\nSeconda moltiplicazione...\n");
		TamaGEMM('T', 'N', Q, AO, QdagaAO, costanti.typed_mone,costanti.typed_zero);
        //tamaClss::genMatrix<T> QAO(QdagaAO.nrows,QdagaAO.ncols);
        //cudaMemcpy(QAO.raw,QdagaAO.raw,QAO.nrows*QAO.ncols*sizeof(T),cudaMemcpyDeviceToHost);
        //QAO.toBinFile("QAO.dat");
		//QdagaAO.toBinFile("QdagaAO.dat");
		////		{int pippo;
		////		cin >> pippo;
		////		}
		//
		//		//Compute QQdagaO
		////		printf("\nTerza moltiplicazione...\n");
		TamaGEMM('N','N',Q,QdagaAO,AO,costanti.typed_one,costanti.typed_one);
		////		{int pippo;
		////		cin >> pippo;
		////		}
		//
		//__DEBUG
		//AO.toBinFile("matAO.dat");
		//		//Now we have AO = A.addOmega and QQdagaO. We must compute the difference: AO-QQdagaO.
		//		//!!!AO is lost
		//		//!!!We operate on the whole raw vector of both matrices!
		//		//AXPY(AO.nrows* AO.ncols, &typed_mone,AO.raw,1,QQdagaO.raw,1);
		//
		max = costanti.real_zero;
//		T appo;
		//__ATTENTION: the comparison term is double precision
		__lpkdoublereal appo2;
		//
		//		//Find the maximal norm of the columns of AO - QQdagaO
		for(int i=0; i<incr; i++){
            appo2=TamaNORM2REAL((AO.nrows),((T *)(AO.raw)) + i*(AO.nrows));
			//printf("\nappo = %.20f\n",appo2);
			if(appo2>max)
				max = appo2;

		}
		//__DEBUG
//		printf("max=%f, Q.ncols = %d",max,Q.ncols);
		//
		if ((max > tolerance) and (Q.ncols <( matA.ncols - incr)) ){
			//We have to include the new omegas in the basis
			//...but only when the number of columns of Q+incr is smaller than the number of columns of A
			//			//Allocate the space for the new Q
			T * tempraw;
            T* pappo;    
#ifdef CUDA_ENABLED

            cudaMalloc((void **) &tempraw, Q.nrows*(Q.ncols+incr)*sizeof(T));
            cudaMemcpy(tempraw,Q.raw,Q.nrows*Q.ncols*sizeof(T),cudaMemcpyDeviceToDevice);
            cudaMemcpy(tempraw+(Q.nrows*Q.ncols),AO2.raw,AO2.nrows*AO2.ncols*sizeof(T),cudaMemcpyDeviceToDevice);
#else
            tempraw= new T[Q.nrows*(Q.ncols + incr)];
            //REPLACE WITH memcpy
			memcpy(tempraw,Q.raw,Q.nrows * Q.ncols * sizeof(T));
			memcpy(tempraw+(Q.nrows*Q.ncols),AO2.raw, AO2.nrows*AO2.ncols*sizeof(T));


#endif

			//__DEBUG
		//	printf("\ntolleranza superata\n");
			//Copy the already existing elements of Q into the tempraw array.

			//point the the old Q.raw
			pappo=Q.raw;
			//set the new Q.raw
			Q.raw = tempraw;
#ifdef CUDA_ENABLED
            cudaFree(pappo);
#else	
            //Delete the old Q.raw
			delete [] pappo;
#endif
			//Update the number of columns
			Q.ncols+=incr;

			//__DEBUG
//			Q.toBinFile("newQno.dat");

			//			if (newQ != NULL)
			//				delete newQ;
			newQ = reorthogonalizeLaPack(Q,'n');
			newQ->nrows = Q.nrows;
			newQ->ncols = Q.ncols;
			//Update Q
			//newQ->mQ->toBinFile("provaOrtho.dat");
			//
		}
		//
		//__DEBUG
		printf("\nmax = %.15f, tolerance = %.15f, Q.ncols = %d matA.ncols-incr = %d\n",max,tolerance,Q.ncols,matA.ncols-incr);
		//newQ->mQ->toBinFile("provaOrtho.dat");

#ifdef CUDA_ENABLED
        //Avoid memory errors on exit when the genMatrix destructor is invoked
        cudaFree(addOmega.raw);
        addOmega.raw=NULL;
        cudaFree(AO.raw);
        AO.raw = NULL;
        cudaFree(QdagaAO.raw);
        QdagaAO.raw = NULL;
        cudaFree(AO2.raw);
        AO2.raw=NULL;
#endif
	}while((max>tolerance) and (Q.ncols <( matA.ncols - incr)));
	//
	printf("\n max= %.15f esco da onlineCheck.\n",max);
	/*****************************/


	return newQ;


}

template <typename T>
QRfactLaPack<T> * RandSubIter( tamaClss::genMatrix<T>  &matA,tamaClss::genMatrix<T>& matO, integer niter,__lpkdoublereal tolerance){

	constants<T> costanti;
	
	assert(matA.ncols == matO.nrows );
    
//__CUDA: here the switch decides where to allocate the matrices
#ifdef CUDA_ENABLED

    //Do not allocate the memory via constructor
    tamaClss::genMatrix<T> *matrY = new tamaClss::genMatrix<T>;

    //Set sizes
    matrY->nrows = matA.nrows;
    matrY->ncols = matO.ncols;
    if(cudaSuccess!= cudaMalloc((void **) &(matrY->raw),matrY->nrows *matrY->ncols*sizeof(T))){
        cout << "\nError in allocation of Ydev\n";
    }
   
#else

    tamaClss::genMatrix<T> *matrY = new tamaClss::genMatrix<T>(matA.nrows,matO.ncols);

#endif
    //tamaClss::genMatrix<T> *matrYtilde = new tamaClss::genMatrix<T>(matO.nrows,matO.ncols);


	//Y=A.O;
	//printf("\nFirst product\n");
	//GEMM(CblasColMajor,CblasNoTrans,CblasNoTrans,matA.nrows,matO.ncols,matO.nrows,/*(const void *)*/&typed_one,matA.raw,matA.nrows,matO.raw,matO.nrows,/*(const void *)*/&typed_zero,matrY->raw,matrY->nrows);
	TamaGEMM('N','N',matA, matO, *matrY, costanti.typed_one, costanti.typed_zero);
   // tamaClss::genMatrix<T> Yprova(matrY->nrows, matrY->ncols);
   // cudaMemcpy(Yprova.raw,matrY->raw,Yprova.nrows*Yprova.ncols*sizeof(T),cudaMemcpyDeviceToHost);
   // Yprova.toBinFile("Yprova.dat");


   	tamaClss::genMatrix<T> *Q;
	tamaClss::genMatrix<T> *Qtilde;

	QRfactLaPack<T> *appoQtildeQR=NULL;
	QRfactLaPack<T> *appoQQR2 = reorthogonalizeLaPack( *matrY,'n');
     //__DEBUG
   // tamaClss::genMatrix<T> Yprovadopo(matrY->nrows, matrY->ncols);
   // cudaMemcpy(Yprovadopo.raw,appoQQR2->mQ->raw,Yprovadopo.nrows*Yprovadopo.ncols*sizeof(T),cudaMemcpyDeviceToHost);
   // Yprovadopo.toBinFile("Yprovadopo.dat");


	//matrY->toBinFile("prima.dat");

	QRfactLaPack<T> *appoQQR = onlineCheck(matA, *(appoQQR2->mQ),tolerance, matO.ncols);
	//delete appoQQR2;
	//appoQQR2=NULL;


	Q = appoQQR->mQ;

#ifdef CUDA_ENABLED

	tamaClss::genMatrix<T> *matrYtilde = new tamaClss::genMatrix<T>;
    matrYtilde->nrows = matA.ncols;
    matrYtilde->ncols = Q->ncols;
    cudaMalloc((void **)&(matrYtilde->raw),matrYtilde->nrows *matrYtilde->ncols*sizeof(T));

#else

	tamaClss::genMatrix<T> *matrYtilde = new tamaClss::genMatrix<T>(matA.ncols,Q->ncols);

#endif

	//We have Y=A.Omega reorthogonalized.


	//Iteration cycle
	for(int j=1;j<=niter;j++){
		TamaGEMM('T','N',matA,*Q,*matrYtilde,costanti.typed_one,costanti.typed_zero);
		//printf("\nSecond product\n");
		delete appoQtildeQR;
		//Ytilde = Qtilde.Rtilde(no Rtilde produced...).
		appoQtildeQR = reorthogonalizeLaPack(*matrYtilde,'n');
		Qtilde = appoQtildeQR->mQ;
		//Replace old Y with Y = A Qtilde
		//printf("\nThird product\n");
		TamaGEMM('N','N',matA,*Qtilde,*matrY,costanti.typed_one,costanti.typed_zero);
		//printf("\nThird product\n");
		//Delete the QR decomposition of Y_{j-1}

		//__DEBUG
		//printf("\nSto per cancellare\n");
		delete appoQQR;
		//printf("\ncancellato!\n");
		//{int check; scanf("%d",&check);}
		appoQQR = reorthogonalizeLaPack( *matrY,'n');
		//{int check; scanf("%d",&check);}
		Q = appoQQR->mQ;
	}
 
#ifdef CUDA_ENABLED
 cudaFree(matrYtilde->raw);
 matrYtilde->raw = NULL;
#endif    

    
    delete matrYtilde;
	
    return appoQQR;

}



/************************************************************************/
template <typename T>
void * SVD(tamaClss::genMatrix<T> &A, tamaClss::genMatrix<T> **U, void **S, tamaClss::genMatrix<T> **VT){
	//The matrices U,S and VT are created here and their addresses are assigned to the pointers ginven in input.
	//__ATTENTION: the code is meant to work with JOBU ==  JOBVT  == 'S'

//	char jobu = 'S';
//	char jobvt = 'S';

	integer leadingDimU = A.nrows;
	integer leadingDimVT = (A.nrows>=A.ncols?A.ncols:A.nrows); //min(M,N)

	constants<T> costanti;


	integer info;

	printf("\n%d, %d\n",leadingDimU,leadingDimVT);

#ifdef CUDA_ENABLED
//CUDA CODE
    T *matrU;
    T *matrVT;
    //T* matrA;
    //T* matrAappo;

    cudaMalloc((void **) &matrU, leadingDimU*leadingDimVT*sizeof(T));

    cudaMalloc((void **)&matrVT, leadingDimVT*(A.ncols)*sizeof(T));

    //__ATTENTION: the following lines must be revised
    //cudaMalloc((void **)&matrA, A.nrows*A.ncols*sizeof(T));

    //Here we copy the matrix A onto device
    //cudaMemcpy(matrA, A.raw, A.nrows*A.ncols*sizeof(T),cudaMemcpyHostToDevice);

    //We store the address of A.raw into matrAappo
    //matrAappo = A.raw;

    //For compatibility with the remaining part of the code
    //set A.raw to the device memory address. We will clean the device memory and 
    //restore A.raw to the original value before leaving the routine.

    //A.raw = matrA;  

	void *matrS;
	if((costanti.typeCheck =='c') or (costanti.typeCheck == 's')){
        cudaMalloc((void **)&matrS, leadingDimVT*sizeof(__lpkreal));
        //	matrS = (void *) new __lpkreal [leadingDimVT];
	}
	else{
        cudaMalloc((void **)&matrS, leadingDimVT*sizeof(__lpkdoublereal));
		//matrS = (void *) new __lpkdoublereal [leadingDimVT];
	}


#else
    T *matrU = new T[leadingDimU*leadingDimVT];

	T *matrVT = new T[leadingDimVT*(A.ncols)];


	void *matrS;
	if((costanti.typeCheck =='c') or (costanti.typeCheck == 's')){
		matrS = (void *) new __lpkreal [leadingDimVT];
	}
	else{
		matrS = (void *) new __lpkdoublereal [leadingDimVT];
	}
#endif

	void * arrayLWork;

	//!!!
	//First call to zgesvd to determine the optimal size of the LWORK array
	//GESVD_(&jobu,&jobvt,&(A.nrows),&(A.ncols),A.raw,&(A.nrows),matrS,matrU,&leadingDimU,matrVT,&leadingDimVT,arrayLWork,&lworkDim,arrayRWork,&info);
    cout << "\nCall tamaGESVD\n";
arrayLWork = 	TamaGESVD(&(A.nrows),&(A.ncols),A.raw,&(A.nrows),matrS,matrU,&leadingDimU,matrVT,&leadingDimVT, &info);
cout <<"\n...done!\n";
  //  {
  //      int pippo;
  //      cin >> pippo;
  //  }
	//
	tamaClss::genMatrix<T> *UObj =new tamaClss::genMatrix<T>;
    UObj->nrows = leadingDimU;
    UObj->ncols = leadingDimVT;
	//delete [] UObj->raw;
	UObj->raw= matrU;


//	tamaClss::genMatrix<T> *SObj =new tamaClss::genMatrix<T>(leadingDimVT,leadingDimVT);
	//tamaClss::genMatrix<T> *SObj =new tamaClss::genMatrix<T>(1,leadingDimVT);
	//delete [] SObj->raw;
	//SObj->raw = matrSmatr;
	//SObj->raw= matrS;


	tamaClss::genMatrix<T> *VTObj =new tamaClss::genMatrix<T>;
    VTObj->nrows = leadingDimVT;
    VTObj->ncols = A.ncols;
	//delete [] VTObj->raw;
	VTObj->raw= matrVT;
//If the code is compiled for cuda, the matrices *U, *VT and *S point to device memory
	*U = UObj;
	//*S = SObj;
	*S = matrS;
	*VT = VTObj;

    //__ATTENTION: this must be removed when the whole code will run on the GPU
#ifdef CUDA_ENABLED
    //Free the memory allocated on the device for A
    //cudaFree(A.raw);
    //Point to the host memory A matrix.
    //A.raw = matrAappo;

#endif

	//delete [] matrS;
	//delete [] arrayLWork;

//	if (info == 0){
//		delete [] arrayRWork;
//		arrayRWork = NULL;
//	}
//
	if (info == 0){
		if ((costanti.typeCheck == 'z') or (costanti.typeCheck =='d')){
			__lpkdoublereal *appo = (__lpkdoublereal *) arrayLWork;
#ifdef CUDA_ENABLED
            cudaFree(appo);
#else
            delete  [] appo;
#endif
		}
		else{
			__lpkreal *appo = (__lpkreal *) arrayLWork;
#ifdef CUDA_ENABLED
            cudaFree(appo);
#else
			delete [] appo;
#endif
		}

		arrayLWork = NULL;
	}

	//return arrayRWork

	return arrayLWork;

}

template <typename T>
void RRSVD(tamaClss::genMatrix<T> &A,integer relevant, integer niter, __lpkdoublereal tolerance, tamaClss::genMatrix<T> **U, /*tamaClss::genMatrix<T>*/void  **S, tamaClss::genMatrix<T> **VT){

	
	integer dimRange = min(2*relevant, A.ncols);

	constants<T> costanti;

	printf("\nThe number of columns of A is %d.\n Number of retained dimensions:%d\n",A.ncols, dimRange);

	

	using namespace std::chrono;
	high_resolution_clock::time_point t1overall = high_resolution_clock::now();
	//	printf("\nThe number of columns of A is %d; choose a number < %d\n",ncolsA, ncolsA);
	//	scanf("%d",&dimRange);

    //We do not allocate the raw here.
    //The allocation depends on CUDA_ENABLED
	tamaClss::genMatrix<T> C2appo;
    C2appo.nrows = A.nrows;
    C2appo.ncols = A.ncols;

    //We do not allocate the raw here.
    //The allocation depends on CUDA_ENABLED
    tamaClss::genMatrix<T> omega;
    omega.nrows = A.ncols;
    omega.ncols = dimRange;
    printf("\nGenerating Omega matrix...");

#ifdef CUDA_ENABLED
    //Initialize the cublas environment
    //We must check whether it this initialization suffices for alla the calls....
    cublasInit();
//CUDA CODE
    //Define the generator
    curandGenerator_t cuGen;
    if(CURAND_STATUS_SUCCESS != curandCreateGenerator(&cuGen,CURAND_RNG_PSEUDO_MTGP32)/*MT19937*/){
        printf("\nProblema createGenerator!\n");
    }
    if(CURAND_STATUS_SUCCESS != curandSetPseudoRandomGeneratorSeed(cuGen,time(NULL))){
        printf("\nProblema setSeed\n");
    }

    if(cudaSuccess!=cudaMalloc((void **)&(C2appo.raw),C2appo.nrows*C2appo.ncols*sizeof(T))){

        printf("\nError in allocation of C2appo\n");
    }
    if(cudaSuccess!=cudaMemcpy(C2appo.raw,A.raw,C2appo.nrows*C2appo.ncols*sizeof(T),cudaMemcpyHostToDevice)){
        printf("\nError in cudaMemcpy C2appo\n");
    }

    //Declare and allocate the omega matrix on the device
    //tamaClss::genMatrix<T> omegaDev;
    //omegaDev.nrows = omega.nrows;
    //omegaDev.ncols = omega.ncols;

    if(cudaSuccess != cudaMalloc((void **)&(omega.raw),omega.nrows*omega.ncols*sizeof(T) ) ){
        cout << "\nError in allocation of omegaDev\n";
    }

    if (costanti.typeCheck == 'z'){
        printf("\nCUDA generates the matrix Omega\n");

        cudaError_t err;
        err = cudaGetLastError();
        //printf("\nRRSVD error check: %d\n",err);
       // {int pippo; cin>>pippo;}

         if(CURAND_STATUS_SUCCESS !=  curandGenerateNormalDouble(cuGen,(__lpkdoublereal *) (omega.raw), 2 * omega.nrows * omega.ncols,0.,1.)){

             cout << endl << "RRSVD:Problem in random number generation"<< endl;
         }

    }
    else if (costanti.typeCheck == 'c'){
         printf("\nCUDA generates the matrix Omega\n");

         if(CURAND_STATUS_SUCCESS != curandGenerateNormal(cuGen,(__lpkreal *)(omega.raw), 2 * omega.nrows * omega.ncols,0.,1.)){

             cout << endl << "Problem in random number generation"<< endl;
         }
    }
    else if(costanti.typeCheck == 'd'){
        printf("\nCUDA generates the matrix Omega\n");
         if(CURAND_STATUS_SUCCESS != curandGenerateNormalDouble(cuGen,(__lpkdoublereal *)(omega.raw),omega.nrows * omega.ncols,0.,1.)){

             cout << endl << "Problem in random number generation"<< endl;
         }
    }
    else{
        printf("\nCUDA generates the matrix Omega\n");

         if(CURAND_STATUS_SUCCESS != curandGenerateNormal(cuGen,(__lpkreal *)(omega.raw),  omega.nrows * omega.ncols,0.,1.)){

             cout << endl << "Problem in random number generation"<< endl;
         }
    }

#else
//NORMAL CODE
    mt19937_64 generator;
	//normal_distribution<lpkreal> norm_distribution(0.,1.);
	//normal_distribution<lpkreal> norm_generator(0.,1.);
	normal_distribution<__lpkdoublereal> norm_distribution_double(0.,1.);
	normal_distribution<__lpkreal> norm_distribution_single(0.,1.);
    generator.seed(time(NULL));
    
    //C2appo.raw IS A.raw. No copy
    C2appo.raw = A.raw;

    //Allocate omega.raw for standard code
    omega.raw = new T[omega.nrows*omega.ncols];

    if (costanti.typeCheck == 'z'){
        for(int i=0; i<2 * omega.nrows * omega.ncols ;i++){
            *(((__lpkdoublereal *) (omega.raw)) +i)=norm_distribution_double(generator);
        }
    }
    else if (costanti.typeCheck == 'c'){
        for(int i=0; i<2 * omega.nrows * omega.ncols ;i++){
            *(((__lpkreal *) (omega.raw)) +i)=norm_distribution_single(generator);
        }
    }
    else if(costanti.typeCheck == 'd'){
        for(int i=0; i< omega.nrows * omega.ncols ;i++){
            *(( (__lpkdoublereal *)(omega.raw)) +i)=norm_distribution_double(generator);
        }
    }
    else{
        for(int i=0; i< omega.nrows * omega.ncols ;i++){
            *(((__lpkreal *) (omega.raw)) +i)=norm_distribution_single(generator);
        }
    }
#endif

    /******************************************************/
    /*Here we have C2appo and omega allocated where needed*/
    /******************************************************/

	printf("...done!\n");

	QRfactLaPack<T> *appoQiterQR;

	//__DEBUG
	//printf("\nmatrice c2 prima di interazioni...\n");

    //__DEBUG
    
   // tamaClss::genMatrix<T> C2prova(C2appo.nrows, C2appo.ncols);
   // cudaMemcpy(C2prova.raw,C2appo.raw, C2appo.nrows*C2appo.ncols*sizeof(T),cudaMemcpyDeviceToHost);

   // 
   // C2prova.toBinFile("C2prova.dat");

   // tamaClss::genMatrix<T> omegaprova(omega.nrows, omega.ncols);
   // cudaMemcpy(omegaprova.raw,omega.raw, omega.nrows*omega.ncols*sizeof(T),cudaMemcpyDeviceToHost);
   // 
   // omegaprova.toBinFile("omegaprova.dat");



    printf("\nRandomized power iteration: q= %d...",niter);
	appoQiterQR = RandSubIter( C2appo,omega,niter,tolerance);
	printf("...done!\n");

	//tamaClss::genMatrix<T> *iterQ = appoQiterQR->mQ;
    //__ATTENTION: Q0???
	tamaClss::genMatrix<T> *Q0;
	
	Q0 = appoQiterQR->mQ;


    printf("\nStarting direct SVD...");

#ifdef CUDA_ENABLED

    tamaClss::genMatrix<T>B;
    B.nrows = Q0->ncols;
    B.ncols = C2appo.ncols;
    cudaMalloc((void **)&(B.raw),B.nrows*B.ncols*sizeof(T));

#else
	tamaClss::genMatrix<T> B(Q0->ncols,C2appo.ncols);
#endif
	//B = Q*.A

	//__DEBUG
	//	printf("\nMatrice B\n");
	//	printf("\nQ0 =[%d,%d] \n",Q0->nrows,Q0->ncols);
	//	printf("\nC4 =[%d,%d] \n",C2.nrows,C2.ncols);
	//	printf("\nB =[%d,%d] \n",B.nrows,B.ncols);
	TamaGEMM('T','N',*Q0,C2appo,B,costanti.typed_one,costanti.typed_zero);


	void * genoutcome;
	tamaClss::genMatrix<T> * outUtilde;

//	tamaClss::genMatrix<T> * outS;
	void * outS;
	tamaClss::genMatrix<T> * outVT;

	using namespace std::chrono;

	high_resolution_clock::time_point t1 = high_resolution_clock::now();

    genoutcome = SVD(B,&outUtilde,&outS,&outVT);

	high_resolution_clock::time_point t2 = high_resolution_clock::now();

	duration<double> time_span = duration_cast<duration<double> > (t2-t1);

	cout <<endl << "time required for DirectSVD:" << time_span.count() << "seconds";



	if(genoutcome !=NULL)
		printf("\nProblema con svd\n");
#ifdef CUDA_ENABLED
    tamaClss::genMatrix<T> *UU = new tamaClss::genMatrix<T>;
    UU->nrows = Q0->nrows;
    UU->ncols = outUtilde->ncols;
    cudaMalloc((void **)&(UU->raw),UU->nrows*UU->ncols*sizeof(T));
#else
	tamaClss::genMatrix<T> *UU = new tamaClss::genMatrix<T>(Q0->nrows,outUtilde->ncols);

#endif
	//GEMM(CblasColMajor,CblasNoTrans,CblasNoTrans,UU->nrows,UU->ncols,Q0->ncols,/*(const double *)*/&typed_one,Q0->raw,Q0->nrows,outUtilde->raw,outUtilde->nrows,/*(const void *)*/&typed_zero,UU->raw,UU->nrows);
	TamaGEMM('N','N',*Q0,*outUtilde,*UU,costanti.typed_one,costanti.typed_zero);

	
	//Cleaning up...
#ifdef CUDA_ENABLED
    //We must copy the device allocated data to host
    //UU
    T * UrawAppo = new T[UU->nrows*UU->ncols];
    cudaMemcpy(UrawAppo,UU->raw,UU->nrows * UU->ncols *sizeof(T),cudaMemcpyDeviceToHost);
    cudaFree(UU->raw);
    UU->raw = UrawAppo;


    //VT
    T* VTrawAppo = new T[outVT->nrows * outVT->ncols];
     cudaMemcpy(VTrawAppo,outVT->raw,outVT->nrows * outVT->ncols *sizeof(T),cudaMemcpyDeviceToHost);
    cudaFree(outVT->raw);
    outVT->raw = VTrawAppo;

    //S
    void * SrawAppo;

    if((costanti.typeCheck == 'z') or (costanti.typeCheck =='d')){
        //Double precision
        SrawAppo = new __lpkdoublereal[UU->ncols];
        cudaMemcpy(SrawAppo,outS,UU->ncols *sizeof(__lpkdoublereal),cudaMemcpyDeviceToHost);
        cudaFree(outS);
        outS = SrawAppo;
    }
    else{
        //Single precision
         SrawAppo = new __lpkreal[UU->ncols];
        cudaMemcpy(SrawAppo,outS,UU->ncols *sizeof(__lpkreal),cudaMemcpyDeviceToHost);
        cudaFree(outS);
        outS = SrawAppo;
    }

    cudaFree(Q0->raw);
    Q0->raw= NULL;
    cudaFree(outUtilde->raw);
    outUtilde->raw=NULL;

#endif
    
	*U = UU;
	*S = outS;
	*VT = outVT;

    high_resolution_clock::time_point t2overall = high_resolution_clock::now();

	duration<double> time_spanoverall = duration_cast<duration<double> > (t2overall-t1overall);

	double tOverall;
	tOverall = time_spanoverall.count();

	cout <<endl << "time overall time required for the randomized SVD:" << tOverall << "seconds";


#ifdef CUDA_ENABLED
    //Delete the device A
    cudaFree(C2appo.raw);
    C2appo.raw = NULL;
    cudaFree(B.raw);
    B.raw = NULL;
    cudaFree(omega.raw);
    omega.raw = NULL;
#else

    //Do not delete the original matrix A on exit!!!
    C2appo.raw=NULL;
#endif

   // {int pippo;
   //     cin >>pippo;
   // }
    delete appoQiterQR;
    delete Q0;
    delete outUtilde;

	    //__DEBUG
    //We shutdown the cula and cublas environments.
    //We reset the device (in case we left some memory allocated)

#ifdef CUDA_ENABLED
    culaShutdown();
    cublasShutdown();
#endif
    //cudaDeviceReset();    
    cout<< "\nend of RRSVD \n";
	//delete outS;
	//delete outVT;

}

//#ifndef REAL
//Double precision complex
void z_rrsvd_(integer *nrowsA, integer *ncolsA, __lpkdoublecomplex **rawA, integer *relevant, integer *niter, __lpkdoublereal *tolerance, integer * nrowsU, integer * ncolsU, __lpkdoublecomplex ** rawU,integer * nrowsS, integer * ncolsS, __lpkdoublereal **rawS,integer * nrowsVT, integer * ncolsVT, __lpkdoublecomplex **rawVT){

	tamaClss::genMatrix<__lpkdoublecomplex> * matA = new tamaClss::genMatrix<__lpkdoublecomplex>;
    matA->nrows = *nrowsA;
    matA->ncols = *ncolsA;
	matA->raw = *rawA;
	
    tamaClss::genMatrix<__lpkdoublecomplex> * U;
//	tamaClss::genMatrix<__lpkdoublereal> * S;
	void * Sraw;
	tamaClss::genMatrix<__lpkdoublecomplex> * VT;

	RRSVD(*matA, *relevant,*niter, *tolerance, &U, &Sraw, &VT);
	*nrowsU = U->nrows;
	*ncolsU = U->ncols;
	*rawU = U->raw;

	*nrowsS = U->ncols;
	*ncolsS = 1;//VT->nrows;
	*rawS = (__lpkdoublereal *)Sraw;

	*nrowsVT = VT->nrows;
	*ncolsVT = VT->ncols;
	*rawVT = VT->raw;

}
//#else


void c_rrsvd_(integer *nrowsA, integer *ncolsA, __lpkcomplex **rawA, integer *relevant, integer *niter, __lpkdoublereal *tolerance, integer * nrowsU, integer * ncolsU, __lpkcomplex ** rawU,integer * nrowsS, integer * ncolsS, __lpkreal **rawS,integer * nrowsVT, integer * ncolsVT, __lpkcomplex **rawVT){

	tamaClss::genMatrix<__lpkcomplex> * matA = new tamaClss::genMatrix<__lpkcomplex>;
    matA->nrows = *nrowsA;
    matA->ncols = *ncolsA;
    matA->raw = *rawA;

	tamaClss::genMatrix<__lpkcomplex> * U;
//	tamaClss::genMatrix<__lpkdoublereal> * S;
	void * Sraw;
	tamaClss::genMatrix<__lpkcomplex> * VT;

	RRSVD(*matA, *relevant,*niter, *tolerance, &U, &Sraw, &VT);
	*nrowsU = U->nrows;
	*ncolsU = U->ncols;
	*rawU = U->raw;

	*nrowsS = U->ncols;
	*ncolsS = 1;//VT->nrows;
	*rawS = (__lpkreal *)Sraw;

	*nrowsVT = VT->nrows;
	*ncolsVT = VT->ncols;
	*rawVT = VT->raw;

}




//Double precision real
void d_rrsvd_(integer *nrowsA, integer *ncolsA, __lpkdoublereal **rawA, integer *relevant, integer *niter, __lpkdoublereal *tolerance, integer * nrowsU, integer * ncolsU, __lpkdoublereal ** rawU,integer * nrowsS, integer * ncolsS, __lpkdoublereal **rawS,integer * nrowsVT, integer * ncolsVT, __lpkdoublereal **rawVT){

	tamaClss::genMatrix<__lpkdoublereal> * matA = new tamaClss::genMatrix<__lpkdoublereal>;
    matA->nrows = *nrowsA;
    matA->ncols = *ncolsA;
    matA->raw = *rawA;

	tamaClss::genMatrix<__lpkdoublereal> * U;
//	tamaClss::genMatrix<lpkreal> * S;
	void * Sraw;
	tamaClss::genMatrix<__lpkdoublereal> * VT;

	RRSVD(*matA, *relevant,*niter, *tolerance, &U, &Sraw, &VT);
	*nrowsU = U->nrows;
	*ncolsU = U->ncols;
	*rawU = U->raw;

	*nrowsS = U->ncols;
	*ncolsS = 1;//VT->nrows;
	//*rawS = S->raw;
	*rawS = (__lpkdoublereal *)Sraw;

	*nrowsVT = VT->nrows;
	*ncolsVT = VT->ncols;
	*rawVT = VT->raw;

}




//Single precision real
void s_rrsvd_(integer *nrowsA, integer *ncolsA, __lpkreal **rawA, integer *relevant, integer *niter, __lpkdoublereal *tolerance, integer * nrowsU, integer * ncolsU, __lpkreal ** rawU,integer * nrowsS, integer * ncolsS, __lpkreal **rawS,integer * nrowsVT, integer * ncolsVT, __lpkreal **rawVT){

	tamaClss::genMatrix<__lpkreal> * matA = new tamaClss::genMatrix<__lpkreal>;
    matA->nrows = *nrowsA;
    matA->ncols = *ncolsA;
	matA->raw = *rawA;

	tamaClss::genMatrix<__lpkreal> * U;
//	tamaClss::genMatrix<lpkreal> * S;
	void * Sraw;
	tamaClss::genMatrix<__lpkreal> * VT;

	RRSVD(*matA, *relevant,*niter, *tolerance, &U, &Sraw, &VT);
	*nrowsU = U->nrows;
	*ncolsU = U->ncols;
	*rawU = U->raw;

	*nrowsS = U->ncols;
	*ncolsS = 1;//VT->nrows;
	//*rawS = S->raw;
	*rawS = (__lpkreal *)Sraw;

	*nrowsVT = VT->nrows;
	*ncolsVT = VT->ncols;
	*rawVT = VT->raw;

}
//#endif


void z_dealloc_(__lpkdoublecomplex ** data){
 delete [] (*data);
}

void c_dealloc_(__lpkcomplex **data){
 delete [] (*data);
}

void d_dealloc_(__lpkdoublereal **data){
delete[] (*data);
}

void s_dealloc_(__lpkreal **data){
delete [] (*data);
}




template <typename T>
__lpkdoublereal frobeniusNorm(tamaClss::genMatrix<T> *B){

//	__lpkdoublereal frobNorm /*= 1e-8*/;

	printf("\nFunction Frobenius_norm to be defined\n");
	return 1e10;
	//T appo;
//	int dimRatio;
//	dimRatio = sizeof(T)/sizeof(lpkreal);
//
//	frobNorm = NORM2REAL(dimRatio* (B->nrows) * (B->ncols), (lpkreal *) B->raw,1);
////#ifdef REAL //if we are dealing with reals...
////			frobNorm = sqrt(DOT((B->nrows* B->ncols),B->raw, 1 , B->raw,1));
////
////
////#else  //If we are dealing with complexes....
////			DOT((B->nrows* B->ncols),B->raw, 1 , B->raw,1,&appo);
////			frobNorm = sqrt(appo.r);
////#endif

		//	return frobNorm;

}


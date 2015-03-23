/*
 * decompositions.h
 *
 *  Created on: Sep 30, 2014
 *      Author: Tama
 */

#ifndef DECOMPOSITIONS_H_
#define DECOMPOSITIONS_H_
#include "genMatrix.h"
#include <random>
#include <ctime>
#include <chrono>
#include <cassert>

#include "ZRRSVD.h"

using namespace std;

template<typename T>
class QRfactLaPack{
public:
	integer nrows;
	integer ncols;
	tamaClss::genMatrix<T> *mQ;
	//doublecomplex *tau;
	tamaClss::genMatrix<T> *mR;


	QRfactLaPack(int,int); //Constructor
	QRfactLaPack(QRfactLaPack&); //Copy constructor
	~QRfactLaPack(); //Desctructor

};

template <typename T>
struct svdTriple{
	tamaClss::genMatrix<T> U;
	tamaClss::genMatrix<T> Sigma;
	tamaClss::genMatrix<T> VT;

};


template <typename T>
QRfactLaPack<T> * reorthogonalizeLaPack( tamaClss::genMatrix<T>  &matObj, char r);

template <typename T>
QRfactLaPack<T> * onlineCheck(tamaClss::genMatrix<T> &matA, tamaClss::genMatrix<T> &Q, __lpkdoublereal tolerance, integer incr, mt19937_64 &generator);

template <typename T>
QRfactLaPack<T> * RandSubIter( tamaClss::genMatrix<T>  &matA,tamaClss::genMatrix<T>& matO, integer niter,__lpkdoublereal tolerance,mt19937_64 &generator);

template <typename T>
void * SVD(tamaClss::genMatrix<T> &A,tamaClss::genMatrix<T> **U, void **S, tamaClss::genMatrix<T> **VT);

template <typename T>
void RRSVD(tamaClss::genMatrix<T> &A,integer relevant, integer niter, __lpkdoublereal tolerance, tamaClss::genMatrix<T> **U, /*tamaClss::genMatrix<T>*/void **S, tamaClss::genMatrix<T> **VT);

template <typename T>
__lpkdoublereal frobeniusNorm(tamaClss::genMatrix<T> *B);

template <typename T>
void TamaGEMM(char opA, char opB, tamaClss::genMatrix<T> &A, tamaClss::genMatrix<T> &B, tamaClss::genMatrix<T> &C, const T c1, const T c2);

template<>
void TamaGEMM<__lpkdoublecomplex>(char opA, char opB, tamaClss::genMatrix<__lpkdoublecomplex> &A, tamaClss::genMatrix<__lpkdoublecomplex> &B, tamaClss::genMatrix<__lpkdoublecomplex> &C, const __lpkdoublecomplex c1, const __lpkdoublecomplex c2);

template<>
void TamaGEMM<__lpkreal>(char opA, char opB, tamaClss::genMatrix<__lpkreal> &A, tamaClss::genMatrix<__lpkreal> &B, tamaClss::genMatrix<__lpkreal> &C, const __lpkreal c1, const __lpkreal c2);

template<>
void TamaGEMM<__lpkcomplex>(char opA, char opB, tamaClss::genMatrix<__lpkcomplex> &A, tamaClss::genMatrix<__lpkcomplex> &B, tamaClss::genMatrix<__lpkcomplex> &C, const __lpkcomplex c1, const __lpkcomplex c2);

template<>
void TamaGEMM<__lpkreal>(char opA, char opB, tamaClss::genMatrix<__lpkreal> &A, tamaClss::genMatrix<__lpkreal> &B, tamaClss::genMatrix<__lpkreal> &C, const __lpkreal c1, const __lpkreal c2);

template<typename T>
void TamaGEQRF(int *m, int *n, T *mat, int *lda, T * tau, T * lwork, int *lworkSize, int *info);

template <>
void TamaGEQRF<__lpkdoublecomplex>(int *m, int *n, __lpkdoublecomplex *mat, int *lda, __lpkdoublecomplex * tau, __lpkdoublecomplex * lwork, int *lworkSize, int *info);

template <>
void TamaGEQRF<__lpkcomplex>(int *m, int *n, __lpkcomplex *mat, int *lda, __lpkcomplex * tau, __lpkcomplex * lwork, int *lworkSize, int *info);

template <>
void TamaGEQRF<__lpkdoublereal>(int *m, int *n, __lpkdoublereal *mat, int *lda, __lpkdoublereal * tau, __lpkdoublereal * lwork, int *lworkSize, int *info);

template <>
void TamaGEQRF<__lpkreal>(int *m, int *n, __lpkreal *mat, int *lda, __lpkreal * tau, __lpkreal * lwork, int *lworkSize, int *info);

template <typename T>
void TamaUNGQR(int *m, int *n, int * kk,  T *mat, int *lda, T * tau, T * lwork, int *lworkSize, int *info);

template <>
void TamaUNGQR<__lpkdoublecomplex>(int *m, int *n, int * kk, __lpkdoublecomplex *mat, int *lda, __lpkdoublecomplex * tau, __lpkdoublecomplex * lwork, int *lworkSize, int *info);

template <>
void TamaUNGQR<__lpkcomplex>(int *m, int *n, int * kk, __lpkcomplex *mat, int *lda, __lpkcomplex * tau, __lpkcomplex * lwork, int *lworkSize, int *info);

template <>
void TamaUNGQR<__lpkdoublereal>(int *m, int *n, int * kk, __lpkdoublereal *mat, int *lda, __lpkdoublereal * tau, __lpkdoublereal * lwork, int *lworkSize, int *info);


template <>
void TamaUNGQR<__lpkreal>(int *m, int *n, int * kk, __lpkreal *mat, int *lda, __lpkreal * tau, __lpkreal * lwork, int *lworkSize, int *info);

template <typename T>
__lpkdoublereal TamaNORM2REAL(integer dim, T *data);

template <>
__lpkdoublereal TamaNORM2REAL<__lpkdoublecomplex>(integer dim, __lpkdoublecomplex * data);

template <>
__lpkdoublereal TamaNORM2REAL<__lpkcomplex>(integer dim, __lpkcomplex * data);

template <>
__lpkdoublereal TamaNORM2REAL<__lpkdoublereal>(integer dim, __lpkdoublereal * data);

template <>
__lpkdoublereal TamaNORM2REAL<__lpkreal>(integer dim, __lpkreal * data);

#endif /* DECOMPOSITIONS_H_ */

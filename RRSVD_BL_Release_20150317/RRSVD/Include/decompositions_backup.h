/*
 * decompositions.h
 *
 *  Created on: Sep 30, 2014
 *      Author: Tama
 */

#ifndef DECOMPOSITIONS_H_
#define DECOMPOSITIONS_H_
#include "genMatrix.h"
//#include "printStuff.h"
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
QRfactLaPack<T> * onlineCheck(tamaClss::genMatrix<T> &matA, tamaClss::genMatrix<T> &Q, __lpkdoublereal tolerance, integer incr);

template <typename T>
QRfactLaPack<T> * RandSubIter( tamaClss::genMatrix<T>  &matA,tamaClss::genMatrix<T>& matO, integer niter,__lpkdoublereal tolerance);

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

#endif /* DECOMPOSITIONS_H_ */

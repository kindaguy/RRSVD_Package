/*
 * ZRRSVD.h
 *
 *  Created on: Oct 28, 2014
 *      Author: Tama
 */

#ifndef ZRRSVD_H_
#define ZRRSVD_H_

#ifdef __cplusplus
extern "C"{
#endif
//#ifndef REAL
//Double precision complex
void z_rrsvd_(integer *nrowsA, integer *ncolsA, __lpkdoublecomplex **rawA, integer *relevant, integer *niter, __lpkdoublereal *tolerance, integer * nrowsU, integer * ncolsU, __lpkdoublecomplex **rawU,integer * nrowsS, integer * ncolsS, __lpkdoublereal **rawS,integer * nrowsVT, integer * ncolsVT, __lpkdoublecomplex **rawVT);

//Single precision complex
void c_rrsvd_(integer *nrowsA, integer *ncolsA, __lpkcomplex **rawA, integer *relevant, integer *niter, __lpkdoublereal *tolerance, integer * nrowsU, integer * ncolsU, __lpkcomplex **rawU,integer * nrowsS, integer * ncolsS, __lpkreal **rawS,integer * nrowsVT, integer * ncolsVT, __lpkcomplex **rawVT);

//Double precision real
void d_rrsvd_(integer *nrowsA, integer *ncolsA, __lpkdoublereal **rawA, integer *relevant, integer *niter, __lpkdoublereal *tolerance, integer * nrowsU, integer * ncolsU, __lpkdoublereal **rawU,integer * nrowsS, integer * ncolsS, __lpkdoublereal **rawS,integer * nrowsVT, integer * ncolsVT, __lpkdoublereal **rawVT);

//Single precision real
void s_rrsvd_(integer *nrowsA, integer *ncolsA, __lpkreal **rawA, integer *relevant, integer *niter, __lpkdoublereal *tolerance, integer * nrowsU, integer * ncolsU, __lpkreal **rawU,integer * nrowsS, integer * ncolsS, __lpkreal **rawS,integer * nrowsVT, integer * ncolsVT, __lpkreal **rawVT);

void z_dealloc_(__lpkdoublecomplex ** data);
void c_dealloc_(__lpkcomplex ** data);
void d_dealloc_(__lpkdoublereal ** data);
void s_dealloc_(__lpkreal ** data);


#ifdef __cplusplus
}
#endif


#endif /* ZRRSVD_H_ */

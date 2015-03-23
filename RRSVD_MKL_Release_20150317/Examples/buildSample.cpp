#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <random>
#include "accelerate_includes.h"
//#include "printStuff.h"
//#include "loadStuff.h"
//#include "genMatrix.h"
#include "decompositions.h"
//Timing
#include <ctime>
#include <ratio>
#include <chrono>

using namespace std;
/*****Input PARAMETERS***/
//M: rows of A
//N: columns of A
//K: number of retained columns of the range base matrix
//NITER: number of iterations of the Random Range Finder
//CHECK: a flag: TRUE(1) means that we want to compute and store the 'standard' SVD of A.
//SETSEED: a flag: TRUE(1) means that we want the seed of the random number generator to be set at each execution
//SVFILE: optional: if present, the program will load the singular values from the specified file. The number of singular values must be equal to N

/***ATTENTION****/
/*For the time being the program works only with tall matrices (M>=N) */


int main(int argc, char **argv) {


	puts("!!!Build Sample!!!");

    constants<__lpkdoublecomplex> costanti;
	//Random number generation tools
	//__ATTENTION: set the seed randomly for real executions
	mt19937_64 generator;

	bool setseed = atoi(argv[3]);

	if(setseed){
		cout <<endl << "Setseed...";
		generator.seed(time(NULL));
		cout << "...done!" << endl;
	}


	//Set two random numbers generators
	//Uniform distribution
	uniform_real_distribution<double> unif_distribution(0.,1.);

	//Normal distribution
	normal_distribution<double> norm_distribution(0.,1.);


//	/****Preparation of input matrices**********/
//	/*******************************************/
//	//Set sizes
	integer nrowsA, ncolsA, nsvA, nrowsU, ncolsU, nrowsS, ncolsS, nrowsV, ncolsV;
////
	printf("\nPrepare an input matrix A.\n");
	nrowsA = atoi(argv[1]);
	ncolsA = atoi(argv[2]);
////	//__DEBUG
////	//	printf("\nNumber of rows and columns: \n");
////	//scanf("%d %d",&nrowsA, &ncolsA);
////

    nrowsU = nrowsA;
	ncolsV = ncolsA;

	nsvA = min(nrowsA,ncolsA);

	ncolsU = nsvA;

	nrowsV = nsvA;

	nrowsS = nsvA;

	ncolsS = nsvA;

	printf("\nAllocating U,V,S\n");
	tamaClss::genMatrix<__lpkdoublecomplex> U(nrowsA,ncolsA);
	tamaClss::genMatrix<__lpkdoublecomplex> V(nrowsV,ncolsV);
	tamaClss::genMatrix<__lpkdoublecomplex> S(nrowsS,ncolsS);

//	//Generate the matrices at random (elements extracted from a unif(0,1) distribution)
//
//	printf("\nGenerating random matrices...");
//
	//U
	for(int i=0; i< U.nrows*U.ncols; i++){
		U.raw[i].real = unif_distribution(generator);
		U.raw[i].imag =/*0.;*/ unif_distribution(generator);
	}
	//V
	for(int i=0; i< nrowsV*ncolsV; i++){
		V.raw[i].real = unif_distribution(generator);
		V.raw[i].imag = unif_distribution(generator);
	}



//
//	printf("...done!\n");
//
//	//__DEBUG
//	//	printf("\nMatrice U prima di riortho\n");
//	//
//	//	printMatrixComplex(U.format,U.raw,U.nrows,U.ncols);
//
	tamaClss::genMatrix<__lpkdoublecomplex> *UQ;
//
	printf("\nReorthogonalizing random matrices....");
//	//Do the renormalization
	QRfactLaPack<__lpkdoublecomplex> *appoUQR = reorthogonalizeLaPack( U,'n');
	UQ = appoUQR->mQ;
//
	printf("\nU riorthogonalized!\n");
//
//	//__DEBUG
//	//printf("\nMatrice U\n");
//	//	printMatrixComplex(UQ->format,UQ->raw,UQ->nrows,UQ->ncols);
//	//
//	//	{
//	//		int check;
//	//		scanf("%d",&check);
//	//	}
//
	UQ->toBinFile("matrixU.dat");
//	printMatrixComplexBinFile("matriceU.dat", UQ->format,UQ->raw, UQ->nrows, UQ->ncols);
//
//
//	//__DEBUG
//	//		printf("\nMatrice V prima di riortogonalizzazione\n");
//	//
//	//		printMatrixComplex(V.format,V.raw,V.nrows,V.ncols);
//
	tamaClss::genMatrix<__lpkdoublecomplex> *VQ;
	QRfactLaPack<__lpkdoublecomplex> *appoVQR = reorthogonalizeLaPack( V,'n');
	VQ = appoVQR->mQ;

//
	printf("\nV riorthogonalized!\n");
	printf("\nRiorthogonalization completed!\n");
//	//delete matRV;
//
//	//__DEBUG
//	//	printf("\nMatrice V\n");
//	//	printMatrixComplex(VQ->format,VQ->raw,VQ->nrows,VQ->ncols);
//
//
	VQ->toBinFile("matrixV.dat");
//	printMatrixComplexBinFile("matriceV.dat",VQ->format, VQ->raw, VQ->nrows, VQ->ncols);
//
//	//Now we generate the matrix S. It elements are not random BUT deterministically computed.
//
	//	//	integer ldS,svNum;
//	//	ldS = (UQ->ncols>VQ->nrows?UQ->ncols:VQ->nrows);
//	//	svNum = (UQ->ncols>VQ->nrows?VQ->nrows:UQ->ncols);
	printf("\nU Generating singular values....");
	nrowsS = UQ->ncols;
	ncolsS = VQ->nrows;
//
	assert(nrowsS==ncolsS);
	__lpkdoublereal appo[nrowsS];
//
//
	if (argc==6){
	int numSV;

	ifstream fileSVIn;
	fileSVIn.open(argv[5],ios::binary);
	if(fileSVIn.fail()){
		printf("\nProblema apertura file %s\n",argv[5]);
	}

	//Read the number of singular values
	fileSVIn.read((char*) &numSV,sizeof(integer));

    if( numSV != min(nrowsA,ncolsA)){
        cout <<"\n The number of singular values in the file do not match the dimensions!!!\n";
        return -1;
    
    }
	//Read the singular values
	for(int i = 0; i<numSV;i++){
		fileSVIn.read((char*) &(appo[i]),sizeof(__lpkdoublereal));
	}

	fileSVIn.close();//Close the file
	}
	else{
		printf("\nGenerating singular values 1/x \n");
		//__TEST CASE
		for(int i=0; i<nrowsS;i++){
			appo[i]=/* exp(-i/10.);*/1./(pow(i+1,1));
		//	printf(" %f ",appo[i]);
			/*exp(-i/100.);*/
		}

	}
//
//	//The positive thing is that LaPack does the same!!!!
//
	__lpkdoublereal norm=TamaNORM2REAL(nrowsS,appo);
//
//	printf("\n\nnorm= %f\n\n",norm);
//
	cblas_dscal(nrowsS,1/norm,appo,1);
//
	for(int i=0; i<ncolsS; i++){
		//__DEBUG
		printf(" %f\n ",appo[i]);
		S.raw[i*S.nrows +i].real =appo[i];
		S.raw[i*S.nrows +i].imag =0.;
	}
////
//	S.toBinFile("matriceS.dat");
//	//printMatrixComplexBinFile("matriceS.dat",S.format, S.raw, S.nrows, S.ncols);
////
//	printf("...done!\n");
////	//__DEBUG
//////		printf("\nMatrice S\n");
//////		printMatrixComplex(S.format,S.raw,S.nrows,S.ncols);
//////		{
//////			int check;
//////			cin >> check;
//////		}
////
//	printf("\nBuilding the case-study matrix A...");
	integer nrowsC1 = nrowsS;
	integer ncolsC1 = nrowsV;
	tamaClss::genMatrix<__lpkdoublecomplex> C1(nrowsC1,ncolsC1);

	TamaGEMM('N','T',S,*VQ,C1, costanti.typed_one,costanti.typed_zero);
//	//__DEBUG
//	//	printf("\nMatrice C1\n");
//	//	printMatrixComplex(C1.format,C1.raw,C1.nrows,C1.ncols);
//
	integer nrowsC2 = UQ->nrows;
	integer ncolsC2 = C1.ncols;
	tamaClss::genMatrix<__lpkdoublecomplex> C2(nrowsC2,ncolsC2);
//
	TamaGEMM('N','N',*UQ,C1,C2,costanti.typed_one,costanti.typed_zero);

    C2.toBinFile(argv[4]);
    return 0;
}

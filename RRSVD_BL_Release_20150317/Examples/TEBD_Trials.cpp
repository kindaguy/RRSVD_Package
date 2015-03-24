//============================================================================
// Name        : TEBD_Trials.cpp
// Author      : Dario Tamascelli
// Version     :
// Copyright   : 
// Description : Hello World in C, Ansi-style
//============================================================================

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <random>
#include "accelerate_includes.h"
#include "genMatrix.h"
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


	puts("!!!Hello World!!!");

	//Random number generation tools
	//__ATTENTION: set the seed randomly for real executions
	mt19937_64 generator;

	bool setseed = atoi(argv[6]);

	if(setseed){
		cout <<endl << "Setseed...";
		generator.seed(time(NULL));
		cout << "...done!" << endl;
	}


	/*****When matrix C2 is already in some file, uncomment the following line*****/
    for(int sample = 0; sample<atoi(argv[7]); sample++){

        tamaClss::genMatrix<__lpkdoublecomplex> C2;

        string sampleFileName;
        sampleFileName = argv[1]/*+to_string(sample)+".dat"*/;

	
	C2.fromBinFile(sampleFileName.c_str()/*argv[7]*/);

//	C2.toBinFile("matriceC2.dat");
//
//	C2.toFile("matriceC2text.dat");
//
//	printf("...done!\n");
//	//C2 is the matrix we want to SVDecompose.

	/**Fine Preparazione della matrice di input**/
	/*******************************************/

	/*****SVD of the whole A matrix******/
	//Now we do the SVD of A (we keep A for debugging)

	bool doCheck = atoi(argv[5]);
	if(doCheck){
		tamaClss::genMatrix<__lpkdoublecomplex> Aappo = C2;

		__lpkdoublereal * genoutcomeA;
		tamaClss::genMatrix<__lpkdoublecomplex> * outUA;
		//tamaClss::genMatrix<lpkcomplex> * outSA;
		tamaClss::genMatrix<__lpkdoublereal> * outSA;
		tamaClss::genMatrix<__lpkdoublecomplex> * outVTA;
		void * appoSA;
		using namespace std::chrono;

		high_resolution_clock::time_point t1A = high_resolution_clock::now();
		genoutcomeA = NULL;

		printf("\nPrima SVD\n");
		SVD(Aappo,&outUA,/*&outSA*/&appoSA,&outVTA);
		printf("\nFine SVD\n");
		//__DEBUG

	//	{
	//		int pippo;
	//		cin >>pippo;
	//	}
		outSA = new tamaClss::genMatrix<__lpkdoublereal>();
		outSA->ncols = 1;
		outSA->nrows = outUA->ncols;
		outSA->raw = (__lpkdoublereal*) appoSA;
		high_resolution_clock::time_point t2A = high_resolution_clock::now();

		duration<double> time_spanA = duration_cast<duration<double> > (t2A-t1A);

		double timeSVDA;
		timeSVDA = time_spanA.count();

		cout <<endl << "time required for SVD of A:" << timeSVDA<< "seconds";

		if(genoutcomeA !=NULL){
			printf("\nProblema con svd A\n");
			return -1;
		}

        sampleFileName="decUA_"+to_string(sample)+".dat";
		outUA->toBinFile(sampleFileName.c_str()/*"decUA.dat"*/);
		//printMatrixComplexBinFile("decUA.dat",outUA->format,outUA->raw,outUA->nrows,outUA->ncols);
        sampleFileName="decSA_"+to_string(sample)+".dat";
		outSA->toBinFile(sampleFileName.c_str()/*"decSA.dat"*/);
		//printMatrixComplexBinFile("decSA.dat",outSA->format,outSA->raw,outSA->nrows,outSA->ncols);
        sampleFileName="decVTA_"+to_string(sample)+".dat";
		outVTA->toBinFile(sampleFileName.c_str()/*"decVTA.dat"*/);
//		printMatrixComplexBinFile("decVTA.dat",outVTA->format,outVTA->raw,outVTA->nrows,outVTA->ncols);

		ofstream fileTempo;
        sampleFileName = "A_SVD_Timing_"+to_string(sample)+".dat";
		fileTempo.open(sampleFileName.c_str()/*"A_SVD_Timing.dat"*/,ios::binary);
		fileTempo.write((char *)&(timeSVDA),sizeof(double));
		fileTempo.close();

		delete [] genoutcomeA;
		delete outUA;
		delete outSA;
		delete outVTA;

	}

	//__DEBUG
	//	printf("\nMatrice C2 dopo SVD\n");
	//	printMatrixComplex(C2.format,C2.raw,C2.nrows,C2.ncols);
	//printMatrixComplex('C',matrU,leadingDimU,leadingDimVT);

	//__DEBUG
	//	{
	//		int pappo;
	//		cin >>pappo;
	//	}




	/*DIMENSIONALITY REDUCTION*/



	//normal_distribution<lpkreal> norm_generator(0.,1.);
	integer dimRange = atoi(argv[2]);
	integer niter = atoi(argv[3]);

	printf("\nThe number of columns of A is %d.\n Number of relevant dimensions:%d\n",C2.ncols, dimRange);

	generator.seed(time(NULL));



	/*for(int ii=0; ii<1; ii++){*/
		using namespace std::chrono;
			//	printf("\nThe number of columns of A is %d; choose a number < %d\n",ncolsA, ncolsA);
		//	scanf("%d",&dimRange);

		tamaClss::genMatrix<__lpkdoublecomplex> C2appo = C2;
		tamaClss::genMatrix<__lpkdoublecomplex> * RRdecU = new tamaClss::genMatrix<__lpkdoublecomplex>;
//		tamaClss::genMatrix<__lpkdoublecomplex> * RRdecS = new tamaClss::genMatrix<__lpkdoublecomplex>;

		tamaClss::genMatrix<__lpkdoublereal> *RRdecS = new tamaClss::genMatrix<__lpkdoublereal>;

		tamaClss::genMatrix<__lpkdoublecomplex> * RRdecVT = new tamaClss::genMatrix<__lpkdoublecomplex>;


		__lpkdoublereal tolerance = atof(argv[4]);
		int nrU,ncU,nrS,ncS,nrVT,ncVT;
		__lpkdoublecomplex *rawU,/**rawS,*/*rawVT;
		__lpkdoublereal *rawSappo;

//void Z_RRSVD(integer *nrowsA, integer *ncolsA, __lpkdoublecomplex *rawA, integer *relevant, integer *niter, lpkreal *tolerance, integer * nrowsU, integer * ncolsU, __lpkdoublecomplex ** rawU,integer * nrowsS, integer * ncolsS, __lpkdoublecomplex **rawS,integer * nrowsVT, integer * ncolsVT, __lpkdoublecomplex **rawVT)
		high_resolution_clock::time_point t1overall = high_resolution_clock::now();

		z_rrsvd_(&(C2appo.nrows),&(C2appo.ncols),&(C2appo.raw),&dimRange,&niter,&tolerance,&nrU,&ncU,&rawU,&nrS,&ncS,&rawSappo,&nrVT,&ncVT,&rawVT);

		RRdecU->nrows=nrU;
		RRdecU->ncols = ncU;
		RRdecU->raw = rawU;

		RRdecS->nrows=nrS;
		RRdecS->ncols = 1;
		RRdecS->raw = (__lpkdoublereal *)rawSappo;

		RRdecVT->nrows=nrVT;
		RRdecVT->ncols = ncVT;
		RRdecVT->raw = rawVT;

		//RRSVD(C2appo,dimRange, niter, 1e-4, &RRdecU, &RRdecS, &RRdecVT);


		high_resolution_clock::time_point t2overall = high_resolution_clock::now();
//
		duration<double> time_spanoverall = duration_cast<duration<double> > (t2overall-t1overall);
//
		double tOverall;
		tOverall = time_spanoverall.count();
//
		cout <<endl << "time overall time required for the randomized SVD:" << tOverall << "seconds";
		ofstream fileTempo;
		string fileName2;
		fileName2 = "RRSVD_Timing_"+to_string(sample)+".dat";
		fileTempo.open(fileName2.c_str(),ios::binary);
		fileTempo.write((char *)&(tOverall),sizeof(double));
		fileTempo.close();
//
//		printf("\nStampo la matrice finale su file...\n");
//
//		//printMatrixComplex(UU.format,UU.raw,UU.nrows,UU.ncols);
//
//		//printMatrixComplexBinFile("matriceUU.dat",UU.format,UU.raw,UU.nrows,UU.ncols);
//
//
//		printf("\nFatto!\n");
//
		string fileName;
		//Write RRSVD onto file.
		fileName= "RRdecU"+to_string(sample)+".dat";
		cout << endl <<fileName << endl;
		RRdecU->toBinFile(fileName.c_str());
//
//
//		//printMatrixComplexBinFile(fileName.c_str(),UU.format,UU.raw,UU.nrows,UU.ncols);
//		//printMatrixComplexBinFile("decU.dat",outUtilde->format,outUtilde->raw,outUtilde->nrows,outUtilde->ncols);
		fileName= "RRdecS"+to_string(sample)+".dat";
		cout << endl <<fileName << endl;
		RRdecS->toBinFile(fileName.c_str());
//
//		//printMatrixComplexBinFile(fileName.c_str(),outS->format,outS->raw,outS->nrows,outS->ncols);
		fileName= "RRdecVT"+to_string(sample)+".dat";
		cout << endl <<fileName << endl;
		RRdecVT->toBinFile(fileName.c_str());
//
//		//printMatrixComplexBinFile(fileName.c_str(),outVT->format,outVT->raw,outVT->nrows,outVT->ncols);
//
//		//Cleaning up...
//		delete appoQiterQR;
//		delete Q0;
//		delete outUtilde;
//		delete outS;
//		delete outVT;
	}


	
	cout << "fine programma" << endl;


	return EXIT_SUCCESS;
}

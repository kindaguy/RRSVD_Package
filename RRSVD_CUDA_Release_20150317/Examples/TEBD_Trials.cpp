//
//  main.c
//  RRSVD_extern
//
//  Created by Dario Tamascelli on 28/10/14.
//  Copyright (c) 2014 Dario Tamascelli. All rights reserved.
//

#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <chrono>
#include <cassert>
#include <random>
#include <ctime>

#include "accelerate_includes.h"
#include "genMatrix.h"
#include "decompositions.h"
#include "ZRRSVD.h"
//#include "decompositions.cpp"

//#include <cuda.h>
//#include <cuda_runtime_api.h>
//#include <cublas.h>
//#include <cula_lapack.h>
//#include <cula_lapack_device.h>
//#include <curand.h>
using namespace std;

__host__ void deviceQuery(){

	int count;
	cudaGetDeviceCount(&count);

	cudaDeviceProp prop;

	for(int i=0;i<count; i++){
		cudaGetDeviceProperties(&prop,i);
		cout << "\n Compute capability: " << prop.major <<"."<<prop.minor << endl;
		cout <<  "\n Number of multiprocessors on the device " << i << ": " << prop.multiProcessorCount << endl;
		cout <<  "\n Number of thread per block: " << prop.maxThreadsPerBlock << endl;
		cout <<  "\n Warp dimension: " << prop.warpSize << endl;
		cout << "\n Number of registers per block " << prop.regsPerBlock << endl;
		cout << "\n Max (balanced) number of registers per thread: "<< prop.regsPerBlock /prop.maxThreadsPerBlock << endl;
		cout << "\n Max thread per multiprocessor:  " << prop.maxThreadsPerMultiProcessor << endl;
		cout << "\n Number shared memory (32 bit words) per block " << prop.sharedMemPerBlock/4 << endl;
		cout << "\n Total global memory " << prop.totalGlobalMem << endl;
		cout << "\n Total constant memory " << prop.totalConstMem << endl;
	}

}

using namespace std::chrono;


int main(int argc, const char * argv[]){

	//    deviceQuery();
	// {int pippo; cin >> pippo;} 
	cudaSetDevice(0);

	/*******************************************/
	/*RANDOM NUMBERS SETUP                     */
	/*******************************************/

	//	mt19937_64 generator;
	//	generator.seed(time(NULL));
	//Set two random numbers generators
	//Uniform distribution
	//	uniform_real_distribution<double> unif_distribution(0.,1.);
	//Normal distribution
	//	normal_distribution<double> norm_distribution(0.,1.);
	//
	//
	/*******************************************/
	/*RANDOM NUMBERS SETUP:CUDA                */
	/*******************************************/
	//  size_t cuSampleSize=100;
	//  size_t i;
	curandGenerator_t cuGen;
	// __lpkdoublereal *datiHost;
	// __lpkdoublereal *datiDevice;

	// if(cudaSuccess != cudaMalloc((void **)&datiDevice, cuSampleSize * sizeof(__lpkdoublereal))){
	//     printf("\nProblema allocazione memoria omega device\n");
	//     return -1;
	// }

	curandCreateGenerator(&cuGen,CURAND_RNG_PSEUDO_MTGP32)/*MT19937*/;

	if(atoi(argv[6])){
		curandSetPseudoRandomGeneratorSeed(cuGen,time(NULL));
	}
	//	/****Preparation of input matrices**********/
	//	/*******************************************/

	//{
	//int pippo;
	//cin >> pippo;
	//}	
	//Set sizes
	integer nrowsA, ncolsA, nsvA, nrowsU, ncolsU, nrowsS, ncolsS, nrowsV, ncolsV;

	printf("\nPrepare an input matrix A.\n");
	nrowsA = 720; //atoi(argv[1]);
	ncolsA = 720; //atoi(argv[2]);
	////	//__DEBUG
	////	//	printf("\nNumber of rows and columns: \n");
	////	//scanf("%d %d",&nrowsA, &ncolsA);
	//
	nrowsU = nrowsA;
	ncolsV = ncolsA;

	nsvA = min(nrowsA,ncolsA);
	ncolsU = nsvA;
	nrowsV = nsvA;
	nrowsS = nsvA;
	ncolsS = nsvA;

	printf("\nAllocating U,V,S\n");
	tamaClss::genMatrix<__lpkdoublereal> U(nrowsA,ncolsA);
	tamaClss::genMatrix<__lpkdoublereal> V(nrowsV,ncolsV);
	tamaClss::genMatrix<__lpkdoublereal> S(nrowsS,ncolsS);


	//__CUDA INIT
	cout << endl << "Initializing CUBLAS"<<endl;
	cublasStatus status;
	cublasInit();

	//CULA INIT
	culaStatus culaStat;
	printf("Initializing CULA\n");
	culaStat = culaInitialize();

	/********************************************/
	/*RANDOM MATRICES GENERATION                */
	/********************************************/

	//Generate the matrices at random (elements extracted from a unif(0,1) distribution)

	printf("\nGenerating random matrices...");


	__lpkdoublereal *deviceRandValU;
	//__lpkdoublereal *hostRandVal;

	if(cudaSuccess!= cudaMalloc((void**) &deviceRandValU, U.nrows*U.ncols*sizeof(__lpkdoublereal))){

		cout << endl <<"error in deviceRandValU allocation"<<endl;

	}

	if(CURAND_STATUS_SUCCESS != curandGenerateNormalDouble(cuGen,deviceRandValU, U.nrows*U.ncols,0.,1.)){

		cout << endl << "Problem in random number generation"<< endl;
	}

	if(cudaSuccess!=cudaMemcpy(U.raw,deviceRandValU,U.nrows*U.ncols*sizeof(__lpkdoublereal),cudaMemcpyDeviceToHost)){

		cout << endl << "error in copy dev->host U"<< endl;

	}
	__lpkdoublereal *deviceRandValV;
	//__lpkdoublereal *hostRandVal;

	if(cudaSuccess!= cudaMalloc((void**) &deviceRandValV, V.nrows*V.ncols*sizeof(__lpkdoublereal))){

		cout << endl <<"error in deviceRandValV allocation"<<endl;

	}

	if(CURAND_STATUS_SUCCESS != curandGenerateNormalDouble(cuGen,deviceRandValV, U.nrows*U.ncols,0.,1.)){

		cout << endl << "Problem in random number generation"<< endl;
	}

	if(cudaSuccess!=cudaMemcpy(V.raw,deviceRandValV,V.nrows*V.ncols*sizeof(__lpkdoublereal),cudaMemcpyDeviceToHost)){

		cout << endl << "error in copy dev->host V"<< endl;

	}


	printf("...done!\n");

	//__DEBUG
	//printf("\nPrinting random U\n");
	//for(int i=0; i<10; i++){
	//    cout <<  endl << U.raw[i] << endl;
	//}
	//{
	//    int pippo;
	//    cin >>pippo;
	//}
	/********************************************/
	/*END RANDOM MATRICES GENERATION            */
	/********************************************/


	tamaClss::genMatrix<__lpkdoublereal> *UQ;

	printf("\nReorthogonalizing random matrices....");

	//Do the renormalization
	//__CULA


	using namespace std::chrono;



	high_resolution_clock::time_point t1cula = high_resolution_clock::now();

	printf("Calling culaDgeqrf\n");
	/**********************************************/
	tamaClss::genMatrix<__lpkdoublereal> UcudaAppo;
	UcudaAppo.nrows = U.nrows;
	UcudaAppo.ncols = U.ncols;
	UcudaAppo.raw = deviceRandValU;

	QRfactLaPack<__lpkdoublereal> * appoQRcuda=reorthogonalizeLaPack(UcudaAppo,'n');

	cudaMemcpy(U.raw,appoQRcuda->mQ->raw,U.nrows*U.ncols*sizeof(__lpkdoublereal),cudaMemcpyDeviceToHost);


	delete appoQRcuda;

	high_resolution_clock::time_point t2cula = high_resolution_clock::now();

	U.toBinFile("culaUreortho.dat");
	duration<double> time_spanculaAlloc = duration_cast<duration<double> > (t2cula-t1cula);

	double timeQRcula;
	timeQRcula = time_spanculaAlloc.count();
	printf("\n Tempo richiesto da cula QR: %f\n",timeQRcula);
	/**********************************************/

	/**********************************************/
	tamaClss::genMatrix<__lpkdoublereal> VcudaAppo;
	VcudaAppo.nrows = V.nrows;
	VcudaAppo.ncols = V.ncols;
	VcudaAppo.raw = deviceRandValV;


	appoQRcuda = reorthogonalizeLaPack(VcudaAppo,'n');
	cudaMemcpy(V.raw,appoQRcuda->mQ->raw, V.nrows*V.ncols*sizeof(__lpkdoublereal),cudaMemcpyDeviceToHost);

	V.toBinFile("culaVreortho.dat");

	delete appoQRcuda;

	/**********************************************/             

	printf("\nRiorthogonalization completed!\n");

	/**********************************************/ 

	tamaClss::genMatrix<__lpkdoublereal> svAppo;
	svAppo.fromBinFile("../Data/decSA_backup.dat");
	cout << endl <<svAppo.nrows <<  " " << svAppo.ncols << endl;

	tamaClss::genMatrix<__lpkdoublereal> svMatr(svAppo.nrows,svAppo.nrows);
	for (int i=0; i< svMatr.nrows ; i++){
		for (int j=0; j<svMatr.ncols; j++){
			if(i==j){
				svMatr.raw[svMatr.nrows *j+i] = (__lpkdoublereal) (svAppo.raw[i]);
				//				svMatr.raw[svMatr.nrows *j+i].r = (__lpkdoublereal) (svAppo.raw[i]);
				//				svMatr.raw[svMatr.nrows *j+i].i = 0.;
			}
			else{
				svMatr.raw[svMatr.nrows *j+i] = 0.;
				//				svMatr.raw[svMatr.nrows *j+i].r = 0;
				//				svMatr.raw[svMatr.nrows *j+i].i = 0;
			}	
		}	
	}

	svMatr.toBinFile("lapackSreortho.dat");

	tamaClss::genMatrix<__lpkdoublereal> ScudaAppo;
	ScudaAppo.nrows = svMatr.nrows;
	ScudaAppo.ncols = svMatr.ncols;
	cudaMalloc((void **)&(ScudaAppo.raw), ScudaAppo.nrows *ScudaAppo.ncols*sizeof(__lpkdoublereal));

	cudaMemcpy(ScudaAppo.raw,svMatr.raw, ScudaAppo.nrows *ScudaAppo.ncols* sizeof(__lpkdoublereal),cudaMemcpyHostToDevice);

	/**********************************************/
	constants<__lpkdoublereal> costanti;

	printf("\nBuilding the case-study matrix A...");
	integer nrowsC1 = nrowsS;
	integer ncolsC1 = nrowsV;
	tamaClss::genMatrix<__lpkdoublereal> C1(nrowsC1,ncolsC1);

	tamaClss::genMatrix<__lpkdoublereal> C1cudaAppo;
	C1cudaAppo.nrows = C1.nrows;
	C1cudaAppo.ncols = C1.ncols;


	status = cublasAlloc(C1.nrows * C1.ncols,sizeof(__lpkdoublereal),(void **)&(C1cudaAppo.raw));
	if(status != CUBLAS_STATUS_SUCCESS){
		printf("\nBad allocation of C1 on device\n");
	}


	TamaGEMM('n','t',ScudaAppo,VcudaAppo,C1cudaAppo,costanti.typed_one,costanti.typed_zero);

	//    cublasDgemm('n','t',svMatr.ncols,VQ->ncols,VQ->nrows, costanti.typed_one,cudaS,svMatr.nrows, cudaVQ, VQ->nrows, costanti.typed_zero,cudaC1,C1.nrows);


	status = cublasGetMatrix(C1.nrows,C1.ncols,sizeof(__lpkdoublereal),C1cudaAppo.raw,C1.nrows, C1.raw, C1.nrows);

	//
	//if (status != CUBLAS_STATUS_SUCCESS){
	//
	//    printf("\nErrore in trasferimento risultato device->host\n");
	//}



	C1.toBinFile("C1cuda.dat");


	integer nrowsC2 = U.nrows;
	integer ncolsC2 = C1.ncols;
	tamaClss::genMatrix<__lpkdoublereal> C2(nrowsC2,ncolsC2);
	tamaClss::genMatrix<__lpkdoublereal> C2cudaAppo;
	C2cudaAppo.nrows = C2.nrows;
	C2cudaAppo.ncols = C2.ncols;


	if(cudaSuccess!=cudaMalloc((void **) &(C2cudaAppo.raw), C2cudaAppo.nrows*C2cudaAppo.ncols*sizeof(__lpkdoublereal))){

		cout << "\nProblem while allocating matC2cuda\n";

	}
	else{

		TamaGEMM('N', 'N', UcudaAppo,C1cudaAppo, C2cudaAppo, costanti.typed_one, costanti.typed_zero);

	}

	if(cudaSuccess!=  cudaMemcpy(C2.raw, C2cudaAppo.raw, C2.nrows *C2.ncols*sizeof(__lpkdoublereal),cudaMemcpyDeviceToHost)){

		cout <<"\nError while copying C2!\n";

	}

	// TamaGEMM(CblasColMajor,CblasNoTrans,CblasNoTrans,C2.nrows,C2.ncols,UQ->ncols,/*(const double *)*/&typed_one,UQ->raw,UQ->nrows,C1.raw,C1.ncols,/*(const void *)*/&typed_zero,C2.raw,C2.nrows);
	//	//
	//C2.toBinFile("matriceC2.dat");
	//	
	//Cleaning up;
	cudaFree(UcudaAppo.raw);
	UcudaAppo.raw = NULL;
	cudaFree(VcudaAppo.raw);
	VcudaAppo.raw = NULL;
	cudaFree(ScudaAppo.raw);
	ScudaAppo.raw=NULL;
	cudaFree(C1cudaAppo.raw);
	C1cudaAppo.raw = NULL;
	cudaFree(C2cudaAppo.raw);
	C2cudaAppo.raw = NULL;

	//	
	__lpkdoublereal tolerance  = atof(argv[4]);
	integer niter = atoi(argv[3]);
	integer dimRange = atoi(argv[2]);
	//
	//	//Carica matrice SilviMat2 come complessa doppia precisione
	//
	//	//Cominciamo dalla funzione reale in singola precisione;
	//
	//	tamaClss::genMatrix<__lpkdoublereal> matDoubleReal = tamaClss::genMatrix<__lpkdoublereal> (C2.nrows,C2.ncols);
	//	for(int i=0; i<matDoubleReal.nrows *matDoubleReal.ncols; i++){
	//		matDoubleReal.raw[i] = C2.raw[i];
	//////		matDoubleReal.raw[i].r = C2.raw[i].r;
	//////		matDoubleReal.raw[i].i = C2.raw[i].i;
	//	}
	////
	//	matDoubleReal.toBinFile("matDoubleReal.dat");
	////
	if(atoi(argv[5])){
		tamaClss::genMatrix<__lpkdoublecomplex> C2Silvi;
		////
		//	C2Silvi.fromBinFile("SilviMat3.dat");
		C2Silvi.fromBinFile(argv[1]);

		cout << "\n Matrix dimension: ";
		cout << C2Silvi.nrows << endl;
		cout << C2Silvi.ncols << endl;

		tamaClss::genMatrix<__lpkdoublecomplex> Aappo;
		Aappo.nrows = C2Silvi.nrows;
		Aappo.ncols = C2Silvi.ncols;
		cudaMalloc((void **)&(Aappo.raw),Aappo.nrows*Aappo.ncols*sizeof(__lpkdoublecomplex));

		cudaMemcpy(Aappo.raw,C2Silvi.raw, Aappo.nrows *Aappo.ncols*sizeof(__lpkdoublecomplex),cudaMemcpyHostToDevice);
		__lpkdoublereal * genoutcomeA;
		tamaClss::genMatrix<__lpkdoublecomplex> * outUA;

		tamaClss::genMatrix<__lpkdoublereal> * outSA;

		tamaClss::genMatrix<__lpkdoublecomplex> * outVTA;



		void * appoSA;

		genoutcomeA = NULL;

		cout<< "\nCall SVD...\n";

		high_resolution_clock::time_point t1SVD = high_resolution_clock::now();
		SVD(Aappo,&outUA,&appoSA,&outVTA);

		cout <<"\n...done!\n";

		tamaClss::genMatrix<__lpkdoublecomplex> Ucheck(C2Silvi.nrows,C2Silvi.ncols);
		cudaMemcpy(Ucheck.raw,outUA->raw,Ucheck.nrows*Ucheck.ncols*sizeof(__lpkdoublecomplex),cudaMemcpyDeviceToHost);


		tamaClss::genMatrix<__lpkdoublecomplex> VTcheck(C2Silvi.nrows,C2Silvi.ncols);
		cudaMemcpy(VTcheck.raw,outVTA->raw, VTcheck.nrows * VTcheck.ncols * sizeof(__lpkdoublecomplex),cudaMemcpyDeviceToHost);


		tamaClss::genMatrix<__lpkdoublereal> Scheck(C2Silvi.ncols,1);
		cudaMemcpy(Scheck.raw,appoSA, Scheck.nrows * Scheck.ncols * sizeof(__lpkdoublereal),cudaMemcpyDeviceToHost);

		high_resolution_clock::time_point t2SVD = high_resolution_clock::now();

		duration<double> time_spanSVD = duration_cast<duration<double> > (t2SVD-t1SVD);

		double timeSVD;
		timeSVD = time_spanSVD.count();
		printf("\n Tempo richiesto da culaSVD: %f\n",timeSVD);


		Ucheck.toBinFile("Ucheck.dat");
		VTcheck.toBinFile("VTcheck.dat");
		Scheck.toBinFile("Scheck.dat");


		cudaFree(outVTA->raw);
		cudaFree(outUA->raw);
		cudaFree(appoSA);
		cudaFree(Aappo.raw);
		Aappo.raw = NULL;
		cudaError_t err;
		err = cudaGetLastError();
		printf("\nMain: error %d\n",err);

	}

	for(int sample = 0; sample <atoi(argv[7]); sample++){

		tamaClss::genMatrix<__lpkdoublecomplex> C2Silvi;

		C2Silvi.fromBinFile(argv[1]);

		cout << "\n Matrix dimension: ";
		cout << C2Silvi.nrows << endl;
		cout << C2Silvi.ncols << endl;

		tamaClss::genMatrix<__lpkdoublecomplex> * RRdecU = new tamaClss::genMatrix<__lpkdoublecomplex>;

		tamaClss::genMatrix<__lpkdoublereal> *RRdecS = new tamaClss::genMatrix<__lpkdoublereal>;

		tamaClss::genMatrix<__lpkdoublecomplex> * RRdecVT = new tamaClss::genMatrix<__lpkdoublecomplex>;


		integer nrU,ncU,nrS,ncS,nrVT,ncVT;
		__lpkdoublecomplex *rawU,*rawVT;
		__lpkdoublereal *rawSappo;

		cout <<"\nCall z_rrsvd\n";
		z_rrsvd_(&(C2Silvi.nrows),&(C2Silvi.ncols),&(C2Silvi.raw),&dimRange,&niter,&tolerance,&nrU,&ncU,&rawU,&nrS,&ncS,&rawSappo,&nrVT,&ncVT,&rawVT);
		cout <<"\nEnd z_rrsvd\n";

		cout<< endl  << nrU << endl;
		cout<< endl  << ncU << endl;
		cout<< endl  << nrS << endl;
		cout<< endl  << ncS << endl;


		RRdecU->nrows=nrU;
		RRdecU->ncols = ncU;
		RRdecU->raw = rawU;
		//
		RRdecS->nrows=nrS;
		RRdecS->ncols = 1;
		RRdecS->raw = (__lpkdoublereal *)rawSappo;
		//
		RRdecVT->nrows=nrVT;
		RRdecVT->ncols = ncVT;
		RRdecVT->raw = rawVT;
		string fileName;

		//Write RRSVD onto file.
		fileName= "RRdecU"+to_string(sample)+".dat";
		cout << endl <<fileName << endl;
		RRdecU->toBinFile(fileName.c_str());

		fileName= "RRdecS"+to_string(sample)+".dat";
		cout << endl <<fileName << endl;
		RRdecS->toBinFile(fileName.c_str());

		fileName= "RRdecVT"+to_string(sample)+".dat";
		cout << endl <<fileName << endl;
		RRdecVT->toBinFile(fileName.c_str());

		//Cleaning up
		delete RRdecU;

		delete RRdecS;

		delete RRdecVT;
	}

 culaShutdown();
 cublasShutdown();
 return 0;
}



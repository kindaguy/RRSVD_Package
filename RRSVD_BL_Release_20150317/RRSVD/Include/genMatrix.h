/*
 * genMatrix.h
 *
 *  Created on: Sep 30, 2014
 *      Author: Tama
 */

#ifndef GENMATRIX_H_
#define GENMATRIX_H_

#include "accelerate_includes.h"
#include <iostream>
#include <fstream>
#include <cstring>
#include "tamaIO.h"


namespace tamaClss {
	
	
	//Class genMatrix declaration
	template<typename T>
	class genMatrix {
	public:
		//Fields:
		int nrows;
		int ncols;
		char format; //'C' Column-Major, 'R Row-Major','D' diagonal
		T *raw;
		
		genMatrix();
		genMatrix(int, int,char f ='C');
		genMatrix(const genMatrix<T> &);
		~genMatrix();
		
		
		//Functionalities:
		//Load matrix from binary file
		void fromBinFile(const char * filename);
		//Store matrix to binary file
		void toBinFile(const char * filename);
		
		//__ATTENTION
		//method fromFile not implemented: we do not load data from text files, yet!
		//void fromFile(const char *filename);
		void toFile(const char *filename);
		
	};
	
	
	
	template <typename T>
	genMatrix<T>::genMatrix(){
		
		nrows = 0;
		ncols = 0;
		format = 'C';
		raw = NULL;
	}
	
	
	template <typename T>
	genMatrix<T>::genMatrix(int r, int c, char f){
		
		nrows = r;
		ncols = c;
		format = f;
		
		raw = new T[nrows*ncols];
		
	}
	
	template <typename T>
	genMatrix<T>::genMatrix(const genMatrix<T> &other){
		nrows = other.nrows;
		ncols = other.ncols;
		format = other.format;
//		if (raw != NULL)
//			delete[]raw;
		raw = new T[nrows*ncols];
		for(int i=0; i<nrows*ncols; i++){
			raw[i]=other.raw[i];
		}
	}
	
	template <typename T>
	genMatrix<T>::~genMatrix(){
		//__DEBUG
		//printf("\nDestructor of complexMatrix invoked\n");
		if (raw !=NULL) {
			delete [] raw;
		}
		
	}
	
	
	template <typename T>
	void genMatrix<T>::toBinFile(const char *filename){
		using namespace std;
		
		ofstream fileOut;
		fileOut.open(filename,ios::binary);
		if(fileOut.fail()){
			printf("\ntoBinFile: Problema apertura file %s\n",filename);
		}
		
		fileOut.write((char*) &(nrows),sizeof(integer));
		fileOut.write((char*) &(ncols),sizeof(integer));
		
		fileOut.write((char*) raw,nrows*ncols*sizeof(T));
		fileOut.close();
	}
	
	template <typename T>
	void genMatrix<T>::fromBinFile(const char *filename){
		
		using namespace std;
		
		ifstream fileIn;
		fileIn.open(filename,ios::binary);
		
		
		//Detect file error
		if (fileIn.fail()){
			cout <<endl<< "fromBinFileProblema apertura file di input: " << filename << endl;
			
		}
		
		
		fileIn.read((char*) &(nrows),sizeof(integer));
		fileIn.read((char*) &(ncols),sizeof(integer));
		
		if(raw==NULL)
			raw=new T[nrows*ncols];
		
		fileIn.read((char*) raw, sizeof(T)*nrows*ncols);
		
		fileIn.close();
	}
	
	template<typename T>
	void genMatrix<T>::toFile(const char *filename){
		using namespace std;
		ofstream fileOut;
		fileOut.open(filename);
		if(fileOut.fail()){
			printf("\ntoFile:problema apertura file %s\n",filename);
		}
		
		
		fileOut << "{";
		for(int i=0; i< nrows; i++){
			fileOut << "{";
			
			for(int j=0; j<ncols; j++){
				fileOut.precision(20);
				tamaIO::TamaWrite(fileOut,raw[j*nrows + i]);
				fileOut << (j<ncols-1?",":"");
			}
			fileOut << "}"<<(i<nrows-1?",":" ");
		}
		
		fileOut <<"}" << endl;
		fileOut.close();	}

	
}//End of namespace


#endif /* GENMATRIX_H_ */

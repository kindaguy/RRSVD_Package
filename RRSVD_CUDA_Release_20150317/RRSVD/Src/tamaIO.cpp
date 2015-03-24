//
//  tamaIO.cpp
//  RRSVD
//
//  Created by Dario Tamascelli on 20/10/14.
//  Copyright (c) 2014 Dario Tamascelli. All rights reserved.
//

#include "tamaIO.h"


void tamaIO::TamaPrint(__lpkdoublereal num){

	printf("%f",num);

}
void tamaIO::TamaPrint(__lpkreal num){
	
	printf("%f",num);
	
}

void tamaIO::TamaPrint(__lpkdoublecomplex num){
	
	printf("%f%cI%f ",num.r,(num.i>=0?'+':'-'),num.i);
	
}

void tamaIO::TamaPrint(__lpkcomplex num){

	printf("%f%cI%f ",num.r,(num.i>=0?'+':'-'),num.i);

}


void tamaIO::TamaWrite(std::ofstream & stream, __lpkdoublereal num){

	stream << num;
}


void tamaIO::TamaWrite(std::ofstream & stream, __lpkreal num){
	
	stream << num;
}

void tamaIO::TamaWrite(std::ofstream & stream, __lpkdoublecomplex num){

	stream << num.r << (num.i >= 0? '+' : '-')  << "I*" << abs(num.i);
}

void tamaIO::TamaWrite(std::ofstream & stream, __lpkcomplex num){
	
	stream << num.r << (num.i >= 0? '+' : '-')  << "I*" << abs(num.i);
}

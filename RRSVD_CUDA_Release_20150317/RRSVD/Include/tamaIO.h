//
//  tamaIO.h
//  RRSVD
//
//  Created by Dario Tamascelli on 20/10/14.
//  Copyright (c) 2014 Dario Tamascelli. All rights reserved.
//

#ifndef RRSVD_tamaIO_h
#define RRSVD_tamaIO_h

#include <iostream>
#include <fstream>
#include "accelerate_includes.h"

namespace tamaIO {
	
	void TamaPrint(__lpkdoublereal);
	void TamaPrint(__lpkdoublecomplex);
	void TamaPrint(__lpkreal);
	void TamaPrint(__lpkcomplex);
	void TamaWrite(std::ofstream &,__lpkdoublereal);
	void TamaWrite(std::ofstream &,__lpkdoublecomplex);
	void TamaWrite(std::ofstream &,__lpkreal);
	void TamaWrite(std::ofstream &,__lpkcomplex);
}

#endif

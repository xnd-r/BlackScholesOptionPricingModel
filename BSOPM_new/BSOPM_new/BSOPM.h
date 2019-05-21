// (c) 2019 Alexander Romanov

#ifndef ____BLACK_SHOLES_OPTION_PRICING_MODELING____
#define ____BLACK_SHOLES_OPTION_PRICING_MODELING____
#include <omp.h>
#include "mkl.h"
#include <math.h>
#include <time.h>
#include <string>
#include <algorithm>
#include <vector>
//#include <stdio.h>
#include <assert.h>
#include <iostream>

#if defined(_WIN64) || defined(_WIN32) 
#define _CRT_SECURE_NO_WARNINGS // using unsafe functions
#pragma warning(disable : 4068) // fopen warning disable
//#pragma warning(disable : 4005) // fopen warning disable
#pragma warning(disable : 4005) // fopen warning disable
#pragma warning(disable : 4996) // REDEFINITION WARNING DISABLED!
#endif

class BSOPM {
public:
	BSOPM() {};
	~BSOPM() {};
/*protected*/
	VSLStreamStatePtr initGen(unsigned int seed, int indexGen/*, int dim = 1*/);
	void wienerProcess(VSLStreamStatePtr stream, int nsteps, float time, float *dw);
	void wAndZProcesses(VSLStreamStatePtr stream, int nsteps, float time, float *buffer);
	std::string getTimestamp(std::string fileName);
};

#endif // !____BLACK_SHOLES_OPTION_PRICING_MODELING____

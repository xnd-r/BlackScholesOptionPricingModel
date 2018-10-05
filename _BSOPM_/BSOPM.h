// Copyright 2018 Romanov Alexander

#ifndef ____BLACK_SHOLES_OPTION_PRICING_MODELING____
#define ____BLACK_SHOLES_OPTION_PRICING_MODELING____
#include <omp.h>
#include <mkl.h>
#include <math.h>
#include <time.h>
#include <string>
#include <algorithm>
#include <vector>
#include <stdio.h>

#if defined(_WIN64) || defined(_WIN32) 
#define _CRT_SECURE_NO_WARNINGS // using unsafe functions
#pragma warning(disable : 4005) // fopen warning disable
#pragma warning(disable : 4996) // REDEFINITION WARNING DISABLED!
#endif

#define __SEED__	20000000l

class BSOPM {
public:
	BSOPM() {};
	~BSOPM() {};
	
	float getStockPrice(float S0, float R, float Sig, float dw, float T); // Analytical S(T)
protected:
	VSLStreamStatePtr initGen();
	void freeGen(VSLStreamStatePtr stream);
	void normalGenerator(float mean, float deviation, int amou, VSLStreamStatePtr stream, float *destArray);
	void wienerProcess(VSLStreamStatePtr stream, int nsteps, float time, float *dw);
};

#endif // !____BLACK_SHOLES_OPTION_PRICING_MODELING____

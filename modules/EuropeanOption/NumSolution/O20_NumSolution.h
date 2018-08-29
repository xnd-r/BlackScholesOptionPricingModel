#ifndef ____NUMERICAL_SOLUTION_OPTION____
#define ____NUMERICAL_SOLUTION_OPTION____
#include "../include/O00_EuropeanOption.h"
#include <assert.h>
#include <vector>
#include <omp.h>
#include <iostream>
#include <algorithm>

class NumSolutionOption : public EuropeanOption {

public:
	const unsigned int bufsize = 5000;
	float tmp1 = (R - SIG * SIG * 0.5f) * TIME;
	float tmp2 = SIG * sqrtf(TIME);
	clock_t start, finish;
	double t;
	float sum = 0.0f; 
	float s; 
	std::vector<double> Times;
	std::vector<double> Prices;
	
	//virtual float GetMCPrice();
};
#endif // !____NUMERICAL_SOLUTION_OPTION____
#ifndef ____NUMERICAL_SOLUTION_OPTION____
#define ____NUMERICAL_SOLUTION_OPTION____
#include "../include/O00_EuropeanOption.h"
#include <assert.h>

class NumSolutionOption : public EuropeanOption {
public:
	float tmp1 = (R - SIG * SIG * 0.5f) * TIME;
	float tmp2 = SIG * sqrtf(TIME);
	clock_t start, finish;
	double t;
	float sum = 0.0f; 
	float s; 
	
	//virtual float GetMCPrice();
};
#endif // !____NUMERICAL_SOLUTION_OPTION____
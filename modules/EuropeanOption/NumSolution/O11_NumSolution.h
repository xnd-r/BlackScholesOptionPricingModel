#ifndef ____NUMERICAL_SOLUTION_OPTION____
#define ____NUMERICAL_SOLUTION_OPTION____
#include "../include/O00_EuropeanOption.h"
#include <assert.h>

class NumSolutionOption : public EuropeanOption {
	float tmp1 = (R - SIG * SIG * 0.5f) * TIME;
	float tmp2 = SIG * sqrtf(TIME);
	const unsigned int bufsize = 1000;
	const unsigned int seed[2] = { __SEED__, __SEED__ };
public:
	clock_t start, finish;
	double t;
	float sum = 0.0f; 
	float s; 
	
	float GetMCPrice(int indexGen);
};
#endif // !____NUMERICAL_SOLUTION_OPTION____
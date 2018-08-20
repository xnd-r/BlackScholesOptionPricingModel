#ifndef ____NUMERICAL_SOLUTION_OPTION____
#define ____NUMERICAL_SOLUTION_OPTION____
#include "../include/O00_EuropeanOption.h"
#include <assert.h>

class NumSolutionOption : public EuropeanOption {
public:
	clock_t start, finish;
	double t;
	float sum = 0.0f; 
	float s; 
	
	virtual float GetMCPrice();
};
#endif // !____NUMERICAL_SOLUTION_OPTION____
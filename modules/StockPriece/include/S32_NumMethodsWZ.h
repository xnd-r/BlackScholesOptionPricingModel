#ifndef ____NUMERICAL_METHODS_W_AND_Z____
#define ____NUMERICAL_METHODS_W_AND_Z____

#include "S30_NumMethods.h"

class  NumMethodWZ : public NumMethods {
public:
	typedef double(NumMethodWZ::*Step) (double, double, double, double);

	NumMethodWZ::Step step_array[2] = { &NumMethodWZ::BurragePlatenStep, &NumMethodWZ::Taylor2Step };
	
	double BurragePlatenStep(double S, double dt, double dw, double dz);
	double Taylor2Step(double S, double dt, double dw, double dz);

	void SimulateStockPrices(Step _step, VSLStreamStatePtr stream, int nPaths, int nSteps, double Time, double *Error);

	void Execute(Step _step);
};

#endif // !____NUMERICAL_METHODS_W_AND_Z____
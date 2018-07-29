#ifndef ____NUMERICAL_METHODS_W____
#define ____NUMERICAL_METHODS_W____

#include "S30_NumMethods.h"


class  NumMethodW : public NumMethods {
public:
	typedef double(NumMethodW::*Step) (double, double, double);

	//Step step_array[] = { EulMarStep, MilsteinStep, RK1Step };

	double EulMarStep(double S, double dt, double dw);
	double MilsteinStep(double S, double dt, double dw);
	double RK1Step(double S, double dt, double dw);

	void SimulateStockPrices(Step _step, VSLStreamStatePtr stream, int nPaths, int nSteps, double Time, double *Error);
};

#endif // !____NUMERICAL_METHODS_W____

#ifndef ____NUMERICAL_METHODS_W____
#define ____NUMERICAL_METHODS_W____

#include "S30_NumMethods.h"

class  NumMethodW : public NumMethods {
public:
	typedef double(NumMethodW::*Step) (double, double, double);

	NumMethodW::Step step_array[3] = { &NumMethodW::EulMarStep, &NumMethodW::MilsteinStep, &NumMethodW::RK1Step };

	//double (NumMethodW::*StepArray[3])(double, double, double) = { &NumMethodW::EulMarStep, &NumMethodW::MilsteinStep, &NumMethodW::RK1Step };
	//Step step_array[3] = { this->EulMarStep, *this->MilsteinStep, *this->RK1Step };

	double EulMarStep(double S, double dt, double dw);
	double MilsteinStep(double S, double dt, double dw);
	double RK1Step(double S, double dt, double dw);

	void SimulateStockPrices(Step _step, VSLStreamStatePtr stream, int nPaths, int nSteps, double Time, double *Error);

	void Execute(Step _step);
};

#endif // !____NUMERICAL_METHODS_W____

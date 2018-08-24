#ifndef ____NUMERICAL_METHODS_W____
#define ____NUMERICAL_METHODS_W____

#include "S30_NumMethods.h"

class  NumMethodW : public NumMethods {
public:
	typedef float(NumMethodW::*Step) (float, float, float);

	NumMethodW::Step StepArray[3] = { &NumMethodW::EulMarStep, &NumMethodW::MilsteinStep, &NumMethodW::RK1Step };

	float EulMarStep(float S, float dt, float dw);
	float MilsteinStep(float S, float dt, float dw);
	float RK1Step(float S, float dt, float dw);

	__declspec(noinline) void SimulateStockPrices(Step _step, VSLStreamStatePtr stream, int nPaths, int nSteps, float Time, float *Error);

	void Execute(Step _step, char* FileName);
};

#endif // !____NUMERICAL_METHODS_W____

#ifndef ____NUMERICAL_METHODS_W_AND_Z____
#define ____NUMERICAL_METHODS_W_AND_Z____

#include "S30_NumMethods.h"

class  NumMethodWZ : public NumMethods {
public:
	typedef float(NumMethodWZ::*Step) (float, float, float, float);

	NumMethodWZ::Step step_array[2] = { &NumMethodWZ::BurragePlatenStep, &NumMethodWZ::Taylor2Step };
	
	float BurragePlatenStep(float S, float dt, float dw, float dz);
	float Taylor2Step(float S, float dt, float dw, float dz);

	__declspec(noinline) void SimulateWandZProcesses(VSLStreamStatePtr stream, int nSteps, float Time, float *buffer);
	__declspec(noinline) void SimulateStockPrices(Step _step, VSLStreamStatePtr stream, int nPaths, int nSteps, float Time, float *Error);

	void Execute(Step _step, char* FileName);
};

#endif // !____NUMERICAL_METHODS_W_AND_Z____
#ifndef ____NUMERICAL_METHODS_W_AND_Z____
#define ____NUMERICAL_METHODS_W_AND_Z____

#include "S30_NumMethods.h"

class  NumMethodWZ : public NumMethods {
public:
	virtual double Step(double S, double dt, double dw, double dz) = 0;
	void SimulateStockPricesNumAn(VSLStreamStatePtr stream, int nPaths, int nSteps, double Time, double *Error);
	void SimulateWandZProcesses(VSLStreamStatePtr stream, int nSteps, double Time, double *buffer);
};

#endif // !____NUMERICAL_METHODS_W_AND_Z____
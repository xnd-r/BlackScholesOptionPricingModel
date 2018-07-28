#ifndef ____NUMERICAL_METHODS_W____
#define ____NUMERICAL_METHODS_W____

#include "S30_NumMethods.h"

class  NumMethodW : public NumMethods {
public:
	virtual double Step(double S, double dt, double dw, double dz) = 0;
	void SimulateStockPricesNumAn(VSLStreamStatePtr stream, int nPaths, int nSteps, double Time, double *Error);
	void SimulateWandZProcesses(VSLStreamStatePtr stream, int nSteps, double Time, double *buffer);
};

#endif // !____NUMERICAL_METHODS_W____
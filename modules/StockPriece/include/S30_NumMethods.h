#ifndef ____NUMERICAL_METHODS_____
#define ____NUMERICAL_METHODS_____

#include "S21_AnExtended.h"
#include "../../../etc/RowTable/RowTable.h"

class  NumMethods : public AnExtended {
#define NPATHS 1000
#define NSTEPS 512
public:
//	virtual void SimulateStockPrices(VSLStreamStatePtr stream, int nPaths, int nSteps, float Time, float *Error);
	void WriteToCsv(float* Errors, int nSteps, int nRows, float Time, int scale, char* FileName);
};

#endif // !____NUMERICAL_METHODS_____
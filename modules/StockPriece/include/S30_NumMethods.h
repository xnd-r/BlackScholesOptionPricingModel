#ifndef ____NUMERICAL_METHODS_____
#define ____NUMERICAL_METHODS_____

#include "S21_AnExtended.h"
#include "../../../etc/RowTable/RowTable.h"

class  NumMethods : public AnExtended {
#define NPATHS 1000
#define NSTEPS 512
public:
	virtual void SimulateStockPrices() = 0;
	void WriteToCsv(double* Errors, int nSteps, int nRows, double Time, int scale, char* FileName);
};

#endif // !____NUMERICAL_METHODS_____
#ifndef ____SDE_ANALITYCAL_EQUATION_EXTENDED____
#define ____SDE_ANALITYCAL_EQUATION_EXTENDED____

#include "S20_AnSimple.h"

class AnExtended : public AnSimple {
#define NPATHS 300
public:
	void SimulateWienerProcess(VSLStreamStatePtr stream, int nSteps, double Time, double *buffer);
	void SimulateStockPrices(VSLStreamStatePtr stream, int nPaths, int nSteps, double Time, double **sBuffer);
	void WriteToCsv(double **buffer, int nRows, int nColumns);
	void Execute();
};

#endif // !____SDE_ANALITYCAL_EQUATION_EXTENDED____
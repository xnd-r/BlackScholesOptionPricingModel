#ifndef ____SDE_ANALITYCAL_EQUATION_EXTENDED____
#define ____SDE_ANALITYCAL_EQUATION_EXTENDED____

#include "S20_AnSimple.h"

class AnExtended : public AnSimple {
#define NPATHS 500
public:
	// InitGen() and FreeGen should be added in main
	void SimulateWienerProcess(VSLStreamStatePtr stream, int nSteps, double Time, double *buffer);
	void SimulateStockPricesExt(VSLStreamStatePtr stream, int nPaths, int nSteps, double Time, double **sBuffer);
	void WriteToCsv(double **buffer, int nRows, int nColumns, int partRows);
};

#endif // !____SDE_ANALITYCAL_EQUATION_EXTENDED____
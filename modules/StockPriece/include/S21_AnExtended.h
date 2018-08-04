#ifndef ____SDE_ANALITYCAL_EQUATION_EXTENDED____
#define ____SDE_ANALITYCAL_EQUATION_EXTENDED____

#include "S20_AnSimple.h"

class AnExtended : public AnSimple {
#define NPATHS 300
public:

	__declspec(noinline) void SimulateWienerProcess(VSLStreamStatePtr stream, int nSteps, float Time, float *buffer);
	__declspec(noinline) void SimulateStockPrices(VSLStreamStatePtr stream, int nPaths, int nSteps, float Time, float **sBuffer);
	void WriteToCsv(float **buffer, int nRows, int nColumns);
	void Execute();
};

#endif // !____SDE_ANALITYCAL_EQUATION_EXTENDED____
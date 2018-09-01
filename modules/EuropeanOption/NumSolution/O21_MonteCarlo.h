#ifndef ____MONTE_CARLO_METHOD____
#define ____MONTE_CARLO_METHOD____
#include "../NumSolution/O20_NumSolution.h"
#include <assert.h>
#include <vector>
#include <omp.h>

class MonteCarlo : public NumSolutionOption {

	const unsigned int seed[2] = { __SEED__, __SEED__ };
public:
	float GetMCPrice(int indexGen);
	float GetMCPrice(int indexGen, int NumThreads);
	void Execute();
};
#endif // !____MONTE_CARLO_METHOD____
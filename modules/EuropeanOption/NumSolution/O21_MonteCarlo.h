#ifndef ____MONTE_CARLO_METHOD____
#define ____MONTE_CARLO_METHOD____
#include "../NumSolution/O20_NumSolution.h"
#include <assert.h>
#include <omp.h>

class MonteCarlo : public NumSolutionOption {
	float tmp1 = (R - SIG * SIG * 0.5f) * TIME;
	float tmp2 = SIG * sqrtf(TIME);
	const unsigned int bufsize = 1000;
	const unsigned int seed[2] = { __SEED__, __SEED__ };
public:

	float GetMCPrice(int indexGen);
	float GetMCPrice(int indexGen, int NumThreads);
};
#endif // !____MONTE_CARLO_METHOD____
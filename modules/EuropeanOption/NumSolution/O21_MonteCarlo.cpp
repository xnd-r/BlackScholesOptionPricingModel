#include "O21_MonteCarlo.h"
#include "../RNG/ER21_Normal.h"

float MonteCarlo::GetMCPrice(int indexGen) {
	assert(N % bufsize == 0);
	start = clock();
	float *gauss = new float[bufsize];

	VSLStreamStatePtr stream;

	if (indexGen == 0) {
		vslNewStreamEx(&stream, VSL_BRNG_MCG59, 2, seed);
	}
	else if (indexGen == 1) {
		vslNewStream(&stream, VSL_BRNG_SOBOL, 1);
	}

	//if (indexGen == 0) {
	//	//VSLStreamStatePtr stream; 
	//	//vslNewStreamEx(&stream, VSL_BRNG_MCG59, 2, seed); 
	//	NormalGen ng(0.0f, 1.0f, __SEED__);
	//}
	//else if (indexGen == 1) {
	//	VSLStreamStatePtr stream;
	//	vslNewStream(&stream, VSL_BRNG_SOBOL, 1);
	//}

	for (unsigned int portion = 0; portion < N / bufsize; portion++) {

		vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, bufsize, gauss, 0.0, 1.0);
		for (int i = 0; i < bufsize; i++) {
			float payoff;
			s = S0 * expf(tmp1 + tmp2 * gauss[i]);
			payoff = s - K;
			if (payoff > 0.0)
				sum = sum + payoff;
		}
	}
	sum = sum / N * exp(-R * TIME);
	//vslDeleteStream(&stream); 
	delete[] gauss;
	finish = clock();
	t = (double)(finish - start) / CLOCKS_PER_SEC;
	return sum;
}

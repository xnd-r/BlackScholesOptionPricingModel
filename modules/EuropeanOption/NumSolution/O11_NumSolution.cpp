#include "O11_NumSolution.h"

float NumSolutionOption::GetMCPrice() {
	assert(N % bufsize == 0); 
	start = clock(); 
	double *gauss = new double[bufsize];
	VSLStreamStatePtr stream; 
	vslNewStreamEx(&stream, VSL_BRNG_MCG59, 2, seed); 
	for (unsigned int portion = 0; portion < N / bufsize; portion++) { 
		vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, bufsize, gauss, 0.0, 1.0);
		for (int i = 0; i < bufsize; i++) { 
			float payoff; 
			s = S0 * exp(tmp1 + tmp2 * gauss[i]); 
			payoff = s - K; 
			if (payoff > 0.0) 
				sum = sum + payoff; 
		} 
	} 
	sum = sum / N * exp(-R * TIME); 
	vslDeleteStream(&stream); 
	delete [] gauss; 
	finish = clock(); 
	t = (double)(finish - start) / CLOCKS_PER_SEC; 
	return sum; 
}
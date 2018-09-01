#include "O21_MonteCarlo.h"
#include "../RNG/ER21_Normal.h"
//#include "../../../etc/RNG/ER21_Normal.h"

// TODO: O21_MonteCarlo.cpp(50): warning C4101: 's': unreferenced local variable


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

	for (unsigned int portion = 0; portion < N / (unsigned int)bufsize; portion++) {

		vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, bufsize, gauss, 0.0, 1.0);
		for (int i = 0; i < bufsize; i++) {
			float payoff;
			s = S0 * expf(tmp1 + tmp2 * gauss[i]);
			payoff = s - K;
			if (payoff > 0.0f)
				sum = sum + payoff;
		}
	}
	sum = sum / N * expf(-R * TIME);
	vslDeleteStream(&stream); 
	delete[] gauss;
	finish = clock();
	t = (double)(finish - start) / CLOCKS_PER_SEC;
	return sum;
}

float MonteCarlo::GetMCPrice(int indexGen, int NumThreads) {
	float s;
	float sum = 0.0f;
	assert(N % bufsize == 0);
	start = clock();
#pragma omp parallel private(s)
	{
		int count = omp_get_num_threads();
		int num = omp_get_thread_num();

		float *gauss = new float[bufsize];
		VSLStreamStatePtr stream;

		if (indexGen == 0) {
			vslNewStreamEx(&stream, VSL_BRNG_MCG59, 2, seed);
		}
		else if (indexGen == 1) {
			vslNewStream(&stream, VSL_BRNG_SOBOL, 1);
		}
		vslSkipAheadStream(stream, N / count * num);

#pragma omp for reduction(+:sum)
		for (int portion = 0; portion < N / bufsize; portion++) {

			vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, bufsize, gauss, 0.0, 1.0);
			for (int i = 0; i < bufsize; i++) {
				float payoff;
				s = S0 * expf(tmp1 + tmp2 * gauss[i]);
				payoff = s - K;
				if (payoff > 0.0)
					sum = sum + payoff;
			}
		}
		vslDeleteStream(&stream);
		delete[] gauss;
	}
	sum = sum / N * exp(-R * TIME);
	finish = clock();
	t = (double)(finish - start) / CLOCKS_PER_SEC;
	return sum;
}

void MonteCarlo::Execute() {
	int tmp = 1;
	double t1, t2;
	for (int k = 1; k < 4; ++k) {
		omp_set_num_threads(tmp);
		for (int j = 0; j < 10; ++j) {
			t1 = omp_get_wtime();
				GetMCPrice(1, tmp);
			t2 = omp_get_wtime();
			Times.push_back(t2 - t1);
			Prices.push_back(GetMCPrice(1, tmp));
		}
		std::sort(Times.begin(), Times.end());
		for (std::vector<double>::const_iterator it = Times.begin(); it != Times.end(); ++it)
			std::cout << *it << std::endl;
		std::cout << " ####################### " << std::endl;
		for (std::vector<double>::const_iterator it = Prices.begin(); it != Prices.end(); ++it)
			std::cout << *it << std::endl;
		std::cout << tmp << " thread(s) done" << std::endl;
		tmp *= 2;
		Times.clear();
		Prices.clear();
	}
}
#include <iostream>
#include <cmath>
#define M_PIF 3.1415926
#include "omp.h"
#include "MDim.h"
#include "NumSolution.h"
#include "BSOPM.h"

float getMDimDensity(int dim, float*x) {
	float denom = 0.f;
	for (int i = 0; i < dim; ++i) {
		denom += x[i] * x[i];
	}
	return expf(-denom / 2.f);
}

void getPointCoords(float *tmpVector, long long ind, int nsamples, int dim, float a, float b)
{
	long long *CountInBlock = new long long[dim];

	CountInBlock[0] = static_cast<long long>(powl(nsamples, dim - 1));
	float h = (b - a) / nsamples;
	for (int i = 1; i < dim; ++i)
		CountInBlock[i] = CountInBlock[i - 1] / static_cast<long long>(nsamples);

	for (int i = 0; i < dim; i++) {
		tmpVector[i] = a + h * (static_cast<float>(static_cast<int>(ind / CountInBlock[i]) % nsamples) + 0.5f);
	}
	delete[] CountInBlock;
}

float getRRNG(int dim, int scale, float a, float b, float r, float time, float* SP, float S0, float K) {

	float tmp1 = (r - 0.04f * .5f) * time;
	float tmp2 = .2f * sqrtf(time);
	long long amo_of_points = static_cast<long long>(powl(scale, dim));
	float* mDimPoint = new float[dim];
	float maxSP;
	long double sum = 0.0;
	double t1 = omp_get_wtime();
	omp_set_num_threads(dim);
	{
		for (long long i = 0; i < amo_of_points; i++) {
			getPointCoords(mDimPoint, i, scale, dim, a, b);
			for (int j = 0; j < dim; ++j) {
				S0 * expf(tmp1 + tmp2 * mDimPoint[j]) - K < 0.f ? SP[j] = 0.f : SP[j] = S0 * expf(tmp1 + tmp2 * mDimPoint[j]) - K;
			}
			maxSP = SP[0];
			//std::cout << mDimPoint[0] << "\t" << SP[0] << std::endl;
			for (int i = 1; i < dim; ++i) {
				if (SP[i] > maxSP)
					maxSP = SP[i];
			}
			maxSP *= getMDimDensity(dim, mDimPoint);
			sum += maxSP;
		}
	}
	float det = 1.f;
	det = det * powf(2.f * M_PIF, dim);
	sum = sum * static_cast<long double>((1 / sqrtf(det)) * (expf(-r * time) * powf((b - a) / scale, dim)) / dim);
	delete[] mDimPoint;
	return static_cast<float>(sum);
}

float getMCPriceParNew(int NumThreads, int StepIndex, int nsteps, int indexGen, int N, int bufsize, unsigned int seed, float K, float R, float Time, float SIG, float pS0, int dim) {
	
	VSLStreamStatePtr stream;
	double t1, t2;
	omp_set_num_threads(NumThreads);
	float sum = .0f, stockPrice;
	t1 = omp_get_wtime();	
	std::cout << "Alive" << std::endl;

	//#pragma omp parallel private(stockPrice)
	{
		//int count = omp_get_num_threads();
		//int num = omp_get_thread_num();
		float payoff, maxSP;
		float dt = Time / nsteps;
		float* wtraject = new float[dim * nsteps];
		float* mean = new float[dim] { 0 };
		float* SP = new float[dim];
		for (int i = 0; i < dim; ++i)
			SP[i] = 100.f;
		float* cov = new float[dim * dim];
		for (int i = 0; i < dim * dim; ++i)
			cov[i] = 0.f;
		for (int i = 0; i < dim; ++i)
			cov[dim * i + i] = sqrtf(dt);

		if (indexGen == 0)
		{
			const unsigned int _seed[2] = { seed, seed };
			vslNewStreamEx(&stream, VSL_BRNG_MCG59, 2, _seed);
		}
		else
		{
			vslNewStream(&stream, VSL_BRNG_SOBOL, dim);
		}

		//vslSkipAheadStream(stream, dim * nsteps / count * num);

		//#pragma omp for reduction(+:sum)
		for (int ind = 0; ind < N; ++ind) {
			vsRngGaussianMV(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, nsteps, wtraject, dim, VSL_MATRIX_STORAGE_FULL, mean, cov);
			for (int k = 0; k < nsteps; ++k) {
				for (int j = 0; j < dim; ++j) {
					SP[j] = SP[j] + SP[j] * (R * dt + SIG * wtraject[dim * k + j]);
					//std::cout << wtraject[dim * k + j] << "\t" << SP[j] << std::endl;

				}
			}
			maxSP = 0.f;
			for (int i = 0; i < dim; ++i) {
				if (SP[i] - 100.f > maxSP)
					maxSP = SP[i] - 100.f;
			}
			sum += maxSP;

			for (int i = 0; i < dim; ++i) {
				SP[i] = 100.f;
			}
		}
		delete[] wtraject;
		delete[] mean;
		delete[] SP;
		delete[] cov;
	}
	t2 = omp_get_wtime();
	//std::cout << "Time = " << t2 - t1 << std::endl;
	sum = sum / N * expf(-R * Time) / dim;
	//std::cout << "MaxCall = " << sum << std::endl;

	vslDeleteStream(&stream);
	return sum;
}
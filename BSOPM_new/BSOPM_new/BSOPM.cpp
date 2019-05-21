#include "BSOPM.h"
#include <iostream>
void BSOPM::wienerProcess(VSLStreamStatePtr stream, int nsteps, float time, float *w) {
	// only for convergence
	float *gaussBuf = new float[nsteps]; // Random values buffer
	float dt = time / nsteps;
	vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, nsteps, gaussBuf, 0.f, sqrtf(dt));
	w[0] = .0f;
	for (int j = 1; j <= nsteps; ++j)
		w[j] = w[j-1] + gaussBuf[j - 1];
	delete[] gaussBuf;
}

void BSOPM::wAndZProcesses(VSLStreamStatePtr stream, int nsteps, float time, float *buffer) {
	// only for convergence
	float *dw = new float[nsteps * 2];
	float dt = time / nsteps;
	float mean[2] = { 0.f, 0.f };
	float hh = dt * sqrtf(dt);
	float cov[3] = { sqrtf(dt) , 0.5f * hh , 1.0f / sqrtf(12.0f) * hh };
	vsRngGaussianMV(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, nsteps, dw, 2, VSL_MATRIX_STORAGE_PACKED, mean, cov);
	buffer[0] = 0; buffer[1] = 0;
	for (int j = 2; j <= nsteps * 2; j += 2)
	{
		buffer[j]	  = buffer[j - 2] + dw[j - 2];
		buffer[j + 1] = buffer[j - 1] + dw[j - 1];
	}
	delete[] dw;
}


VSLStreamStatePtr BSOPM::initGen(unsigned int seed, int indexGen/*, int dim = 1*/) {
	VSLStreamStatePtr stream;
	if (indexGen == 0) {
		const unsigned int _seed[2] = { seed, seed };
		vslNewStreamEx(&stream, VSL_BRNG_MCG59, 2, _seed);
	}
	else if (indexGen == 1) {
		vslNewStream(&stream, VSL_BRNG_SOBOL, 1);
	}
	return stream;
}


std::string BSOPM::getTimestamp(std::string fileName) {
	time_t timestamp;
	time(&timestamp);
	std::string date = fileName.append(asctime(localtime(&timestamp)));
	date.pop_back();
	for (std::string::iterator it = date.begin(); it<date.end(); ++it) {
		if (*it == ':') {
			date.erase(it);
		}
		std::replace(date.begin(), date.end(), ' ', '_');
	}
	date.append(".csv");
	return date;
}

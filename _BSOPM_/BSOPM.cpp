#include "BSOPM.h"
#include <iostream>
void BSOPM::wienerProcess(VSLStreamStatePtr stream, int nsteps, float time, float *w) {
	float *gaussBuf = new float[nsteps]; // Random values buffer
	float dt = time / nsteps;
	normalGenerator(.0f, sqrtf(dt), nsteps, stream, gaussBuf);
	w[0] = .0f;
	for (int j = 1; j <= nsteps; ++j)
		w[j] = /*w[j-1] + */gaussBuf[j - 1];
	delete[] gaussBuf;
}

float BSOPM::getStockPrice(float s0, float r, float sig, float dw, float t) {
	return s0 * expf((r - sig * sig / 2.f) * t + sig * dw);
}

void BSOPM::normalGenerator(float mean, float deviation, int amou, VSLStreamStatePtr stream, float *destArray) {
	//Getting amount random values and writing them to the destination array
	vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, amou, destArray, mean, deviation);
}

VSLStreamStatePtr BSOPM::initGen(unsigned int seed, int indexGen) {
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

void BSOPM::freeGen(VSLStreamStatePtr stream) {
	vslDeleteStream(&stream);
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

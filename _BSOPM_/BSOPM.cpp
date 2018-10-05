#include "BSOPM.h"

void BSOPM::wienerProcess(VSLStreamStatePtr stream, int nsteps, float time, float *dw) {
	float *gaussBuf = new float[nsteps]; // Random values buffer
	float dt = time / nsteps; 
	// getting nsteps random values with N(0, dt):
	normalGenerator(.0f, sqrtf(dt), nsteps, stream, gaussBuf);
	// getting Wiener differentials
	dw[0] = .0f;
	for (int j = 1; j <= nsteps; j++)
		dw[j] = dw[j - 1] + gaussBuf[j - 1];
	delete[] gaussBuf;
}

float BSOPM::getStockPrice(float s0, float r, float sig, float dw, float t) {
		return s0 * expf((r - sig * sig / 2.f) * t + sig * dw);
}

void BSOPM::normalGenerator(float mean, float deviation, int amou, VSLStreamStatePtr stream, float *destArray) {
	//Getting amount random values and writing them to the destination array
	vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, amou, destArray, mean, deviation);
}

VSLStreamStatePtr BSOPM::initGen() {
	VSLStreamStatePtr stream;
	const unsigned int seed[2] = { __SEED__, __SEED__ };
	vslNewStreamEx(&stream, VSL_BRNG_MCG59, 2, seed); // base RNG 
	// TODO: add another generator
	return stream;
}

void BSOPM::freeGen(VSLStreamStatePtr stream) { 
	vslDeleteStream(&stream);
}

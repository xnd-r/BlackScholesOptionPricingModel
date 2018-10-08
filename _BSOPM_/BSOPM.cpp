#include "BSOPM.h"
#include <iostream>
void BSOPM::wienerProcess(VSLStreamStatePtr stream, int nsteps, float time, float *dw) {
	float *gaussBuf = new float[nsteps]; // Random values buffer
	float dt = time / nsteps; 
	// getting nsteps random values with N(0, dt):
	normalGenerator(.0f, sqrtf(dt), nsteps, stream, gaussBuf);
	//for (int i = 0; i < nsteps; ++i) {
	//	std::cout << gaussBuf[i] << "\n";
	//}
	// getting Wiener differentials
	dw[0] = .0f;

	//std::string date = "_Distribution.csv";
	//FILE *f = fopen(date.c_str(), "w");
	for (int j = 1; j <= nsteps; j++) {
		dw[j] = dw[j - 1] + gaussBuf[j - 1];
		//fprintf(f, "%lf;\n", gaussBuf[j - 1]);
	}
	//fclose(f);
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

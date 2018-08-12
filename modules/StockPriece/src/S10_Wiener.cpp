#include "../include/S10_Wiener.h"

VSLStreamStatePtr StockPrice::InitGen() {
	VSLStreamStatePtr stream;
	const unsigned int seed[2] = { __SEED__, __SEED__ };
	vslNewStreamEx(&stream, VSL_BRNG_MCG59, 2, seed); // base RNG 
	return stream;
}

void StockPrice::FreeGen(VSLStreamStatePtr stream) { // deleting of datastructure of gen
	vslDeleteStream(&stream); // is it really needed to write one-row func?
}

__declspec(noinline) void StockPrice::GenerateGauss(float expect_val, float deviation, int amou,
	VSLStreamStatePtr stream, float *dest_array) {
	//Getting amount random values and writing them to the destination array
	vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, amou, dest_array, expect_val, deviation);
}


__declspec(noinline) void StockPrice::SimulateWienerProcess(int nPaths, int nSteps, float Time, float **buffer) {

	VSLStreamStatePtr stream = InitGen();

	float *wiener_diff = new float[nSteps]; // Random values buffer
	float h = Time / (float)nSteps; // step

	#if defined(__INTEL_COMPILER) 
		#pragma simd
		#pragma vector always	
	#endif

	for (int i = 0; i < nPaths; i++) {
		// getting nSteps random values with N(0, h)
		GenerateGauss(0, sqrtf(h), nSteps, stream, wiener_diff);
		// TODO: Add N(0, __STEP__) correcthness
		buffer[i][0] = 0;
		for (int j = 1; j <= nSteps; j++) {
			buffer[i][j] = buffer[i][j - 1] + wiener_diff[j - 1];
		}
	}
	delete[] wiener_diff;
	FreeGen(stream); // Generator data sturcture deleting
}

void StockPrice::WriteToCsv(float **buffer, int nRows, int nColumns) {

	time_t rawtime;
	time(&rawtime);
	std::string date = asctime(localtime(&rawtime));
	date.pop_back();
	date.append("_Wiener_Trajectories_Results.csv");

	for (std::string::iterator it = date.begin(); it<date.end(); ++it) {
		if (*it == ':') {
			date.erase(it);
		}
		std::replace(date.begin(), date.end(), ' ', '_');
	}

	FILE *f = fopen(date.c_str(), "w");
	for (int i = 0; i < nRows; i++)
	{
		for (int j = 0; j < nColumns; j++) {
			//std::string tmp_cell = std::to_string(buffer[i][j]);

			//for (std::string::iterator it = tmp_cell.begin(); it<tmp_cell.end(); ++it) {
			//	std::replace(tmp_cell.begin(), tmp_cell.end(), '.', ',');
			//}
			fprintf(f, "%lf;"/*"%s;"*/, buffer[i][j]/*tmp_cell.c_str()*/);
		}
		fprintf(f, "\n");
	}
	fclose(f);
}
void StockPrice::Execute() {
	float **buffer = new float*[NPATHS];
	for (int i = 0; i < NPATHS; i++)
		buffer[i] = new float[NSTEPS + 1];
	SimulateWienerProcess(NPATHS, NSTEPS, TIME, buffer);
	WriteToCsv(buffer, NPATHS, NSTEPS + 1);
	for (int i = 0; i < NPATHS; i++)
		delete[] buffer[i];
	delete[] buffer;
}
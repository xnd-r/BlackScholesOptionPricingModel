#include "../include/S21_AnExtended.h"

__declspec(noinline) void AnExtended::SimulateWienerProcess(VSLStreamStatePtr stream, int nSteps, float Time, float *buffer) {

	float *wiener_diff = new float[nSteps]; // Random values buffer
	float h = Time / (float)nSteps; // step
	GenerateGauss(0, sqrt(h), nSteps, stream, wiener_diff); 
	buffer[0] = 0.0f;
	// getting nSteps random values with N(0, __STEP__)
	for (int j = 1; j <= nSteps; j++) {
		buffer[j] = buffer[j - 1] + wiener_diff[j - 1];
	}

	delete[] wiener_diff;
}

__declspec(noinline) void AnExtended::SimulateStockPrices(VSLStreamStatePtr stream, int nPaths, int nSteps, float Time, float **sBuffer) {
	// simulating stock priese according to an intermediate values of segment 
	float *wiener_diff = new float[nSteps + 1];
	float h = Time / (float)nSteps; // step
	for (int i = 0; i < nPaths; i++) {
		// Wiener process trajectory simulating
		SimulateWienerProcess(stream, nSteps, Time, wiener_diff);
		sBuffer[i][0] = S0;
		for (int j = 1; j <= nSteps; j++) {
			sBuffer[i][j] = GetStockPrice(wiener_diff[j], h * j);
		}
	}
}

void AnExtended::WriteToCsv(float **buffer, int nRows, int nColumns) {

	time_t rawtime;
	time(&rawtime);
	std::string date = asctime(localtime(&rawtime));
	date.pop_back();
	date.append("_Analitycal_Extended_Results.csv.csv");

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
	float tmp_avg;

	#if defined(__INTEL_COMPILER) 
		#pragma simd
		#pragma vector always	
	#endif

	for (int i = 0; i < nColumns; i++) {
		tmp_avg = 0.0;
		for (int j = 0; j < nRows; j++) {
			tmp_avg += buffer[j][i];
		}
		tmp_avg /= nRows;
		std::string tmp_cell = std::to_string(tmp_avg);
		for (std::string::iterator it = tmp_cell.begin(); it < tmp_cell.end(); ++it) {
			std::replace(tmp_cell.begin(), tmp_cell.end(), '.', ',');
		}
		fprintf(f, "%s;", tmp_cell.c_str());
	}
	fprintf(f, "\n");
	fclose(f);
}

void AnExtended::Execute() {
	VSLStreamStatePtr stream = InitGen();
	float **sBuffer = new float*[NPATHS];
	for (int i = 0; i < NPATHS; i++)
		sBuffer[i] = new float[NSTEPS + 1];
	SimulateStockPrices(stream, NPATHS, NSTEPS, TIME,
		sBuffer);
	WriteToCsv(sBuffer, NPATHS, NSTEPS + 1);
	for (int i = 0; i < NPATHS; i++)
		delete[] sBuffer[i];
	delete[] sBuffer;
	FreeGen(stream);
}

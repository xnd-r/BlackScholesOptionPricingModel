#include "../include/S10_Wiener.h"

VSLStreamStatePtr StockPriece::InitGen() {
	VSLStreamStatePtr stream;
	const unsigned int seed[2] = { __SEED__, __SEED__ };
	vslNewStreamEx(&stream, VSL_BRNG_MCG59, 2, seed); // base RNG 
	return stream;
}

void StockPriece::FreeGen(VSLStreamStatePtr stream) { // deleting of datastructure of gen
	vslDeleteStream(&stream); // is it really needed to write one-row func?
}

void StockPriece::GenerateGauss(double expect_val, double deviation, int amou,
	VSLStreamStatePtr stream, double *dest_array) {
	//Getting amount random values and writing them to the destination array
	vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, amou, dest_array, expect_val, deviation);
}


void StockPriece::SimulateWienerProcess(int nPaths, int nSteps, double Time, double **buffer) {

	VSLStreamStatePtr stream = InitGen();

	double *wiener_diff = new double[nSteps]; // Random values buffer
	double h = Time / (double)nSteps; // step

	for (int i = 0; i < nPaths; i++) {
		// getting nSteps random values with N(0, h)
		GenerateGauss(0, sqrt(h), nSteps, stream, wiener_diff);
		buffer[i][0] = 0;
		for (int j = 1; j <= nSteps; j++) {
			buffer[i][j] = buffer[i][j - 1] + wiener_diff[j - 1];
		}
	}
	delete[] wiener_diff;
	FreeGen(stream); // Generator data sturcture deleting
}

void StockPriece::WriteToCsv(double **buffer, int nRows, int nColumns) {

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


void StockPriece::PrintToFile(double **buffer, int nRows, int nColumns,
	char * FileName)
{
	FILE *f = fopen(FileName, "w");
	for (int i = 0; i < nRows; i++)
	{
		for (int j = 0; j < nColumns; j++)
			fprintf(f, "%lf;", buffer[i][j]);
		fprintf(f, "\n");
	}
	fclose(f);
}

#include "../include/S20_AnSimple.h"

double AnSimple::SimulateStockPrices(int nPaths, double Time, double *sBuffer) {
	VSLStreamStatePtr stream = InitGen();

	double *wiener_diff = new double[nPaths]; // random values with N(0, Time) buffer
	GenerateGauss(0.0, sqrt(Time), nPaths, stream, wiener_diff);
	// As W(0)=0 and W(Time) = dW =>
	// W == dW
	//result can be calculated as nPath vaues of W(Time)
	double sum = 0.0;
	for (int i = 0; i < nPaths; i++) {
		sBuffer[i] = GetStockPrice(wiener_diff[i], Time);
		sum += sBuffer[i];
	}
	sum = sum / (double)nPaths;
	FreeGen(stream);
	delete[] wiener_diff;
	return sum;
}

void AnSimple::WriteToCsv(double *buffer, int nRows, int Time, double avg) {

	time_t rawtime;
	time(&rawtime);
	std::string date = asctime(localtime(&rawtime));
	date.pop_back();
	date.append("_GBM_Analitycal_Simple.csv");

	for (std::string::iterator it = date.begin(); it<date.end(); ++it) {
		if (*it == ':') {
			date.erase(it);
		}
		std::replace(date.begin(), date.end(), ' ', '_');
	}

	FILE *f = fopen(date.c_str(), "w");
	for (int j = 0; j < nRows; j++) {
		std::string tmp_cell = std::to_string(buffer[j]);

		for (std::string::iterator it = tmp_cell.begin(); it<tmp_cell.end(); ++it) {
			std::replace(tmp_cell.begin(), tmp_cell.end(), '.', ',');
		}
		fprintf(f, /*"%lf;\n"*/"%s;", /*buffer[j]*/tmp_cell.c_str());
	}
	fprintf(f, "\n");
	fprintf(f, "Average price = %lf;\n", avg);
	fclose(f);
}

void AnSimple::Execute() {
	double *sBuffer = new double[NPATHS];
	double avg = SimulateStockPrices(NPATHS, TIME, sBuffer);
	WriteToCsv(sBuffer, NPATHS, NSTEPS, avg);
	printf("Average price = %lf\n", avg);
	delete[] sBuffer;
}

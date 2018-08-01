#include "../include/S32_NumMethodsWZ.h"

double NumMethodWZ::BurragePlatenStep(double S, double dt, double dw, double dz) {
	return S * (1.0 + R * dt + SIG * dw) + 0.5 * S * SIG * SIG * (dw * dw - dt) +
		R * SIG * S * dz + 0.5 * R * R * S * dt * dt +
		R * SIG * S * (dw * dt - dz) + 0.5 * SIG * SIG * SIG * S *
		(1.0 / 3.0 * dw * dw - dt) * dw;
}

double NumMethodWZ::Taylor2Step(double S, double dt, double dw, double dz) {
	return S * (1.0 + R * dt + SIG * dw + 0.5 * SIG * SIG * (dw * dw - dt) +
		R * SIG * dz + 0.5 * R * R * dt * dt + R * SIG * (dw * dt - dz));
}

void NumMethodWZ::SimulateWandZProcesses(VSLStreamStatePtr stream, int nSteps, double Time, double *buffer) {
	double *dw = new double[nSteps * 2];

	double h = Time / (double)nSteps;

	// Getting nSteps random vectors, dim == 2
	double mean[2]; mean[0] = 0.0; mean[1] = 0.0;
	double hh = h * sqrt(h);
	double cov[3]; cov[0] = sqrt(h); cov[1] = 0.5 * hh;
	cov[2] = 1 / sqrt(12.0) * hh;

	vdRngGaussianMV(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, nSteps, dw, 2, VSL_MATRIX_STORAGE_PACKED, mean, cov);

	buffer[0] = 0; buffer[1] = 0;

	for (int j = 2; j <= nSteps * 2; j += 2)
	{
		buffer[j] = buffer[j - 2] + dw[j - 2];
		buffer[j + 1] = buffer[j - 1] + dw[j - 1];
	}

	delete[] dw;
}

void NumMethodWZ::SimulateStockPrices(Step _step, VSLStreamStatePtr stream, int nPaths, int nSteps, double Time, double *Error) {
	double *w_traject = new double[(nSteps + 1) * 2];

	for (int i = 0; i < nPaths; i++) {
		SimulateWandZProcesses(stream, nSteps, Time, w_traject);

		double S_An = GetStockPrice(w_traject[nSteps * 2], Time);

		double S_Num[8]; // array of S(t) for different scaling 
		int scale = 128;

		for (int j = 0; j < 8; j++) {
			S_Num[j] = S0;
			int numMethodSteps = nSteps / scale;
			double h = Time / (double)numMethodSteps;
			double t = 0;
			int index = 0;
			for (int k = 0; k < numMethodSteps; k++) {
				t = t + h;
				index = index + scale;
				double dw = w_traject[index * 2] - w_traject[(index - scale) * 2]; 
				double dz = w_traject[index * 2 + 1] - w_traject[(index - scale) * 2 + 1];
				S_Num[j] = (this->*_step)(S_Num[j], h, dw, dz);
			}

			Error[j] = Error[j] + fabs(S_An - S_Num[j]);
			scale = scale / 2;
		}
	}

	for (int j = 0; j < 8; j++)
		Error[j] = Error[j] / nPaths;
	delete[] w_traject;
}

void NumMethodWZ::Execute(Step _step, char* FileName) {
	VSLStreamStatePtr stream = InitGen();
	double *Error = new double[8];
	for (int i = 0; i < 8; i++)
		Error[i] = 0.0;

	SimulateStockPrices(_step, stream, NPATHS, NSTEPS, TIME, Error);
	WriteToCsv(Error, NSTEPS, 8, TIME, 128, FileName);
	for (int i = 0; i < 8; i++)
		printf("Error %d = %lf\n", i + 1, Error[i]);
	delete[] Error;
	FreeGen(stream);
}
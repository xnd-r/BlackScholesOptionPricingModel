#include "../include/S31_NumMethodsW.h"

//NumMethodW::Step step_array[3] = { NumMethodW::EulMarStep, MilsteinStep, RK1Step };

double NumMethodW::EulMarStep(double S, double dt, double dw) {
	return S + S * (R * dt + SIG * dw);
}

double NumMethodW::MilsteinStep(double S, double dt, double dw) {
	return S + S * (R * dt + SIG * dw) +
		0.5 * S * SIG * SIG * (dw * dw - dt);
}

double NumMethodW::RK1Step(double S, double dt, double dw) {
	return S + S * (R * dt + SIG * dw) + 0.5 * (SIG * SIG * sqrt(dt)) * (dw * dw - dt);
}

void NumMethodW::SimulateStockPrices(Step _step, VSLStreamStatePtr stream, int nPaths, int nSteps, double Time, double *Error)
{
	double *w_traject = new double[nSteps + 1];

	for (int i = 0; i < nPaths; i++) {
		SimulateWienerProcess(stream, nSteps, Time, w_traject);

		double S_An = GetStockPrice(w_traject[nSteps], Time);

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
				S_Num[j] = _step(S_Num[j], h, w_traject[index] - w_traject[index - scale]);
			}

			Error[j] = Error[j] + fabs(S_An - S_Num[j]);
			scale = scale / 2;
		}
	}

	for (int j = 0; j < 8; j++)
		Error[j] = Error[j] / nPaths;
	delete[] w_traject;
}
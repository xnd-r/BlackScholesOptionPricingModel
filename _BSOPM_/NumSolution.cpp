#include "NumSolution.h"
#include "AnSolution.h"
#include <iostream>

float NumSolution::EulMarStep(float S, float dt, float dw, float r, float sig) {
	return S + S * (r * dt + sig * dw);
}

float NumSolution::MilsteinStep(float S, float dt, float dw, float r, float sig) {
	return S + S * (r * dt + sig * dw) + 0.5f * S * sig * sig * (dw * dw - dt);
}

float NumSolution::RK1Step(float S, float dt, float dw, float r, float sig) {
	return S * (1.0f + r * dt + sig * dw + 0.5f * sig * (r * dt + sig * sqrtf(dt)) * (dw * dw - dt) / sqrtf(dt));
}

float NumSolution::BurragePlatenStep(float S, float dt, float dw, float dz, float r, float sig) {
	return S * (1.f + r * dt + sig * dw) + .5f * S * sig * sig * (dw * dw - dt) +
		r * sig * S * dz + .5f * r * r * S * dt * dt +
		r * sig * S * (dw * dt - dz) + .5f * sig * sig * sig * S *
		(1.f / 3.f * dw * dw - dt) * dw;
}

void NumSolution::wAndZProcesses(VSLStreamStatePtr stream, int nsteps, float time, float *buffer) {
	float *w = new float[nsteps * 2];
	float dt = time / nsteps;
	float mean[2] = { 0.f, 0.f };
	float hh = dt * sqrtf(dt);
	float cov[3] = { sqrtf(dt) , 0.5f * hh , 1.0f / sqrtf(12.0f) * hh };

	vsRngGaussianMV(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, nsteps, w, 2, VSL_MATRIX_STORAGE_PACKED, mean, cov);
	buffer[0] = 0; buffer[1] = 0;
	for (int j = 2; j <= nsteps * 2; j += 2)
	{
		buffer[j]		= w[j - 2];
		buffer[j + 1]	= w[j - 1];
	}
	delete[] w;
}

void NumSolution::wAndZProcessesLarge(VSLStreamStatePtr stream, int nsteps, int nsamples, float time, float *w) {
	float *gaussBuf = new float[(nsteps + 1) * nsamples * 2];
	float dt = time / nsteps;
	float mean[2] = { 0.f, 0.f };
	float hh = dt * sqrtf(dt);
	float cov[3] = { sqrtf(dt) , 0.5f * hh , 1.0f / sqrtf(12.0f) * hh };

	vsRngGaussianMV(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, (nsteps + 1) * nsamples, gaussBuf, 2, VSL_MATRIX_STORAGE_PACKED, mean, cov);
	for (int j = 0; j <= (nsteps + 1) * nsamples * 2; j += 2)
	{
		w[j] = gaussBuf[j - 2];
		w[j + 1] = gaussBuf[j - 1];
	}
	for (int j = 0; j < (nsteps + 1) * nsamples * 2; j += (nsamples + 1)) {
		w[j] = 0.f;
		w[j+1] = 0.f;
	}
	delete[] gaussBuf;
}

void NumSolution::wienerProcessLarge(VSLStreamStatePtr stream, int nsteps, int nsamples, float time, float *w) {

	float *gaussBuf = new float[(nsteps + 1) * nsamples]; // Random values buffer
	float dt = time / nsteps;
	normalGenerator(.0f, sqrtf(dt), (nsteps + 1) * nsamples, stream, gaussBuf);
	for (int j = 0; j < (nsteps + 1) * nsamples; ++j) {
		w[j] = gaussBuf[j];
	}
	for (int j = 0; j < (nsteps + 1) * nsamples; j += (nsamples + 1)) {
		w[j] = 0.f;
	}
	delete[] gaussBuf;
}

float* NumSolution::memoryWrajectAlloc(VSLStreamStatePtr stream, int StepIndex, int npaths, int nsteps, float time, float* wtraject) {
	if (StepIndex != 3) {
		wtraject = new float[npaths * (nsteps + 1)];
		wienerProcessLarge(stream, nsteps, npaths, time, wtraject); // npaths random numbers became zero
	}
	else {
		wtraject = new float[npaths * (nsteps + 1) * 2];
		wAndZProcessesLarge(stream, nsteps, npaths, time, wtraject);
	}
	return wtraject;
}

bool NumSolution::checkConvergence(int StepIndex, VSLStreamStatePtr stream, int npaths, int nsteps, float pS0, float pR, float pSig, float time, float *Error, unsigned int seed, int indexGen)
{
// NOT WORKS!!
	for (int i = 0; i < 8; i++)
		Error[i] = 0.0f;
	float* wtraject, *sAnMem;
	if (StepIndex != 3) {
		wtraject = new float[nsteps + 1];
		sAnMem = new float[nsteps];
	}
	else {
		wtraject = new float[(nsteps + 1) * 2];
		sAnMem = new float[nsteps * 2];
	}
	float S_An;
	AnSolution as;
	for (int i = 0; i < npaths; i++) {

		if (StepIndex != 3) {
			S_An = as.simulateStockPriceAn(nsteps, pS0, pR, pSig, time, wtraject, seed, indexGen); //getStockPrice(pS0, pR, pSig, wtraject[nsteps], time); // one or scope of trajectories?
			wienerProcess(stream, nsteps, time, wtraject);
		}
		else {
			S_An = as.simulateStockPriceAn(nsteps/** 2*/, pS0, pR, pSig, time, wtraject, seed, indexGen); // getStockPrice(pS0, pR, pSig, wtraject[nsteps * 2], time);
			wAndZProcesses(stream, nsteps, time, wtraject);
		}
		float S_Num[8]; // array of S(t) for different scaling 
		int scale = 128;

		for (int j = 0; j < 8; j++) {
			S_Num[j] = pS0;
			int numMethodSteps = nsteps / scale;
			float dt = time / numMethodSteps;
			int index = 0;
			switch (StepIndex) {
			case 0: // EulMarStep
				for (int k = 0; k < numMethodSteps; k++) {
					index += scale;
					S_Num[j] = EulMarStep(S_Num[j], dt, wtraject[index]/* - wtraject[index - scale]*/, pR, pSig);
				}
				break;
			case 1: // MilsteinStep
				for (int k = 0; k < numMethodSteps; k++) {
					index += scale;
					S_Num[j] = MilsteinStep(S_Num[j], dt, wtraject[index]/* - wtraject[index - scale]*/, pR, pSig);
				}
				break;
			case 2: // RK1Step
				for (int k = 0; k < numMethodSteps; k++) {			
					index += scale;
					S_Num[j] = RK1Step(S_Num[j], dt, wtraject[index]/* - wtraject[index - scale]*/, pR, pSig);
				}
				break;
			case 3: // BurragePlatenStep
				float dw, dz;
				for (int k = 0; k < numMethodSteps; k++) {
					index += scale;
					dw = wtraject[index * 2];
					dz = wtraject[index * 2 + 1];
					S_Num[j] = BurragePlatenStep(S_Num[j], dt, dw, dz, pR, pSig);
				}
				break;
			}
			Error[j] += fabs(S_An - S_Num[j]);
			scale = scale / 2;
		}
	}

	for (int j = 0; j < 8; j++)
		Error[j] = Error[j] / npaths;
	delete[] wtraject;
	return Error[7] < 0.01 ? true : false;
}

float NumSolution::stockPricesIntegrator(VSLStreamStatePtr stream, float* wtraject, int StepIndex, int nsteps, float pS0, float pR, float pSig, float time)
{
	float dt = time / nsteps;
	float stockPrice = pS0;

	switch (StepIndex) {
		case 0: // EulMarStep
			for (int k = 0; k < nsteps; k++) {
				stockPrice = EulMarStep(stockPrice, dt, wtraject[k], pR, pSig);
			}
			break;
		case 1: // MilsteinStep
			for (int k = 0; k < nsteps; k++) {
				stockPrice = MilsteinStep(stockPrice, dt, wtraject[k], pR, pSig);
			}
			break;
		case 2: // RK1Step
			for (int k = 0; k < nsteps; k++) {
				stockPrice = RK1Step(stockPrice, dt, wtraject[k], pR, pSig);
			}
			break;
		case 3: // BurragePlatenStep
			float dw, dz;
			for (int k = 0; k < nsteps; k++) {
				dw = wtraject[k * 2];
				dz = wtraject[k * 2 + 1];
				stockPrice = BurragePlatenStep(stockPrice, dt, dw, dz, pR, pSig);
			}
			break;
		}
	return stockPrice;
}

float NumSolution::stockPricesIntegratorVol(VSLStreamStatePtr stream, float* wtraject, int StepIndex, int nsteps, float pS0, float* pR, float* pSig, float time)
{
	float dt = time / nsteps;
	float stockPrice = pS0;

	switch (StepIndex) {
	case 0: // EulMarStep
		for (int k = 0; k < nsteps; k++) {
			stockPrice = EulMarStep(stockPrice, dt, wtraject[k], pR[k], pSig[k]);
		}
		break;
	case 1: // MilsteinStep
		for (int k = 0; k < nsteps; k++) {
			stockPrice = MilsteinStep(stockPrice, dt, wtraject[k], pR[k], pSig[k]);
		}
		break;
	case 2: // RK1Step
		for (int k = 0; k < nsteps; k++) {
			stockPrice = RK1Step(stockPrice, dt, wtraject[k], pR[k], pSig[k]);
		}
		break;
	case 3: // BurragePlatenStep
		float dw, dz;
		for (int k = 0; k < nsteps; k++) {
			dw = wtraject[k * 2];
			dz = wtraject[k * 2 + 1];
			stockPrice = BurragePlatenStep(stockPrice, dt, dw, dz, pR[k], pSig[k]);
		}
		break;
	}
	return stockPrice;
}

float NumSolution::SimulateStockPrices(int StepIndex, int indexGen, int npaths, int nsteps, float pS0, float pR, float pSig, float time, unsigned int seed)
{
	VSLStreamStatePtr stream = initGen(seed, indexGen);
	float* wtraject = memoryWrajectAlloc(stream, StepIndex, npaths, nsteps, time, wtraject);

	float stockPrice = 0.f;
	float dt = time / nsteps;
	for (int i = 0; i < npaths; i++) {
		stockPrice += stockPricesIntegrator(stream, &wtraject[i * (nsteps + 1)], StepIndex, nsteps, pS0, pR, pSig, time);
	}
	delete[] wtraject;
	freeGen(stream);
	return stockPrice / npaths;
}

float NumSolution::SimulateStockPricesVol(int StepIndex, int npaths, int nsteps, float pS0, float* pR, float* pSig, float time, unsigned int seed, int indexGen)
{
	VSLStreamStatePtr stream = initGen(seed, indexGen);
	float* wtraject = memoryWrajectAlloc(stream, StepIndex, npaths, nsteps, time, wtraject);

	float stockPrice = 0.f;
	float dt = time / nsteps;
	for (int i = 0; i < npaths; i++) {
		stockPrice += stockPricesIntegratorVol(stream, &wtraject[i * (nsteps + 1)], StepIndex, nsteps, pS0, pR, pSig, time);
	}
	delete[] wtraject;
	freeGen(stream);
	return stockPrice / npaths;
}

float NumSolution::getMCPrice(int StepIndex, int nsteps, int indexGen, int N, unsigned int seed, float K, float R, float Time, float SIG, float pS0) {
	VSLStreamStatePtr stream = initGen(seed, indexGen);
	float* wtraject = memoryWrajectAlloc(stream, StepIndex, N, nsteps, Time, wtraject);

	float payoff, sum = .0f, stockPrice;
	for (unsigned int i = 0; i < N; i++) {
			stockPrice = stockPricesIntegrator(stream, &wtraject[i * (nsteps + 1) ], StepIndex, nsteps, pS0, R, SIG, Time);
			payoff = stockPrice - K;
			if (payoff > 0.f)
				sum += payoff;
	}
	sum = sum / N * expf(-R * Time);
	vslDeleteStream(&stream);
	delete[] wtraject;
	return sum;
}

float NumSolution::getMCPricePar(int NumThreads, int StepIndex, int nsteps, int indexGen, int N, unsigned int seed, float K, float R, float Time, float SIG, float pS0, double& workTime) {
	VSLStreamStatePtr stream = initGen(seed, indexGen);
	float* wtraject = memoryWrajectAlloc(stream, StepIndex, N, nsteps, Time, wtraject);
	double t1, t2;
	omp_set_num_threads(NumThreads);
	float sum = .0f, stockPrice;
	t1 = omp_get_wtime();
#pragma omp parallel private(stockPrice)
	{
		int count = omp_get_num_threads();
		int num = omp_get_thread_num();
		vslSkipAheadStream(stream, N * (nsteps + 1) / count * num); // if index == 3?
		float payoff;
#pragma omp for private(payoff), reduction(+:sum)
		for (int i = 0; i < N; i++) {
			stockPrice = stockPricesIntegrator(stream, &wtraject[i * (nsteps + 1)], StepIndex, nsteps, pS0, R, SIG, Time);
			payoff = stockPrice - K;
			if (payoff > 0.f) 
				sum += payoff;
		}
	}
	t2 = omp_get_wtime();
	workTime = t2 - t1;
	sum = sum / N * expf(-R * Time);
	vslDeleteStream(&stream);
	delete[] wtraject;
	return sum;
}

void NumSolution::MCExecute(int StepIndex, int nsteps, int indexGen, int N, unsigned int seed, float K, float R, float Time, float SIG, float pS0) {
	int tmp = 1;
	double workTime = 0.0;
	float tmp_price;
	std::vector<double> Times;
	std::vector<double> Prices;
	for (int k = 1; k < 5; ++k) {
		for (int j = 0; j < 7; ++j) {

			tmp_price = getMCPricePar(tmp, 2, nsteps, indexGen, N, seed, K, R, Time, SIG, pS0, workTime);
			Times.push_back(workTime);
			Prices.push_back(tmp_price);
			std::cout << j << " Writing Done " << std::endl;
		}
		std::sort(Times.begin(), Times.end());
		for (std::vector<double>::const_iterator it = Times.begin(); it != Times.end(); ++it)
			std::cout << *it << std::endl;
		std::cout << " ####################### " << std::endl;
		for (std::vector<double>::const_iterator it = Prices.begin(); it != Prices.end(); ++it)
			std::cout << *it << std::endl;
		std::cout << tmp << " thread(stockPrice) done" << std::endl;
		tmp *= 2;
		Times.clear();
		Prices.clear();
	}

}
void NumSolution::WriteToCsv(float* Errors, int nsteps, int nrows, float Time, int scale, int stepIndex) {
	row_table rt;
	time_t timestamp;
	time(&timestamp);
	std::string date = asctime(localtime(&timestamp));
	date.pop_back();
	std::string fileName;
	switch (stepIndex) {
	case 0:
		fileName = "_EulerMar.csv";
		break;
	case 1:
		fileName = "_Milstein.csv";
		break;
	case 2:
		fileName = "_RK1.csv";
		break;
	case 3:
		fileName = "_BurragePlaten.csv";
		break;

	}
	date.append(fileName);

	for (std::string::iterator it = date.begin(); it<date.end(); ++it) {
		if (*it == ':') {
			date.erase(it);
		}
		std::replace(date.begin(), date.end(), ' ', '_');
	}


	FILE *f = fopen(date.c_str(), "w");
	fprintf(f, "step;e;log(step);log(e);\n");
	for (int i = 0; i < nrows; ++i) {
		rt.setValues(Errors[i], nsteps, static_cast<int>(scale / pow(2, i)), nrows, Time);
		std::string tmp_string = rt.createString();
		for (std::string::iterator it = tmp_string.begin(); it<tmp_string.end(); ++it) {
			std::replace(tmp_string.begin(), tmp_string.end(), '.', ',');
		}
		fprintf(f, "%stockPrice\n", tmp_string.c_str());
	}
	fclose(f);
}

void NumSolution::getErrors(int nsteps, int indexGen, int N, unsigned int seed, float K, float R, float Time, float SIG, float pS0, int sampleStep) {
	time_t timestamp;
	time(&timestamp);
	std::string date = asctime(localtime(&timestamp));
	date.pop_back();
	std::string fileName;
	date.append(std::to_string(N));
	date.append("_Samples_");
	date.append(std::to_string(sampleStep));
	date.append("_Step_");
	fileName = "_Errors.csv";
	date.append(fileName);


	for (std::string::iterator it = date.begin(); it<date.end(); ++it) {
		if (*it == ':') {
			date.erase(it);
		}
		std::replace(date.begin(), date.end(), ' ', '_');
	}


	FILE *f = fopen(date.c_str(), "w");

	fprintf(f, "%stockPrice;%stockPrice;%stockPrice;%stockPrice\n;", "Euler", "Mils", "RK1", "Platen");
	for (int i = 0; i < N; i += sampleStep) {
		fprintf(f, "%lf;%lf;%lf;%lf\n", 
			(getMCPrice(0, nsteps, indexGen, i, seed, K, R, Time, SIG, pS0) - 20.9244),
			(getMCPrice(1, nsteps, indexGen, i, seed, K, R, Time, SIG, pS0) - 20.9244),
			(getMCPrice(2, nsteps, indexGen, i, seed, K, R, Time, SIG, pS0) - 20.9244), 
			(getMCPrice(3, nsteps, indexGen, i, seed, K, R, Time, SIG, pS0) - 20.9244));
		}
	
	fclose(f);

}

void NumSolution::Execute(int StepIndex, int indexGen, int npaths, int nsteps, float pS0, float pR, float pSig, float time, unsigned int seed) {
	VSLStreamStatePtr stream = initGen(seed, indexGen);
	float *Error = new float[8];
	checkConvergence(StepIndex, stream, npaths, nsteps, pS0, pR, pSig, time, Error, seed, indexGen);
	WriteToCsv(Error, nsteps, 8, time, 128, StepIndex);
	for (int i = 0; i < 8; i++)
		printf("Error %d = %lf\n", i + 1, Error[i]);
	freeGen(stream);
	delete[] Error;
}

//float NumSolution::Integrand(float z) {
//	float payoff;	
//	int scale = 2000;
//	stockPrice = pS0 * expf(tmp1 + tmp2 * z);
//	payoff = stockPrice - K;
//	payoff > 0.0f ? payoff *= expf(-z * z / 2.0f) : payoff = 0.0f;
//	return payoff;
//}
//
//float NumSolution::GetRPrice(float a, float b) {
//	int scale = 2000;
//	float h = (b - a) / scale;
//	start = clock();
//	for (int i = 0; i < scale; ++i)
//		sum += Integrand(a + h * (i + 0.5f));
//	sum = sum * expf(-R * TIME) / sqrtf(2.0f * M_PIF) * h;
//
//	finish = clock();
//	t = (double)(finish - start) / CLOCKS_PER_SEC;
//
//	return sum;
//}
//
//float NumSolution::GetTPrice(float a, float b) {
//	int scale = 2000;
//	float h = (b - a) / scale;
//	start = clock();
//
//	sum += Integrand(a) * 0.5f * h + Integrand(b) * 0.5f * (b - a);
//	for (int i = 1; i < scale - 1; ++i)
//		sum += Integrand(a + h * i);
//	sum = sum * expf(-R * TIME) / sqrtf(2.0f * M_PIF) * h;
//
//	finish = clock();
//	t = (double)(finish - start) / CLOCKS_PER_SEC;
//
//	return sum;
//}
//
//float NumSolution::GetSPrice(float a, float b) {
//	int scale = 2000;
//	float h = (b - a) / scale;
//	float sum2 = 0.0f;
//
//	start = clock();
//	float sum4 = Integrand(a + h);
//	sum += Integrand(a) + Integrand(b);
//
//	for (int i = 1; i < scale - 2; i += 2) {
//		sum4 += Integrand(a + (i + 1)*h);
//		sum2 += Integrand(a + i * h);
//	}
//
//	finish = clock();
//	t = (double)(finish - start) / CLOCKS_PER_SEC;
//	sum = (sum + 4 * sum4 + 2 * sum2) * expf(-R * TIME) / sqrtf(2.0f * M_PIF) * (h / 3.0f);
//
//	return sum;
//}
//
//float NumSolution::Get3_8Price(float a, float b) {
//	int scale = 2000;
//	float h = (b - a) / scale;
//	float sum2 = 0.0f;
//	float sum3 = Integrand(a + h);
//
//	start = clock();
//	sum += Integrand(a) + Integrand(b);
//
//	for (int i = 1; i < scale - 2; i += 3) {
//		sum3 += Integrand(a + i * h);
//	}
//
//	for (int i = 2; i < scale - 1; i += 3) {
//		sum3 += Integrand(a + i * h);
//	}
//
//	for (int i = 3; i < scale; i += 3) {
//		sum2 += Integrand(a + i * h);
//	}
//
//	finish = clock();
//	t = (double)(finish - start) / CLOCKS_PER_SEC;
//	sum = (sum + 3 * sum3 + 2 * sum2) * expf(-R * TIME) / sqrtf(2.0f * M_PIF) * (3.0f * h / 8.0f);
//
//	return sum;
//}
//
//void NumSolution::SetS(int amo) {
//	int scale = 2000;
//	s_array = new float[amo];
//	exp_array = new float[amo];
//#if defined(__INTEL_COMPILER) 
//#pragma ivdep
//#pragma vector always	
//#endif
//	int j = 0;
//#pragma omp parallel for private(j) 
//	for (j = 0; j < amo; ++j) {
//		s_array[j] = pS0 * expf(tmp1 + tmp2 * (a + h * ((float)j + 0.5f)));
//		exp_array[j] = expf(-(a + h * ((float)j + 0.5f)) * (a + h * ((float)j + 0.5f)) / 2.0f);
//	}
//}
//
//float NumSolution::GetRPrice(float a, float b, int NumThreads/*, float* s_array, float* expf_array*/) {
//	float sum = 0.0f;
//	int scale = 2000;
//	float h = (b - a) / scale;
//	float payoff;
//	int i, ind;
//	start = clock();
//	omp_set_num_threads(NumThreads);
//#if defined(__INTEL_COMPILER) 
//#pragma ivdep
//#pragma vector always	
//#endif
//#pragma omp parallel for private(i, payoff, ind) reduction(+:sum)
//	for (ind = 0; ind < N; ind++) {
//		for (i = 0; i < scale; ++i) {
//			payoff = s_array[i] - K;
//			payoff > 0.0f ? payoff *= exp_array[i] : payoff = 0.0f;
//			sum = sum + payoff;
//		}
//	}
//	sum = (sum * expf(-R * TIME) * h) / sqrtf(2.0f * M_PIF);
//	sum /= N;
//	finish = clock();
//	t = (double)(finish - start) / CLOCKS_PER_SEC;
//
//	return sum;
//}
//
//float NumSolution::GetTPrice(float a, float b, int NumThreads) {
//	float sum;
//	int scale = 2000;
//	float h = (b - a) / scale;
//	float payoff;
//	int i, ind;
//	start = clock();
//	omp_set_num_threads(NumThreads);
//#if defined(__INTEL_COMPILER) 
//#pragma ivdep
//#pragma vector always	
//#endif
//
//	sum = (Integrand(a) + Integrand(b)) * 0.5f * h * (b - a);
//#pragma omp parallel for private(i, payoff, ind) reduction(+:sum)
//	for (ind = 0; ind < N; ind++) {
//		for (int i = 1; i < scale - 1; ++i) {
//			payoff = s_array[i] - K;
//			payoff > 0.0f ? payoff *= exp_array[i] : payoff = 0.0f;
//			sum = sum + payoff;
//		}
//	}
//	sum = sum * expf(-R * TIME) / sqrtf(2.0f * M_PIF) * h;
//	sum /= N;
//	finish = clock();
//	t = (double)(finish - start) / CLOCKS_PER_SEC;
//
//	return sum;
//}
//
//float NumSolution::GetSPrice(float a, float b, int NumThreads) {
//	float sum;
//	int scale = 2000;
//	float h = (b - a) / scale;
//	float payoff;
//	int i, ind;
//	float sum2 = 0.0f;
//	start = clock();
//	omp_set_num_threads(NumThreads);
//#if defined(__INTEL_COMPILER) 
//#pragma ivdep
//#pragma vector always	
//#endif
//	float sum4 = Integrand(a + h);
//	sum = Integrand(a) + Integrand(b);
//
//#pragma omp parallel for private(i, payoff, ind) reduction(+:sum2, sum4)
//	for (ind = 0; ind < N; ind++) {
//		for (int i = 1; i < scale - 2; i += 2) {
//			payoff = s_array[i] - K;
//			payoff > 0.0f ? payoff *= exp_array[i] : payoff = 0.0f;
//			sum2 = sum2 + payoff;
//
//			payoff = s_array[i + 1] - K;
//			payoff > 0.0f ? payoff *= exp_array[i + 1] : payoff = 0.0f;
//			sum4 = sum4 + payoff;
//		}
//	}
//	finish = clock();
//	t = (double)(finish - start) / CLOCKS_PER_SEC;
//	sum = (sum + 4 * sum4 + 2 * sum2) * expf(-R * TIME) / sqrtf(2.0f * M_PIF) * (h / 3.0f);
//	sum /= N;
//	return sum;
//}
//
//
//float NumSolution::Get3_8Price(float a, float b, int NumThreads) {
//	float sum;
//	int scale = 2000;
//	float h = (b - a) / scale;
//	float payoff;
//	int i, ind;
//	start = clock();
//	omp_set_num_threads(NumThreads);
//#if defined(__INTEL_COMPILER) 
//#pragma ivdep
//#pragma vector always	
//#endif
//	float sum2 = 0.0f;
//	float sum3 = Integrand(a + h);
//	sum += Integrand(a) + Integrand(b);
//
//#pragma omp parallel for private(i, payoff, ind) reduction(+:sum2, sum3)
//	for (ind = 0; ind < N; ind++) {
//		for (int i = 1; i < scale - 2; i += 3) {
//			payoff = s_array[i] - K;
//			payoff > 0.0f ? payoff *= exp_array[i] : payoff = 0.0f;
//			sum3 = sum3 + payoff;
//		}
//
//		for (int i = 2; i < scale - 1; i += 3) {
//			payoff = s_array[i] - K;
//			payoff > 0.0f ? payoff *= exp_array[i] : payoff = 0.0f;
//			sum3 = sum3 + payoff;
//		}
//
//		for (int i = 3; i < scale; i += 3) {
//			payoff = s_array[i] - K;
//			payoff > 0.0f ? payoff *= exp_array[i] : payoff = 0.0f;
//			sum2 = sum2 + payoff;
//		}
//	}
//	finish = clock();
//	t = (double)(finish - start) / CLOCKS_PER_SEC;
//	sum = (sum + 3 * sum3 + 2 * sum2) * expf(-R * TIME) / sqrtf(2.0f * M_PIF) * (3.0f * h / 8.0f);
//	sum /= N;
//	return sum;
//}
//

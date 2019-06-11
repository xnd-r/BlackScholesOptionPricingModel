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


void NumSolution::wAndZArrayLarge(VSLStreamStatePtr stream, int nsteps, int nsamples, float time, float *w) {
	float dt = time / nsteps;
	float mean[2] = { 0.f, 0.f };
	float hh = dt * sqrtf(dt);
	float cov[3] = { sqrtf(dt) , 0.5f * hh , 1.0f / sqrtf(12.0f) * hh };
	vsRngGaussianMV(VSL_RNG_METHOD_GAUSSIANMV_ICDF, stream, nsteps * nsamples, w, 2, VSL_MATRIX_STORAGE_PACKED, mean, cov);
}

void BPMemAlloc(VSLStreamStatePtr stream, int nsteps, float time, float *w) {
	float dt = time / nsteps;
	float mean[2] = { 0.f, 0.f };
	float hh = dt * sqrtf(dt);
	float cov[3] = { sqrtf(dt) , 0.5f * hh , 1.0f / sqrtf(12.0f) * hh };
	vsRngGaussianMV(VSL_RNG_METHOD_GAUSSIANMV_ICDF, stream, nsteps, w, 2, VSL_MATRIX_STORAGE_PACKED, mean, cov);
}


float* NumSolution::memoryWrajectAlloc(VSLStreamStatePtr stream, int StepIndex, int npaths, int nsteps, float time, float* wtraject) {
// legacy
	if (StepIndex != 3) {
		wtraject = new float[npaths * nsteps];
		vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, nsteps * npaths, wtraject, 0.f, sqrtf(time/nsteps));	}
	else {
		wtraject = new float[npaths * (nsteps) * 2];
		vsRngGaussian(VSL_RNG_METHOD_GAUSSIANMV_ICDF, stream, nsteps * npaths, wtraject, 0.f, sqrtf(time/nsteps));
	}
	return wtraject;
}

void NumSolution::checkConvergence(int StepIndex, VSLStreamStatePtr stream, int npaths, int nsteps, float pS0, float pR, float pSig, float time, float *Error, unsigned int seed, int indexGen)
{
	for (int i = 0; i < 12; i++)
		Error[i] = 0.0f;
	float* wtraject;
	if (StepIndex != 3) {
		wtraject = new float[nsteps + 1];
	}
	else {
		wtraject = new float[(nsteps + 1) * 2];
	}
	float S_An;
	for (int i = 0; i < npaths; i++) {

		if (StepIndex != 3) {
			wienerProcess(stream, nsteps, time, wtraject);
			S_An = pS0 * expf((pR - pSig * pSig / 2.f) * time + pSig * wtraject[nsteps]);

		}
		else {
			wAndZProcesses(stream, nsteps, time, wtraject);
			S_An = pS0 * expf((pR - pSig * pSig / 2.f) * time + pSig * wtraject[nsteps * 2]);
		}
		float S_Num[12]; // array of S(t) for different scaling 
		int scale = 2048;

		for (int j = 0; j < 12; j++) {
			S_Num[j] = pS0;
			int numMethodSteps = nsteps / scale;
			float dt = time / numMethodSteps;
			int index = 0;
			switch (StepIndex) {
			case 0: // EulMarStep
				for (int k = 0; k < numMethodSteps; k++) {
					index += scale;
					S_Num[j] = EulMarStep(S_Num[j], dt, wtraject[index] - wtraject[index - scale], pR, pSig);
				}
				//std::cout << i << "\t" << j << "\t" << numMethodSteps << "\t" << scale <<"\n";
				break;
			case 1: // MilsteinStep
				for (int k = 0; k < numMethodSteps; k++) {
					index += scale;
					S_Num[j] = MilsteinStep(S_Num[j], dt, wtraject[index] - wtraject[index - scale], pR, pSig);
				}
				break;
			case 2: // RK1Step
				for (int k = 0; k < numMethodSteps; k++) {			
					index += scale;
					S_Num[j] = RK1Step(S_Num[j], dt, wtraject[index] - wtraject[index - scale], pR, pSig);
				}
				break;
			case 3: // BurragePlatenStep
				float dw, dz;
				for (int k = 0; k < numMethodSteps; k++) {
					index += scale;
					dw = wtraject[index * 2] - wtraject[(index - scale) * 2];
					dz = wtraject[index * 2 + 1] - wtraject[(index - scale) * 2 + 1];
					S_Num[j] = BurragePlatenStep(S_Num[j], dt, dw, dz, pR, pSig);
				}
				break;
			}
			Error[j] += fabs(S_An - S_Num[j]);
			scale = scale / 2;
		}
	}

	for (int j = 0; j < 12; j++)
		Error[j] = Error[j] / npaths;
	delete[] wtraject;
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
	return 0.f;
}

float NumSolution::SimulateStockPrices(int StepIndex, int indexGen, int npaths, int nsteps, float pS0, float pR, float pSig, float time, unsigned int seed)
{
	VSLStreamStatePtr stream = initGen(seed, indexGen);
	float* wtraject;//= memoryWrajectAlloc(stream, StepIndex, npaths, nsteps, time, wtraject);
	
		//if (StepIndex != 3) {
	//	wtraject = new float[npaths * nsteps];
	//	vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, nsteps * npaths, wtraject, 0.f, sqrtf(time/nsteps));
	//}
	//else {
	//	wtraject = new float[npaths * nsteps * 2];
	//	wAndZArrayLarge(stream, nsteps, npaths, time, wtraject);
	//}
	float stockPrice = 0.f;
	float dt = time / nsteps;
	double t1, t2, sum_time = 0.0;
	StepIndex != 3 ? wtraject = new float[nsteps] : wtraject = new float[2 * nsteps];
	for (int ind = 0; ind < npaths; ++ind) {
		StepIndex != 3 ? vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, nsteps, wtraject, 0.f, sqrtf(time / nsteps)) :
				BPMemAlloc(stream, nsteps, time, wtraject);
		t1 = omp_get_wtime();
			stockPrice += stockPricesIntegrator(stream, wtraject, StepIndex, nsteps, pS0, pR, pSig, time);
		t2 = omp_get_wtime();
		sum_time += (t2 - t1);
	}
	delete[] wtraject; // ???????????
	vslDeleteStream(&stream);
	std::cout << "Stock Price via " << StepIndex << " method is " << stockPrice / npaths << std::endl;
	std::cout << sum_time << std::endl;
	return stockPrice / npaths;
}

float NumSolution::SimulateStockPricesVol(int StepIndex, int npaths, int nsteps, float pS0, float* pR, float* pSig, float time, unsigned int seed, int indexGen)
{
	VSLStreamStatePtr stream = initGen(seed, indexGen);
	float* wtraject;// = memoryWrajectAlloc(stream, StepIndex, npaths, nsteps, time, wtraject);
	if (StepIndex != 3) {
		wtraject = new float[npaths * (nsteps)];
		vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, nsteps * npaths, wtraject, 0.f, sqrtf(time/nsteps));
	}
	else {
		wtraject = new float[npaths * (nsteps) * 2];
		wAndZArrayLarge(stream, nsteps, npaths, time, wtraject);
	}
	float stockPrice = 0.f;
	float dt = time / nsteps;
	for (int i = 0; i < npaths; i++) {
		stockPrice += stockPricesIntegratorVol(stream, &wtraject[i * (nsteps)], StepIndex, nsteps, pS0, pR, pSig, time);
	}
	delete[] wtraject;
	vslDeleteStream(&stream);
	return stockPrice / npaths;
}

float NumSolution::MC(int StepIndex, int nsteps, int indexGen, int N, unsigned int seed, float K, float R, float Time, float SIG, float pS0)
{
	VSLStreamStatePtr stream = initGen(seed, indexGen);
	float* wtraject;//= memoryWrajectAlloc(stream, StepIndex, npaths, nsteps, time, wtraject);

		//if (StepIndex != 3) {
	//	wtraject = new float[npaths * nsteps];
	//	vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, nsteps * npaths, wtraject, 0.f, sqrtf(time/nsteps));
	//}
	//else {
	//	wtraject = new float[npaths * nsteps * 2];
	//	wAndZArrayLarge(stream, nsteps, npaths, time, wtraject);
	//}
	float summ = 0.f, payoff, SP = 0.f;
	float dt = Time / nsteps;
	double t1, t2, sum_time = 0.0;
	StepIndex != 3 ? wtraject = new float[nsteps] : wtraject = new float[2 * nsteps];
	for (int ind = 0; ind < N; ++ind) {
		StepIndex != 3 ? vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, nsteps, wtraject, 0.f, sqrtf(Time / nsteps)) :
			BPMemAlloc(stream, nsteps, Time, wtraject);
		t1 = omp_get_wtime();
		SP = stockPricesIntegrator(stream, wtraject, StepIndex, nsteps, pS0, R, SIG, Time);
		payoff = SP - K;
		if (payoff > 0.f)
			summ += payoff;
		t2 = omp_get_wtime();
		sum_time += (t2 - t1);
	}
	summ = summ / N * expf(-R * Time);
	delete[] wtraject; // ???????????
	vslDeleteStream(&stream);
	std::cout << "Stock Price via " << StepIndex << " method is " << summ / N << std::endl;
	//std::cout << sum_time << std::endl;
	return summ / N;
}


float NumSolution::getMCPrice(int StepIndex, int nsteps, int indexGen, int N, unsigned int seed, float K, float R, float Time, float SIG, float pS0) {
	VSLStreamStatePtr stream = initGen(seed, indexGen);
	float* wtraject;// = memoryWrajectAlloc(stream, StepIndex, N, nsteps, Time, wtraject);
	if (StepIndex != 3) {
		wtraject = new float[N * (nsteps)];
		vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, nsteps * N, wtraject, 0.f, sqrtf(Time/nsteps));
	}
	else {
		wtraject = new float[N * (nsteps) * 2];
		wAndZArrayLarge(stream, nsteps, N, Time, wtraject);
	}
	float payoff, sum = .0f, stockPrice;
	double t1 = omp_get_wtime();
	for (int i = 0; i < N; i++) {
		StepIndex != 3 ?
			stockPrice = 116.177 :// stockPricesIntegrator(stream, &wtraject[i * (nsteps)], StepIndex, nsteps, pS0, R, SIG, Time) :
			stockPrice = stockPricesIntegrator(stream, &wtraject[i * (nsteps) * 2], StepIndex, nsteps, pS0, R, SIG, Time);
		payoff = stockPrice - K;
		if (payoff > 0.f)
			sum += payoff;
	}
	sum = sum / N * expf(-R * Time);
	double t2 = omp_get_wtime();
	vslDeleteStream(&stream);
	delete[] wtraject;
	//std::cout << "MC EU_OP " << sum << std::endl;
	//std::cout << "Time is " << t2 - t1 << std::endl;
	return sum;
}

float NumSolution::getMCPricePar(int NumThreads, int StepIndex, int nsteps, int indexGen, int N, unsigned int seed, float K, float R, float Time, float SIG, float pS0, double& workTime) {
	VSLStreamStatePtr stream = initGen(seed, indexGen);
	float* wtraject;
		if (StepIndex != 3) {
		wtraject = new float[N * (nsteps)];
		vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, nsteps * N, wtraject, 0.f, sqrtf(Time/nsteps));	}
	else {
		wtraject = new float[N * (nsteps) * 2];
		wAndZArrayLarge(stream, nsteps, N, Time, wtraject);
	}
	double t1, t2;
	omp_set_num_threads(NumThreads);
	float sum = .0f, stockPrice;
	t1 = omp_get_wtime();
#pragma omp parallel private(stockPrice)
	{
		int count = omp_get_num_threads();
		int num = omp_get_thread_num();
		StepIndex != 3 ? vslSkipAheadStream(stream, N * (nsteps) / count * num) : vslSkipAheadStream(stream, N * 2 * (nsteps) / count * num);
		float payoff;
#pragma omp for private(payoff), reduction(+:sum)
		for (int i = 0; i < N; i++) {
			StepIndex != 3 ? 
				stockPrice = stockPricesIntegrator(stream, &wtraject[i * (nsteps)	 ], StepIndex, nsteps, pS0, R, SIG, Time) : 
				stockPrice = stockPricesIntegrator(stream, &wtraject[i * (nsteps) * 2], StepIndex, nsteps, pS0, R, SIG, Time);
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

void NumSolution::MCParExecute(int StepIndex, int nsteps, int indexGen, int N, unsigned int seed, float K, float R, float Time, float SIG, float pS0) {
	int tmp = 1;
	double workTime = 0.0;
	float tmp_price;
	std::vector<double> Times;
	std::vector<double> Prices;
	for (int k = 1; k < 6; ++k) {
		for (int j = 0; j < 5; ++j) {

			tmp_price = getMCPricePar(tmp, StepIndex, nsteps, indexGen, N, seed, K, R, Time, SIG, pS0, workTime);
			Times.push_back(workTime);
			Prices.push_back(tmp_price);
			std::cout << j << " Writing Done " << std::endl;
		}
		std::sort(Times.begin(), Times.end());
		for (std::vector<double>::const_iterator it = Times.begin(); it != Times.end(); ++it)
			std::cout << *it << std::endl;
		std::cout << "####################### " << std::endl;
		for (std::vector<double>::const_iterator it = Prices.begin(); it != Prices.end(); ++it)
			std::cout << *it << std::endl;
		std::cout << tmp << " thread (stockPrice MCPar) Done" << std::endl;
		tmp *= 2;
		Times.clear();
		Prices.clear();
		std::cout << "MCPar Dode \n";
	}

}


std::string createString(double Error, int nsteps, int scale, double time) {
	std::string res;

	res.append(std::to_string(time / static_cast<float>(nsteps / scale)));
	res.append(";");
	res.append(std::to_string(Error));
	res.append(";");
	res.append(std::to_string(log10(time / static_cast<float>(nsteps / scale))));
	res.append(";");
	res.append(std::to_string(log10(Error)));
	res.append(";");
	
	return res;
}



void NumSolution::WriteMethodErrors(float* Errors, int nsteps, int npaths, int nrows, float Time, int scale, int stepIndex) {

	std::string fileName;
	switch (stepIndex) {
	case 0:
		fileName = "EulerMar_";
		break;
	case 1:
		fileName = "Milstein_";
		break;
	case 2:
		fileName = "RK1_";
		break;
	case 3:
		fileName = "BurragePlaten_";
		break;
	}
	std::string date = std::to_string(nsteps).append("_Steps_").append(std::to_string(npaths)).append("_Paths_").append(std::to_string(Time)).append("_Years_").append(getTimestamp(fileName));

	FILE *f = fopen(date.c_str(), "w");
	fprintf(f, "step;e;log(step);log(e);\n");
	for (int i = 0; i < nrows; ++i) {
		std::string tmp_string = createString(Errors[i], nsteps, static_cast<int>(scale / pow(2, i)), Time);
		for (std::string::iterator it = tmp_string.begin(); it<tmp_string.end(); ++it) {
			std::replace(tmp_string.begin(), tmp_string.end(), '.', ',');
		}
		fprintf(f, "%s\n", tmp_string.c_str());
	}
	fclose(f);
}

void NumSolution::getErrors(int nsteps, int indexGen, int N, unsigned int seed, float K, float R, float Time, float SIG, float pS0, int sampleStep, float fair) {

	std::string date = std::to_string(N).append("_Samples_").append(std::to_string(sampleStep)).append("_Step_").append(getTimestamp("Errors_"));
	FILE *f = fopen(date.c_str(), "w");
	fprintf(f, "%sstockPrice;%sstockPrice;%sstockPrice;%sstockPrice\n;", "Euler ", "Mils ", "RK1 ", "Platen ");
	for (int i = 0; i <= N; i += sampleStep) {
		fprintf(f, "%lf;%lf;%lf;%lf\n",
			(getMCPrice(0, nsteps, indexGen, i, seed, K, R, Time, SIG, pS0) - fair),
			(getMCPrice(1, nsteps, indexGen, i, seed, K, R, Time, SIG, pS0) - fair),
			(getMCPrice(2, nsteps, indexGen, i, seed, K, R, Time, SIG, pS0) - fair),
			(getMCPrice(3, nsteps, indexGen, i, seed, K, R, Time, SIG, pS0) - fair));
	}
	fclose(f);

}

void NumSolution::writeMethodConvergence(int StepIndex, int indexGen, int npaths, int nsteps, float pS0, float pR, float pSig, float time, unsigned int seed) {
	VSLStreamStatePtr stream = initGen(seed, indexGen);
	float *Error = new float[12];
	checkConvergence(StepIndex, stream, npaths, nsteps, pS0, pR, pSig, time, Error, seed, indexGen);
	WriteMethodErrors(Error, nsteps, npaths, 12, time, 2048, StepIndex);
	for (int i = 0; i < 12; i++)
		printf("Error %d = %lf\n", i + 1, Error[i]);
	vslDeleteStream(&stream);
	delete[] Error;
}

float NumSolution::Integrand(float z, float pS0, float K, float tmp1, float tmp2, int scale) {
	float payoff;
	float stockPrice = pS0 * expf(tmp1 + tmp2 * z);
	payoff = stockPrice - K;
	payoff > 0.0f ? payoff *= expf(-z * z / 2.0f) : payoff = 0.0f;
	return payoff;
}

void NumSolution::GetRPricePar(float a, float b, int scale, int NumThreads, int N, float r, float sig, float* pT, float* pK, float* pS0, float* pC) {
	float h = (b - a) / scale;
	float sum = 0.f;
	omp_set_nested(1); // does it works?
	omp_set_num_threads(NumThreads);
	float tmp1, tmp2;
#if defined(__INTEL_COMPILER) 
#pragma ivdep
#pragma vector always	
#endif
#pragma omp parallel for private(tmp1, tmp2, sum)
	for (int j = 0; j < N; ++j) {
		tmp1 = (r - sig * sig * 0.5f) * pT[j];
		tmp2 = sig * sqrtf(pT[j]);
		#pragma omp parallel for reduction(+:sum)
		for (int i = 0; i < scale; ++i) {
			sum += Integrand(a + h * (i + 0.5f), pS0[j], pK[j], tmp1, tmp2, scale);
		}
		sum = sum * expf(-r * pT[j]) / sqrtf(2.0f * M_PIF) * h;
		pC[j] = sum;
	}
}

// void NumSolution::GetTPricePar(float a, float b, int scale, int NumThreads, int N, float r, float sig, float* pT, float* pK, float* pS0, float* pC) {
// 	float h = (b - a) / scale;
// 	float sum = 0.f;
// //	omp_set_nested(1); // does it works?
// 	omp_set_num_threads(NumThreads);
// 	float tmp1, tmp2;
// #if defined(__INTEL_COMPILER) 
// #pragma ivdep
// #pragma vector always	
// #endif
// //#pragma omp parallel for private(tmp1, tmp2, sum)
// 	for (int j = 0; j < N; ++j) {
// 		tmp1 = (r - sig * sig * 0.5f) * pT[j];
// 		tmp2 = sig * sqrtf(pT[j]);
// 		sum += Integrand(a, pS0[j], pK[j], tmp1, tmp2, scale) * 0.5f * h + Integrand(b, pS0[j], pK[j], tmp1, tmp2, scale) * 0.5f * (b - a);
// #pragma omp parallel for reduction(+:sum)
// 	for (int i = 1; i < scale; ++i) {
// 		sum += Integrand(a + h * i, pS0[j], pK[j], tmp1, tmp2, scale);
// 	}
// 	sum = sum * expf(-r * pT[j]) / sqrtf(2.0f * M_PIF) * h;
// 	pC[j] = sum;
// 	}
// }

void NumSolution::GetSPricePar(float a, float b, int scale, int NumThreads, int N, float r, float sig, float* pT, float* pK, float* pS0, float* pC) {
	float h = (b - a) / scale;
	float sum2 = 0.f, sum = 0.f, sum4 = 0.f;
//	omp_set_nested(1); // almost always is incorrect!!
	omp_set_num_threads(NumThreads);
	float tmp1, tmp2;
#if defined(__INTEL_COMPILER) 
#pragma ivdep
#pragma vector always	
#endif
//#pragma omp parallel for private(tmp1, tmp2, sum, sum2, sum4)
	for (int j = 0; j < N; ++j) {
		tmp1 = (r - sig * sig * 0.5f) * pT[j];
		tmp2 = sig * sqrtf(pT[j]);
		sum4 += Integrand(a + h, pS0[j], pK[j], tmp1, tmp2, scale);
		sum += Integrand(a, pS0[j], pK[j], tmp1, tmp2, scale) + Integrand(b, pS0[j], pK[j], tmp1, tmp2, scale);
#pragma omp parallel for reduction(+:sum2, sum4)
	for (int i = 1; i < scale - 2; i += 2) {
		sum4 += Integrand(a + (i - 1) * h, pS0[j], pK[j], tmp1, tmp2, scale);
		sum2 += Integrand(a + i * h, pS0[j], pK[j], tmp1, tmp2, scale);
	}

	sum = (sum + 4 * sum4 + 2 * sum2) * expf(-r * pT[j]) / sqrtf(2.0f * M_PIF) * (h / 3.0f);
	pC[j] = sum;
	}
}

// void NumSolution::GetSPrice(float a, float b, int scale, int N, float r, float sig, float* pT, float* pK, float* pS0, float* pC) {
// 	float h = (b - a) / scale;
// 	float sum2 = 0.f, sum = 0.f, sum4 = 0.f;
// 	float tmp1, tmp2;
// #if defined(__INTEL_COMPILER) 
// #pragma ivdep
// #pragma vector always	
// #endif
// 	for (int j = 0; j < N; ++j) {
// 		tmp1 = (r - sig * sig * 0.5f) * pT[j];
// 		tmp2 = sig * sqrtf(pT[j]);
// 		sum4 += Integrand(a + h, pS0[j], pK[j], tmp1, tmp2, scale);
// 		sum += Integrand(a, pS0[j], pK[j], tmp1, tmp2, scale) + Integrand(b, pS0[j], pK[j], tmp1, tmp2, scale);
// 		for (int i = 1; i < scale - 2; i += 2) {
// 			sum4 += Integrand(a + (i - 1)*h, pS0[j], pK[j], tmp1, tmp2, scale);
// 			sum2 += Integrand(a + i * h, pS0[j], pK[j], tmp1, tmp2, scale);
// 		}
// 		sum = (sum + 4 * sum4 + 2 * sum2) * expf(-r * pT[j]) / sqrtf(2.0f * M_PIF) * (h / 3.0f);
// 		pC[j] = sum;
// 	}
// }

void NumSolution::Get3_8PricePar(float a, float b, int scale, int NumThreads, int N, float r, float sig, float* pT, float* pK, float* pS0, float* pC) {
// legacy
	assert(scale / 3 == 0);
	float h = (b - a) / scale;
	float sum3 = 0.f, sum2 = 0.f, sum = 0.f;
	omp_set_num_threads(NumThreads);
	float tmp1, tmp2;
#if defined(__INTEL_COMPILER) 
#pragma ivdep
#pragma vector always	
#endif
	for (int j = 0; j < N; ++j) {
		tmp1 = (r - sig * sig * 0.5f) * pT[j];
		tmp2 = sig * sqrtf(pT[j]);
		sum3 = Integrand(a + h, pS0[j], pK[j], tmp1, tmp2, scale);
		sum += Integrand(a, pS0[j], pK[j], tmp1, tmp2, scale) + Integrand(b, pS0[j], pK[j], tmp1, tmp2, scale);
#pragma omp parallel for reduction(+:sum2, sum3)
		for (int i = 3; i < scale; i += 3) {
			sum2 += Integrand(a + i * h, pS0[j], pK[j], tmp1, tmp2, scale);
			sum3 += Integrand(a + (i - 2) * h, pS0[j], pK[j], tmp1, tmp2, scale);
			sum3 += Integrand(a + (i - 1) * h, pS0[j], pK[j], tmp1, tmp2, scale);
		}
	sum = (sum + 3 * sum3 + 2 * sum2) * expf(-r * pT[j]) / sqrtf(2.0f * M_PIF) * (3.0f * h / 8.0f);
	pC[j] = sum;
	}
}
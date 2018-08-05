#include "EM10_MathStat.h"

RVCharacteristics::RVCharacteristics(float* wd, int _len, float _h) {
	QuickSort(wd, 0, _len - 1);
	len = _len;
	VarTh = _h;
	wd_sorted = wd;
	MeanTh = 0.0;

	MeanCh = this->GetMean();
	VarCh = GetVariance();
	MedCh = GetMediane();
	MedTh = 0.0f;
	VarDiff = abs(VarCh - VarTh);
}

float RVCharacteristics::GetMean() {
	float res = 0.0;
	for (int i = 0; i < len; ++i)
		res += wd_sorted[i];
	return res / (float)len;
}

bool RVCharacteristics::IsExpValCorrect(float _eps) {
	return GetMean() < _eps ? true : false;
}

float RVCharacteristics::GetVariance() {
	float res = 0.0;
	for (int i = 0; i < len; ++i) {
		res += pow(wd_sorted[i] - MeanCh, 2);
	}
	return res / (float)len;
}

bool RVCharacteristics::IsVArianceCorrect(float var, float _eps) {
	return abs(var - GetVariance()) < _eps ? true : false;
}

float RVCharacteristics::GetMediane() {
	
	return len % 2 == 0 ? (wd_sorted[len / 2] + wd_sorted[(len / 2) - 1]) : wd_sorted[((len + 1)/ 2) - 1];
}


void RVCharacteristics::QuickSort(float* array, int first, int last) {

	int i = first;
	int j = last;
	srand(time(0));
	auto random = rand;
	int randElem = random();
	float base = array[(randElem % (last - first)) + first];

	do {
		while (array[i] < base) i++;
		while (array[j] > base) j--;

		if (i <= j) {
			if (array[i] > array[j])
				std::swap(array[i], array[j]);
			i++;
			j--;
		}
	} while (i <= j);

	if (i < last)
		QuickSort(array, i, last);
	if (first < j)
		QuickSort(array, first, j);
}

void RVCharacteristics::WriteToCsv(float eps) {
	time_t rawtime;
	time(&rawtime);
	std::string date = asctime(localtime(&rawtime));
	date.pop_back();
	date.append("_Random_Generator_Statistics.csv");

	for (std::string::iterator it = date.begin(); it<date.end(); ++it) {
		if (*it == ':') {
			date.erase(it);
		}
		std::replace(date.begin(), date.end(), ' ', '_');
	}

	FILE *f = fopen(date.c_str(), "w");
	for (int j = 0; j < len; j++) {
		std::string tmp_cell = std::to_string(wd_sorted[j]);

		for (std::string::iterator it = tmp_cell.begin(); it<tmp_cell.end(); ++it) {
			std::replace(tmp_cell.begin(), tmp_cell.end(), '.', ',');
		}
		fprintf(f, /*"%lf;\n"*/"%s;", /*buffer[j]*/tmp_cell.c_str());
		fprintf(f, "\n");
	}
	fprintf(f, "\n");
	fprintf(f, "ExpectedVal: %lf;\n", MeanCh);
	fprintf(f, "Variance: %lf;\n", VarCh);
	fprintf(f, "Step: %lf;\n", VarTh);
	fprintf(f, "Variance Difference: %lf;\n", VarDiff);
	fprintf(f, "Mediane: %lf;\n", MedCh);
	fprintf(f, "Expected Val correcthness: %i;\n", IsExpValCorrect(_EPSILON_));
	fprintf(f, "Variance correcthness: %i;\n", IsVArianceCorrect(VarTh, _EPSILON_));
	fclose(f);
}

void RVCharacteristics::Execute() {
	//StockPrice sp;
	//VSLStreamStatePtr stream = sp.InitGen();

	//float *wiener_diff = new float[NSTEPS]; // Random values buffer
	//float h = TIME / (float)NSTEPS; // step

	//sp.GenerateGauss(0, sqrt(h), NSTEPS, stream, wiener_diff);
	//WriteToCsv(_EPSILON_);
	//delete[] wiener_diff;
	//sp.FreeGen(stream);
}
#include "MathStat.h"

RVCharacteristics::RVCharacteristics(double* wd, int len) {
	QuickSort(wd, 0, len - 1);
}

double RVCharacteristics::GetExpectedVal(double* wd, int len) {
	double res = 0.0;
	for (int i = 0; i < len; ++i)
		res += wd[i];
	return res / (double)len;
}

bool RVCharacteristics::IsExpValCorrect(double ev, double _eps) {
	return ev < _eps ? true : false;
}

double RVCharacteristics::GetVariance(double* wd, double eval, int len) {
	double res = 0.0;
	for (int i = 0; i < len; ++i) {
		res += pow(wd[i] - eval, 2);
	}
	return res / (double)len;
}

bool RVCharacteristics::IsVArianceCorrect(double var, double var_stat, double _eps) {
	return abs(var - var_stat) < _eps ? true : false;
}

void RVCharacteristics::QuickSort(double* array, int first, int last) {

	int i = first;
	int j = last;
	srand(time(0));
	auto random = rand;
	int randElem = random();
	double base = array[(randElem % (last - first)) + first];

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
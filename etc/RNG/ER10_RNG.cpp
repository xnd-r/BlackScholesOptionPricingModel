#include "ER10_RNG.h"

MCG59::MCG59(long int seed) : a(302875106592253), seed(GetLongSeed(seed)) {
	seed % 2 == 0 ? m = 576460752303423488 : m = 144115188075855872;
};

// need test
long long int MCG59::GetLongSeed(long int seed) {
	return (uint64_t) seed << 32 | seed;
}

void MCG59::RandomArray(float* _array, int len) {
	for (int i = 0; i < len; ++i) {
		seed = (a * seed) % m;
		if (seed < 0)
			seed += m;
		_array[i] = (float)seed / m;
	}
}

float MCG59::GetFloat() {
	seed = (a * seed) % m;
	if (seed < 0)
		seed += m;
	return (float)seed / m;
}

float MCG59::GetFloatFromRange(float min, float max) {
	return floorf(min + GetFloat() * (max - min + 1));
}

void MCG59::WriteToCsv(float* arr, int len) {
	time_t rawtime;
	time(&rawtime);
	std::string date = asctime(localtime(&rawtime));
	date.pop_back();
	date.append("_MCG59_ResultValues.csv");

	for (std::string::iterator it = date.begin(); it<date.end(); ++it) {
		if (*it == ':') {
			date.erase(it);
		}
		std::replace(date.begin(), date.end(), ' ', '_');
	}

	FILE *f = fopen(date.c_str(), "w");
	for (int j = 0; j < len; j++) {
		std::string tmp_cell = std::to_string(arr[j]);

		for (std::string::iterator it = tmp_cell.begin(); it<tmp_cell.end(); ++it)
			std::replace(tmp_cell.begin(), tmp_cell.end(), '.', ',');

		fprintf(f, /*"%lf;\n"*/"%s;\n", /*buffer[j]*/tmp_cell.c_str());
	}
	fprintf(f, "\n");
	fclose(f);
}
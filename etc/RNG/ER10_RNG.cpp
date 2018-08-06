#include "ER10_RNG.h"

MCG59::MCG59(long long int seed) : a(302875106592253), seed(seed) {
	seed % 2 == 0 ? m = 576460752303423488 : m = 144115188075855872;
};

void MCG59::RandomArray(float* _array, int len) {
	for (int i = 0; i < len; ++i) {
		seed = abs((a * seed) % m);
		_array[i] = (float)seed / m;
	}
}

float MCG59::GetFloat() {
	seed = abs((a * seed) % m);
	return (float)seed / m;
}

float MCG59::GetFloatFromRange(float min, float max) {
	return floorf(min + GetFloat() * (max - min + 1));
}

#include "ER20_Fibonacci.h"

FibonacciGen::FibonacciGen(const int _len, long int _seed) : MCG59(_seed) { len = _len; }

void FibonacciGen::RandomArray(float* dest_arr) {
	MCG59::RandomArray(resArray, maxAB);
	for (int i = maxAB + 1; i < len; ++i) {
		GetFloat(i);
	}
	// Why should we use queue based on array?
	//std::queue<float> tmpArray;
	//for (int i = 0; i < maxAB; ++i)
	//	tmpArray.push(resArray[i]);
}

float FibonacciGen::GetFloat(int i) {
	return resArray[i - a] >= resArray[i - b] ? resArray[i - a] - resArray[i - b] : resArray[i - a] - resArray[i - b] + 1.0f;
}
FibonacciGen::~FibonacciGen() {
	delete[] resArray;
}
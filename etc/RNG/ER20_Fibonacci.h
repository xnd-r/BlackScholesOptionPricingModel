#ifndef ____FIBONACCI_GENERATOR____
#define ____FIBONACCI_GENERATOR____

#include <math.h>
#include <queue>
#include "ER10_RNG.h"

class  FibonacciGen : public MCG59 {
	const int a = 55; // lag 
	const int b = 24; // lag
	const int maxAB = a ^ ((a ^ b) & - (a < b)); // max(a, b);
	int len;
public:
	FibonacciGen(int _len, long int _seed);
	float* resArray = new float[len];

	void RandomArray(float* dest_arr);
	float GetFloat(int i);
	~FibonacciGen();
};
#endif // !____FIBONACCI_GENERATOR____

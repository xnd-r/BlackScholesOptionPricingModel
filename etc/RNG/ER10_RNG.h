#ifndef ____RANDOM_NUMBER_GENERATOR____
#define ____RANDOM_NUMBER_GENERATOR____

#include <math.h>

class  MCG59 {
	long long int a;
	long long int m;
	long long int seed;

public:
	MCG59() {};
	MCG59(long long int seed);
	~MCG59() {};
	virtual void RandomArray(float* dest_arr, int len);
	virtual float GetFloat();
	virtual float GetFloatFromRange(float min_float, float max_float); // returns random value from [min_float, max_float]
};

#endif // !____RANDOM_NUMBER_GENERATOR____
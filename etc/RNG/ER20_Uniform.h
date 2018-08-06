#ifndef ____UNIFORM_GENERATOR____
#define ____UNIFORM_GENERATOR____

#include <math.h>
#include "ER10_RNG.h"

class  UnifromGen : public MCG59 {
public:
	void RandomArray(float* dest_arr, int len);
};

#endif // !____UNIFORM_GENERATOR____
#ifndef ____ROW_TABLE_____
#define ____ROW_TABLE_____

#include <string>
#include <math.h>

struct row_table {
	float error;
	int nSteps;
	int scale;
	int nRows;
	float time;

	void SetValues(float e, int nst, int sc, int nr, float t);
	std::string CreateString();

};

#endif // !____ROW_TABLE_____
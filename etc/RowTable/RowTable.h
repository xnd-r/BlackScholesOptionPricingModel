#ifndef ____ROW_TABLE_____
#define ____ROW_TABLE_____

#include <string>
#include <math.h>

struct row_table {
	double error;
	int nSteps;
	int scale;
	int nRows;
	double time;

	void SetValues(double e, int nst, int sc, int nr, double t);
	std::string CreateString();

};

#endif // !____ROW_TABLE_____
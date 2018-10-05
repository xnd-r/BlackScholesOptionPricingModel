#ifndef ____ROW_TABLE_____
#define ____ROW_TABLE_____

#include <string>
#include <math.h>

struct row_table {
	row_table() {};

	double Error;
	int nsteps;
	int scale;
	int nrows;
	float time;

	void setValues(double e, int nst, int sc, int nr, float t);
	std::string createString();
};

#endif // !____ROW_TABLE_____
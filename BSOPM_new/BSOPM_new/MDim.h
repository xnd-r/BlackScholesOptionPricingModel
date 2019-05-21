float getMDimDensity(int dim, float*x);
void getPointCoords(float *tmpVector, long long ind, int nsamples, int dim, float a, float b);
float getRRNG(int dim, int scale, float a, float b, float r, float time, float* SP, float S0, float K);
float getMCPriceParNew(int NumThreads, int StepIndex, int nsteps, int indexGen, int N, int bufsize, unsigned int seed, float K, float R, float Time, float SIG, float pS0, int dim);

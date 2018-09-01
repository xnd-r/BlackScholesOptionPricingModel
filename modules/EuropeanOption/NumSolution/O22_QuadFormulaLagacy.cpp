	float sum = 0.0f;
	float s;
	h = (b - a) / scale;
	float payoff;
	int i, k, ind;
	//int points_amo = (int)(b - a) / h + 1;
	int buflen = 128;// points_amo / 20;
	assert(N % bufsize == 0);
	//assert(points_amo % buflen == 0);
	float* s_array = new float[buflen];
	start = clock();
//#pragma omp parallel private(s, portion, i)
//	{
//#pragma omp parallel for private(s, i) reduction(+:sum)
	for (ind = 0; ind < N; ind++) {
		for (int k = 0; k < scale; k += buflen) {
			for (i = 0; i < buflen; ++i) {
				s_array[i] = S0 * expf(tmp1 + tmp2 * (a + h * (k + i + 0.5f)));
				payoff = s_array[i] - K;
				payoff > 0.0f ? payoff *= expf(-(a + h * (k + i + 0.5f)) * (a + h * (k + i + 0.5f)) / 2.0f) : payoff = 0.0f;
				sum += payoff;
			}
		}
	}
//	}
	sum = sum * expf(-R * TIME) / sqrtf(2.0f * M_PIF) * h / N ;
	//std::cout << k / buflen << " ########### " << std::endl;

	finish = clock();
	t = (double)(finish - start) / CLOCKS_PER_SEC;
	delete[] s_array;
	return sum;
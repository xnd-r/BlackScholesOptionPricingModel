// Copyright 2018 Romanov Alexander

#include "gtest.h"

#include "../modules/EuropeanOption/AnSolution/O11_CallOption.h"

using ::testing::internal::RE;

class Romanov_Alexander_Call_option : public ::testing::Test {
protected:
	float pT, pK, pS0;
	float pC;
	CallOption co;

	void set_values() {
		pT = 3.0f;
		pK = 100.0f;
		pS0 = 100.0f;
		#define N 1
		//pC = 20.9244f;
	}

};

TEST_F(Romanov_Alexander_Call_option, test_base_version_0)
{
	// Arrange
	set_values();

	// Act
	(this->co->*version_array[0])(pT, pK, pS0, pC);
	 
	// Assert
	ASSERT_EQ(pC, 20.9244);
}

TEST_F(Romanov_Alexander_Call_option, test_floating_point_1)
{
	// Arrange
	set_values();

	// Act
	version_array[1](&pT, &pK, &pS0, &pC);

	// Assert
	ASSERT_EQ(pC, 20.9244);
}

TEST_F(Romanov_Alexander_Call_option, test_erf_2)
{
	// Arrange
	set_values();

	// Act
	version_array[2](&pT, &pK, &pS0, &pC);

	// Assert
	ASSERT_EQ(pC, 20.9244);
}

#if defined(__INTEL_COMPILER) 
TEST_F(Romanov_Alexander_Call_option, test_restrict_3)
{
	// Arrange
	set_values();

	// Act
	option_array[3](&pT, &pK, &pS0, &pC);

	// Assert
	ASSERT_EQ(pC, 20.9244);
}
#endif

TEST_F(Romanov_Alexander_Call_option, test_fp_simd_vector_always_4)
{
	// Arrange
	set_values();

	// Act
	version_array[4](&pT, &pK, &pS0, &pC);

	// Assert
	ASSERT_EQ(pC, 20.9244);
}

TEST_F(Romanov_Alexander_Call_option, test_fp_erf_simd_vector_always_5)
{
	// Arrange
	set_values();

	// Act
	version_array[5](&pT, &pK, &pS0, &pC);

	// Assert
	ASSERT_EQ(pC, 20.9244);
}

TEST_F(Romanov_Alexander_Call_option, test_defined_invsqrt2_6)
{
	// Arrange
	set_values();

	// Act
	version_array[6](&pT, &pK, &pS0, &pC);

	// Assert
	ASSERT_EQ(pC, 20.9244);
}

TEST_F(Romanov_Alexander_Call_option, test_omp_8)
{
	// Arrange
	set_values();

	// Act
	version_array[7](&pT, &pK, &pS0, &pC);

	// Assert
	ASSERT_EQ(pC, 20.9244);
}

TEST_F(Romanov_Alexander_Call_option, test_omp_nontemporal_9)
{
	// Arrange
	set_values();

	// Act
	version_array[8](&pT, &pK, &pS0, &pC);

	// Assert
	ASSERT_EQ(pC, 20.9244);
}
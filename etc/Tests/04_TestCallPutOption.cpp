// Copyright 2018 Romanov Alexander

#include "gtest.h"
#include "../modules/EuropeanOption/AnSolution/O12_CallPutOption.h"
//#include "../../modules/EuropeanOption/AnSolution/O12_CallPutOption.h"
//#include "../../3rdparty/gtest/gtest.h"
using ::testing::internal::RE;

class Romanov_Alexander_Call_put_option : public ::testing::Test {
protected:
	float pT, pK, pS0;
	float pC, pP;
	CallPutOption cpo;

	void set_values() {
		pT = 3.0f;
		pK = 100.0f;
		pS0 = 100.0f;
		#define N 1
		//pC = 20.9244f;
		//pP = 8.66673f;
	}

};

//TEST_F(Romanov_Alexander_Call_put_option, test_base_version_0)
//{
//	// Arrange
//	set_values();
//
//	// Act
//	version_array[0](&pT, &pK, &pS0, &pC, &pP);
//
//	// Assert
//	ASSERT_EQ(pC, 20.9244);
//	ASSERT_EQ(pP, 8.66673);
//}
//
//TEST_F(Romanov_Alexander_Call_put_option, test_floating_point_1)
//{
//	// Arrange
//	set_values();
//
//	// Act
//	version_array[1](&pT, &pK, &pS0, &pC, &pP);
//
//	// Assert
//	ASSERT_EQ(pC, 20.9244);
//	ASSERT_EQ(pP, 8.66673);
//}
//
//TEST_F(Romanov_Alexander_Call_put_option, test_erf_2)
//{
//	// Arrange
//	set_values();
//
//	// Act
//	version_array[2](&pT, &pK, &pS0, &pC, &pP);
//
//	// Assert
//	ASSERT_EQ(pC, 20.9244);
//	ASSERT_EQ(pP, 8.66673);
//}
//
//TEST_F(Romanov_Alexander_Call_put_option, test_restrict_3)
//{
//	// Arrange
//	set_values();
//
//	// Act
//	version_array[3](&pT, &pK, &pS0, &pC, &pP);
//
//	// Assert
//	ASSERT_EQ(pC, 20.9244);
//	ASSERT_EQ(pP, 8.66673);
//}
//
//TEST_F(Romanov_Alexander_Call_put_option, test_fp_simd_vector_always_4)
//{
//	// Arrange
//	set_values();
//
//	// Act
//	version_array[4](&pT, &pK, &pS0, &pC, &pP);
//
//	// Assert
//	ASSERT_EQ(pC, 20.9244);
//	ASSERT_EQ(pP, 8.66673);
//}
//
//TEST_F(Romanov_Alexander_Call_put_option, test_fp_erf_simd_vector_always_5)
//{
//	// Arrange
//	set_values();
//
//	// Act
//	version_array[5](&pT, &pK, &pS0, &pC, &pP);
//
//	// Assert
//	ASSERT_EQ(pC, 20.9244);
//	ASSERT_EQ(pP, 8.66673);
//}
//
//TEST_F(Romanov_Alexander_Call_put_option, test_defined_invsqrt2_6)
//{
//	// Arrange
//	set_values();
//
//	// Act
//	version_array[6](&pT, &pK, &pS0, &pC, &pP);
//
//	// Assert
//	ASSERT_EQ(pC, 20.9244);
//	ASSERT_EQ(pP, 8.66673);
//}
//
//TEST_F(Romanov_Alexander_Call_put_option, test_omp_8)
//{
//	// Arrange
//	set_values();
//
//	// Act
//	version_array[7](&pT, &pK, &pS0, &pC, &pP);
//
//	// Assert
//	ASSERT_EQ(pC, 20.9244);
//	ASSERT_EQ(pP, 8.66673);
//}
//
//TEST_F(Romanov_Alexander_Call_put_option, test_omp_nontemporal_9)
//{
//	// Arrange
//	set_values();
//
//	// Act
//	version_array[8](&pT, &pK, &pS0, &pC, &pP);
//
//	// Assert
//	ASSERT_EQ(pC, 20.9244);
//	ASSERT_EQ(pP, 8.66673);
//}
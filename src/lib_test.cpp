#include <iostream>
#include <Eigen/Dense>
#include "../include/modern_robotics.h"
#include "gtest/gtest.h"

# define M_PI           3.14159265358979323846  /* pi */

TEST(MRTest, VecToSO3Test) {
	mr::Vector3x vec(1, 2, 3);
	mr::Matrix3x result(3, 3);
	result << 0, -3, 2, 3, 0, -1, -2, 1, 0;
	EXPECT_EQ(result, mr::VecToso3(vec));
}

TEST(MRTest, JacobianSpaceTest) {
	mr::MatrixXx s_list(6, 3);
	s_list << 0, 0, 0,
		0, 1, -1,
		1, 0, 0,
		0, -0.0711, 0.0711,
		0, 0, 0,
		0, 0, -0.2795;
	mr::VectorXx theta(3);
	theta << 1.0472, 1.0472, 1.0472;
	mr::MatrixXx result(6, 3);
	result << 0, -0.866, 0.866,
		0, 0.5, -0.5,
		1, 0, 0,
		0, -0.0355, -0.0855,
		0, -0.0615, -0.1481,
		0, 0, -0.1398;
	mr::MatrixXx tmp_result = mr::JacobianSpace(s_list, theta);
	// std::cout << tmp_result << std::endl;
	ASSERT_TRUE(mr::JacobianSpace(s_list, theta).isApprox(result, 4));
}


TEST(MRTest, JacobianBodyTest) {
	mr::MatrixXx b_list(6, 3);
	b_list << 0, 0, 0,
		0, 1, -1,
		1, 0, 0,
		0.0425, 0, 0,
		0.5515, 0, 0,
		0, -0.5515, 0.2720;
	mr::VectorXx theta(3);
	theta << 0, 0, 1.5708;
	mr::MatrixXx result(6, 3);
	result << 1, 0, 0,
		0, 1, -1,
		0, 0, 0,
		0, -0.2795, 0,
		0.2795, 0, 0,
		-0.0425, -0.2720, 0.2720;
	mr::MatrixXx tmp_result = mr::JacobianBody(b_list, theta);
	// std::cout << tmp_result << std::endl;
	ASSERT_TRUE(mr::JacobianBody(b_list, theta).isApprox(result, 4));
}

TEST(MRTest, adTest) {
	mr::VectorXx V(6);
	V << 1, 2, 3, 4, 5, 6;

	mr::MatrixXx result(6, 6);
	result << 0, -3, 2, 0, 0, 0,
		3, 0, -1, 0, 0, 0,
		-2, 1, 0, 0, 0, 0,
		0, -6, 5, 0, -3, 2,
		6, 0, -4, 3, 0, -1,
		-5, 4, 0, -2, 1, 0;

	ASSERT_TRUE(mr::ad(V).isApprox(result, 4));
}

TEST(MRTest, TransInvTest) {
	mr::MatrixXx input(4, 4);
	input << 1, 0, 0, 0,
		0, 0, -1, 0,
		0, 1, 0, 3,
		0, 0, 0, 1;
	mr::MatrixXx result(4, 4);
	result << 1, 0, 0, 0,
		0, 0, 1, -3,
		0, -1, 0, 0,
		0, 0, 0, 1;

	auto inv = mr::TransInv(input);
	ASSERT_TRUE(inv.isApprox(result, 4));
}

TEST(MRTest, RotInvTest) {
	mr::MatrixXx input(3, 3);
	input << 0, 0, 1,
		1, 0, 0,
		0, 1, 0;
	mr::MatrixXx result(3, 3);
	result << 0, 1, 0,
		0, 0, 1,
		1, 0, 0;

	auto inv = mr::RotInv(input);
	ASSERT_TRUE(inv.isApprox(result, 4));
}

TEST(MRTest, ScrewToAxisTest) {
	mr::Vector3x q, s;
	q << 3, 0, 1;
	s << 0, 0, 1;
	double h = 2;

	mr::VectorXx axis = mr::ScrewToAxis(q, s, h);
	mr::VectorXx result(6);
	result << 0, 0, 1, 0, -3, 2;

	ASSERT_TRUE(axis.isApprox(result, 4));
}

TEST(MRTest, FKInBodyTest) {
	mr::MatrixXx M(4, 4);
	M << -1, 0, 0, 0,
		0, 1, 0, 6,
		0, 0, -1, 2,
		0, 0, 0, 1;
	mr::MatrixXx Blist(6, 3);
	Blist << 0, 0, 0,
		0, 0, 0,
		-1, 0, 1,
		2, 0, 0,
		0, 1, 0,
		0, 0, 0.1;
	mr::VectorXx thetaList(3);
	thetaList << M_PI / 2.0, 3, M_PI;

	mr::MatrixXx result(4, 4);
	result << 0, 1, 0, -5,
		1, 0, 0, 4,
		0, 0, -1, 1.68584073,
		0, 0, 0, 1;
	mr::MatrixXx FKCal = mr::FKinBody(M, Blist, thetaList);

	ASSERT_TRUE(FKCal.isApprox(result, 4));
}

TEST(MRTest, FKInSpaceTest) {
	mr::MatrixXx M(4, 4);
	M << -1, 0, 0, 0,
		0, 1, 0, 6,
		0, 0, -1, 2,
		0, 0, 0, 1;
	mr::MatrixXx Slist(6, 3);
	Slist << 0, 0, 0,
		0, 0, 0,
		1, 0, -1,
		4, 0, -6,
		0, 1, 0,
		0, 0, -0.1;
	mr::VectorXx thetaList(3);
	thetaList << M_PI / 2.0, 3, M_PI;

	mr::MatrixXx result(4, 4);
	result << 0, 1, 0, -5,
		1, 0, 0, 4,
		0, 0, -1, 1.68584073,
		0, 0, 0, 1;
	mr::MatrixXx FKCal = mr::FKinBody(M, Slist, thetaList);

	ASSERT_TRUE(FKCal.isApprox(result, 4));
}

TEST(MRTest, AxisAng6Test) {
	mr::VectorXx input(6);
	mr::VectorXx result(7);
	input << 1.0, 0.0, 0.0, 1.0, 2.0, 3.0;
	result << 1.0, 0.0, 0.0, 1.0, 2.0, 3.0, 1.0;

	mr::VectorXx output = mr::AxisAng6(input);
	ASSERT_TRUE(output.isApprox(result, 4));
}

TEST(MRTest, MatrixLog6Test) {
	mr::MatrixXx Tinput(4, 4);
	mr::MatrixXx result(4, 4);
	Tinput << 1, 0, 0, 0,
		0, 0, -1, 0,
		0, 1, 0, 3,
		0, 0, 0, 1;

	result << 0, 0, 0, 0,
		0, 0, -1.57079633, 2.35619449,
		0, 1.57079633, 0, 2.35619449,
		0, 0, 0, 0;

	mr::MatrixXx Toutput = mr::MatrixLog6(Tinput);
	ASSERT_TRUE(Toutput.isApprox(result, 4));
}

TEST(MRTest, DistanceToSO3Test) {
	mr::Matrix3x input;
	double result = 0.088353;
	input << 1.0, 0.0, 0.0,
		0.0, 0.1, -0.95,
		0.0, 1.0, 0.1;
	EXPECT_NEAR(result, mr::DistanceToSO3(input), 3);
}

TEST(MRTest, DistanceToSE3Test) {
	mr::Matrix4x input;
	double result = 0.134931;
	input << 1.0, 0.0, 0.0, 1.2,
		0.0, 0.1, -0.95, 1.5,
		0.0, 1.0, 0.1, -0.9,
		0.0, 0.0, 0.1, 0.98;
	EXPECT_NEAR(result, mr::DistanceToSE3(input), 3);
}

TEST(MRTest, TestIfSO3Test) {
	mr::Matrix3x input;
	bool result = false;
	input << 1.0, 0.0, 0.0,
		0.0, 0.1, -0.95,
		0.0, 1.0, 0.1;
	ASSERT_EQ(result, mr::TestIfSO3(input));
}

TEST(MRTest, TestIfSE3Test) {
	mr::Matrix4x input;
	bool result = false;
	input << 1.0, 0.0, 0.0, 1.2,
		0.0, 0.1, -0.95, 1.5,
		0.0, 1.0, 0.1, -0.9,
		0.0, 0.0, 0.1, 0.98;
	ASSERT_EQ(result, mr::TestIfSE3(input));
}

TEST(MRTest, IKinBodyTest) {
	mr::MatrixXx BlistT(3, 6);
	BlistT << 0, 0, -1, 2, 0, 0,
		0, 0, 0, 0, 1, 0,
		0, 0, 1, 0, 0, 0.1;
	mr::MatrixXx Blist = BlistT.transpose();
	mr::Matrix4x M;
	M << -1, 0, 0, 0,
		0, 1, 0, 6,
		0, 0, -1, 2,
		0, 0, 0, 1;
	mr::Matrix4x T;
	T << 0, 1, 0, -5,
		1, 0, 0, 4,
		0, 0, -1, 1.6858,
		0, 0, 0, 1;
	mr::VectorXx thetalist(3);
	thetalist << 1.5, 2.5, 3;
	mr::numericX eomg = 0.01;
	mr::numericX ev = 0.001;
	bool b_result = true;
	mr::VectorXx theta_result(3);
	theta_result << 1.57073819, 2.999667, 3.14153913;
	bool iRet = mr::IKinBody(Blist, M, T, thetalist, eomg, ev);
	ASSERT_EQ(b_result, iRet);
	ASSERT_TRUE(thetalist.isApprox(theta_result, 4));
}

TEST(MRTest, IKinSpaceTest) {
	mr::MatrixXx SlistT(3, 6);
	SlistT << 0, 0, 1, 4, 0, 0,
		0, 0, 0, 0, 1, 0,
		0, 0, -1, -6, 0, -0.1;
	mr::MatrixXx Slist = SlistT.transpose();
	mr::Matrix4x M;
	M << -1, 0, 0, 0,
		0, 1, 0, 6,
		0, 0, -1, 2,
		0, 0, 0, 1;
	mr::Matrix4x T;
	T << 0, 1, 0, -5,
		1, 0, 0, 4,
		0, 0, -1, 1.6858,
		0, 0, 0, 1;
	mr::VectorXx thetalist(3);
	thetalist << 1.5, 2.5, 3;
	mr::numericX eomg = 0.01;
	mr::numericX ev = 0.001;
	bool b_result = true;
	mr::VectorXx theta_result(3);
	theta_result << 1.57073783, 2.99966384, 3.1415342;
	bool iRet = mr::IKinSpace(Slist, M, T, thetalist, eomg, ev);
	ASSERT_EQ(b_result, iRet);
	ASSERT_TRUE(thetalist.isApprox(theta_result, 4));
}

TEST(MRTest, AdjointTest) {
	mr::Matrix4x T;
	T << 1, 0, 0, 0,
		0, 0, -1, 0,
		0, 1, 0, 3,
		0, 0, 0, 1;
	mr::MatrixXx result(6, 6);
	result <<
		1, 0, 0, 0, 0, 0,
		0, 0, -1, 0, 0, 0,
		0, 1, 0, 0, 0, 0,
		0, 0, 3, 1, 0, 0,
		3, 0, 0, 0, 0, -1,
		0, 0, 0, 0, 1, 0;

	ASSERT_TRUE(mr::Adjoint(T).isApprox(result, 4));
}

TEST(MRTest, InverseDynamicsTest) {
	mr::VectorXx thetalist(3);
	thetalist << 0.1, 0.1, 0.1;
	mr::VectorXx dthetalist(3);
	dthetalist << 0.1, 0.2, 0.3;
	mr::VectorXx ddthetalist(3);
	ddthetalist << 2, 1.5, 1;
	mr::VectorXx g(3);
	g << 0, 0, -9.8;
	mr::VectorXx Ftip(6);
	Ftip << 1, 1, 1, 1, 1, 1;

	std::vector<mr::MatrixXx> Mlist;
	std::vector<mr::MatrixXx> Glist;

	mr::Matrix4x M01;
	M01 << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0.089159,
		0, 0, 0, 1;
	mr::Matrix4x M12;
	M12 << 0, 0, 1, 0.28,
		0, 1, 0, 0.13585,
		-1, 0, 0, 0,
		0, 0, 0, 1;
	mr::Matrix4x M23;
	M23 << 1, 0, 0, 0,
		0, 1, 0, -0.1197,
		0, 0, 1, 0.395,
		0, 0, 0, 1;
	mr::Matrix4x M34;
	M34 << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0.14225,
		0, 0, 0, 1;

	Mlist.push_back(M01);
	Mlist.push_back(M12);
	Mlist.push_back(M23);
	Mlist.push_back(M34);

	mr::VectorXx G1(6);
	G1 << 0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7;
	mr::VectorXx G2(6);
	G2 << 0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393;
	mr::VectorXx G3(6);
	G3 << 0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275;

	Glist.push_back(G1.asDiagonal());
	Glist.push_back(G2.asDiagonal());
	Glist.push_back(G3.asDiagonal());

	mr::MatrixXx SlistT(3, 6);
	SlistT << 1, 0, 1, 0, 1, 0,
		0, 1, 0, -0.089, 0, 0,
		0, 1, 0, -0.089, 0, 0.425;
	mr::MatrixXx Slist = SlistT.transpose();

	mr::VectorXx taulist = mr::InverseDynamics(thetalist, dthetalist, ddthetalist, g,
		Ftip, Mlist, Glist, Slist);

	mr::VectorXx result(3);
	result << 74.6962, -33.0677, -3.23057;

	ASSERT_TRUE(taulist.isApprox(result, 4));
}

TEST(MRTest, GravityForcesTest) {
	mr::VectorXx thetalist(3);
	thetalist << 0.1, 0.1, 0.1;
	mr::VectorXx g(3);
	g << 0, 0, -9.8;

	std::vector<mr::MatrixXx> Mlist;
	std::vector<mr::MatrixXx> Glist;

	mr::Matrix4x M01;
	M01 << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0.089159,
		0, 0, 0, 1;
	mr::Matrix4x M12;
	M12 << 0, 0, 1, 0.28,
		0, 1, 0, 0.13585,
		-1, 0, 0, 0,
		0, 0, 0, 1;
	mr::Matrix4x M23;
	M23 << 1, 0, 0, 0,
		0, 1, 0, -0.1197,
		0, 0, 1, 0.395,
		0, 0, 0, 1;
	mr::Matrix4x M34;
	M34 << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0.14225,
		0, 0, 0, 1;

	Mlist.push_back(M01);
	Mlist.push_back(M12);
	Mlist.push_back(M23);
	Mlist.push_back(M34);

	mr::VectorXx G1(6);
	G1 << 0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7;
	mr::VectorXx G2(6);
	G2 << 0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393;
	mr::VectorXx G3(6);
	G3 << 0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275;

	Glist.push_back(G1.asDiagonal());
	Glist.push_back(G2.asDiagonal());
	Glist.push_back(G3.asDiagonal());

	mr::MatrixXx SlistT(3, 6);
	SlistT << 1, 0, 1, 0, 1, 0,
		0, 1, 0, -0.089, 0, 0,
		0, 1, 0, -0.089, 0, 0.425;
	mr::MatrixXx Slist = SlistT.transpose();

	mr::VectorXx grav = mr::GravityForces(thetalist, g, Mlist, Glist, Slist);

	mr::VectorXx result(3);
	result << 28.4033, -37.6409, -5.4416;

	ASSERT_TRUE(grav.isApprox(result, 4));
}

TEST(MRTest, MassMatrixTest) {
	mr::VectorXx thetalist(3);
	thetalist << 0.1, 0.1, 0.1;

	std::vector<mr::MatrixXx> Mlist;
	std::vector<mr::MatrixXx> Glist;

	mr::Matrix4x M01;
	M01 << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0.089159,
		0, 0, 0, 1;
	mr::Matrix4x M12;
	M12 << 0, 0, 1, 0.28,
		0, 1, 0, 0.13585,
		-1, 0, 0, 0,
		0, 0, 0, 1;
	mr::Matrix4x M23;
	M23 << 1, 0, 0, 0,
		0, 1, 0, -0.1197,
		0, 0, 1, 0.395,
		0, 0, 0, 1;
	mr::Matrix4x M34;
	M34 << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0.14225,
		0, 0, 0, 1;

	Mlist.push_back(M01);
	Mlist.push_back(M12);
	Mlist.push_back(M23);
	Mlist.push_back(M34);

	mr::VectorXx G1(6);
	G1 << 0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7;
	mr::VectorXx G2(6);
	G2 << 0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393;
	mr::VectorXx G3(6);
	G3 << 0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275;

	Glist.push_back(G1.asDiagonal());
	Glist.push_back(G2.asDiagonal());
	Glist.push_back(G3.asDiagonal());

	mr::MatrixXx SlistT(3, 6);
	SlistT << 1, 0, 1, 0, 1, 0,
		0, 1, 0, -0.089, 0, 0,
		0, 1, 0, -0.089, 0, 0.425;
	mr::MatrixXx Slist = SlistT.transpose();

	mr::MatrixXx M = mr::MassMatrix(thetalist, Mlist, Glist, Slist);

	mr::MatrixXx result(3, 3);
	result << 22.5433, -0.3071, -0.0072,
		-0.3071, 1.9685, 0.4322,
		-0.0072, 0.4322, 0.1916;

	ASSERT_TRUE(M.isApprox(result, 4));
}

TEST(MRTest, VelQuadraticForcesTest) {
	mr::VectorXx thetalist(3);
	thetalist << 0.1, 0.1, 0.1;
	mr::VectorXx dthetalist(3);
	dthetalist << 0.1, 0.2, 0.3;

	std::vector<mr::MatrixXx> Mlist;
	std::vector<mr::MatrixXx> Glist;

	mr::Matrix4x M01;
	M01 << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0.089159,
		0, 0, 0, 1;
	mr::Matrix4x M12;
	M12 << 0, 0, 1, 0.28,
		0, 1, 0, 0.13585,
		-1, 0, 0, 0,
		0, 0, 0, 1;
	mr::Matrix4x M23;
	M23 << 1, 0, 0, 0,
		0, 1, 0, -0.1197,
		0, 0, 1, 0.395,
		0, 0, 0, 1;
	mr::Matrix4x M34;
	M34 << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0.14225,
		0, 0, 0, 1;

	Mlist.push_back(M01);
	Mlist.push_back(M12);
	Mlist.push_back(M23);
	Mlist.push_back(M34);

	mr::VectorXx G1(6);
	G1 << 0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7;
	mr::VectorXx G2(6);
	G2 << 0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393;
	mr::VectorXx G3(6);
	G3 << 0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275;

	Glist.push_back(G1.asDiagonal());
	Glist.push_back(G2.asDiagonal());
	Glist.push_back(G3.asDiagonal());

	mr::MatrixXx SlistT(3, 6);
	SlistT << 1, 0, 1, 0, 1, 0,
		0, 1, 0, -0.089, 0, 0,
		0, 1, 0, -0.089, 0, 0.425;
	mr::MatrixXx Slist = SlistT.transpose();

	mr::VectorXx c = mr::VelQuadraticForces(thetalist, dthetalist, Mlist, Glist, Slist);

	mr::VectorXx result(3);
	result << 0.2645, -0.0551, -0.0069;

	ASSERT_TRUE(c.isApprox(result, 4));
}

TEST(MRTest, EndEffectorForcesTest) {
	mr::VectorXx thetalist(3);
	thetalist << 0.1, 0.1, 0.1;
	mr::VectorXx Ftip(6);
	Ftip << 1, 1, 1, 1, 1, 1;

	std::vector<mr::MatrixXx> Mlist;
	std::vector<mr::MatrixXx> Glist;

	mr::Matrix4x M01;
	M01 << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0.089159,
		0, 0, 0, 1;
	mr::Matrix4x M12;
	M12 << 0, 0, 1, 0.28,
		0, 1, 0, 0.13585,
		-1, 0, 0, 0,
		0, 0, 0, 1;
	mr::Matrix4x M23;
	M23 << 1, 0, 0, 0,
		0, 1, 0, -0.1197,
		0, 0, 1, 0.395,
		0, 0, 0, 1;
	mr::Matrix4x M34;
	M34 << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0.14225,
		0, 0, 0, 1;

	Mlist.push_back(M01);
	Mlist.push_back(M12);
	Mlist.push_back(M23);
	Mlist.push_back(M34);

	mr::VectorXx G1(6);
	G1 << 0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7;
	mr::VectorXx G2(6);
	G2 << 0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393;
	mr::VectorXx G3(6);
	G3 << 0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275;

	Glist.push_back(G1.asDiagonal());
	Glist.push_back(G2.asDiagonal());
	Glist.push_back(G3.asDiagonal());

	mr::MatrixXx SlistT(3, 6);
	SlistT << 1, 0, 1, 0, 1, 0,
		0, 1, 0, -0.089, 0, 0,
		0, 1, 0, -0.089, 0, 0.425;
	mr::MatrixXx Slist = SlistT.transpose();

	mr::VectorXx JTFtip = mr::EndEffectorForces(thetalist, Ftip, Mlist, Glist, Slist);

	mr::VectorXx result(3);
	result << 1.4095, 1.8577, 1.3924;

	ASSERT_TRUE(JTFtip.isApprox(result, 4));
}


TEST(MRTest, ForwardDynamicsTest) {
	mr::VectorXx thetalist(3);
	thetalist << 0.1, 0.1, 0.1;
	mr::VectorXx dthetalist(3);
	dthetalist << 0.1, 0.2, 0.3;
	mr::VectorXx taulist(3);
	taulist << 0.5, 0.6, 0.7;
	mr::VectorXx g(3);
	g << 0, 0, -9.8;
	mr::VectorXx Ftip(6);
	Ftip << 1, 1, 1, 1, 1, 1;

	std::vector<mr::MatrixXx> Mlist;
	std::vector<mr::MatrixXx> Glist;

	mr::Matrix4x M01;
	M01 << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0.089159,
		0, 0, 0, 1;
	mr::Matrix4x M12;
	M12 << 0, 0, 1, 0.28,
		0, 1, 0, 0.13585,
		-1, 0, 0, 0,
		0, 0, 0, 1;
	mr::Matrix4x M23;
	M23 << 1, 0, 0, 0,
		0, 1, 0, -0.1197,
		0, 0, 1, 0.395,
		0, 0, 0, 1;
	mr::Matrix4x M34;
	M34 << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0.14225,
		0, 0, 0, 1;

	Mlist.push_back(M01);
	Mlist.push_back(M12);
	Mlist.push_back(M23);
	Mlist.push_back(M34);

	mr::VectorXx G1(6);
	G1 << 0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7;
	mr::VectorXx G2(6);
	G2 << 0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393;
	mr::VectorXx G3(6);
	G3 << 0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275;

	Glist.push_back(G1.asDiagonal());
	Glist.push_back(G2.asDiagonal());
	Glist.push_back(G3.asDiagonal());

	mr::MatrixXx SlistT(3, 6);
	SlistT << 1, 0, 1, 0, 1, 0,
		0, 1, 0, -0.089, 0, 0,
		0, 1, 0, -0.089, 0, 0.425;
	mr::MatrixXx Slist = SlistT.transpose();

	mr::VectorXx ddthetalist = mr::ForwardDynamics(thetalist, dthetalist, taulist, g,
		Ftip, Mlist, Glist, Slist);

	mr::VectorXx result(3);
	result << -0.9739, 25.5847, -32.9150;

	ASSERT_TRUE(ddthetalist.isApprox(result, 4));
}

TEST(MRTest, EulerStepTest) {
	mr::VectorXx thetalist(3);
	thetalist << 0.1, 0.1, 0.1;
	mr::VectorXx dthetalist(3);
	dthetalist << 0.1, 0.2, 0.3;
	mr::VectorXx ddthetalist(3);
	ddthetalist << 2, 1.5, 1;
	mr::numericX dt = 0.1;

	mr::EulerStep(thetalist, dthetalist, ddthetalist, dt);

	mr::VectorXx result_thetalistNext(3);
	result_thetalistNext << 0.11, 0.12, 0.13;
	mr::VectorXx result_dthetalistNext(3);
	result_dthetalistNext << 0.3, 0.35, 0.4;

	ASSERT_TRUE(thetalist.isApprox(result_thetalistNext, 4));
	ASSERT_TRUE(dthetalist.isApprox(result_dthetalistNext, 4));
}

TEST(MRTest, ComputedTorqueTest) {
	mr::VectorXx thetalist(3);
	thetalist << 0.1, 0.1, 0.1;
	mr::VectorXx dthetalist(3);
	dthetalist << 0.1, 0.2, 0.3;
	mr::VectorXx eint(3);
	eint << 0.2, 0.2, 0.2;
	mr::VectorXx g(3);
	g << 0, 0, -9.8;

	std::vector<mr::MatrixXx> Mlist;
	std::vector<mr::MatrixXx> Glist;

	mr::Matrix4x M01;
	M01 << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0.089159,
		0, 0, 0, 1;
	mr::Matrix4x M12;
	M12 << 0, 0, 1, 0.28,
		0, 1, 0, 0.13585,
		-1, 0, 0, 0,
		0, 0, 0, 1;
	mr::Matrix4x M23;
	M23 << 1, 0, 0, 0,
		0, 1, 0, -0.1197,
		0, 0, 1, 0.395,
		0, 0, 0, 1;
	mr::Matrix4x M34;
	M34 << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0.14225,
		0, 0, 0, 1;

	Mlist.push_back(M01);
	Mlist.push_back(M12);
	Mlist.push_back(M23);
	Mlist.push_back(M34);

	mr::VectorXx G1(6);
	G1 << 0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7;
	mr::VectorXx G2(6);
	G2 << 0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393;
	mr::VectorXx G3(6);
	G3 << 0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275;

	Glist.push_back(G1.asDiagonal());
	Glist.push_back(G2.asDiagonal());
	Glist.push_back(G3.asDiagonal());

	mr::MatrixXx SlistT(3, 6);
	SlistT << 1, 0, 1, 0, 1, 0,
		0, 1, 0, -0.089, 0, 0,
		0, 1, 0, -0.089, 0, 0.425;
	mr::MatrixXx Slist = SlistT.transpose();

	mr::VectorXx thetalistd(3);
	thetalistd << 1.0, 1.0, 1.0;
	mr::VectorXx dthetalistd(3);
	dthetalistd << 2, 1.2, 2;
	mr::VectorXx ddthetalistd(3);
	ddthetalistd << 0.1, 0.1, 0.1;
	mr::numericX Kp = 1.3;
	mr::numericX Ki = 1.2;
	mr::numericX Kd = 1.1;

	mr::VectorXx taulist = mr::ComputedTorque(thetalist, dthetalist, eint, g,
		Mlist, Glist, Slist, thetalistd, dthetalistd, ddthetalistd, Kp, Ki, Kd);

	mr::VectorXx result(3);
	result << 133.00525246, -29.94223324, -3.03276856;

	ASSERT_TRUE(taulist.isApprox(result, 4));
}

TEST(MRTest, CubicTimeScalingTest) {
	mr::numericX Tf = 2.0;
	mr::numericX t = 0.6;
	mr::numericX result = 0.216;

	EXPECT_NEAR(result, mr::CubicTimeScaling(Tf, t), 3);
}

TEST(MRTest, QuinticTimeScalingTest) {
	mr::numericX Tf = 2.0;
	mr::numericX t = 0.6;
	mr::numericX result = 0.16308;

	EXPECT_NEAR(result, mr::QuinticTimeScaling(Tf, t), 3);
}

TEST(MRTest, JointTrajectoryTest) {
	int dof = 8;
	mr::VectorXx thetastart(dof);
	thetastart << 1, 0, 0, 1, 1, 0.2, 0, 1;
	mr::VectorXx thetaend(dof);
	thetaend << 1.2, 0.5, 0.6, 1.1, 2, 2, 0.9, 1;
	mr::numericX Tf = 4.0;
	int N = 6;
	int method = 3;

	mr::MatrixXx result(N, dof);
	result << 1, 0, 0, 1, 1, 0.2, 0, 1,
		1.0208, 0.052, 0.0624, 1.0104, 1.104, 0.3872, 0.0936, 1,
		1.0704, 0.176, 0.2112, 1.0352, 1.352, 0.8336, 0.3168, 1,
		1.1296, 0.324, 0.3888, 1.0648, 1.648, 1.3664, 0.5832, 1,
		1.1792, 0.448, 0.5376, 1.0896, 1.896, 1.8128, 0.8064, 1,
		1.2, 0.5, 0.6, 1.1, 2, 2, 0.9, 1;

	mr::MatrixXx traj = mr::JointTrajectory(thetastart, thetaend, Tf, N, method);
	ASSERT_TRUE(traj.isApprox(result, 4));
}

TEST(MRTest, ScrewTrajectoryTest) {
	mr::MatrixXx Xstart(4, 4);
	Xstart << 1, 0, 0, 1,
		0, 1, 0, 0,
		0, 0, 1, 1,
		0, 0, 0, 1;
	mr::MatrixXx Xend(4, 4);
	Xend << 0, 0, 1, 0.1,
		1, 0, 0, 0,
		0, 1, 0, 4.1,
		0, 0, 0, 1;
	mr::numericX Tf = 5.0;
	int N = 4;
	int method = 3;

	std::vector<mr::MatrixXx> result(N);
	result[0] = Xstart;
	mr::Matrix4x X12;
	X12 << 0.904, -0.25, 0.346, 0.441,
		0.346, 0.904, -0.25, 0.529,
		-0.25, 0.346, 0.904, 1.601,
		0, 0, 0, 1;
	mr::Matrix4x X23;
	X23 << 0.346, -0.25, 0.904, -0.117,
		0.904, 0.346, -0.25, 0.473,
		-0.25, 0.904, 0.346, 3.274,
		0, 0, 0, 1;
	result[1] = X12;
	result[2] = X23;
	result[3] = Xend;

	std::vector<mr::MatrixXx> traj = mr::ScrewTrajectory(Xstart, Xend, Tf, N, method);

	for (int i = 0; i < N; ++i) {
		ASSERT_TRUE(traj[i].isApprox(result[i], 4));
	}
}

TEST(MRTest, CartesianTrajectoryTest) {
	mr::MatrixXx Xstart(4, 4);
	Xstart << 1, 0, 0, 1,
		0, 1, 0, 0,
		0, 0, 1, 1,
		0, 0, 0, 1;
	mr::MatrixXx Xend(4, 4);
	Xend << 0, 0, 1, 0.1,
		1, 0, 0, 0,
		0, 1, 0, 4.1,
		0, 0, 0, 1;
	mr::numericX Tf = 5.0;
	int N = 4;
	int method = 5;

	std::vector<mr::MatrixXx> result(N);
	result[0] = Xstart;
	mr::Matrix4x X12;
	X12 << 0.937, -0.214, 0.277, 0.811,
		0.277, 0.937, -0.214, 0,
		-0.214, 0.277, 0.937, 1.651,
		0, 0, 0, 1;
	mr::Matrix4x X23;
	X23 << 0.277, -0.214, 0.937, 0.289,
		0.937, 0.277, -0.214, 0,
		-0.214, 0.937, 0.277, 3.449,
		0, 0, 0, 1;
	result[1] = X12;
	result[2] = X23;
	result[3] = Xend;

	std::vector<mr::MatrixXx> traj = mr::CartesianTrajectory(Xstart, Xend, Tf, N, method);

	for (int i = 0; i < N; ++i) {
		ASSERT_TRUE(traj[i].isApprox(result[i], 4));
	}
}

TEST(MRTest, InverseDynamicsTrajectoryTest) {
	int dof = 3;
	mr::VectorXx thetastart(dof);
	thetastart << 0, 0, 0;
	mr::VectorXx thetaend(dof);
	thetaend << M_PI / 2, M_PI / 2, M_PI / 2;
	mr::numericX Tf = 3.0;
	int N = 1000;
	int method = 5;

	mr::MatrixXx traj = mr::JointTrajectory(thetastart, thetaend, Tf, N, method);
	mr::MatrixXx thetamat = traj;
	mr::MatrixXx dthetamat = mr::MatrixXx::Zero(N, dof);
	mr::MatrixXx ddthetamat = mr::MatrixXx::Zero(N, dof);
	mr::numericX dt = Tf / (N - 1.0);
	for (int i = 0; i < N - 1; ++i) {
		dthetamat.row(i + 1) = (thetamat.row(i + 1) - thetamat.row(i)) / dt;
		ddthetamat.row(i + 1) = (dthetamat.row(i + 1) - dthetamat.row(i)) / dt;
	}
	mr::VectorXx g(3);
	g << 0, 0, -9.8;
	mr::MatrixXx Ftipmat = mr::MatrixXx::Zero(N, 6);

	std::vector<mr::MatrixXx> Mlist;
	std::vector<mr::MatrixXx> Glist;
	mr::Matrix4x M01;
	M01 << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0.089159,
		0, 0, 0, 1;
	mr::Matrix4x M12;
	M12 << 0, 0, 1, 0.28,
		0, 1, 0, 0.13585,
		-1, 0, 0, 0,
		0, 0, 0, 1;
	mr::Matrix4x M23;
	M23 << 1, 0, 0, 0,
		0, 1, 0, -0.1197,
		0, 0, 1, 0.395,
		0, 0, 0, 1;
	mr::Matrix4x M34;
	M34 << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0.14225,
		0, 0, 0, 1;
	Mlist.push_back(M01);
	Mlist.push_back(M12);
	Mlist.push_back(M23);
	Mlist.push_back(M34);

	mr::VectorXx G1(6);
	G1 << 0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7;
	mr::VectorXx G2(6);
	G2 << 0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393;
	mr::VectorXx G3(6);
	G3 << 0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275;
	Glist.push_back(G1.asDiagonal());
	Glist.push_back(G2.asDiagonal());
	Glist.push_back(G3.asDiagonal());

	mr::MatrixXx SlistT(3, 6);
	SlistT << 1, 0, 1, 0, 1, 0,
		0, 1, 0, -0.089, 0, 0,
		0, 1, 0, -0.089, 0, 0.425;
	mr::MatrixXx Slist = SlistT.transpose();

	int numTest = 3;
	mr::MatrixXx result(numTest, dof);
	mr::VectorXx tau_timestep_beg(3);
	tau_timestep_beg << 13.22970794, -36.262108, -4.181341;
	mr::VectorXx tau_timestep_mid(3);
	tau_timestep_mid << 115.55863434, -22.05129215, 1.00916115;
	mr::VectorXx tau_timestep_end(3);
	tau_timestep_end << 81.12700926, -23.20753925, 2.48432708;
	result << tau_timestep_beg.transpose(),
		tau_timestep_mid.transpose(),
		tau_timestep_end.transpose();

	mr::MatrixXx taumat = mr::InverseDynamicsTrajectory(thetamat, dthetamat, ddthetamat, g, Ftipmat, Mlist, Glist, Slist);
	mr::MatrixXx taumat_timestep(numTest, dof);
	taumat_timestep << taumat.row(0),
		taumat.row(int(N / 2) - 1),
		taumat.row(N - 1);
	ASSERT_TRUE(taumat_timestep.isApprox(result, 4));
}

TEST(MRTest, ForwardDynamicsTrajectoryTest) {
	mr::VectorXx thetalist(3);
	thetalist << 0.1, 0.1, 0.1;
	mr::VectorXx dthetalist(3);
	dthetalist << 0.1, 0.2, 0.3;
	int N = 10, dof = 3;
	mr::MatrixXx taumat(N, 3);
	taumat << 3.63, -6.58, -5.57,
		3.74, -5.55, -5.5,
		4.31, -0.68, -5.19,
		5.18, 5.63, -4.31,
		5.85, 8.17, -2.59,
		5.78, 2.79, -1.7,
		4.99, -5.3, -1.19,
		4.08, -9.41, 0.07,
		3.56, -10.1, 0.97,
		3.49, -9.41, 1.23;
	mr::VectorXx g(3);
	g << 0, 0, -9.8;
	mr::MatrixXx Ftipmat = mr::MatrixXx::Zero(N, 6);

	std::vector<mr::MatrixXx> Mlist;
	std::vector<mr::MatrixXx> Glist;
	mr::Matrix4x M01;
	M01 << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0.089159,
		0, 0, 0, 1;
	mr::Matrix4x M12;
	M12 << 0, 0, 1, 0.28,
		0, 1, 0, 0.13585,
		-1, 0, 0, 0,
		0, 0, 0, 1;
	mr::Matrix4x M23;
	M23 << 1, 0, 0, 0,
		0, 1, 0, -0.1197,
		0, 0, 1, 0.395,
		0, 0, 0, 1;
	mr::Matrix4x M34;
	M34 << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0.14225,
		0, 0, 0, 1;
	Mlist.push_back(M01);
	Mlist.push_back(M12);
	Mlist.push_back(M23);
	Mlist.push_back(M34);

	mr::VectorXx G1(6);
	G1 << 0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7;
	mr::VectorXx G2(6);
	G2 << 0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393;
	mr::VectorXx G3(6);
	G3 << 0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275;
	Glist.push_back(G1.asDiagonal());
	Glist.push_back(G2.asDiagonal());
	Glist.push_back(G3.asDiagonal());

	mr::MatrixXx SlistT(3, 6);
	SlistT << 1, 0, 1, 0, 1, 0,
		0, 1, 0, -0.089, 0, 0,
		0, 1, 0, -0.089, 0, 0.425;
	mr::MatrixXx Slist = SlistT.transpose();
	mr::numericX dt = 0.1;
	int intRes = 8;

	mr::MatrixXx result_thetamat(N, dof);
	mr::MatrixXx result_dthetamat(N, dof);
	result_thetamat << 0.1, 0.1, 0.1,
		0.10643138, 0.2625997, -0.22664947,
		0.10197954, 0.71581297, -1.22521632,
		0.0801044, 1.33930884, -2.28074132,
		0.0282165, 2.11957376, -3.07544297,
		-0.07123855, 2.87726666, -3.83289684,
		-0.20136466, 3.397858, -4.83821609,
		-0.32380092, 3.73338535, -5.98695747,
		-0.41523262, 3.85883317, -7.01130559,
		-0.4638099, 3.63178793, -7.63190052;
	result_dthetamat << 0.1, 0.2, 0.3,
		0.01212502, 3.42975773, -7.74792602,
		-0.13052771, 5.55997471, -11.22722784,
		-0.35521041, 7.11775879, -9.18173035,
		-0.77358795, 8.17307573, -7.05744594,
		-1.2350231, 6.35907497, -8.99784746,
		-1.31426299, 4.07685875, -11.18480509,
		-1.06794821, 2.49227786, -11.69748583,
		-0.70264871, -0.55925705, -8.16067131,
		-0.1455669, -4.57149985, -3.43135114;

	std::vector<mr::MatrixXx> traj = mr::ForwardDynamicsTrajectory(thetalist, dthetalist, taumat, g, Ftipmat, Mlist, Glist, Slist, dt, intRes);
	mr::MatrixXx traj_theta = traj.at(0);
	mr::MatrixXx traj_dtheta = traj.at(1);

	ASSERT_TRUE(traj_theta.isApprox(result_thetamat, 4));
	ASSERT_TRUE(traj_dtheta.isApprox(result_dthetamat, 4));
}

TEST(MRTest, SimulateControlTest) {
	mr::VectorXx thetalist(3);
	thetalist << 0.1, 0.1, 0.1;
	mr::VectorXx dthetalist(3);
	dthetalist << 0.1, 0.2, 0.3;
	mr::VectorXx g(3);
	g << 0, 0, -9.8;

	std::vector<mr::MatrixXx> Mlist;
	std::vector<mr::MatrixXx> Glist;
	mr::Matrix4x M01;
	M01 << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0.089159,
		0, 0, 0, 1;
	mr::Matrix4x M12;
	M12 << 0, 0, 1, 0.28,
		0, 1, 0, 0.13585,
		-1, 0, 0, 0,
		0, 0, 0, 1;
	mr::Matrix4x M23;
	M23 << 1, 0, 0, 0,
		0, 1, 0, -0.1197,
		0, 0, 1, 0.395,
		0, 0, 0, 1;
	mr::Matrix4x M34;
	M34 << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0.14225,
		0, 0, 0, 1;
	Mlist.push_back(M01);
	Mlist.push_back(M12);
	Mlist.push_back(M23);
	Mlist.push_back(M34);

	mr::VectorXx G1(6);
	G1 << 0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7;
	mr::VectorXx G2(6);
	G2 << 0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393;
	mr::VectorXx G3(6);
	G3 << 0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275;
	Glist.push_back(G1.asDiagonal());
	Glist.push_back(G2.asDiagonal());
	Glist.push_back(G3.asDiagonal());

	mr::MatrixXx SlistT(3, 6);
	SlistT << 1, 0, 1, 0, 1, 0,
		0, 1, 0, -0.089, 0, 0,
		0, 1, 0, -0.089, 0, 0.425;
	mr::MatrixXx Slist = SlistT.transpose();
	mr::numericX dt = 0.01;
	mr::VectorXx thetaend(3);
	thetaend << M_PI / 2, M_PI / 2, M_PI / 2;
	mr::numericX Tf = 1.0;
	int N = int(1.0*Tf / dt);
	int method = 5;

	mr::MatrixXx traj = mr::JointTrajectory(thetalist, thetaend, Tf, N, method);
	mr::MatrixXx thetamatd = traj;
	mr::MatrixXx dthetamatd = mr::MatrixXx::Zero(N, 3);
	mr::MatrixXx ddthetamatd = mr::MatrixXx::Zero(N, 3);
	dt = Tf / (N - 1.0);
	for (int i = 0; i < N - 1; ++i) {
		dthetamatd.row(i + 1) = (thetamatd.row(i + 1) - thetamatd.row(i)) / dt;
		ddthetamatd.row(i + 1) = (dthetamatd.row(i + 1) - dthetamatd.row(i)) / dt;
	}

	mr::VectorXx gtilde(3);
	gtilde << 0.8, 0.2, -8.8;

	std::vector<mr::MatrixXx> Mtildelist;
	std::vector<mr::MatrixXx> Gtildelist;
	mr::Matrix4x Mhat01;
	Mhat01 << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0.1,
		0, 0, 0, 1;
	mr::Matrix4x Mhat12;
	Mhat12 << 0, 0, 1, 0.3,
		0, 1, 0, 0.2,
		-1, 0, 0, 0,
		0, 0, 0, 1;
	mr::Matrix4x Mhat23;
	Mhat23 << 1, 0, 0, 0,
		0, 1, 0, -0.2,
		0, 0, 1, 0.4,
		0, 0, 0, 1;
	mr::Matrix4x Mhat34;
	Mhat34 << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0.2,
		0, 0, 0, 1;
	Mtildelist.push_back(Mhat01);
	Mtildelist.push_back(Mhat12);
	Mtildelist.push_back(Mhat23);
	Mtildelist.push_back(Mhat34);

	mr::VectorXx Ghat1(6);
	Ghat1 << 0.1, 0.1, 0.1, 4, 4, 4;
	mr::VectorXx Ghat2(6);
	Ghat2 << 0.3, 0.3, 0.1, 9, 9, 9;
	mr::VectorXx Ghat3(6);
	Ghat3 << 0.1, 0.1, 0.1, 3, 3, 3;
	Gtildelist.push_back(Ghat1.asDiagonal());
	Gtildelist.push_back(Ghat2.asDiagonal());
	Gtildelist.push_back(Ghat3.asDiagonal());
	mr::MatrixXx Ftipmat = mr::MatrixXx::Ones(N, 6);
	mr::numericX Kp = 20.0;
	mr::numericX Ki = 10.0;
	mr::numericX Kd = 18.0;
	int intRes = 8;

	int numTest = 3;  // test 0, N/2-1, N-1 indices of results
	mr::MatrixXx result_taumat(numTest, 3);
	mr::MatrixXx result_thetamat(numTest, 3);

	mr::VectorXx tau_timestep_beg(3);
	tau_timestep_beg << -14.2640765, -54.06797429, -11.265448;
	mr::VectorXx tau_timestep_mid(3);
	tau_timestep_mid << 31.98269367, 9.89625811, 1.47810165;
	mr::VectorXx tau_timestep_end(3);
	tau_timestep_end << 57.04391384, 4.75360586, -1.66561523;
	result_taumat << tau_timestep_beg.transpose(),
		tau_timestep_mid.transpose(),
		tau_timestep_end.transpose();

	mr::VectorXx theta_timestep_beg(3);
	theta_timestep_beg << 0.10092029, 0.10190511, 0.10160667;
	mr::VectorXx theta_timestep_mid(3);
	theta_timestep_mid << 0.85794085, 1.55124503, 2.80130978;
	mr::VectorXx theta_timestep_end(3);
	theta_timestep_end << 1.56344023, 3.07994906, 4.52269971;
	result_thetamat << theta_timestep_beg.transpose(),
		theta_timestep_mid.transpose(),
		theta_timestep_end.transpose();

	std::vector<mr::MatrixXx> controlTraj = mr::SimulateControl(thetalist, dthetalist, g, Ftipmat, Mlist, Glist, Slist, thetamatd, dthetamatd,
		ddthetamatd, gtilde, Mtildelist, Gtildelist, Kp, Ki, Kd, dt, intRes);
	mr::MatrixXx traj_tau = controlTraj.at(0);
	mr::MatrixXx traj_theta = controlTraj.at(1);
	mr::MatrixXx traj_tau_timestep(numTest, 3);
	traj_tau_timestep << traj_tau.row(0),
		traj_tau.row(int(N / 2) - 1),
		traj_tau.row(N - 1);
	mr::MatrixXx traj_theta_timestep(numTest, 3);
	traj_theta_timestep << traj_theta.row(0),
		traj_theta.row(int(N / 2) - 1),
		traj_theta.row(N - 1);

	ASSERT_TRUE(traj_tau_timestep.isApprox(result_taumat, 4));
	ASSERT_TRUE(traj_theta_timestep.isApprox(result_thetamat, 4));
}
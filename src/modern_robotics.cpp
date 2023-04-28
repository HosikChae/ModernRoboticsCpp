#include "../include/modern_robotics.h"

/*
 * modernRobotics.cpp
 * Adapted from modern_robotics.py provided by modernrobotics.org
 * Provides useful Jacobian and frame representation functions
 */
#include <Eigen/Dense>
#include <cmath>
#include <vector>

# define M_PI           3.14159265358979323846  /* pi */

namespace mr {

	/* Function: Find if the value is negligible enough to consider 0
	 * Inputs: value to be checked as a double
	 * Returns: Boolean of true-ignore or false-can't ignore
	 */
	bool NearZero(const numericX val) {
		return (std::abs(val) < .000001);
	}

	/*
	 * Function: Calculate the 6x6 matrix [adV] of the given 6-vector
	 * Input: VectorXx (6x1)
	 * Output: MatrixXx (6x6)
	 * Note: Can be used to calculate the Lie bracket [V1, V2] = [adV1]V2
	 */
	MatrixXx ad(VectorXx V) {
		Matrix3x omgmat = VecToso3(Vector3x(V(0), V(1), V(2)));

		MatrixXx result(6, 6);
		result.topLeftCorner<3, 3>() = omgmat;
		result.topRightCorner<3, 3>() = Matrix3x::Zero(3, 3);
		result.bottomLeftCorner<3, 3>() = VecToso3(Vector3x(V(3), V(4), V(5)));
		result.bottomRightCorner<3, 3>() = omgmat;
		return result;
	}

	/* Function: Returns a normalized version of the input vector
	 * Input: MatrixXx
	 * Output: MatrixXx
	 * Note: MatrixXx is used instead of VectorXx for the case of row vectors
	 * 		Requires a copy
	 *		Useful because of the MatrixXx casting
	 */
	MatrixXx Normalize(MatrixXx V) {
		V.normalize();
		return V;
	}


	/* Function: Returns the skew symmetric matrix representation of an angular velocity vector
	 * Input: Vector3x 3x1 angular velocity vector
	 * Returns: MatrixXx 3x3 skew symmetric matrix
	 */
	Matrix3x VecToso3(const Vector3x& omg) {
		Matrix3x m_ret;
		m_ret << 0, -omg(2), omg(1),
			omg(2), 0, -omg(0),
			-omg(1), omg(0), 0;
		return m_ret;
	}


	/* Function: Returns angular velocity vector represented by the skew symmetric matrix
	 * Inputs: MatrixXx 3x3 skew symmetric matrix
	 * Returns: Vector3x 3x1 angular velocity
	 */
	Vector3x so3ToVec(const MatrixXx& so3mat) {
		Vector3x v_ret;
		v_ret << so3mat(2, 1), so3mat(0, 2), so3mat(1, 0);
		return v_ret;
	}


	/* Function: Translates an exponential rotation into it's individual components
	 * Inputs: Exponential rotation (rotation matrix in terms of a rotation axis
	 *				and the angle of rotation)
	 * Returns: The axis and angle of rotation as [x, y, z, theta]
	 */
	Vector4x AxisAng3(const Vector3x& expc3) {
		Vector4x v_ret;
		v_ret << Normalize(expc3), expc3.norm();
		return v_ret;
	}


	/* Function: Translates an exponential rotation into a rotation matrix
	 * Inputs: exponenential representation of a rotation
	 * Returns: Rotation matrix
	 */
	Matrix3x MatrixExp3(const Matrix3x& so3mat) {
		Vector3x omgtheta = so3ToVec(so3mat);

		Matrix3x m_ret = Matrix3x::Identity();
		if (NearZero(so3mat.norm())) {
			return m_ret;
		}
		else {
			numericX theta = (AxisAng3(omgtheta))(3);
			Matrix3x omgmat = so3mat * (1 / theta);
			return m_ret + std::sin(theta) * omgmat + ((1 - std::cos(theta)) * (omgmat * omgmat));
		}
	}


	/* Function: Computes the matrix logarithm of a rotation matrix
	 * Inputs: Rotation matrix
	 * Returns: matrix logarithm of a rotation
	 */
	Matrix3x MatrixLog3(const Matrix3x& R) {
		numericX acosinput = (R.trace() - 1) / 2.0;
		MatrixXx m_ret = MatrixXx::Zero(3, 3);
		if (acosinput >= 1)
			return m_ret;
		else if (acosinput <= -1) {
			Vector3x omg;
			if (!NearZero(1 + R(2, 2)))
				omg = (1.0 / std::sqrt(2 * (1 + R(2, 2))))*Vector3x(R(0, 2), R(1, 2), 1 + R(2, 2));
			else if (!NearZero(1 + R(1, 1)))
				omg = (1.0 / std::sqrt(2 * (1 + R(1, 1))))*Vector3x(R(0, 1), 1 + R(1, 1), R(2, 1));
			else
				omg = (1.0 / std::sqrt(2 * (1 + R(0, 0))))*Vector3x(1 + R(0, 0), R(1, 0), R(2, 0));
			m_ret = VecToso3(M_PI * omg);
			return m_ret;
		}
		else {
			numericX theta = std::acos(acosinput);
			m_ret = theta / 2.0 / sin(theta)*(R - R.transpose());
			return m_ret;
		}
	}

	/* Function: Combines a rotation matrix and position vector into a single
	 * 				Special Euclidian Group (SE3) homogeneous transformation matrix
	 * Inputs: Rotation Matrix (R), Position Vector (p)
	 * Returns: Matrix of T = [ [R, p],
	 *						    [0, 1] ]
	 */
	MatrixXx RpToTrans(const Matrix3x& R, const Vector3x& p) {
		MatrixXx m_ret(4, 4);
		m_ret << R, p,
			0, 0, 0, 1;
		return m_ret;
	}


	/* Function: Separates the rotation matrix and position vector from
	 *				the transfomation matrix representation
	 * Inputs: Homogeneous transformation matrix
	 * Returns: std::vector of [rotation matrix, position vector]
	 */
	std::vector<MatrixXx> TransToRp(const MatrixXx& T) {
		std::vector<MatrixXx> Rp_ret;
		Matrix3x R_ret;
		// Get top left 3x3 corner
		R_ret = T.block<3, 3>(0, 0);

		Vector3x p_ret(T(0, 3), T(1, 3), T(2, 3));

		Rp_ret.push_back(R_ret);
		Rp_ret.push_back(p_ret);

		return Rp_ret;
	}


	/* Function: Translates a spatial velocity vector into a transformation matrix
	 * Inputs: Spatial velocity vector [angular velocity, linear velocity]
	 * Returns: Transformation matrix
	 */
	MatrixXx VecTose3(const VectorXx& V) {
		// Separate angular (exponential representation) and linear velocities
		Vector3x exp(V(0), V(1), V(2));
		Vector3x linear(V(3), V(4), V(5));

		// Fill in values to the appropriate parts of the transformation matrix
		MatrixXx m_ret(4, 4);
		m_ret << VecToso3(exp), linear,
			0, 0, 0, 0;

		return m_ret;
	}


	/* Function: Translates a transformation matrix into a spatial velocity vector
	 * Inputs: Transformation matrix
	 * Returns: Spatial velocity vector [angular velocity, linear velocity]
	 */
	VectorXx se3ToVec(const MatrixXx& T) {
		VectorXx m_ret(6);
		m_ret << T(2, 1), T(0, 2), T(1, 0), T(0, 3), T(1, 3), T(2, 3);

		return m_ret;
	}


	/* Function: Provides the adjoint representation of a transformation matrix
	 *			 Used to change the frame of reference for spatial velocity vectors
	 * Inputs: 4x4 Transformation matrix SE(3)
	 * Returns: 6x6 Adjoint Representation of the matrix
	 */
	MatrixXx Adjoint(const MatrixXx& T) {
		std::vector<MatrixXx> R = TransToRp(T);
		MatrixXx ad_ret(6, 6);
		ad_ret = MatrixXx::Zero(6, 6);
		MatrixXx zeroes = MatrixXx::Zero(3, 3);
		ad_ret << R[0], zeroes,
			VecToso3(R[1]) * R[0], R[0];
		return ad_ret;
	}


	/* Function: Rotation expanded for screw axis
	 * Inputs: se3 matrix representation of exponential coordinates (transformation matrix)
	 * Returns: 6x6 Matrix representing the rotation
	 */
	MatrixXx MatrixExp6(const MatrixXx& se3mat) {
		// Extract the angular velocity vector from the transformation matrix
		Matrix3x se3mat_cut = se3mat.block<3, 3>(0, 0);
		Vector3x omgtheta = so3ToVec(se3mat_cut);

		MatrixXx m_ret(4, 4);

		// If negligible rotation, m_Ret = [[Identity, angular velocty ]]
		//									[	0	 ,		1		   ]]
		if (NearZero(omgtheta.norm())) {
			// Reuse previous variables that have our required size
			se3mat_cut = MatrixXx::Identity(3, 3);
			omgtheta << se3mat(0, 3), se3mat(1, 3), se3mat(2, 3);
			m_ret << se3mat_cut, omgtheta,
				0, 0, 0, 1;
			return m_ret;
		}
		// If not negligible, MR page 105
		else {
			numericX theta = (AxisAng3(omgtheta))(3);
			Matrix3x omgmat = se3mat.block<3, 3>(0, 0) / theta;
			Matrix3x expExpand = MatrixXx::Identity(3, 3) * theta + (1 - std::cos(theta)) * omgmat + ((theta - std::sin(theta)) * (omgmat * omgmat));
			Vector3x linear(se3mat(0, 3), se3mat(1, 3), se3mat(2, 3));
			Vector3x GThetaV = (expExpand*linear) / theta;
			m_ret << MatrixExp3(se3mat_cut), GThetaV,
				0, 0, 0, 1;
			return m_ret;
		}

	}

	MatrixXx MatrixLog6(const MatrixXx& T) {
		MatrixXx m_ret(4, 4);
		auto rp = mr::TransToRp(T);
		Matrix3x omgmat = MatrixLog3(rp.at(0));
		Matrix3x zeros3x = Matrix3x::Zero(3, 3);
		if (NearZero(omgmat.norm())) {
			m_ret << zeros3x, rp.at(1),
				0, 0, 0, 0;
		}
		else {
			numericX theta = std::acos((rp.at(0).trace() - 1) / 2.0);
			Matrix3x logExpand1 = MatrixXx::Identity(3, 3) - omgmat / 2.0;
			Matrix3x logExpand2 = (1.0 / theta - 1.0 / std::tan(theta / 2.0) / 2)*omgmat*omgmat / theta;
			Matrix3x logExpand = logExpand1 + logExpand2;
			m_ret << omgmat, logExpand*rp.at(1),
				0, 0, 0, 0;
		}
		return m_ret;
	}


	/* Function: Compute end effector frame (used for current spatial position calculation)
	 * Inputs: Home configuration (position and orientation) of end-effector
	 *		   The joint screw axes in the space frame when the manipulator
	 *             is at the home position
	 * 		   A list of joint coordinates.
	 * Returns: Transfomation matrix representing the end-effector frame when the joints are
	 *				at the specified coordinates
	 * Notes: FK means Forward Kinematics
	 */
	MatrixXx FKinSpace(const MatrixXx& M, const MatrixXx& Slist, const VectorXx& thetaList) {
		MatrixXx T = M;
		for (int i = (thetaList.size() - 1); i > -1; i--) {
			T = MatrixExp6(VecTose3(Slist.col(i)*thetaList(i))) * T;
		}
		return T;
	}

	/*
	 * Function: Compute end effector frame (used for current body position calculation)
	 * Inputs: Home configuration (position and orientation) of end-effector
	 *		   The joint screw axes in the body frame when the manipulator
	 *             is at the home position
	 * 		   A list of joint coordinates.
	 * Returns: Transfomation matrix representing the end-effector frame when the joints are
	 *				at the specified coordinates
	 * Notes: FK means Forward Kinematics
	 */
	MatrixXx FKinBody(const MatrixXx& M, const MatrixXx& Blist, const VectorXx& thetaList) {
		MatrixXx T = M;
		for (int i = 0; i < thetaList.size(); i++) {
			T = T * MatrixExp6(VecTose3(Blist.col(i)*thetaList(i)));
		}
		return T;
	}


	/* Function: Gives the space Jacobian
	 * Inputs: Screw axis in home position, joint configuration
	 * Returns: 6xn Spatial Jacobian
	 */
	MatrixXx JacobianSpace(const MatrixXx& Slist, const MatrixXx& thetaList) {
		MatrixXx Js = Slist;
		MatrixXx T = MatrixXx::Identity(4, 4);
		VectorXx sListTemp(Slist.col(0).size());
		for (int i = 1; i < thetaList.size(); i++) {
			sListTemp << Slist.col(i - 1) * thetaList(i - 1);
			T = T * MatrixExp6(VecTose3(sListTemp));
			// std::cout << "array: " << sListTemp << std::endl;
			Js.col(i) = Adjoint(T) * Slist.col(i);
		}

		return Js;
	}

	/*
	 * Function: Gives the body Jacobian
	 * Inputs: Screw axis in BODY position, joint configuration
	 * Returns: 6xn Bobdy Jacobian
	 */
	MatrixXx JacobianBody(const MatrixXx& Blist, const MatrixXx& thetaList) {
		MatrixXx Jb = Blist;
		MatrixXx T = MatrixXx::Identity(4, 4);
		VectorXx bListTemp(Blist.col(0).size());
		for (int i = thetaList.size() - 2; i >= 0; i--) {
			bListTemp << Blist.col(i + 1) * thetaList(i + 1);
			T = T * MatrixExp6(VecTose3(-1 * bListTemp));
			// std::cout << "array: " << sListTemp << std::endl;
			Jb.col(i) = Adjoint(T) * Blist.col(i);
		}
		return Jb;
	}

	MatrixXx TransInv(const MatrixXx& transform) {
		auto rp = mr::TransToRp(transform);
		auto Rt = rp.at(0).transpose();
		auto t = -(Rt * rp.at(1));
		MatrixXx inv(4, 4);
		inv = MatrixXx::Zero(4,4);
		inv.block(0, 0, 3, 3) = Rt;
		inv.block(0, 3, 3, 1) = t;
		inv(3, 3) = 1;
		return inv;
	}

	MatrixXx RotInv(const MatrixXx& rotMatrix) {
		return rotMatrix.transpose();
	}

	VectorXx ScrewToAxis(Vector3x q, Vector3x s, numericX h) {
		VectorXx axis(6);
		axis.segment(0, 3) = s;
		axis.segment(3, 3) = q.cross(s) + (h * s);
		return axis;
	}

	VectorXx AxisAng6(const VectorXx& expc6) {
		VectorXx v_ret(7);
		numericX theta = Vector3x(expc6(0), expc6(1), expc6(2)).norm();
		if (NearZero(theta))
			theta = Vector3x(expc6(3), expc6(4), expc6(5)).norm();
		v_ret << expc6 / theta, theta;
		return v_ret;
	}

	MatrixXx ProjectToSO3(const MatrixXx& M) {
		Eigen::JacobiSVD<MatrixXx> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
		MatrixXx R = svd.matrixU() * svd.matrixV().transpose();
		if (R.determinant() < 0)
			// In this case the result may be far from M; reverse sign of 3rd column
			R.col(2) *= -1;
		return R;
	}

	MatrixXx ProjectToSE3(const MatrixXx& M) {
		Matrix3x R = M.block<3, 3>(0, 0);
		Vector3x t = M.block<3, 1>(0, 3);
		MatrixXx T = RpToTrans(ProjectToSO3(R), t);
		return T;
	}

	numericX DistanceToSO3(const Matrix3x& M) {
		if (M.determinant() > 0)
			return (M.transpose() * M - Matrix3x::Identity()).norm();
		else
			return 1.0e9;
	}

	numericX DistanceToSE3(const Matrix4x& T) {
		Matrix3x matR = T.block<3, 3>(0, 0);
		if (matR.determinant() > 0) {
			Matrix4x m_ret;
			m_ret << matR.transpose()*matR, Vector3x::Zero(3),
				T.row(3);
			m_ret = m_ret - Matrix4x::Identity();
			return m_ret.norm();
		}
		else
			return 1.0e9;
	}

	bool TestIfSO3(const Matrix3x& M) {
		return std::abs(DistanceToSO3(M)) < 1e-3;
	}

	bool TestIfSE3(const Matrix4x& T) {
		return std::abs(DistanceToSE3(T)) < 1e-3;
	}
	bool IKinBody(const MatrixXx& Blist, const MatrixXx& M, const MatrixXx& T,
		VectorXx& thetalist, numericX eomg, numericX ev) {
		int i = 0;
		int maxiterations = 20;
		MatrixXx Tfk = FKinBody(M, Blist, thetalist);
		MatrixXx Tdiff = TransInv(Tfk)*T;
		VectorXx Vb = se3ToVec(MatrixLog6(Tdiff));
		Vector3x angular(Vb(0), Vb(1), Vb(2));
		Vector3x linear(Vb(3), Vb(4), Vb(5));

		bool err = (angular.norm() > eomg || linear.norm() > ev);
		MatrixXx Jb;
		while (err && i < maxiterations) {
			Jb = JacobianBody(Blist, thetalist);
			thetalist += Jb.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(Vb);
			i += 1;
			// iterate
			Tfk = FKinBody(M, Blist, thetalist);
			Tdiff = TransInv(Tfk)*T;
			Vb = se3ToVec(MatrixLog6(Tdiff));
			angular = Vector3x(Vb(0), Vb(1), Vb(2));
			linear = Vector3x(Vb(3), Vb(4), Vb(5));
			err = (angular.norm() > eomg || linear.norm() > ev);
		}
		return !err;
	}

	bool IKinSpace(const MatrixXx& Slist, const MatrixXx& M, const MatrixXx& T,
		VectorXx& thetalist, numericX eomg, numericX ev) {
		int i = 0;
		int maxiterations = 20;
		MatrixXx Tfk = FKinSpace(M, Slist, thetalist);
		MatrixXx Tdiff = TransInv(Tfk)*T;
		VectorXx Vs = Adjoint(Tfk)*se3ToVec(MatrixLog6(Tdiff));
		Vector3x angular(Vs(0), Vs(1), Vs(2));
		Vector3x linear(Vs(3), Vs(4), Vs(5));

		bool err = (angular.norm() > eomg || linear.norm() > ev);
		MatrixXx Js;
		while (err && i < maxiterations) {
			Js = JacobianSpace(Slist, thetalist);
			thetalist += Js.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(Vs);
			i += 1;
			// iterate
			Tfk = FKinSpace(M, Slist, thetalist);
			Tdiff = TransInv(Tfk)*T;
			Vs = Adjoint(Tfk)*se3ToVec(MatrixLog6(Tdiff));
			angular = Vector3x(Vs(0), Vs(1), Vs(2));
			linear = Vector3x(Vs(3), Vs(4), Vs(5));
			err = (angular.norm() > eomg || linear.norm() > ev);
		}
		return !err;
	}

	/*
	* Function: This function uses forward-backward Newton-Euler iterations to solve the
	* equation:
	* taulist = Mlist(thetalist) * ddthetalist + c(thetalist, dthetalist) ...
	*           + g(thetalist) + Jtr(thetalist) * Ftip
	* Inputs:
	*  thetalist: n-vector of joint variables
	*  dthetalist: n-vector of joint rates
	*  ddthetalist: n-vector of joint accelerations
	*  g: Gravity vector g
	*  Ftip: Spatial force applied by the end-effector expressed in frame {n+1}
	*  Mlist: List of link frames {i} relative to {i-1} at the home position
	*  Glist: Spatial inertia matrices Gi of the links
	*  Slist: Screw axes Si of the joints in a space frame, in the format
	*         of a matrix with the screw axes as the columns.
	*
	* Outputs:
	*  taulist: The n-vector of required joint forces/torques
	*
	*/
	VectorXx InverseDynamics(const VectorXx& thetalist, const VectorXx& dthetalist, const VectorXx& ddthetalist,
									const VectorXx& g, const VectorXx& Ftip, const std::vector<MatrixXx>& Mlist,
									const std::vector<MatrixXx>& Glist, const MatrixXx& Slist) {
	    // the size of the lists
		int n = thetalist.size();

		MatrixXx Mi = MatrixXx::Identity(4, 4);
		MatrixXx Ai = MatrixXx::Zero(6,n);
		std::vector<MatrixXx> AdTi;
		for (int i = 0; i < n+1; i++) {
			AdTi.push_back(MatrixXx::Zero(6,6));
		}
		MatrixXx Vi = MatrixXx::Zero(6,n+1);    // velocity
		MatrixXx Vdi = MatrixXx::Zero(6,n+1);   // acceleration

		Vdi.block(3, 0, 3, 1) = - g;
		AdTi[n] = mr::Adjoint(mr::TransInv(Mlist[n]));
		VectorXx Fi = Ftip;

		VectorXx taulist = VectorXx::Zero(n);

		// forward pass
		for (int i = 0; i < n; i++) {
			Mi = Mi * Mlist[i];
			Ai.col(i) = mr::Adjoint(mr::TransInv(Mi))*Slist.col(i);

			AdTi[i] = mr::Adjoint(mr::MatrixExp6(mr::VecTose3(Ai.col(i)*-thetalist(i)))
			          * mr::TransInv(Mlist[i]));

			Vi.col(i+1) = AdTi[i] * Vi.col(i) + Ai.col(i) * dthetalist(i);
			Vdi.col(i+1) = AdTi[i] * Vdi.col(i) + Ai.col(i) * ddthetalist(i)
						   + ad(Vi.col(i+1)) * Ai.col(i) * dthetalist(i); // this index is different from book!
		}

		// backward pass
		for (int i = n-1; i >= 0; i--) {
			Fi = AdTi[i+1].transpose() * Fi + Glist[i] * Vdi.col(i+1)
			     - ad(Vi.col(i+1)).transpose() * (Glist[i] * Vi.col(i+1));
			taulist(i) = Fi.transpose() * Ai.col(i);
		}
		return taulist;
	}

	/*
	 * Function: This function calls InverseDynamics with Ftip = 0, dthetalist = 0, and
	 *   ddthetalist = 0. The purpose is to calculate one important term in the dynamics equation
	 * Inputs:
	 *  thetalist: n-vector of joint variables
	 *  g: Gravity vector g
	 *  Mlist: List of link frames {i} relative to {i-1} at the home position
	 *  Glist: Spatial inertia matrices Gi of the links
	 *  Slist: Screw axes Si of the joints in a space frame, in the format
	 *         of a matrix with the screw axes as the columns.
	 *
	 * Outputs:
	 *  grav: The 3-vector showing the effect force of gravity to the dynamics
	 *
	 */
	VectorXx GravityForces(const VectorXx& thetalist, const VectorXx& g,
									const std::vector<MatrixXx>& Mlist, const std::vector<MatrixXx>& Glist, const MatrixXx& Slist) {
	    int n = thetalist.size();
		VectorXx dummylist = VectorXx::Zero(n);
		VectorXx dummyForce = VectorXx::Zero(6);
		VectorXx grav = mr::InverseDynamics(thetalist, dummylist, dummylist, g,
                                                dummyForce, Mlist, Glist, Slist);
		return grav;
	}

	/*
  	 * Function: This function calls InverseDynamics n times, each time passing a
	 * ddthetalist vector with a single element equal to one and all other
	 * inputs set to zero. Each call of InverseDynamics generates a single
	 * column, and these columns are assembled to create the inertia matrix.
	 *
	 * Inputs:
	 *  thetalist: n-vector of joint variables
	 *  Mlist: List of link frames {i} relative to {i-1} at the home position
	 *  Glist: Spatial inertia matrices Gi of the links
	 *  Slist: Screw axes Si of the joints in a space frame, in the format
	 *         of a matrix with the screw axes as the columns.
	 *
	 * Outputs:
	 *  M: The numerical inertia matrix M(thetalist) of an n-joint serial
	 *     chain at the given configuration thetalist.
	 */
	MatrixXx MassMatrix(const VectorXx& thetalist,
                                const std::vector<MatrixXx>& Mlist, const std::vector<MatrixXx>& Glist, const MatrixXx& Slist) {
		int n = thetalist.size();
		VectorXx dummylist = VectorXx::Zero(n);
		VectorXx dummyg = VectorXx::Zero(3);
		VectorXx dummyforce = VectorXx::Zero(6);
		MatrixXx M = MatrixXx::Zero(n,n);
		for (int i = 0; i < n; i++) {
			VectorXx ddthetalist = VectorXx::Zero(n);
			ddthetalist(i) = 1;
			M.col(i) = mr::InverseDynamics(thetalist, dummylist, ddthetalist,
                             dummyg, dummyforce, Mlist, Glist, Slist);
		}
		return M;
	}

	/*
  	 * Function: This function calls InverseDynamics with g = 0, Ftip = 0, and
     * ddthetalist = 0.
	 *
	 * Inputs:
	 *  thetalist: n-vector of joint variables
	 *  dthetalist: A list of joint rates
	 *  Mlist: List of link frames {i} relative to {i-1} at the home position
	 *  Glist: Spatial inertia matrices Gi of the links
	 *  Slist: Screw axes Si of the joints in a space frame, in the format
	 *         of a matrix with the screw axes as the columns.
	 *
	 * Outputs:
	 *  c: The vector c(thetalist,dthetalist) of Coriolis and centripetal
	 *     terms for a given thetalist and dthetalist.
	 */
	VectorXx VelQuadraticForces(const VectorXx& thetalist, const VectorXx& dthetalist,
                                const std::vector<MatrixXx>& Mlist, const std::vector<MatrixXx>& Glist, const MatrixXx& Slist) {
		int n = thetalist.size();
		VectorXx dummylist = VectorXx::Zero(n);
		VectorXx dummyg = VectorXx::Zero(3);
		VectorXx dummyforce = VectorXx::Zero(6);
		VectorXx c = mr::InverseDynamics(thetalist, dthetalist, dummylist,
                             dummyg, dummyforce, Mlist, Glist, Slist);
		return c;
	}

	/*
  	 * Function: This function calls InverseDynamics with g = 0, dthetalist = 0, and
     * ddthetalist = 0.
	 *
	 * Inputs:
	 *  thetalist: n-vector of joint variables
	 *  Ftip: Spatial force applied by the end-effector expressed in frame {n+1}
	 *  Mlist: List of link frames {i} relative to {i-1} at the home position
	 *  Glist: Spatial inertia matrices Gi of the links
	 *  Slist: Screw axes Si of the joints in a space frame, in the format
	 *         of a matrix with the screw axes as the columns.
	 *
	 * Outputs:
	 *  JTFtip: The joint forces and torques required only to create the
	 *     end-effector force Ftip.
	 */
	VectorXx EndEffectorForces(const VectorXx& thetalist, const VectorXx& Ftip,
								const std::vector<MatrixXx>& Mlist, const std::vector<MatrixXx>& Glist, const MatrixXx& Slist) {
		int n = thetalist.size();
		VectorXx dummylist = VectorXx::Zero(n);
		VectorXx dummyg = VectorXx::Zero(3);

		VectorXx JTFtip = mr::InverseDynamics(thetalist, dummylist, dummylist,
                             dummyg, Ftip, Mlist, Glist, Slist);
		return JTFtip;
	}

	/*
	 * Function: This function computes ddthetalist by solving:
	 * Mlist(thetalist) * ddthetalist = taulist - c(thetalist,dthetalist)
	 *                                  - g(thetalist) - Jtr(thetalist) * Ftip
	 * Inputs:
	 *  thetalist: n-vector of joint variables
	 *  dthetalist: n-vector of joint rates
	 *  taulist: An n-vector of joint forces/torques
	 *  g: Gravity vector g
	 *  Ftip: Spatial force applied by the end-effector expressed in frame {n+1}
	 *  Mlist: List of link frames {i} relative to {i-1} at the home position
	 *  Glist: Spatial inertia matrices Gi of the links
	 *  Slist: Screw axes Si of the joints in a space frame, in the format
	 *         of a matrix with the screw axes as the columns.
	 *
	 * Outputs:
	 *  ddthetalist: The resulting joint accelerations
	 *
	 */
	VectorXx ForwardDynamics(const VectorXx& thetalist, const VectorXx& dthetalist, const VectorXx& taulist,
									const VectorXx& g, const VectorXx& Ftip, const std::vector<MatrixXx>& Mlist,
									const std::vector<MatrixXx>& Glist, const MatrixXx& Slist) {

		VectorXx totalForce = taulist - mr::VelQuadraticForces(thetalist, dthetalist, Mlist, Glist, Slist)
                 							 - mr::GravityForces(thetalist, g, Mlist, Glist, Slist)
                                             - mr::EndEffectorForces(thetalist, Ftip, Mlist, Glist, Slist);

		MatrixXx M = mr::MassMatrix(thetalist, Mlist, Glist, Slist);

		// Use LDLT since M is positive definite
        VectorXx ddthetalist = M.ldlt().solve(totalForce);

		return ddthetalist;
	}

	void EulerStep(VectorXx& thetalist, VectorXx& dthetalist, const VectorXx& ddthetalist, numericX dt) {
		thetalist += dthetalist * dt;
		dthetalist += ddthetalist * dt;
		return;
	}

	MatrixXx InverseDynamicsTrajectory(const MatrixXx& thetamat, const MatrixXx& dthetamat, const MatrixXx& ddthetamat,
		const VectorXx& g, const MatrixXx& Ftipmat, const std::vector<MatrixXx>& Mlist, const std::vector<MatrixXx>& Glist,
		const MatrixXx& Slist) {
		MatrixXx thetamatT = thetamat.transpose();
		MatrixXx dthetamatT = dthetamat.transpose();
		MatrixXx ddthetamatT = ddthetamat.transpose();
		MatrixXx FtipmatT = Ftipmat.transpose();

		int N = thetamat.rows();  // trajectory points
		int dof = thetamat.cols();
		MatrixXx taumatT = MatrixXx::Zero(dof, N);
		for (int i = 0; i < N; ++i) {
			taumatT.col(i) = InverseDynamics(thetamatT.col(i), dthetamatT.col(i), ddthetamatT.col(i), g, FtipmatT.col(i), Mlist, Glist, Slist);
		}
		MatrixXx taumat = taumatT.transpose();
		return taumat;
	}

	std::vector<MatrixXx> ForwardDynamicsTrajectory(const VectorXx& thetalist, const VectorXx& dthetalist, const MatrixXx& taumat,
		const VectorXx& g, const MatrixXx& Ftipmat, const std::vector<MatrixXx>& Mlist, const std::vector<MatrixXx>& Glist,
		const MatrixXx& Slist, numericX dt, int intRes) {
		MatrixXx taumatT = taumat.transpose();
		MatrixXx FtipmatT = Ftipmat.transpose();
		int N = taumat.rows();  // force/torque points
		int dof = taumat.cols();
		MatrixXx thetamatT = MatrixXx::Zero(dof, N);
		MatrixXx dthetamatT = MatrixXx::Zero(dof, N);
		thetamatT.col(0) = thetalist;
		dthetamatT.col(0) = dthetalist;
		VectorXx thetacurrent = thetalist;
		VectorXx dthetacurrent = dthetalist;
		VectorXx ddthetalist;
		for (int i = 0; i < N - 1; ++i) {
			for (int j = 0; j < intRes; ++j) {
				ddthetalist = ForwardDynamics(thetacurrent, dthetacurrent, taumatT.col(i), g, FtipmatT.col(i), Mlist, Glist, Slist);
				EulerStep(thetacurrent, dthetacurrent, ddthetalist, 1.0*dt / intRes);
			}
			thetamatT.col(i + 1) = thetacurrent;
			dthetamatT.col(i + 1) = dthetacurrent;
		}
		std::vector<MatrixXx> JointTraj_ret;
		JointTraj_ret.push_back(thetamatT.transpose());
		JointTraj_ret.push_back(dthetamatT.transpose());
		return JointTraj_ret;
	}

	VectorXx ComputedTorque(const VectorXx& thetalist, const VectorXx& dthetalist, const VectorXx& eint,
		const VectorXx& g, const std::vector<MatrixXx>& Mlist, const std::vector<MatrixXx>& Glist,
		const MatrixXx& Slist, const VectorXx& thetalistd, const VectorXx& dthetalistd, const VectorXx& ddthetalistd,
		numericX Kp, numericX Ki, numericX Kd) {

		VectorXx e = thetalistd - thetalist;  // position err
		VectorXx tau_feedforward = MassMatrix(thetalist, Mlist, Glist, Slist)*(Kp*e + Ki * (eint + e) + Kd * (dthetalistd - dthetalist));

		VectorXx Ftip = VectorXx::Zero(6);
		VectorXx tau_inversedyn = InverseDynamics(thetalist, dthetalist, ddthetalistd, g, Ftip, Mlist, Glist, Slist);

		VectorXx tau_computed = tau_feedforward + tau_inversedyn;
		return tau_computed;
	}

	numericX CubicTimeScaling(numericX Tf, numericX t) {
		numericX timeratio = 1.0*t / Tf;
		numericX st = 3 * pow(timeratio, 2) - 2 * pow(timeratio, 3);
		return st;
	}

	numericX QuinticTimeScaling(numericX Tf, numericX t) {
		numericX timeratio = 1.0*t / Tf;
		numericX st = 10 * pow(timeratio, 3) - 15 * pow(timeratio, 4) + 6 * pow(timeratio, 5);
		return st;
	}

	MatrixXx JointTrajectory(const VectorXx& thetastart, const VectorXx& thetaend, numericX Tf, int N, int method) {
		numericX timegap = Tf / (N - 1);
		MatrixXx trajT = MatrixXx::Zero(thetastart.size(), N);
		numericX st;
		for (int i = 0; i < N; ++i) {
			if (method == 3)
				st = CubicTimeScaling(Tf, timegap*i);
			else
				st = QuinticTimeScaling(Tf, timegap*i);
			trajT.col(i) = st * thetaend + (1 - st)*thetastart;
		}
		MatrixXx traj = trajT.transpose();
		return traj;
	}
	std::vector<MatrixXx> ScrewTrajectory(const MatrixXx& Xstart, const MatrixXx& Xend, numericX Tf, int N, int method) {
		numericX timegap = Tf / (N - 1);
		std::vector<MatrixXx> traj(N);
		numericX st;
		for (int i = 0; i < N; ++i) {
			if (method == 3)
				st = CubicTimeScaling(Tf, timegap*i);
			else
				st = QuinticTimeScaling(Tf, timegap*i);
			MatrixXx Ttemp = MatrixLog6(TransInv(Xstart)*Xend);
			traj.at(i) = Xstart * MatrixExp6(Ttemp*st);
		}
		return traj;
	}

	std::vector<MatrixXx> CartesianTrajectory(const MatrixXx& Xstart, const MatrixXx& Xend, numericX Tf, int N, int method) {
		numericX timegap = Tf / (N - 1);
		std::vector<MatrixXx> traj(N);
		std::vector<MatrixXx> Rpstart = TransToRp(Xstart);
		std::vector<MatrixXx> Rpend = TransToRp(Xend);
		Matrix3x Rstart = Rpstart[0]; Vector3x pstart = Rpstart[1];
		Matrix3x Rend = Rpend[0]; Vector3x pend = Rpend[1];
		numericX st;
		for (int i = 0; i < N; ++i) {
			if (method == 3)
				st = CubicTimeScaling(Tf, timegap*i);
			else
				st = QuinticTimeScaling(Tf, timegap*i);
			Matrix3x Ri = Rstart * MatrixExp3(MatrixLog3(Rstart.transpose() * Rend)*st);
			Vector3x pi = st*pend + (1 - st)*pstart;
			MatrixXx traji(4, 4);
			traji << Ri, pi,
				0, 0, 0, 1;
			traj.at(i) = traji;
		}
		return traj;
	}
	std::vector<MatrixXx> SimulateControl(const VectorXx& thetalist, const VectorXx& dthetalist, const VectorXx& g,
		const MatrixXx& Ftipmat, const std::vector<MatrixXx>& Mlist, const std::vector<MatrixXx>& Glist,
		const MatrixXx& Slist, const MatrixXx& thetamatd, const MatrixXx& dthetamatd, const MatrixXx& ddthetamatd,
		const VectorXx& gtilde, const std::vector<MatrixXx>& Mtildelist, const std::vector<MatrixXx>& Gtildelist,
		numericX Kp, numericX Ki, numericX Kd, numericX dt, int intRes) {
		MatrixXx FtipmatT = Ftipmat.transpose();
		MatrixXx thetamatdT = thetamatd.transpose();
		MatrixXx dthetamatdT = dthetamatd.transpose();
		MatrixXx ddthetamatdT = ddthetamatd.transpose();
		int m = thetamatdT.rows(); int n = thetamatdT.cols();
		VectorXx thetacurrent = thetalist;
		VectorXx dthetacurrent = dthetalist;
		VectorXx eint = VectorXx::Zero(m);
		MatrixXx taumatT = MatrixXx::Zero(m, n);
		MatrixXx thetamatT = MatrixXx::Zero(m, n);
		VectorXx taulist;
		VectorXx ddthetalist;
		for (int i = 0; i < n; ++i) {
			taulist = ComputedTorque(thetacurrent, dthetacurrent, eint, gtilde, Mtildelist, Gtildelist, Slist, thetamatdT.col(i),
				dthetamatdT.col(i), ddthetamatdT.col(i), Kp, Ki, Kd);
			for (int j = 0; j < intRes; ++j) {
				ddthetalist = ForwardDynamics(thetacurrent, dthetacurrent, taulist, g, FtipmatT.col(i), Mlist, Glist, Slist);
				EulerStep(thetacurrent, dthetacurrent, ddthetalist, dt / intRes);
			}
			taumatT.col(i) = taulist;
			thetamatT.col(i) = thetacurrent;
			eint += dt * (thetamatdT.col(i) - thetacurrent);
		}
		std::vector<MatrixXx> ControlTauTraj_ret;
		ControlTauTraj_ret.push_back(taumatT.transpose());
		ControlTauTraj_ret.push_back(thetamatT.transpose());
		return ControlTauTraj_ret;
	}
}

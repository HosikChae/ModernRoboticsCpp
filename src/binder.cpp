//
// Created by pong0923 on 4/27/23.
//

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "modern_robotics.h"

namespace py = pybind11;
using namespace mr;

PYBIND11_MODULE(modern_robotics, m) {
    m.attr("__name__") = "ModernRobotics";
    m.doc() = "Modern Robotics: Mechanics, Planning, and Control C++ Library --- The primary purpose of the provided software is to be easy to read and educational, reinforcing the concepts in the book. The code is optimized neither for efficiency nor robustness. http://modernrobotics.org/";
    m.def("NearZero", &NearZero, "Find if the value is negligible enough to consider 0");
    m.def("ad", &ad, "Calculate the 6x6 matrix [adV] of the given 6-vector");
    m.def("Normalize", &Normalize, "Returns a normalized version of the input vector");
    m.def("VecToso3", &VecToso3, "Returns the skew symmetric matrix representation of an angular velocity vector");
    m.def("so3ToVec", &so3ToVec, "Returns angular velocity vector represented by the skew symmetric matrix");
    m.def("AxisAng3", &AxisAng3, "Tranlates an exponential rotation into it's individual components");
    m.def("MatrixExp3", &MatrixExp3, "Translates an exponential rotation into a rotation matrix");
    m.def("MatrixLog3", &MatrixLog3, "Computes the matrix logarithm of a rotation matrix");
    m.def("RpToTrans", &RpToTrans, "Combines a rotation matrix and position vector into a single Special Euclidian Group (SE3) homogeneous transformation matrix");
    m.def("TransToRp", &TransToRp, "Separates the rotation matrix and position vector from the transfomation matrix representation");
    m.def("VecTose3", &VecTose3, "Translates a spatial velocity vector into a transformation matrix");
    m.def("se3ToVec", &se3ToVec, "Translates a transformation matrix into a spatial velocity vector");
    m.def("Adjoint", &Adjoint, "Provides the adjoint representation of a transformation matrix Used to change the frame of reference for spatial velocity vectors");
    m.def("MatrixExp6", &MatrixExp6, "Rotation expanded for screw axis");
    m.def("MatrixLog6", &MatrixLog6, "Computes the matrix logarithm of a homogeneous transformation matrix");
    m.def("FKinSpace", &FKinSpace, "Compute end effector frame (used for current spatial position calculation)");
    m.def("FKinBody", &FKinBody, "Compute end effector frame (used for current body position calculation)");
    m.def("JacobianSpace", &JacobianSpace, "Gives the space Jacobian");
    m.def("JacobianBody", &JacobianBody, "Gives the body Jacobian");
    m.def("TransInv", &TransInv, "Inverts a homogeneous transformation matrix");
    m.def("RotInv", &RotInv, "Inverts a rotation matrix");
    m.def("ScrewToAxis", &ScrewToAxis, "Takes a parametric description of a screw axis and converts it to a normalized screw axis");
    m.def("AxisAng6", &AxisAng6, "Translates a 6-vector of exponential coordinates into screw axis-angle form");
    m.def("ProjectToSO3", &ProjectToSO3, "Returns projection of one matrix into SO(3)");
    m.def("ProjectToSE3", &ProjectToSE3, "Returns projection of one matrix into SE(3)");
    m.def("DistanceToSO3", &DistanceToSO3, "Returns the Frobenius norm to describe the distance of M from the SO(3) manifold");
    m.def("DistanceToSE3", &DistanceToSE3, "Returns the Frobenius norm to describe the distance of mat from the SE(3) manifold");
    m.def("TestIfSO3", &TestIfSO3, "Returns true if M is close to or on the manifold SO(3)");
    m.def("TestIfSE3", &TestIfSE3, "Returns true if T is close to or on the manifold SE(3)");
    m.def("IKinBody", &IKinBody, "Computes inverse kinematics in the body frame for an open chain robot");
    m.def("IKinSpace", &IKinSpace, "Computes inverse kinematics in the space frame for an open chain robot");
    m.def("InverseDynamics", &InverseDynamics, "This function uses forward-backward Newton-Euler iterations to solve the manipulator equation equation:");
    m.def("GravityForces", &GravityForces, "This function calls InverseDynamics with Ftip = 0, dthetalist = 0, and ddthetalist = 0. The purpose is to calculate one important term in the dynamics equation       ");
    m.def("MassMatrix", &MassMatrix, "This function calls InverseDynamics n times, each time passing a ddthetalist vector with a single element equal to one and all other inputs set to zero. Each call of InverseDynamics generates a single column, and these columns are assembled to create the inertia matrix.       ");
    m.def("VelQuadraticForces", &VelQuadraticForces, "This function calls InverseDynamics with g = 0, Ftip = 0, and ddthetalist = 0.      ");
    m.def("EndEffectorForces", &EndEffectorForces, "This function calls InverseDynamics with g = 0, dthetalist = 0, and ddthetalist = 0.  ");
    m.def("ForwardDynamics", &ForwardDynamics, "This function computes ddthetalist by solving: Mlist(thetalist) * ddthetalist = taulist - c(thetalist,dthetalist) - g(thetalist) - Jtr(thetalist) * Ftip");
    m.def("EulerStep", &EulerStep, "Compute the joint angles and velocities at the next timestep using first order Euler integration");
    m.def("InverseDynamicsTrajectory", &InverseDynamicsTrajectory, "Compute the joint forces/torques required to move the serial chain along the given trajectory using inverse dynamics");
    m.def("ForwardDynamicsTrajectory", &ForwardDynamicsTrajectory, "Compute the motion of a serial chain given an open-loop history of joint forces/torques");
    m.def("ComputedTorque", &ComputedTorque, "Compute the joint control torques at a particular time instant");
    m.def("CubicTimeScaling", &CubicTimeScaling, "Compute s(t) for a cubic time scaling");
    m.def("QuinticTimeScaling", &QuinticTimeScaling, "Compute s(t) for a quintic time scaling");
    m.def("JointTrajectory", &JointTrajectory, "Compute a straight-line trajectory in joint space");
    m.def("ScrewTrajectory", &ScrewTrajectory, "Compute a trajectory as a list of N SE(3) matrices corresponding to the screw motion about a space screw axis");
    m.def("CartesianTrajectory", &CartesianTrajectory, "Compute a trajectory as a list of N SE(3) matrices corresponding to the origin of the end-effector frame following a straight line");
    m.def("SimulateControl", &SimulateControl, "Compute the motion of a serial chain given an open-loop history of joint forces/torques");
}
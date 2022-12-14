// Copyright 2022 VorteX-co
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "../include/plant.hpp"
// ================================================================
void Plant::initialize(
  const double & m, const double & volume,
  const Vector6d & Ib, const Vector3d & r_cob,
  const Vector3d & r_cog, const Vector6d & Ma,
  const Vector6d & Dlinear, const Vector6d & Dquad)
{
  /********************************************************************************
   * Paramteres initilization
   * r_cob_,    vector from center of Buoyancy to the origin
   * r_cog_,    vector from center of Gravity to the origin
   * Ib_,          Inertia vector [ixx, iyy, izz, ixy, ixz, iyz]
   * Mr_,         Rigid body mass matrix
   * Ma_,        Added mass matrix
   * DLinear_, Linear damping coefficients vector
   * Dquad_ ,  Quadratic damping coefficients vector
   ********************************************************************************/
  mass_ = m;
  volume_ = volume;
  r_cob_ = r_cob.cast<CppAD::AD<double>>();
  r_cog_ = r_cog.cast<CppAD::AD<double>>();
  Dlinear_ = Dlinear.cast<CppAD::AD<double>>();
  Dquad_ = Dquad.cast<CppAD::AD<double>>();
  // Rigid-Body System Inertia Matrix, (equation 3.44) Fossen 2011 book
  // Rigid-body mass matrix
  Matrix6d Mr;
  Mr.block<3, 3>(0, 0) = m * Matrix3d::Identity();
  Matrix3d inertia_matrix;
  inertia_matrix << Ib(0), Ib(3), Ib(4), Ib(3), Ib(1), Ib(5), Ib(4), Ib(5),
    Ib(2);
  Mr.setZero();
  Mr.block<3, 3>(3, 3) = inertia_matrix;
  Mr.block<3, 3>(0, 3) = -m * skew(r_cog);
  Mr.block<3, 3>(3, 0) = m * skew(r_cog);
  // Total system Mass Matrix
  Matrix6d _Ma = Ma.asDiagonal();
  M_ = Mr + _Ma;
  M_ad_ = M_.cast<CppAD::AD<double>>();
}
// ================================================================
Vector12d Plant::nonlinear_update(const Vector12d & x, const Vector6d & u)
{
  // Nonlinear motion model update
  // Used for simulation mainly
  // ??? = f (x, u)
  // Since all the equations in this class is developed using CppAD based data types
  // We convert the input interface with <double> based  to CppAD<double>
  Vector12dAD _x = x.cast<CppAD::AD<double>>();
  Vector6dAD _u = u.cast<CppAD::AD<double>>();
  const Vector6dAD & nu = _x.tail<6>();
  const Vector3dAD & euler = _x.segment<3>(3);
  Vector12dAD dx;
  dx.segment<6>(0) = nonlinear_kinematics(euler, nu);
  Vector6dAD fd = nonlinear_damping_forces(nu);
  Vector6dAD fc = nonlinear_coriolis_forces(nu);
  Vector6dAD fg = nonlinear_restoring_forces(euler);
  dx.segment<6>(6) = M_ad_.inverse() * (_u - fd - fc - fg);
  Vector12d return_val;
  using CppAD::Value;
  // Converting from CppAD<double> based datatype
  // to it's base type <double>
  for (size_t i = 0; i < 12; i++) {
    return_val(i) = Value(dx[i]);
  }
  return return_val;
}
// ================================================================
void Plant::linearize(
  const Vector12d & x, Matrix12d & A, Eigen::Matrix<double, 12, 6> & B)
{
  // Computing LTV system parameters A, B
  // See, H2 and H??? Designs for Diving and Course Control of an Autonomous
  // Underwater Vehicle in Presence of Waves, L??cia Moreira. Equation (15) for
  // the notations.
  Matrix6d J, J_star, C, D, G;
  kinematics_jacobian(x.segment<6>(0), x.segment<6>(6), J, J_star);
  G = restoring_forces_jacobian(x.segment<6>(0));
  D = damping_forces_jacobian(x.segment<6>(6));
  C = coriolis_forces_jacobian(x.segment<6>(6));
  A.block<6, 6>(0, 0) = J_star;
  A.block<6, 6>(0, 6) = J;
  A.block<6, 6>(6, 0) = -M_.inverse() * G;
  A.block<6, 6>(6, 6) = -M_.inverse() * (C + D);
  B.setZero();
  B.block<6, 6>(6, 0) = M_.inverse();
}
// ================================================================
Vector6dAD Plant::nonlinear_kinematics(
  const Vector3dAD & euler,
  const Vector6dAD & nu)
{
  /* kinematics transformation matrix (J) 6x6 matrix for transforming  6DOF
   * velocity from Body to Inertial frame Reference: equation (2.40) Fossen 2011
   */
  Matrix6dAD J = Matrix6dAD::Zero();
  J.block<3, 3>(0, 0) = RBtoI(euler);
  J.block<3, 3>(3, 3) = TBtoI(euler);
  Vector6dAD eta_dot = J * nu;
  return eta_dot;
}
// ================================================================
Vector6dAD Plant::nonlinear_damping_forces(const Vector6dAD & nu)
{
  /* Calculating the total damping matrix D
   * D(v) = D_linear + D_nonlinear * |nu|)
   * [nu] is the body velocity vector
   */
  Matrix6dAD D = -1 * Dlinear_.asDiagonal();
  D += -1 * (Dquad_.cwiseProduct(nu.cwiseAbs())).asDiagonal();
  Vector6dAD damping_forces = D * nu;
  return damping_forces;
}
// ================================================================
Vector6dAD Plant::nonlinear_coriolis_forces(const Vector6dAD & nu)
{
  /* Coriolis???Centripetal forces fc.
   * [nu] is the body velocity vector
   * fc = C(nu) * nu
   * C(nu) is the  Coriolis???Centripetal Matrix calculated from the System mass
   * matrix M Reference: equation (3.46) Fossen 2011
   */
  Matrix6dAD C = Matrix6dAD::Zero();
  C.block<3, 3>(0, 3) = -skew_ad(
    M_ad_.block<3, 3>(0, 0) * nu.head<3>() +
    M_ad_.block<3, 3>(0, 3) * nu.tail<3>());
  C.block<3, 3>(3, 0) = -skew_ad(
    M_ad_.block<3, 3>(0, 0) * nu.head<3>() +
    M_ad_.block<3, 3>(0, 3) * nu.tail<3>());
  C.block<3, 3>(3, 3) = -skew_ad(
    M_ad_.block<3, 3>(3, 0) * nu.head<3>() +
    M_ad_.block<3, 3>(3, 3) * nu.tail<3>());
  Vector6dAD coriolis_forces = C * nu;
  return coriolis_forces;
}
// ================================================================
Vector6dAD Plant::nonlinear_restoring_forces(const Vector3dAD & euler)
{
      /* Restoring forces and moments
       * Reference: equation (4.5) Fossen 2011
       */
        Vector6dAD g = Vector6dAD::Zero();
        g << 0.0, 0.0, 0.0, 0.20 * mass_ * 9.806 *  CppAD::cos(euler(1)) * CppAD::sin(euler(0)), 0.2 * mass_ * 9.806 * CppAD::sin(euler(1)), 0.0;
        return g;
}
// ================================================================
Matrix6d Plant::damping_forces_jacobian(const Vector6d & nu)
{
  // Dynamic vector based on CppAD::AD<double>
  VectorXdAD _nu(6);
  _nu = nu.cast<CppAD::AD<double>>();
  // Declare independent variables and starting recording
  CppAD::Independent(_nu);
  // Nonlinear damping forces
  VectorXdAD fd = nonlinear_damping_forces(_nu);
  // Create f: _nu -> fd and stop tape recording
  auto f = CppAD::ADFun<double>(_nu, fd);
  // Equilibrium state / operating point for linearization
  VectorXd nu_eq(6);
  nu_eq << 0.4, 0.0, 0.3, 0.0, 0.0, 0.0;
  Eigen::Matrix<double, 6, 6, Eigen::RowMajor> D;
  Eigen::Map<VectorXd> J_map(D.data(), D.size());
  // Compute the derivative fd'(nu_eq)
  J_map << f.Jacobian(nu_eq);
  return D;
}

// =========================================================================
void Plant::kinematics_jacobian(
  const Vector6d & eta, const Vector6d & nu,
  Matrix6d & J, Matrix6d & J_star)
{
  // Linearizing the 6DOF kinematics  ???? = J(??) v.
  // See H2 and H??? Designs, L??cia Moreira. Equation (13)
  // ???? ??? J  ????? + J* ?????
  // J = J(??0) * ?????
 // J* = (???J(??)/?????) ??0
 // Stacking the pose ?? and the velocity ??
  VectorXdAD input(12);
  input.segment<6>(0) = eta.cast<CppAD::AD<double>>();
  input.segment<6>(6) = nu.cast<CppAD::AD<double>>();
  // Declare independent variables [domain] and starting recording
  CppAD::Independent(input);
  const Vector3dAD & euler_ad = input.segment<3>(3);
  const Vector6dAD & nu_ad = input.segment<6>(6);
  // Dependent variable [range]
  //  ???? = J(??) v
  Vector6dAD eta_dot = nonlinear_kinematics(euler_ad, nu_ad);
  auto f = CppAD::ADFun<double>(input, VectorXdAD(eta_dot));
  VectorXd input_eq(12);
  // Compute the derivative at this domain values
  // Equilibrium points for linearization
  input_eq.setZero();
 input_eq.segment<3>(6) = Vector3d(0.4, 0.0, 0.3);
  Eigen::Matrix<double, 6, 12, Eigen::RowMajor> Jac;
  Eigen::Map<VectorXd> J_map(Jac.data(), Jac.size());
  J_map << f.Jacobian(VectorXd(input_eq));
  J_star = Jac.block<6, 6>(0, 0);
  J = Jac.block<6, 6>(0, 6);
}
// ================================================================
Matrix6d Plant::restoring_forces_jacobian(const Vector6d & eta)
{
  // Linearizing the restoring forces g(??)
  // g(??) ??? G * ?????
  // G = ???g(??) / ?????
  // Converting input from <double> based
  // to CppAD<double>
  VectorXdAD _eta(6);
  _eta = eta.cast<CppAD::AD<double>>();
  // Declare independent variables [domain] and starting recording
  CppAD::Independent(_eta);
  const Vector3dAD & euler = _eta.tail<3>();
  // Dependent variable [range]
  VectorXdAD g(6);
  g = nonlinear_restoring_forces(euler);
  auto f = CppAD::ADFun<double>(_eta, g);
  // Compute the derivative at this domain values
  // Equilibrium points for linearization
  VectorXd eta_eq(6);
  eta_eq.setZero();
  // Forming the jacobian square matrix. [converting (36,1) vector to (6,6) matrix]
  // The normal output from CppAD is a stacked vector
  Eigen::Matrix<double, 6, 6, Eigen::RowMajor> G;
  Eigen::Map<VectorXd> J_map(G.data(), G.size());
  J_map << f.Jacobian(eta_eq);
  return G;
}
// =========================================================================
Matrix6d Plant::coriolis_forces_jacobian(const Vector6d & nu)
{
  // Declare independent variables [domain] and starting recording
  VectorXdAD _nu(6);
  _nu = nu.cast<CppAD::AD<double>>();
  CppAD::Independent(_nu);
  // Dependent variable [range]
  VectorXdAD fc(6);
  fc = nonlinear_coriolis_forces(_nu);
  auto f = CppAD::ADFun<double>(_nu, fc);
  // Compute the derivative at this domain values
  // Equilibrium points for linearization
  VectorXd nu_eq(6);
  nu_eq << 0.4, 0.0, 0.3, 0.0, 0.0, 0.0;
  // nu_eq.segment<3>(0) = nu.head<3>();
  // Forming the jacobian square matrix. [converting (36,1) vector to (6,6) matrix]
  // The normal output from CppAD is a stacked vector
  Eigen::Matrix<double, 6, 6, Eigen::RowMajor> C;
  Eigen::Map<VectorXd> J_map(C.data(), C.size());
  J_map << f.Jacobian(nu_eq);
  return C;
}
// ================================================================
Matrix3dAD Plant::RBtoI(const Vector3dAD & euler)
{
  /* Rotation matrix from Body frame to Inertial frame
   * Reference: equation (2.18) Fossen 2011
   * euler -> rpy w.r.t  inertial NED
   */
  CppAD::AD<double> cphi = CppAD::cos(euler(0));
  CppAD::AD<double> sphi = CppAD::sin(euler(0));
  CppAD::AD<double> cth = CppAD::cos(euler(1));
  CppAD::AD<double> sth = CppAD::sin(euler(1));
  CppAD::AD<double> cpsi = CppAD::cos(euler(2));
  CppAD::AD<double> spsi = CppAD::sin(euler(2));
  Matrix3dAD R;
  R << cpsi * cth, -spsi * cphi + cpsi * sth * sphi,
    spsi * sphi + cpsi * cphi * sth, spsi * cth,
    cpsi * cphi + sphi * sth * spsi, -cpsi * sphi + sth * spsi * cphi, -sth,
    cth * sphi, cth * cphi;
  return R;
}
// ================================================================
Matrix3dAD Plant::TBtoI(const Vector3dAD & euler)
{
  /* Transformation matrix from Body frame to Inertial frame
   * Reference: equation (2.28.b) Fossen 2011
   */
  CppAD::AD<double> cphi = CppAD::cos(euler(0));
  CppAD::AD<double> sphi = CppAD::sin(euler(0));
  CppAD::AD<double> cth = CppAD::cos(euler(1));
  // cth = CppAD::sign(cth) * CppAD::min(0.01745, CppAD::abs(cth));
  CppAD::AD<double> sth = CppAD::sin(euler(1));
  Matrix3dAD T;
  T << 1.0, sphi * sth / cth, cphi * sth / cth, 0.0, cphi, -sphi, 0.0,
    sphi / cth, cphi / cth;
  return T;
}
// ================================================================
Matrix3dAD Plant::TItoB(const Vector3dAD & euler)
{
  /* Transformation matrix from Inertial frame to body frame
   * Reference: equation (2.28.a) Fossen 2011
   */
  CppAD::AD<double> cphi = CppAD::cos(euler(0));
  CppAD::AD<double> sphi = CppAD::sin(euler(0));
  CppAD::AD<double> cth = CppAD::cos(euler(1));
  CppAD::AD<double> sth = CppAD::sin(euler(1));
  Matrix3dAD T;
  T << 1.0, -sth, 0.0, 0.0, cphi, cth * sphi, 0.0, -sphi, cth * cphi;
  return T;
}
// ================================================================
Matrix3dAD Plant::skew_ad(const Vector3dAD & v)
{
  // Skew symmetric matrix
  // Matrix for cross-product operation of a vector
  Matrix3dAD S;
  S << 0.0, -v(2), v(1), v(2), 0.0, -v(0), -v(1), v(0), 0.0;
  return S;
}
// ================================================================
Matrix3d Plant::skew(const Vector3d & v)
{
  Matrix3d S;
  S << 0.0, -v(2), v(1), v(2), 0.0, -v(0), -v(1), v(0), 0.0;
  return S;
}
// ================================================================

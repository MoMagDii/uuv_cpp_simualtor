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

#include "../include/lqr.hpp"
#include <eigen3/unsupported/Eigen/src/MatrixFunctions/MatrixExponential.h>

// =========================================================================
void LQR::set_params(
  const double & m, const double & volume, const Vector6d & Ib,
  const Vector3d & r_cob, const Vector3d & r_cog,
  const Vector6d & Ma, const Vector6d & Dlinear,
  const Vector6d & Dquad, const Vector12d & Q,
  const Vector6d & R, const Vector6d & tau_max,
  const Vector12d & error_max)
{
  vehicle_.initialize(m, volume, Ib, r_cob, r_cog, Ma, Dlinear, Dquad);
  tau_max_ = tau_max;
  error_max_ = error_max;
  Q_ = Q.asDiagonal();
  R_ = R.asDiagonal();
  error_integral_ = Vector12d::Zero();
}
// =========================================================================
void LQR::saturate_control(Vector6d & tau)
{
  /* Saturating the control forces and moments at ± tau_max [N, N.m]
   */
  for (int i = 0; i <= 5; i++) {
    if (tau(i) > tau_max_(i)) {
      tau(i) = tau_max_(i);
    } else if (tau(i) < -tau_max_(i)) {
      tau(i) = -tau_max_(i);
    }
  }
}
// =========================================================================
void LQR::saturate_error(Vector12d & delta)
{
  /* Saturating the state errors ± error_max
   */
  for (int i = 0; i <= 11; i++) {
    if (delta(i) > error_max_(i)) {
      delta(i) = error_max_(i);
    } else if (delta(i) < -error_max_(i)) {
      delta(i) = -error_max_(i);
    }
  }
}
// =========================================================================
void LQR::to_SNAME(Vector12d & x)
{
  /* from ENU to NED.
   */
  x(1) *= -1;
  x(2) *= -1;
  x(4) *= -1;
  x(5) *= -1;
  x(7) *= -1;
  x(8) *= -1;
  x(10) *= -1;
  x(11) *= -1;
}
//=========================================================================
Matrix6d LQR::ItoB_transformation(Vector3d & euler)
{
    double cphi = cos(euler(0));
    double sphi = sin(euler(0));
    double cth = cos(euler(1));
    double sth = sin(euler(1));
    double cpsi = cos(euler(2));
    double spsi = sin(euler(2));
    Matrix3d R;
    R << cpsi * cth, -spsi * cphi + cpsi * sth * sphi,
      spsi * sphi + cpsi * cphi * sth, spsi * cth,
      cpsi * cphi + sphi * sth * spsi, -cpsi * sphi + sth * spsi * cphi, -sth,
      cth * sphi, cth * cphi;

    Matrix3d T;
    T << 1.0, -sth, 0.0, 0.0, cphi, cth * sphi, 0.0, -sphi, cth * cphi;
    Matrix6d J_inv;
    J_inv << R.transpose(), Matrix3d::Zero(), T, Matrix3d::Zero();
    return J_inv;
}
// =========================================================================
Vector6d LQR::action(Vector12d x, Vector12d xd, Vector6d feedforward_acc)
{
//  to_SNAME(x);
//  to_SNAME(xd);
  // MatrixXd A, B;
  Matrix12d A;
  Eigen::Matrix<double, 12, 6> B;
  // Calculating the Jacobian based linearization of the dynamics
  vehicle_.linearize(x, A, B);
  // Solving the Continous Algebraic Riccati Equation (CARE)
  Matrix12d P;
   care_solver.solve(P, A, B, Q_, R_);
//  bool done = solve_riccati(A,B,Q_,R_,P, 0.00001, 10000);
  // The optimal feedback gain matrix
  Eigen::Matrix<double, 6, 12> K = -R_.inverse() * B.transpose() * P;
  // LQR control law
  Vector12d error = (x - xd);
  saturate_error(error);
  Vector6d tau = K * error;
  saturate_control(tau);
  display_log(x,xd,tau);
  error_integral_ = simulation_rk4(error, tau);
  std::cout << " ************ e_int *****************" << std::endl;
  std::cout << error_integral_ << std::endl;
  return tau;
}
// =========================================================================
Vector12d LQR::simulation_rk4(const Vector12d & x0, const Vector6d & u)
{
  /*
   * Simulating the control forces on nonlinear vehicle model
   *  Runge-Kutta (4th-order) method for integrating the equations of motion..
   */
  Vector12d x;
  Vector12d k1 = 0.05 * vehicle_.nonlinear_update(x0, u);
  x = x0 + 0.5 * k1;
  Vector12d k2 = 0.05 * vehicle_.nonlinear_update(x, u);
  x = x0 + 0.5 * k2;
  Vector12d k3 = 0.05 * vehicle_.nonlinear_update(x, u);
  x = x0 + k3;
  Vector12d k4 = 0.05 * vehicle_.nonlinear_update(x, u);
  Vector12d xnext = x0 + (k1 + 2 * (k2 + k3) + k4) / 6;
  return xnext;
}
// =========================================================================
Vector12d LQR::simulation_euler(const Vector12d & x0, const Vector6d & u)
{
  /*
   * Simulating the control forces on nonlinear vehicle model
   *  Euler method for integrating the equations of motion.
   */
  Vector12d xnext = x0 + 0.05 * vehicle_.nonlinear_update(x0, u);
  return xnext;
}
// =========================================================================
Vector12d LQR::simulation_rk4_linearized(const Vector12d & x0, const Vector6d & u)
{
  /*
   * Simulating the control forces on nonlinear vehicle model
   *  Runge-Kutta (4th-order) method for integrating the equations of motion..
   */
    Vector12d x;
    Vector12d dx;
    Matrix12d A;
    Vector12d x_eq;
    x_eq << 0.0 , 0.0 ,0.0, 0.0, 0.0, 0.0, 0.4, 0.0, 0.3, 0.0, 0.0, 0.0;
    Eigen::Matrix<double, 12, 6> B;
    // Calculating the Jacobian based linearization of the dynamics
  vehicle_.linearize(x0, A, B);
  dx = A * x0 + B * u;
  Vector12d k1 = 0.05 * dx;
  x = x0 + 0.5 * k1;

  vehicle_.linearize(x, A, B);
  dx = A * x + B * u ;
  Vector12d k2 = 0.05 * dx;
  x = x0 + 0.5 * k2;

  vehicle_.linearize(x, A, B);
  dx = A * x + B * u ;
  Vector12d k3 = 0.05 * dx;
  x = x0 + k3;

  vehicle_.linearize(x, A, B);
  dx = A * x + B * u ;
  Vector12d k4 = 0.05 * dx;

  Vector12d xnext = x0 + (k1 + 2 * (k2 + k3) + k4) / 6;
  return xnext;
}
// =========================================================================
void LQR::display_log(const Vector12d & x, const Vector12d & xd, const Vector6d & u)
{
    std::cout << " ###### Logging controller status ##### "  << std::endl;
    std::cout << "$ x:  " << x(0) << " <> " << "xd: " << xd(0) << " <> " << "e: " << x(0) - xd(0) << std::endl;
    std::cout << "$ y:  " << x(1) << " <> " << "yd: " << xd(1) << " <> " << "e: " << x(1) - xd(1) << std::endl;
    std::cout << "$ z:  " << x(2) << " <> " << "zd: " << xd(2) << " <> " << "e: " << x(2) - xd(2) << std::endl;
    std::cout << "$ φ:  " << x(3) << " <> " << "φd: " << xd(3) << " <> " << "e: " << x(3) - xd(3) << std::endl;
    std::cout << "$ θ:  " << x(4) << " <> " << "θd: " << xd(4) << " <> " << "e: " << x(4) - xd(4) << std::endl;
    std::cout << "$ ψ:  " << x(5) << " <> " << "ψd: " << xd(5) << " <> " << "e: " << x(5) - xd(5) << std::endl;
    std::cout << "$ u:  " << x(6) << " <> " << "ud: " << xd(6) << " <> " << "e: " << x(6) - xd(6) << std::endl;
    std::cout << "$ v:  " << x(7) << " <> " << "vd: " << xd(7) << " <> " << "e: " << x(7) - xd(7) << std::endl;
    std::cout << "$ w:  " << x(8) << " <> " << "wd: " << xd(8) << " <> " << "e: " << x(8) - xd(8) << std::endl;
    std::cout << "$ p:  " << x(9) << " <> " << "pd: " << xd(9) << " <> " << "e: " << x(9) - xd(9) << std::endl;
    std::cout << "$ q:  " << x(10) << " <> " << "qd: " << xd(10) << " <> " << "e: " << x(10) - xd(10) << std::endl;
    std::cout << "$ r:  " << x(11) << " <> " << "rd: " << xd(11) << " <> " << "e: " << x(11) - xd(11) << std::endl;
    std::cout << " ##############################"  << std::endl;
     std::cout << u  << std::endl;
}


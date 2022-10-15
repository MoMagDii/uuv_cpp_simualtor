#include <iostream>
#include <cmath>
#include <matplot/matplot.h>
#include <set>
#include <thread>
#include <vector>
#include "../include/lqr.hpp"
#include  "../include/optimal_trajectory.hpp"
#include "matplotlibcpp.h"
 using namespace matplot;
namespace plt = matplotlibcpp;
using namespace std;

int main() {

    /**************************************************************************
     * Simulation Variables
     ***************************************************************************/
    LQR lqr;
    Trajectory trajectory_generator;    

    std::vector<double> TimeStamp;

    std::vector<double> x_vals;
    std::vector<double> y_vals;
    std::vector<double> z_vals;
    std::vector<double> roll_vals;
    std::vector<double> pitch_vals;
    std::vector<double> yaw_vals;
    std::vector<double> u_vals;
    std::vector<double> v_vals;
    std::vector<double> w_vals;
    std::vector<double> r_vals;
    std::vector<double> p_vals;
    std::vector<double> q_vals;


    std::vector<double> x_linearized;
    std::vector<double> y_linearized;
    std::vector<double> z_linearized;
    std::vector<double> roll_linearized;
    std::vector<double> pitch_linearized;
    std::vector<double> yaw_linearized;
    std::vector<double> u_linearized;
    std::vector<double> v_linearized;
    std::vector<double> w_linearized;
    std::vector<double> r_linearized;
    std::vector<double> p_linearized;
    std::vector<double> q_linearized;


    std::vector<double> x_d;
    std::vector<double> y_d;
    std::vector<double> z_d;
    std::vector<double> roll_d;
    std::vector<double> pitch_d;
    std::vector<double> yaw_d;
    std::vector<double> u_d;
    std::vector<double> v_d;
    std::vector<double> w_d;
    std::vector<double> r_d;
    std::vector<double> p_d;
    std::vector<double> q_d;

    std::vector<double> xe_vals;
    std::vector<double> ye_vals;
    std::vector<double> ze_vals;
    std::vector<double> rolle_vals;
    std::vector<double> pitche_vals;
    std::vector<double> yawe_vals;
    std::vector<double> ue_vals;
    std::vector<double> ve_vals;
    std::vector<double> we_vals;
    std::vector<double> re_vals;
    std::vector<double> pe_vals;
    std::vector<double> qe_vals;

    std::vector<double> xe_int_vals;
    std::vector<double> ye_int_vals;
    std::vector<double> ze_int_vals;
    std::vector<double> rolle_int_vals;
    std::vector<double> pitche_int_vals;
    std::vector<double> yawe_int_vals;
    std::vector<double> ue_int_vals;
    std::vector<double> ve_int_vals;
    std::vector<double> we_int_vals;
    std::vector<double> re_int_vals;
    std::vector<double> pe_int_vals;
    std::vector<double> qe_int_vals;

    std::vector<double> x_u;
    std::vector<double> y_u;
    std::vector<double> z_u;
    std::vector<double> roll_u;
    std::vector<double> pitch_u;
    std::vector<double> yaw_u;

    std::vector<double> J_vals;
    std::vector<double> s_vals;
    std::vector<double> s_dot_vals;
    std::vector<double> curvature_vals;
    std::vector<double> e_vals;
    std::vector<double> e_int_vals;


    /**************************************************************************
     * Variables Initialization
     ***************************************************************************/
    /******************************
     *  AUV parameters
     * ****************************/
    double mass = 36.7;
    double volume = 0.036;
    Eigen::VectorXd Ib(6);
    Ib << 0.50, 0.49, 0.8599, -0.0059, 0.0005, -0.0113;
    Eigen::VectorXd Dlinear(6);
    Dlinear <<  -22.8, -22.95, -50.26, -13.05, -13.2, -4.1559;
    Eigen::VectorXd Dquad(6);
    Dquad << -7.43, -7.98, -26.5, -0.0, -0.0, -0.00;
    Eigen::VectorXd Ma(6);
    Ma << 13.77, 20.77, 50.01, 17.02, 18.02, 12.7;
    Eigen::VectorXd r_cog(3);
    r_cog << 0.0, 0.0, 0.0;
    Eigen::VectorXd r_cob(3);
    r_cob << 0.0, 0.0, 0.20;
    /****************************
     * Controller parameters
     * ***************************/
    Eigen::VectorXd Q(12);
    Q << 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0;
    Eigen::VectorXd R1(6);
    R1 << 0.01, 0.025, 0.01, 0.025, 0.025, 0.015;
    Eigen::VectorXd tau_max(6);
    tau_max << 45.0, 45.0, 45.0, 45.0, 45.0, 45.0;
    Eigen::VectorXd error_max(12);
    error_max << 1.0, 1.0, 1.0, 1.0, 1., 1., 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;
    // Setting the Controller parameters
    lqr.set_params(
      mass, volume, Ib, r_cob, r_cog, Ma, Dlinear, Dquad, Q,
      R1, tau_max, error_max);
    /****************************
     * Controller parameters
     * ***************************/

    Eigen::VectorXd translation_constraints(3);
    translation_constraints << 0.6, 0.3, 0.0001;
    Eigen::VectorXd rotation_constraints(3);
    rotation_constraints << 0.6, 0.3, 0.0001;
    trajectory_generator.set_params(translation_constraints, rotation_constraints);


    // Starting variables
    Vector12d state, state_linearized, error, error_int;
     Vector12d state_d;
    state = Vector12d::Zero();
    state_linearized = Vector12d::Zero();
    state_d = Vector12d::Zero();
   error = Vector12d::Zero();
    error_int = Vector12d::Zero();
    Vector12d reference;
    reference  << 1.8, 0.5, 0.5, 0.0, 0.0, 0.30754131 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    Vector12d desired_state;
    Vector12d desired_inertial_state;
    Vector6d desired_pose, desired_velocity, desired_acceleration;

    /**************************************************************************
     * Generating a trajectory
     ***************************************************************************/
    trajectory_generator.generate_trajectory(state, reference);

    /**************************************************************************
     * Simulation loop
     ***************************************************************************/
    double t = 0.0;
    while (t <= 75 ) {

        trajectory_generator.evaluate_generated_trajectory(t, desired_pose, desired_velocity, desired_acceleration);

        Vector3d euler = desired_pose.tail<3>();
        // desired_state << desired_pose, lqr.ItoB_transformation(euler) * desired_velocity ;
        double xdot = desired_velocity(0);
        double ydot = desired_velocity(1);
        double zdot = desired_velocity(2);
        double xdot2 = desired_acceleration(0);
        double ydot2 = desired_acceleration(1);
        double s_dot = std::sqrt(xdot * xdot + ydot * ydot + zdot * zdot);

        // state_d = lqr.simulation_rk4(state, tau);

        double kappa = 1 / sqrt(21 * 21 + 14 * 14);

        double chi = atan2(ydot,xdot);
        double gamma = -atan2(zdot, sqrt(xdot * xdot + ydot * ydot));
        Vector6d virtual_traget_velocity;
        virtual_traget_velocity << cos(chi - 0.30754131) * cos(gamma) * s_dot, sin(chi - 0.30754131) * cos(gamma) * s_dot , -sin(gamma) * s_dot , 0.0, 0.0, -kappa * sin(chi - 0.30754131) * cos(gamma) * s_dot;
        Vector6d virtual_traget_pose;
        virtual_traget_pose = desired_pose;
        double u = state(6);
        double v = state(7);
        double w = state(8);
        double alpha = atan2(w,u);
        double theta = state(4);
        double Uv = sqrt(u * u + w * w);
        double beta_generalized = atan2(v, Uv * cos(theta - alpha));
        desired_state << virtual_traget_pose, virtual_traget_velocity;


        Vector6d tau = lqr.action(state, desired_state, desired_acceleration);
        state = lqr.simulation_rk4(state, tau);
        state_linearized = lqr.simulation_rk4_linearized(state_linearized, tau);
        error =  state - desired_state;
        error_int += lqr.simulation_rk4(error, tau);
        double J = (state - desired_state).cwiseProduct(Q).squaredNorm();
        J += tau.cwiseProduct(R1).squaredNorm();
        TimeStamp.push_back(t);

        x_vals.push_back(state(0));
        y_vals.push_back(state(1));
        z_vals.push_back(state(2));
        roll_vals.push_back(state(3));
        pitch_vals.push_back(state(4));
        yaw_vals.push_back(state(5));
        u_vals.push_back(state(6));
        v_vals.push_back(state(7));
        w_vals.push_back(state(8));
        r_vals.push_back(state(9));
        p_vals.push_back(state(10));
        q_vals.push_back(state(11));

        x_linearized.push_back(state_linearized(0));
        y_linearized.push_back(state_linearized(1));
        z_linearized.push_back(state_linearized(2));
        roll_linearized.push_back(state_linearized(3));
        pitch_linearized.push_back(state_linearized(4));
        yaw_linearized.push_back(state_linearized(5));
        u_linearized.push_back(state_linearized(6));
        v_linearized.push_back(state_linearized(7));
        w_linearized.push_back(state_linearized(8));
        r_linearized.push_back(state_linearized(9));
        p_linearized.push_back(state_linearized(10));
        q_linearized.push_back(state_linearized(11));


        x_d.push_back(desired_state(0));
        y_d.push_back(desired_state(1));
        z_d.push_back(desired_state(2));
        roll_d.push_back(desired_state(3));
        pitch_d.push_back(desired_state(4));
        yaw_d.push_back(desired_state(5));
        u_d.push_back(desired_state(6));
        v_d.push_back(desired_state(7));
        w_d.push_back(desired_state(8));
        r_d.push_back(desired_state(9));
        p_d.push_back(desired_state(10));
        q_d.push_back(desired_state(11));

        xe_vals.push_back(error(0));
        ye_vals.push_back(error(1));
        ze_vals.push_back(error(2));
        rolle_vals.push_back(error(3));
        pitche_vals.push_back(error(4));
        yawe_vals.push_back(error(5));
        ue_vals.push_back(error(6));
        ve_vals.push_back(error(7));
        we_vals.push_back(error(8));
        re_vals.push_back(error(9));
        pe_vals.push_back(error(10));
        qe_vals.push_back(error(11));

        xe_int_vals.push_back(error_int(0));
        ye_int_vals.push_back(error_int(1));
        ze_int_vals.push_back(error_int(2));
        rolle_int_vals.push_back(error_int(3));
        pitche_int_vals.push_back(error_int(4));
        yawe_int_vals.push_back(error_int(5));
        ue_int_vals.push_back(error_int(6));
        ve_int_vals.push_back(error_int(7));
        we_int_vals.push_back(error_int(8));
        re_int_vals.push_back(error_int(9));
        pe_int_vals.push_back(error_int(10));
        qe_int_vals.push_back(error_int(11));

        x_u.push_back(tau(0));
        y_u.push_back(tau(1));
        z_u.push_back(tau(2));
        roll_u.push_back(tau(3));
        pitch_u.push_back(tau(4));
        yaw_u.push_back(tau(5));

        J_vals.push_back(J);
        s_vals.push_back(0.0);
        s_dot_vals.push_back(s_dot);
        curvature_vals.push_back(kappa);

        t += 0.05; // 20 Hz
    }
    /**************************************************************************
     * Plotting
     ***************************************************************************/




     // Plot values
     // NOTE: feel free to play around with this.
     // It's useful for debugging!
     plt::subplot(3, 3, 1);
     plt::title("Position Tracking");
     plt::plot(x_vals, "k");
     plt::plot(x_d, "--r");
     plt::plot(y_vals, "g");
     plt::plot(y_d, "--b");
     plt::plot(z_vals, "m");
     plt::plot(z_d, "--c");


     plt::subplot(3, 3, 2);
     plt::title("Orientation Tracking");
     plt::plot(roll_vals,"k");
     plt::plot(roll_d, "--r");
     plt::plot(pitch_vals,"g");
     plt::plot(pitch_d, "--b");
     plt::plot(yaw_vals,"m");
     plt::plot(yaw_d, "--c");

     plt::subplot(3, 3, 3);
     plt::title("Linear Velocity  Tracking");
     plt::plot(u_vals,"k");
     plt::plot(u_d, "--r");
     plt::plot(v_vals,"g");
     plt::plot(v_d, "--b");
     plt::plot(w_vals,"m");
     plt::plot(w_d, "--c");

     plt::subplot(3, 3, 4);
     plt::title("Angular Velocity  Tracking");
     plt::plot(p_vals,"k");
     plt::plot(p_d, "--r");
     plt::plot(q_vals,"g");
     plt::plot(q_d, "--b");
     plt::plot(r_vals,"m");
     plt::plot(r_d, "--y");

     plt::subplot(3, 3, 5);
     plt::title("Control forces");
     plt::plot(x_u,"k");
     plt::plot(y_u, "r");
     plt::plot(z_u,"m");
     plt::plot(roll_u, "b");
     plt::plot(pitch_u,"g");
     plt::plot(yaw_d, "y");

     plt::subplot(3, 3, 6);
     plt::title("Linearized Pose");
     plt::plot(x_vals,"k");
     plt::plot(x_linearized, "--r");
     plt::plot(y_vals,"b");
     plt::plot(y_linearized, "--m");
     plt::plot(z_vals,"y");
     plt::plot(z_linearized, "--c");
     plt::plot(roll_vals,"y");
     plt::plot(roll_linearized, "--r");
     plt::plot(pitch_vals,"k");
     plt::plot(pitch_linearized, "--c");
     plt::plot(yaw_vals,"k");
     plt::plot(yaw_linearized, "--m");

     plt::subplot(3, 3, 7);
     plt::title("Linearized Pose");
     plt::plot(u_vals,"k");
     plt::plot(u_linearized, "--r");
     plt::plot(v_vals,"m");
     plt::plot(v_linearized, "--c");
     plt::plot(w_vals,"y");
     plt::plot(w_linearized, "--g");
     plt::plot(p_vals,"k");
     plt::plot(p_linearized, "--r");
     plt::plot(q_vals,"k");
     plt::plot(q_linearized, "--r");
     plt::plot(r_vals,"k");
     plt::plot(r_linearized, "--r");


     plt::subplot(3, 3, 8);
     plt::title("Error");
     plt::plot(xe_vals,"k");
     plt::plot(ye_vals, "--r");
     plt::plot(ze_vals,"k");
     plt::plot(rolle_vals, "--r");
     plt::plot(pitche_vals,"k");
     plt::plot(yawe_vals, "--r");
//     plt::plot(ue_vals,"k");
//     plt::plot(ve_vals, "--r");
//     plt::plot(we_vals,"k");
//     plt::plot(pe_vals, "--r");
//     plt::plot(qe_vals,"k");
//     plt::plot(r_vals, "--r");

     plt::subplot(3, 3, 9);
     plt::title("Error Integral");
     plt::plot(xe_int_vals,"k");
     plt::plot(ye_int_vals, "--r");
     plt::plot(ze_int_vals,"m");
     plt::plot(rolle_int_vals, "--c");
     plt::plot(pitche_int_vals,"g");
     plt::plot(yawe_int_vals, "--y");
//     plt::plot(ue_int_vals,"k");
//     plt::plot(ve_int_vals, "--r");
//     plt::plot(we_int_vals,"k");
//     plt::plot(pe_int_vals, "--r");
//     plt::plot(qe_int_vals,"k");
//     plt::plot(re_int_vals, "--r");

    // Set the "super title"




     plt::show();
     return 0;
}

#include "MPC.h"
#include "prm.h"

#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// Instantiate global parameters
Prm prm;

// The solver takes all the state variables and actuator variables in a singular vector. 
// Thus, we need to establish when one variable starts and another ends to make our life easier.
const size_t x_start = 0;
const size_t y_start = x_start + prm.N;
const size_t psi_start = y_start + prm.N;
const size_t v_start = psi_start + prm.N;
const size_t cte_start = v_start + prm.N;
const size_t epsi_start = cte_start + prm.N;
const size_t delta_start = epsi_start + prm.N;
const size_t a_start = delta_start + prm.N - 1;

class FG_eval {
    public:
        // Fitted polynomial coefficients
        Eigen::VectorXd coeffs;
        FG_eval(Eigen::VectorXd coeffs) { this -> coeffs = coeffs; }
        
        typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
        // fg is a vector containing the cost and constraints.
        // vars is a vector containing the variable values (state & actuators).
        void operator()(ADvector& fg, const ADvector& vars) {
            // The cost is stored is the first element of fg.
            fg[0] = 0;
            
            // Reference State Cost: the part of the cost based on the reference state.
            for (int t = 0; t < prm.N; t++) {
                fg[0] += prm.w_cte * CppAD::pow(vars[cte_start + t] - prm.ref_cte, 2);
                fg[0] += prm.w_epsi * CppAD::pow(vars[epsi_start + t] - prm.ref_epsi, 2);
                fg[0] += prm.w_v * CppAD::pow(vars[v_start + t] - prm.ref_v, 2);
            }
            
            // Minimize actuations (to smooth turns and velocity)
            for (int t = 0; t < prm.N - 1; t++) {
                fg[0] += prm.w_delta * CppAD::pow(vars[delta_start + t], 2);
                fg[0] += prm.w_a * CppAD::pow(vars[a_start + t], 2);
                fg[0] += prm.w_delta_v * CppAD::pow(vars[delta_start + t] * vars[v_start + t], 2);
            }
        
            // Minimize gap between sequential actuations (to make decisions consistent over time)
            for (int t = 0; t < prm.N - 2; t++) {
                fg[0] += prm.w_ddelta * CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2);
                fg[0] += prm.w_da * CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
            }
            
            
            // Initial constraints
            // We add 1 to each of the starting indices due to cost being located at index 0 of `fg`.
            // This bumps up the position of all the other values.
            fg[1 + x_start] = vars[x_start];
            fg[1 + y_start] = vars[y_start];
            fg[1 + psi_start] = vars[psi_start];
            fg[1 + v_start] = vars[v_start];
            fg[1 + cte_start] = vars[cte_start];
            fg[1 + epsi_start] = vars[epsi_start];
            
            // Other constraints
            //std::cout << "prm.dt = " << prm.dt << std::endl;
            for (int t = 1; t < prm.N; t++) {
                // State at time t.
                AD<double> x1 = vars[x_start + t];
                AD<double> y1 = vars[y_start + t];
                AD<double> psi1 = vars[psi_start + t];
                AD<double> v1 = vars[v_start + t];
                AD<double> cte1 = vars[cte_start + t];
                AD<double> epsi1 = vars[epsi_start + t];
            
                // State at time t-1.
                AD<double> x0 = vars[x_start + t - 1];
                AD<double> y0 = vars[y_start + t - 1];
                AD<double> psi0 = vars[psi_start + t - 1];
                AD<double> v0 = vars[v_start + t - 1];
                AD<double> cte0 = vars[cte_start + t - 1];
                AD<double> epsi0 = vars[epsi_start + t - 1];
                
                // Actuation at time t-1.
                AD<double> delta0 = vars[delta_start + t - 1];
                AD<double> a0 = vars[a_start + t - 1];
                
                // Predicted trajectory at time t: polynomial
                AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * CppAD::pow(x0, 2) + coeffs[3] * CppAD::pow(x0, 3);
                
                // Predicted heading at time t: arctan of tangent to polynomial
                AD<double> psi_des0 = CppAD::atan(coeffs[1] + 2 * coeffs[2] * x0 + 3 * coeffs[3] * CppAD::pow(x0, 2));
                
                // Kinematic model
                fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * prm.dt);
                fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * prm.dt);
                fg[1 + psi_start + t] = psi1 - (psi0 + v0 * delta0 / prm.Lf * prm.dt);
                fg[1 + v_start + t] = v1 - (v0 + a0 * prm.dt);
                fg[1 + cte_start + t] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * prm.dt));
                fg[1 + epsi_start + t] = epsi1 - ((psi0 - psi_des0) + v0 * delta0 / prm.Lf * prm.dt);
            }
        }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<vector<double>> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
    // Boolean value for solver
    bool ok = true;
    //size_t i;
    typedef CPPAD_TESTVECTOR(double) Dvector;
    
    // State variables for easier use
    const double x = state[0];
    const double y = state[1];
    const double psi = state[2];
    const double v = state[3];
    const double cte = state[4];
    const double epsi = state[5];
    
    // Number of independent variables: N timesteps == N - 1 actuations
    const size_t n_vars = prm.N * 6 + (prm.N - 1) * 2;
    
    // Number of constraints
    const size_t n_constraints = prm.N * 6;
    
    // Initial value of independent variables: should be 0 besides initial state.
    Dvector vars(n_vars);
    for (int i = 0; i < n_vars; i++) {
        vars[i] = 0.0;
    }
    // Variables for lower and upper limits
    Dvector vars_lowerbound(n_vars);
    Dvector vars_upperbound(n_vars);
    
    // Lower and upper limits for state variables
    for (int i = 0; i < delta_start; i++) {
        vars_lowerbound[i] = - prm.max_var;
        vars_upperbound[i] = prm.max_var;
    }
    
    // Lower and upper limits for delta
    for (int i = delta_start; i < a_start; i++) {
        vars_lowerbound[i] = - prm.max_delta;
        vars_upperbound[i] = prm.max_delta;
    }
    
    // Lower and upper limits for acceleration / decceleration
    for (int i = a_start; i < n_vars; i++) {
        vars_lowerbound[i] = prm.min_thr;
        vars_upperbound[i] = prm.max_thr;
    }
    
    
    // Variable for lower and upper limits on constraints
    Dvector constraints_lowerbound(n_constraints);
    Dvector constraints_upperbound(n_constraints);
    
    // Lower and upper limits for constraints: should be 0 besides initial state.
    for (int i = 0; i < n_constraints; i++) {
        constraints_lowerbound[i] = 0.0;
        constraints_upperbound[i] = 0.0;
    }
    
    // Lower limits for state
    constraints_lowerbound[x_start] = x;
    constraints_lowerbound[y_start] = y;
    constraints_lowerbound[psi_start] = psi;
    constraints_lowerbound[v_start] = v;
    constraints_lowerbound[cte_start] = cte;
    constraints_lowerbound[epsi_start] = epsi;
    
    // Upper limits for state
    constraints_upperbound[x_start] = x;
    constraints_upperbound[y_start] = y;
    constraints_upperbound[psi_start] = psi;
    constraints_upperbound[v_start] = v;
    constraints_upperbound[cte_start] = cte;
    constraints_upperbound[epsi_start] = epsi;

    // Object that computes objective and constraints
    FG_eval fg_eval(coeffs);
    
    // Options for IPOPT solver
    std::string options;
    options += "Integer print_level  0\n";
    // NOTE: Setting sparse to true allows the solver to take advantage of sparse routines, this makes the computation MUCH FASTER. If you can uncomment 1 of these and see if it makes a difference or not but if you uncomment both the computation time should go up in orders of magnitude.
    options += "Sparse  true        forward\n";
    options += "Sparse  true        reverse\n";
    // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
    options += "Numeric max_cpu_time          0.5\n";
    
    // Placeholder to return solution
    CppAD::ipopt::solve_result<Dvector> solution;
    //std::cout << "solution.capacity() = " << solution.capacity() << std::endl;
    
    // Solve the problem
    CppAD::ipopt::solve<Dvector, FG_eval>(
        options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
        constraints_upperbound, fg_eval, solution);
    
    // Check some of solution values
    ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

    // Cost
    auto cost = solution.obj_value;
    //std::cout << "solution.status = " << ok << " ; Cost = " << cost << std::endl;
    
    // Actuators and waypoints
    vector<double> mpc_delta(prm.N - 1);
    vector<double> mpc_a(prm.N - 1);
    vector<double> mpc_pts_x(prm.N);
    vector<double> mpc_pts_y(prm.N);
    for (int i = 0; i < prm.N; i++) {
        if (i < prm.N - 1) {
            mpc_delta[i] = solution.x[delta_start + i];
            mpc_a[i] = solution.x[a_start + i];
        }
        mpc_pts_x[i] = solution.x[x_start + i];
        mpc_pts_y[i] = solution.x[y_start + i];
    }
    
    return {mpc_delta, mpc_a, mpc_pts_x, mpc_pts_y};
}
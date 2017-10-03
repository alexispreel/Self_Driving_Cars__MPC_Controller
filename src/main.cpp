#include "MPC.h"
#include "prm.h"

#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <ctime>
#include <unistd.h>
#include <sys/time.h>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"

// For convenience
using json = nlohmann::json;
using namespace std;
using namespace Eigen;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned, else the empty string "" will be returned.
string hasData(string s) {
    auto found_null = s.find("null");
    auto b1 = s.find_first_of("[");
    auto b2 = s.rfind("}]");
    if (found_null != string::npos) {
        return "";
    } else if (b1 != string::npos && b2 != string::npos) {
        return s.substr(b1, b2 - b1 + 2);
    }
    return "";
}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
    double result = 0.0;
    for (int i = 0; i < coeffs.size(); i++) {
        result += coeffs[i] * pow(x, i);
    }
    return result;
}

// Fit a polynomial.
// Adapted from https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals, int order) {
    assert(xvals.size() == yvals.size());
    assert(order >= 1 && order <= xvals.size() - 1);
    Eigen::MatrixXd A(xvals.size(), order + 1);
    
    for (int i = 0; i < xvals.size(); i++) {
        A(i, 0) = 1.0;
    }
    for (int j = 0; j < xvals.size(); j++) {
        for (int i = 0; i < order; i++) {
        A(j, i + 1) = A(j, i) * xvals(j);
        }
    }
    
    auto Q = A.householderQr();
    auto result = Q.solve(yvals);
    return result;
}

int main() {
    uWS::Hub h;
    
    // Time variables
    //double t_cur;
    //double t_prv;
    //t_prv = clock();
    
    // Time variables
    //timeval time_cur;
    //timeval time_prv;
    //gettimeofday(&time_prv, 0);
    
    // Instantiate global parameters
    Prm prm;
    
    // MPC is initialized here!
    MPC mpc;
    
    h.onMessage([&prm, &mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode) {
    //h.onMessage([&prm, &mpc, &t_cur, &t_prv, &time_cur, &time_prv](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode) {
        // "42" at the start of the message means there's a websocket message event.
        // The 4 signifies a websocket message
        // The 2 signifies a websocket event
        string sdata = string(data).substr(0, length);
        using namespace std;
        //std::cout << sdata << std::endl;
        //std::cout << "************************" << std::endl;
        if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
            string s = hasData(sdata);
            if (s != "") {
                auto j = json::parse(s);
                string event = j[0].get<string>();
                if (event == "telemetry") {
                    // j[1] is the data JSON object
                    
                    // Get time delta
                    //t_cur = clock();
                    //std::cout << "time delta = " << (t_cur - t_prv) / CLOCKS_PER_SEC << std::endl;
                    //std::cout << "time delta = " << (t_cur - t_prv) / 1000 << std::endl;
                    //t_prv = t_cur;
                    
                    // Get time delta
                    //gettimeofday(&time_cur, 0);
                    //std::cout << "time delta = " << (time_cur.tv_sec - time_prv.tv_sec) << std::endl;
                    //time_prv = time_cur;
                    
                    // Waypoints
                    vector<double> ptsx = j[1]["ptsx"];
                    vector<double> ptsy = j[1]["ptsy"];
                    //std::cout << "ptsx.size() = " << ptsx.size() << std::endl;
                    
                    // State: position, heading, velocity
                    const double px = j[1]["x"];
                    const double py = j[1]["y"];
                    const double psi = j[1]["psi"];
                    const double speed = j[1]["speed"];
                    const double v = speed * prm.mph2mps; // converted from mph to m/s
                    //std::cout << "v = " << v << std::endl;
                    
                    // Current steering and throttle
                    const double steering_angle = j[1]["steering_angle"];
                    const double throttle = j[1]["throttle"];
                    
                    
                    // Convert to vehicle coordinate system
                    const double cos_psi = cos(psi);
                    const double sin_psi = sin(psi);
                    Eigen::VectorXd ptsx_car(ptsx.size());
                    Eigen::VectorXd ptsy_car(ptsy.size());
                    for(int i = prm.wpts_shift; i < ptsx.size(); i++) {
                        const double d_x = ptsx[i] - px;
                        const double d_y = ptsy[i] - py;
                        ptsx_car[i] = d_x * cos_psi + d_y * sin_psi;
                        ptsy_car[i] = - d_x * sin_psi + d_y * cos_psi;
                    }
                    
                    
                    // Fit polynomial to waypoints
                    auto coeffs = polyfit(ptsx_car, ptsy_car, 3);
                    
                    
                    // State in car coordinates
                    double px_car = 0.0;
                    double py_car = 0.0;
                    double psi_car = 0.0;
                    double v_car = v;
                    
                    // Calculate CTE and heading error
                    double cte = coeffs[0]; // i.e. polyeval(coeffs, 0) - 0
                    double epsi = - atan(coeffs[1]); // i.e. arctan(f'(0))
                    
                    // Account for latency
                    if (prm.latency > 0) {
                        px_car += v * cos_psi * prm.latency;
                        py_car += v * sin_psi * prm.latency;
                        psi_car += v * steering_angle / prm.Lf * prm.latency;
                        v_car += throttle * prm.latency;
                        cte += v * sin_psi * prm.latency;
                        epsi += v * steering_angle / prm.Lf * prm.latency;
                    }
                    
                    // Calculate CTE
                    //const double cte = polyeval(coeffs, px_car) - py_car;
                    
                    // Calculate heading error 
                    //const double epsi = psi_car - atan(coeffs[1] + 2 * coeffs[2] * px_car + 3 * coeffs[3] * pow(px_car, 2));
                    //const double epsi = psi_car - atan(coeffs[1]);
                    //std::cout << "cte = " << cte << " ; epsi = " << epsi << std::endl;
                    
                    
                    // State vector
                    Eigen::VectorXd state(6);
                    state << px_car, py_car, psi_car, v_car, cte, epsi;
                    
                    // Solver to find optimal actuators and return predicted trajectory
                    vector<vector<double>> mpc_output = mpc.Solve(state, coeffs);
                    
                    
                    // Steering and throttle values, averaged
                    double steer_value = accumulate(mpc_output[0].begin(), mpc_output[0].begin() + prm.n_smooth, 0.0) / (double) prm.n_smooth;
                    double throttle_value = accumulate(mpc_output[1].begin(), mpc_output[1].begin() + prm.n_smooth, 0.0) / (double) prm.n_smooth;
                    std::cout << "steer_value = " << steer_value << " ; throttle_value = " << throttle_value << std::endl;
                    
                    // Pass normalized steering and throttle values
                    json msgJson;
                    msgJson["steering_angle"] = - steer_value / deg2rad(25);
                    msgJson["throttle"] = throttle_value;
                    
                    
                    //Display MPC predicted trajectory
                    //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
                    // the points in the simulator are connected by a Green line
                    const vector<double> mpc_x_vals(mpc_output[2].begin() + prm.wpts_shift, mpc_output[2].end());
                    const vector<double> mpc_y_vals(mpc_output[3].begin() + prm.wpts_shift, mpc_output[3].end());
                    msgJson["mpc_x"] = mpc_x_vals;
                    msgJson["mpc_y"] = mpc_y_vals;
                    
                    //Display waypoints / reference line
                    //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
                    // the points in the simulator are connected by a Yellow line
                    vector<double> next_x_vals;
                    vector<double> next_y_vals;
                    for(int i = prm.wpts_shift; i < ptsx.size();i++){
                        next_x_vals.push_back(ptsx_car[i]);
                        next_y_vals.push_back(ptsy_car[i]);
                    }
                    msgJson["next_x"] = next_x_vals;
                    msgJson["next_y"] = next_y_vals;
                    
                    
                    auto msg = "42[\"steer\"," + msgJson.dump() + "]";
                    //std::cout << msg << std::endl;
                    
                    
                    // Latency: the purpose is to mimic real driving conditions where the car does actuate the commands instantly.
                    // Feel free to play around with this value but should be to drive around the track with 100ms latency.
                    // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE SUBMITTING.
                    if (prm.latency > 0) {
                        this_thread::sleep_for(chrono::milliseconds(int(prm.latency * 1000)));
                    }
                    ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
                }
            } else {
                // Manual driving
                std::string msg = "42[\"manual\",{}]";
                ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
            }
        }
    });
    
    // We don't need this since we're not using HTTP but if it's removed the
    // program doesn't compile :-(
    h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data, size_t, size_t) {
        const std::string s = "<h1>Hello world!</h1>";
        if (req.getUrl().valueLength == 1) {
            res->end(s.data(), s.length());
        } else {
        // i guess this should be done more gracefully?
        res->end(nullptr, 0);
        }
    });
    
    h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
        std::cout << "Connected!!!" << std::endl;
    });
    
    h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code, char *message, size_t length) {
        ws.close();
        std::cout << "Disconnected" << std::endl;
    });
    
    int port = 4567;
    if (h.listen(port)) {
        std::cout << "Listening to port " << port << std::endl;
    } else {
        std::cerr << "Failed to listen to port" << std::endl;
        return -1;
    }
    h.run();
}
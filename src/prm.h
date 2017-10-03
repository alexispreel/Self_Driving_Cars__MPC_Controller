#ifndef prm_H
#define prm_H

#include <cstddef>

class Prm {
    public:
        // Constructor
        Prm();
        
        // Destructor
        virtual ~Prm();
        
        // Timestep length and duration: to be reduced as reference speed rises
        const size_t N = 25;
        const double dt = 0.125;
        
        // Length from front to center of gravity
        const double Lf = 2.67;
        
        // Conversion from mph to m/s
        const double mph2mps = 0.447;
        
        // Reference CTE, heading error and velocity
        const double ref_cte = 0.0;
        const double ref_epsi = 0.0;
        const double ref_v = 45.0 * mph2mps; // converted from mph to m/s
        
        // Weights in total cost
        const double w_cte = 15.0;
        const double w_epsi = 2.0;
        const double w_v = 1.0;
        
        const double w_delta = 1.0;
        const double w_a = 1.0;
        const double w_delta_v = 300.0;
        
        const double w_ddelta = 20.0;
        const double w_da = 5.0;
        
        // Lower and upper limits for variables: delta, throttle factor and other variables
        const double max_delta = 0.436332; // 25 deg in rad
        const double max_thr = 1.0;
        const double min_thr = - 0.2;
        const double max_var = 1.0e19;
        
        // Latency factor (in seconds)
        const double latency = 0.1;
        
        // Waypoints to be skipped
        const int wpts_shift = 0;
        
        // Number of forward steps for smoothing (max N-1)
        const int n_smooth = 4;
};

#endif /* prm_H */
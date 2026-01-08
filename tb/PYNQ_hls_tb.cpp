#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "../src/PYNQ_hls_settings.h"

// Simulation Parameters (Matched with Matlab parameters)
struct SimParam {
    double R = 2.5;                // Stator resistance
    double Ld = 18e-3;             // D-axis inductance
    double Lq = 32e-3;             // Q-axis inductance
    double lambda_m = 0.30175;     // Flux linkage
    double p = 6.0;                // Pole pairs
    double J = 0.00782;            // Inertia
    double B = 0.02277;            // Friction coefficient
    double Ib = 10.0;              // Base current
    double Vb = 178.82;            // Base voltage
    double wb = 628.0;             // Base speed (rad/s)
};

// Motor Plant Dynamics
void motor_plant(double state[3], DataType u_pu[2], double dt, const SimParam& p) {
    // State conversion from Per-Unit to SI units
    double we = state[0] * p.wb;
    double iq = state[1] * p.Ib;
    double id = state[2] * p.Ib;
    
    // Voltage inputs conversion
    double vq = (double)u_pu[0] * p.Vb;
    double vd = (double)u_pu[1] * p.Vb;

    // Current derivatives (Electrical Equations)
    double diq = (vq - p.R * iq - we * p.Ld * id - we * p.lambda_m) / p.Lq;
    double did = (vd - p.R * id + we * p.Lq * iq) / p.Ld;
    
    // Torque and Speed derivative (Mechanical Equations)
    double Te = 1.5 * p.p * (p.lambda_m * iq + (p.Ld - p.Lq) * id * iq);
    double dwe = (Te - p.B * (we / p.p)) * (p.p / p.J);

    // State Integration (Euler Method) and normalization back to Per-Unit
    state[0] += (dwe / p.wb) * dt;
    state[1] += (diq / p.Ib) * dt;
    state[2] += (did / p.Ib) * dt;
}

int main() {
    SimParam sp;
    DataType inputArr[N2] = {0};
    DataType outputArr[N2] = {0};
    
    double state[3] = {0, 0, 0}; // Initial states: [omega_pu, iq_pu, id_pu]
    double t_end = 0.1;          // Simulation end time
    double dt = 1e-4;            // Step size
    double ramp_time = 0.0065;    // Speed ramp-up duration
    double final_rpm = 900.0;    // Target speed

    std::ofstream outfile_local("../../../../../csv/sim_results.csv");
    std::ofstream outfile("../../../../../../../../data_and_result/sim_result_csv/sim_results.csv");
    outfile_local << "Time,Ref_RPM,Act_RPM,Vq_pu,Vd_pu,Iq_pu,Id_pu\n";
    outfile << "Time,Ref_RPM,Act_RPM,Vq_pu,Vd_pu,Iq_pu,Id_pu\n";

    std::cout << "Starting SDA HLS Testbench..." << std::endl;

    for (double t = 0; t <= t_end; t += dt) {
        // 1. Reference Generator (Ramp function)
        double ref_rpm = (t < ramp_time) ? (final_rpm / ramp_time) * t : final_rpm;
        double omega_ref_pu = (ref_rpm * M_PI * sp.p / 60.0) / sp.wb;

        // 2. Prepare HLS Inputs
        // Map inputs according to PYNQ_hls.cpp definition
        inputArr[0] = (DataType)0;            // External Torque (Set to 0)
        inputArr[1] = (DataType)omega_ref_pu; // Speed Reference
        inputArr[2] = (DataType)state[0];     // Current Speed (Feedback)
        inputArr[3] = (DataType)state[1];     // Iq Current (Feedback)
        inputArr[4] = (DataType)state[2];     // Id Current (Feedback)

        // 2.1 Debug/Safety Check: Avoid division by zero in the Cubic Equation (q parameter)
        // If q is derived from inputArr[1..4], ensure they are not all zero at t=0
        if (t == 0 && inputArr[1] == (DataType)0) {
             // Optional: Initialize with a very small value to prevent CSim crash
             // inputArr[1] = (DataType)0.0001;
        }

        // 3. Execute HLS Core (SDA Algorithm)
        sda_main(inputArr, outputArr);

        // 4. Update Motor Plant (Numerical Simulation)
        motor_plant(state, outputArr, dt, sp);

        // 5. Data Logging and Console Output
        double act_rpm = state[0] * sp.wb * 60.0 / (M_PI * sp.p);

        // Log data every 10 steps to reduce file size
        if ((int)(t/dt) % 10 == 0) {
        	double iq = state[1];
        	double id = state[2];

        	outfile_local << t << "," << ref_rpm << "," << act_rpm << ","
        	        	            << (float)outputArr[0] << "," << (float)outputArr[1] << ","
        	        	            << iq << "," << id << "\n";

        	outfile << t << "," << ref_rpm << "," << act_rpm << ","
        	            << (float)outputArr[0] << "," << (float)outputArr[1] << ","
        	            << iq << "," << id << "\n";
            
            std::cout << "Time: " << std::fixed << std::setprecision(4) << t 
                      << " | Ref: " << (int)ref_rpm << " | Act: " << (int)act_rpm << std::endl;
        }
    }

    outfile_local.close();
    outfile.close();
    std::cout << "Simulation complete. Results saved to sim_results.csv" << std::endl;

    return 0;
}

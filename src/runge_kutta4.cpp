#include "runge_kutta4.h"
#include <iostream>


std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>> Runge_Kutta::solve(
		double u0, double t0, double eps, double h, const Parameters& parameters, long max_count)
{
	std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>> results;

	double dt = t0, t = t0, n_t, xi = 0;
	double u = u0;
	double prev_u;
	long int count = 0;

	while (std::abs(parameters.u_f - u) >= eps && count < max_count)
	{
		++count;

		results.emplace_back(std::make_pair(dt, xi), std::make_pair(u, t));

		double f_i = calculate_phi(parameters, u);
		auto velocity_function = [&f_i, &parameters](double t, double u) -> double
		{
			return f_i * (parameters.u_f - u);
		};

		double k1 = h * velocity_function(t, u);
		double k2 = h * velocity_function(t + h / 2, u + k1 / 2);
		double k3 = h * velocity_function(t + h / 2, u + k2 / 2);
		double k4 = h * velocity_function(t + h, u + k3);

		double new_u = u + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;

		dt += h;
		prev_u = u;
		u = new_u;
		if ((int)u == 999)
		{
			std::cout << "";
 		}
		n_t = get_r_t(parameters, dt);


		xi += 0.5 * (n_t - t) * (u + new_u);


		t = n_t;
		double delta_u = std::abs(new_u - prev_u);
		if (delta_u > 2)
		{
			h *= 0.8;
		}
		else if (delta_u < 1)
		{
			h *= 1.1;
		}
	}

	return results;
}



void Runge_Kutta::save_to_file(
		const std::string& filename,
		const std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>>& data,
		const Parameters& parameters)
{
	std::ofstream file(filename);
	if (!file.is_open())
	{
		throw std::runtime_error("Failed to open file for writing: " + filename);
	}

	file << "# Constants and Parameters\n\n";
	file << "# Final flow velocity (u_f): " << parameters.u_f << " (m/s) - The final velocity of the flow\n";
	file << "# Initial particle velocity (u_0): " << parameters.u_0 << " (m/s) - The initial velocity of the particle\n";
	file << "# Particle radius (r_is): " << parameters.r_is << " (m) - The radius of the particle\n";
	file << "# Nozzle radius (r_n): " << parameters.r_n << " (m) - The radius of the nozzle\n";
	file << "# Particle density (q): " << parameters.q << " (kg/m^3) - The density of the particle\n";
	file << "# Gas density (q_g): " << parameters.q_g << " (kg/m^3) - The density of the gas\n";
	file << "# Gas constant (R): " << parameters.R << " (J/(kg*K)) - The specific gas constant\n";
	file << "# Gamma: " << parameters.gamma;
	file << "# Molar mass of air (M): " << parameters.M << " (kg/mol) - The molar mass of air\n";
	file << "# Gas temperature (T): " << parameters.T << " (K) - The temperature of the gas\n";
	file << "# Reference temperature (T_0): " << parameters.T_0 << " (K) - The reference temperature for gas\n";
	file << "# Reference viscosity (eta_0): " << parameters.eta_0 << " (Pa*s) - The reference viscosity of the gas\n";
	file << "# Specific heat capacity of the gas (Cp): " << parameters.Cp << " (J/(kg*K)) - The specific heat capacity of the gas\n";
	file << "# Speed of sound (Sound speed): " << parameters.a << " (m/s) - The speed of sound in the gas\n";
	file << "# Dynamic viscosity (Viscosity): " << parameters.viscosity << " (Pa*s) - The dynamic viscosity of the gas\n";
	file << "\n\n\n";



	file << "\n\n\n# iter real_time (s) velocity (m/s) distance (m)\n";

	double max_dt = data[data.size() - 1].first.first;
	for (long i = 0; i < data.size(); ++i)
	{
		double xi = data[i].first.second;
		double u = data[i].second.first;
		double t = data[i].second.second;
		file << i << " " << t / max_dt << " " << u << " " << xi / max_dt << "\n";
	}

	file << "\n\n# Total time: " << parameters.time << " s";
	file << "\n# Total distance: " << data[data.size() - 1].first.second / max_dt << " m";
	std::cout << "\nTotal time: " << parameters.time << " s\n";

	file.close();
}










double get_r_t (const Parameters& p, const double& t)
{
	double r_t = (p.time * t);
	return r_t;
}

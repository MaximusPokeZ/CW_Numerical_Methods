#include "runge_kutta4.h"
#include <iostream>



std::vector<std::pair<double, double>> Runge_Kutta::solve(
		double u0, double t0, double eps, double h, const Parameters& parameters)
{
	std::vector<std::pair<double, double>> results;

	double t = t0;
	double u = u0;
	double prev_u = u0;

	while (std::abs(parameters.u_f - u) >= eps)
	{
		results.emplace_back(t, u);

		double phi = calculate_phi(parameters, u);
		auto velocity_function = [&phi, &parameters](double t, double u) -> double
		{
			return phi * (parameters.u_f - u);
		};

		double k1 = h * velocity_function(t, u);
		double k2 = h * velocity_function(t + h / 2, u + k1 / 2);
		double k3 = h * velocity_function(t + h / 2, u + k2 / 2);
		double k4 = h * velocity_function(t + h, u + k3);

		double new_u = u + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;

		double delta_u = std::abs(new_u - prev_u);

		if (delta_u > 2)
		{
			h *= 0.8;
		}
		else if (delta_u < 1)
		{
			h *= 1.2;
		}

		prev_u = u;
		u = new_u;
		t += h;
	}

	return results;
}


void Runge_Kutta::save_to_file(
		const std::string& filename,
		const std::vector<std::pair<double, double>>& data,
		const Parameters& parameters)
{
	std::ofstream file(filename);
	if (!file.is_open())
	{
		throw std::runtime_error("Failed to open file for writing: " + filename);
	}

	file << "# t real_time u xi\n";

	double xi = 0;
	double tau = (parameters.q * parameters.Cp * parameters.r_is * parameters.r_is * 6) / (3 * calculate_h(parameters));
	double real_time = tau * std::log(std::abs(parameters.u_f - parameters.u_0) / (0.01 * parameters.u_f));

	std::cout << "\nTotal time: " << real_time << "\n";

	for (size_t i = 0; i < data.size(); ++i)
	{
		double t = data[i].first;
		double u = data[i].second;
		double r_t = (real_time * t) / data[data.size()-1].first;
		double r_t_1 = (real_time * data[i + 1].first) / data[data.size()-1].first;
		file << t << " " << r_t << " " << u << " " << xi << "\n";

		if (i < data.size() - 1)
		{
			xi += 0.5 * (r_t_1 - r_t) * (data[i].second + data[i + 1].second);
		}
	}

	file.close();
}

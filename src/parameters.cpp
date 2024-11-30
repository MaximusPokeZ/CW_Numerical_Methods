#include "parameters.h"

Parameters Parameters::load_from_file(const std::string& filename)
{
	double u_f, u_0, r_is, r_n, q, q_g, R, gamma, M, T, T_0, eta_0, Cp;

	std::ifstream input_file(filename);
	if (!input_file.is_open())
	{
		throw std::runtime_error("Failed to open file: " + filename);
	}

	input_file >> u_f >> u_0 >> r_is >> r_n >> q >> q_g >> R >> gamma >> M >> T >> T_0 >> eta_0 >> Cp;
	input_file.close();

	// Проверка диапазонов значений
	if (u_f < 100.0 || u_f > 2000.0)
	{
		throw std::invalid_argument("u_f (Flow velocity) must be between 100 and 2000 m/s");
	}
	if (u_0 < 0.0 || u_0 > 1000.0)
	{
		throw std::invalid_argument("u_0 (Particle velocity) must be between 0 and 1000 m/s");
	}
	if (r_is < 1e-6 || r_is > 1e-4)
	{
		throw std::invalid_argument("r_is (Particle radius) must be between 1e-6 and 1e-4 meters");
	}

	return { u_f, u_0, r_is, r_n, q, q_g, R, gamma, M, T, T_0, eta_0, Cp };
}


double calculate_phi(const Parameters& params, const double& u)
{
	double mach = calculate_Mach_number(params.u_f, params.sound_speed);
	double reynolds = calculate_Reynolds_number(params.q_g, params.r_is, u, params.u_f, params.viscosity);
	double S = calculate_S(mach, params.gamma);
	double C_i = calculate_C_i(mach, reynolds, S);
	double St = calculate_stocks(params.q, params.r_is, params.sound_speed, params.viscosity, params.r_n);
	return C_i / St;
}

double calculate_stocks (double q, double r_is, double a_s, double eta, double r_n)
{
	return (2 * q * r_is * r_is * a_s) / (9 * eta * r_n);
}

double calculate_C_i (double M, double Re, double S)
{
	return (0.9 + 0.34 / (M * M) + 1.86 * std::sqrt(M / Re) * (2 + 2 / (S * S) + 1.058 / S - 1 / pow(S, 4))) /
			(1 + 1.86 * (std::sqrt(M / Re)));
}

double calculate_Mach_number (double u_f, double a_s)
{
	return u_f / a_s;
}

double calculate_Reynolds_number ( double q,  double r_is,  double u,  double u_f,  double eta)
{
	return (2 * q * std::abs(u - u_f) * r_is) / eta;
}

double calculate_S ( double M,  double gamma)
{
	return M * std::sqrt(gamma * 0.5);
}

double get_velocity(const Parameters& p)
{
	return (2.0 / 9.0) * (p.q_g * p.r_is * p.r_is * 9.81 / p.viscosity);
}

double get_Re(const Parameters& parameters)
{
	return parameters.q_g * 2 * parameters.r_is * get_velocity(parameters) / parameters.viscosity;
}

double calculate_h (const Parameters& param)
{
	return 0.023 * std::pow(get_Re(param), 0.8) *
		std::pow(1030 * param.viscosity / 0.023, 0.4) * 0.023 / (2 * param.r_is);
}
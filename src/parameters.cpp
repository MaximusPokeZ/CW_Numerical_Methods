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
	double mach = calculate_Mach_number(params.u_f, u, params.a); mach = (!(int)mach) ? params.u_f / params.a : mach;
	double reynolds = calculate_Reynolds_number(params.q_g, params.r_is, u, params.u_f, params.viscosity);
	double S = calculate_S(mach, params.gamma);
	double C_i = calculate_C_i(mach, reynolds, S, params.u_f, params.a);
	double St = calculate_stocks(params.q, params.r_is, params.a, params.viscosity, params.r_n);
	return C_i / St;
}

double calculate_stocks (double q, double r_is, double a_s, double eta, double r_n)
{
	double temp = (2 * q * r_is * r_is * a_s) / (9 * eta * r_n);
	return temp;
}

double calculate_C_i (double M, double Re, double S, double u, double a)
{
	double f2 = (0.9 + 0.34 / (M * M) + 1.86 * std::sqrt(M / Re) * (2 + 2 / (S * S) + 1.058 / S - 1 / pow(S, 4))) /
				(1 + 1.86 * (std::sqrt(M / Re)));

	double f1 = 24.0 / (Re + S * (4.33 + ((3.65 - 1.53) / (1 + 0.353)) * exp(-0.247 * (Re / S))))
			+ exp((-0.5 * M) / sqrt(Re)) * ((4.5 + 0.38 * (0.03 * Re + 0.48 * sqrt(Re))) / (1 + 0.03 * Re + 0.48 * sqrt(Re))
					+ 0.1 * M * M + 0.2 * pow(M, 8)) + (1 - exp(-M / Re)) * 0.6 * S;

	double f;
	if (M >= 1 && M <= 1.75)
	{
		f = f1 + (4.0 / 3.0) * (M - 1) * (f2 - f1);
	}
	else if (M < 1)
	{
		f = f1;
	}
	else
	{
		f = f2;
	}
	return f;
}

double calculate_Mach_number (double u_f, double u, double a_s)
{
	double tmp = std::abs(u_f - u) / a_s;
	return tmp;
}

double calculate_Reynolds_number ( double q,  double r_is,  double u,  double u_f,  double eta)
{
	double tmp = (2 * q * std::abs(u - u_f) * r_is) / eta;
	return tmp;
}

double calculate_S ( double M,  double gamma)
{
	double tmp = M * std::sqrt(gamma * 0.5);
	return tmp;
}

double get_velocity(const Parameters& p)
{
	double tmp = (2.0 / 9.0) * (p.q_g * p.r_is * p.r_is * 9.81 / p.viscosity);
	return tmp;
}

double get_Re(const Parameters& parameters)
{
	double tmp = parameters.q_g * 2 * parameters.r_is * get_velocity(parameters) / parameters.viscosity;
	return tmp;
}

double calculate_h (const Parameters& param)
{
	double tmp = 0.023 * std::pow(get_Re(param), 0.8) * std::pow(1030 * param.viscosity / 0.023, 0.4) * 0.023 / (2 * param.r_is);
	return tmp;
}

double r_t(const Parameters& parameters)
{
	double tau = (parameters.q * parameters.Cp * parameters.r_is * parameters.r_is * 6) / (3 * calculate_h(parameters));
	double real_time = tau * std::log(std::abs(parameters.u_f - parameters.u_0) / (0.01 * parameters.u_f));
	return real_time;
}

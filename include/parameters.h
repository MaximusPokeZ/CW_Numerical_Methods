#ifndef CW_NM_PARAMETERS_H
#define CW_NM_PARAMETERS_H

#include <string>
#include <fstream>

struct Parameters
{

	const double u_f;      // Скорость потока
	const double u_0;      // Начальная скорость частицы
	const double t_0 = 0.0;  // Начальное время
	const double r_is;     // Радиус частицы
	const double r_n;     // Радиус сопла
	const double q;        // Плотность частицы
	const double q_g;        // Плотность газа
	const double R;        // Газовая постоянная
	const double gamma;    // Отношение молярных теплоемкостей
	const double M;        // Молярная масса воздуха
	const double T;        // Температура газа
	const double T_0;      // Контрольная температура газа
	const double eta_0;    // Контрольная вязкость газа
	const double Cp; // Удельная теплоемкость вещества

	const double a = calculate_sound_speed(R, gamma, M, T); // Скорость звука
	const double viscosity = calculate_viscosity(eta_0, T_0, T);    // Вязкость

	static Parameters load_from_file(const std::string& filename);

	Parameters(double u_f, double u_0, double r_is, double r_n,
			   double q, double q_g, double R, double gamma, double M, double T, double T_0, double eta_0, double Cp)

			: u_f(u_f), u_0(u_0), r_is(r_is), r_n(r_n),
			  q(q), q_g(q_g), R(R), gamma(gamma), M(M),
			  T(T), T_0(T_0), eta_0(eta_0), Cp(Cp),
			  a(calculate_sound_speed(R, gamma, M, T)),
			  viscosity(calculate_viscosity(eta_0, T_0, T)) {}

private:

	static double calculate_sound_speed(double R, double gamma, double M, double T)
	{
		return std::sqrt((gamma * R * T) / M);
	}

	static double calculate_viscosity(double eta_0, double T_0, double T)
	{
		return eta_0 * std::pow(T / T_0, 1.5);
	}
};

double calculate_phi(const Parameters& params, const double& u);

double calculate_stocks (double q, double r_is, double a_s, double eta, double r_n);

double calculate_C_i (double M, double Re, double S, double u, double a, const Parameters& p);

double calculate_Mach_number (double u_f, double u, double a_s);

double calculate_Reynolds_number (double q, double r_is, double u, double u_f, double eta);

double calculate_S (double M, double gamma);


#endif //CW_NM_PARAMETERS_H

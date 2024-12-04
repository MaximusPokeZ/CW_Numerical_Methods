#include <iostream>
#include "runge_kutta4.h"


int main()
{
	try
	{
		Parameters params = Parameters::load_from_file("./resources/input.txt");

		double eps = 0.001, h = 0.001;

		auto result = Runge_Kutta::solve(params.u_0, params.t_0, eps, h, params);
		Runge_Kutta::save_to_file("result/result.txt", result, params);

	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << '\n';
		return 1;
	}


	return 0;
}

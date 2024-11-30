#ifndef CW_NM_RUNGE_KUTTA4_H
#define CW_NM_RUNGE_KUTTA4_H

#include <vector>
#include <functional>
#include <utility>
#include <fstream>
#include "parameters.h"

class Runge_Kutta
{

public:

	static std::vector<std::pair<double, double>> solve (double u0, double t0, double eps, double h, const Parameters& parameters);


	static void save_to_file(const std::string& filename, const std::vector<std::pair<double, double>>& data,
							 const Parameters& parameters);

};



#endif //CW_NM_RUNGE_KUTTA4_H

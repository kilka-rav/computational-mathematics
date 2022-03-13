#include "integrate.hpp"


int main() {
	Functions one(-1, 1, 0.005);
	std::cout << "Rectangle = " << one.rectangle() << " Error = " << fabs(one.rectangle() - M_PI ) << std::endl;
	std::cout << "Trapezoid = " << one.trapezoid() << " Error = " << fabs(one.trapezoid() - M_PI )<< std::endl;
	std::cout << "Sympson = " << one.sympson() << " Error = " << fabs(one.sympson() - M_PI ) << std::endl;
	std::cout << "Rect(eps 0.01) = " << one.rectangle(0.01) << " Error = " << fabs(one.rectangle(0.01) - M_PI ) << " Runge_modif = " << one.runge_modif << std::endl;
        std::cout << "Trap(eps 0.01) = " << one.trapezoid(0.01) << " Error = " << fabs(one.trapezoid(0.01) - M_PI ) << " Runge_modif = " << one.runge_modif << std::endl;
        std::cout << "Symp(eps 0.01) = " << one.sympson(0.01) << " Error = " << fabs(one.sympson(0.01) - M_PI ) << " Runge_modif = " << one.runge_modif << std::endl;  
	Functions two(0, 0.0005, 0.01, true);
	std::cout << "un_own = " << two.unown() << " Error = " << fabs(two.unown() - sqrt(M_PI)) << std::endl;
	return 0;
}

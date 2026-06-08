#include "Particule_flu.h"
#include <iostream>

Particule_flu::Particule_flu(double x, double y, double vx, double vy, double rho, double dx, double dy)
    : x(x), y(y), vx(vx), vy(vy), rho(rho), dx(dx), dy(dy) {}

void Particule_flu::afficher() const {
    std::cout << "Particule fluide" << std::endl;
    std::cout << "  Position : (" << x << ", " << y << ") m" << std::endl;
    std::cout << "  Vitesse  : (" << vx << ", " << vy << ") m/s" << std::endl;
    std::cout << "  Densité  : " << rho << " g/m^2" << std::endl;
    std::cout << "  Taille   : (" << dx << ", " << dy << ") m" << std::endl;
}
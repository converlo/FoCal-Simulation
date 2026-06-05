#include "Particule.h"
#include <iostream>

Particule::Particule(double x, double y, double z, double energie, std::string type)
    : x(x), y(y), z(z), energie(energie), type(type) {}

void Particule::afficher() const {
    std::cout << "Particule : " << type << std::endl;
    std::cout << "  Position : (" << x << ", " << y << ", " << z << ") cm" << std::endl;
    std::cout << "  Energie  : " << energie << " GeV" << std::endl;
}
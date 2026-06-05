#pragma once
#include <string>

class Particule {
public:
    double x, y, z;       // position (cm)
    double energie;        // énergie (GeV)
    std::string type;      // "photon", "electron", "pion"

    Particule(double x, double y, double z, double energie, std::string type);
    void afficher() const;
};
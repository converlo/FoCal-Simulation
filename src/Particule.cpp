#include "Particule_flu.h"
#include <iostream>

Particule_flu::Particule_flu(double x, double y, double masse, double h)
    : x(x), y(y), vx(0), vy(0), rho(0), P(0), masse(masse), h(h) {}
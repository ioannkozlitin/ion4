#ifndef RAIZER_H
#define RAIZER_H

#include <vector>
#include <string>

void calculatorMixRaizer(const std::vector<unsigned int> &Z, const std::vector<double> &x, double rCoeff, double lgRhoMin, double lgRhoMax, double lgRhoStep, double lgTMin, double lgTMax, double lgTStep, std::string filename);
void raizerTest(int argc, char *argv[]);

#endif // RAIZER_H

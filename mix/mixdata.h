#ifndef MIXDATA_H
#define MIXDATA_H

#include <vector>
#include "../saha/src/elements.h"
#include "../saha/src/sahasolver.h"

class MixData
{
public:
    MixData(const std::vector<unsigned int> &Z, const std::vector<double> &x, double rCoeff, bool correctV0, bool newVolumes, double T, double V, bool TVaeFlag);
    MixData(const std::vector<unsigned int> &Z, const std::vector<double> &x, double rCoeff, bool correctV0, bool newVolumes, double TeV, double Rho);
    std::vector<TElement> elements;
    std::vector<double> x;
    std::vector<std::vector<double>> xx;
    double T;
    double maxZ;
    double meanA;

    double SetTeVRho(double TeV, double Rho);
    double SetTVae(double T, double V);
    double GetFullV() const {return V;}
    double xe();
    double vfree();
    double p();
    double e();
    double s();
    double getZd(double xe);
    double getRd(double V);

    SahaPoint GetSahaPoint();

private:
    double V;
    void initZX(const std::vector<unsigned int> &Z, const std::vector<double> &x, double rCoeff, bool correctV0, bool newVolumes);
};

#endif

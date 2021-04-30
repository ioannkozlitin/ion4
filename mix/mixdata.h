#ifndef MIXDATA_H
#define MIXDATA_H

#include <vector>
#include "../saha/src/elements.h"
#include "../saha/src/sahasolver.h"

class MixData
{
public:
    MixData(const std::vector<unsigned int> &Z, const std::vector<double> &x, double rCoeff, bool correctV0, bool newVolumes, double T, double t, double V, bool TVaeFlag);
    MixData(const std::vector<unsigned int> &Z, const std::vector<double> &x, double rCoeff, bool correctV0, bool newVolumes, double TeV, double teV, double Rho);
    std::vector<TElement> elements;
    std::vector<double> x;
    std::vector<std::vector<double>> xx;
    double T;
    double t;
    double maxZ;
    double meanA;

    enum OutputFormat {ae, GPaKJ, TPaMJ};

    double SetTeVRho(double TeV, double teV, double Rho);
    double SetTVae(double T, double t, double V);
    double GetFullV() const {return V;}
    double xe();
    double vfree();
    double p();
    double pi();
    double pe();
    double e();
    double ei();
    double ee();
    double s();
    double si();
    double se();
    double getZd(double xe);
    double getRd(double V);

    SahaPoint GetSahaPoint();
    static void SetOutputFormat(OutputFormat value);

    double formatE(double x);
    double formatP(double x);
    double formatS(double x);
    double formatT(double x);

private:
    static OutputFormat outputFormat;
    double V;
    void initZX(const std::vector<unsigned int> &Z, const std::vector<double> &x, double rCoeff, bool correctV0, bool newVolumes);
};

#endif

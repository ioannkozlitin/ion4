#include "sahaleft.h"
#include <cmath>

double SahaLeft::operator()(const MixData &data, double vFree, double xe)
{
    const std::vector<std::vector<double>> &bb = b(data, vFree, xe);

    double s = 0;
    for(int i = 0; i < bb.size(); i++)
    {
        const std::vector<double> &bi = bb[i];
        double biMax = bi[0];
        for(int j = 1; j < bi.size(); j++) if(biMax < bi[j]) biMax = bi[j];

        double siA = 0, siB = 0;
        for(int j = 0; j < bi.size(); j++)
        {
            siA += j * exp(bi[j]-biMax);
            siB += exp(bi[j]-biMax);
        }

        s += siA/siB * data.x[i];
    }

    return s - xe;
}

double SahaLeft::GetVIon(const MixData &data, double vFree, double xe)
{
    const std::vector<std::vector<double>> &bb = b(data, vFree, xe);

    double vion = 0;
    for(int i = 0; i < bb.size(); i++)
    {
        double biMax = bb[i][0];
        for(int j = 1; j < bb[i].size(); j++) if(biMax < bb[i][j]) biMax = bb[i][j];

        double siBfull = 0;
        for(int j = 0; j < bb[i].size(); j++)
        {
            siBfull += exp(bb[i][j]-biMax);
        }

        double logX0 = -log(siBfull)-biMax;

        for(int j = 0; j < bb[i].size(); j++)
        {
            //double xx_ij = data.x[i] / siBfull * exp(bb[i][j]-biMax);
            double xx_ij = data.x[i] * exp(logX0 + bb[i][j]);
            vion += data.elements[i].v[j] * xx_ij;
        }
    }

    return vion;
}

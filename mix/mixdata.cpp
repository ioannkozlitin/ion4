#include "mixdata.h"
#include <stdio.h>
#include <cmath>
#include "../saha/src/atom_ed.h"
#include "../saha/src/mu_fi.h"

MixData::OutputFormat MixData::outputFormat = MixData::ae;

MixData::MixData(const std::vector<unsigned int> &Z, const std::vector<double> &x, double rCoeff, bool correctV0, bool newVolumes, double T, double t, double V, bool TVaeFlag)
{
    if(TVaeFlag)
    {
        this->T = T;
        this->V = V;
        this->t = t;
    }
    else
    {
        printf("TVaeFlag must be true!\n");
        return;
    }

    initZX(Z, x, rCoeff, correctV0, newVolumes);
}

MixData::MixData(const std::vector<unsigned int> &Z, const std::vector<double> &x, double rCoeff, bool correctV0, bool newVolumes, double TeV, double teV, double Rho)
{
    initZX(Z, x, rCoeff, correctV0, newVolumes);

    this->V = meanA * eRo / Rho;
    this->T = TeV / eFi;
    this->t = teV / eFi;
}

double MixData::SetTeVRho(double TeV, double teV, double Rho)
{
    V = meanA * eRo / Rho;
    T = TeV / eFi;
    t = teV / eFi;
    return 0;
}

double MixData::SetTVae(double T, double t, double V)
{
    this->T = T;
    this->t = t;
    this->V = V;
    return 0;
}

double MixData::xe()
{
    double r = 0;
    for(int i = 0; i < elements.size(); i++)
    {
        for(int j = 1; j <= elements[i].Z; j++)
        {
            r += xx[i][j] * j;
        }
    }

    return r;
}

double MixData::vfree()
{
    double vion = 0;
    for(int i = 0; i < elements.size(); i++)
    {
        for(int j = 0; j <= elements[i].Z; j++)
        {
            vion += xx[i][j] * elements[i].v[j];
        }
    }

    return std::max(V - vion, 0.0);
}

double MixData::p()
{
    //return 2*sqrt(2.0)/(3*M_PI*M_PI) * pow(t,2.5) * I15mu_d_t(t,vfree(),xe()) + T / vfree();
    return pi() + pe();
}

double MixData::pi()
{
    return T / vfree();
}

double MixData::pe()
{
    return 2 * sqrt(2.0)/(3*M_PI*M_PI) * pow(t,2.5) * I15mu_d_t(t,vfree(),xe());
}

double MixData::e()
{
    /*double vFree = vfree();
    double Ee = sqrt(2.0)/(M_PI*M_PI) * pow(t,2.5) * vFree * I15mu_d_t(t,vFree,xe());

    double Efi = 0;
    for(unsigned int i = 0; i < elements.size(); i++)
    {
        for(int j = 1; j <= elements[i].Z; j++)
        {
            Efi += elements[i].cumFi[j-1] * xx[i][j];
        }
    }
    return 1.5 * T + Ee + Efi;*/
    return ei() + ee();
}

double MixData::ei()
{
    double Efi = 0;
    for(unsigned int i = 0; i < elements.size(); i++)
    {
        for(int j = 1; j <= elements[i].Z; j++)
        {
            Efi += elements[i].cumFi[j-1] * xx[i][j];
        }
    }

    return 1.5 * T + Efi;
}

double MixData::ee()
{
    double vFree = vfree();
    return sqrt(2.0)/(M_PI*M_PI) * pow(t,2.5) * vFree * I15mu_d_t(t,vFree,xe());
}

double MixData::s()
{
    return si() + se();
}

double MixData::si()
{
    double vFree = vfree();
    double Si = 2.5;
    for(unsigned int i = 0; i < elements.size(); i++)
    {
        const double M = 1822.887 * elements[i].A;
        for (unsigned int j = 0; j <= elements[i].Z; j++)
        {
            if (xx[i][j] > 0)
            {
                Si += xx[i][j] * (elements[i].logG[j] + log(vFree) - log(xx[i][j]) + 1.5 * log(M * T / 2.0 / M_PI));
            }
        }
    }

    return Si;
}

double MixData::se()
{
    double vFree = vfree();
    double _xe = xe();

    double Se = 0;
    if(_xe > 0)
    {
        Se = sqrt(2.0) / (M_PI*M_PI) * pow(t, 1.5) * vFree * (5.0 / 3.0 * I15mu_d_t(t, vFree, _xe) - mu(t,vFree,_xe) / t * I05mu_d_t(t, vFree, _xe));
    }

    return Se;
}

double MixData::getZd(double xe)
{
    double z2 = xe;
    double allIonsX = 0;

    for(int i = 0; i < elements.size(); i++)
    {
        for(int j = 1; j <= elements[i].Z; j++)
        {
            z2 += xx[i][j] * j * j;
            allIonsX += xx[i][j];
        }
    }

    double z2final = (allIonsX > 1e-16) ? z2 / allIonsX : 1;

    return sqrt(std::max(z2final, 0.0));
}

double MixData::getRd(double V)
{
    double R3 = 3 * V / (4 * M_PI);
    double allIonsX = 0;

    for(int i = 0; i < elements.size(); i++)
    {
        for(int j = 1; j <= elements[i].Z; j++)
        {
            allIonsX += xx[i][j];
        }
    }

    return pow(std::max(R3 / allIonsX, 0.0), 1 / 3.0);
}

SahaPoint MixData::GetSahaPoint()
{
    SahaPoint result;

    double vFree = vfree();
    double _si = si(), _se = se(), _pi = pi(), _pe = pe(), _ei = ei(), _ee = ee();

    result.Z = maxZ;
    result.T = T;
    result.t = t;
    result.V = V;
    result.E = formatE(_ei + _ee);
    result.Ei = formatE(_ei);
    result.Ee = formatE(_ee);
    result.P = formatP(_pi + _pe);
    result.Pi = formatP(_pi);
    result.Pe = formatP(_pe);
    result.S = formatS(_si + _se);
    result.Si = formatS(_si);
    result.Se = formatS(_se);
    result.Xe = xe();
    result.M = formatT(mu(t, vFree, result.Xe));
    result.F = formatE(result.E - t * _se - T * _si);
    result.vFactor = vFree / V;
    result.IMu = I05mu_d_t(t,vFree,result.Xe);
    result.K = pow(result.IMu,1.5);
    result.vError = 0;
    result.auxIt = 0;
    result.zd = getZd(result.Xe);
    result.rd = getRd(V);

    return result;
}

void MixData::SetOutputFormat(MixData::OutputFormat value)
{
    outputFormat = value;
}

double MixData::formatE(double x)
{
    switch(outputFormat)
    {
        case GPaKJ:
            return eE / meanA * x;
        case TPaMJ:
            return eE / meanA * x * 1e-3;
        default:
            return x;
    }
}

double MixData::formatP(double x)
{
    switch(outputFormat)
    {
        case GPaKJ:
            return x * eP;
        case TPaMJ:
            return 1e-3 * x * eP;
        default:
            return x;
    }
}

double MixData::formatS(double x)
{
    switch(outputFormat)
    {
        case GPaKJ:
            return eE/eFi * x / meanA;
        case TPaMJ:
            return eE/eFi * x / meanA;
    default:
        return x;
    }
}

double MixData::formatT(double x)
{
    switch(outputFormat)
    {
        case GPaKJ:
            return eFi * x;
        case TPaMJ:
            return eFi * x * 1e-3;
        default:
            return x;
    }
}

void MixData::initZX(const std::vector<unsigned int> &Z, const std::vector<double> &x, double rCoeff, bool correctV0, bool newVolumes)
{
    if(x.size() != Z.size())
    {
        printf("Error: size(Z) != size(x)\n");
        return;
    }

    for(auto &zi : Z)
    {
        elements.push_back(TElement(zi, rCoeff, correctV0, newVolumes));
        xx.push_back(std::vector<double>(zi+1));
    }

    this->x = x;

    maxZ = 0;for(int i = 0; i < x.size(); i++) maxZ += Z[i] * x[i];
    meanA = 0;for(int i = 0; i < Z.size(); i++) meanA += x[i] * elements::GetA(Z[i]);

    //Soft Ions
    //for(TElement &elem : elements) elem.softIon(V);
}

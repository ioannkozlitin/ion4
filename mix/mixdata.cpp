#include "mixdata.h"
#include <stdio.h>
#include <cmath>
#include "../saha/src/atom_ed.h"
#include "../saha/src/mu_fi.h"

MixData::MixData(const std::vector<unsigned int> &Z, const std::vector<double> &x, double rCoeff, bool correctV0, bool newVolumes, double T, double V, bool TVaeFlag)
{
    if(TVaeFlag)
    {
        this->T = T;
        this->V = V;
    }
    else
    {
        printf("TVaeFlag must be true!\n");
        return;
    }

    initZX(Z, x, rCoeff, correctV0, newVolumes);
}

MixData::MixData(const std::vector<unsigned int> &Z, const std::vector<double> &x, double rCoeff, bool correctV0, bool newVolumes, double TeV, double Rho)
{
    initZX(Z, x, rCoeff, correctV0, newVolumes);

    this->V = meanA * eRo / Rho;
    this->T = TeV / eFi;
}

double MixData::SetTeVRho(double TeV, double Rho)
{
    V = meanA * eRo / Rho;
    T = TeV / eFi;
}

double MixData::SetTVae(double T, double V)
{
    this->T = T;
    this->V = V;
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
    return 2*sqrt(2.0)/(3*M_PI*M_PI) * pow(T,2.5) * I15mu_d_t(T,vfree(),xe()) + T / vfree();
}

double MixData::e()
{
    double vFree = vfree();
    double Ee = sqrt(2.0)/(M_PI*M_PI) * pow(T,2.5) * vFree * I15mu_d_t(T,vFree,xe());

    double Efi = 0;
    for(unsigned int i = 0; i < elements.size(); i++)
    {
        for(int j = 1; j <= elements[i].Z; j++)
        {
            Efi += elements[i].cumFi[j-1] * xx[i][j];
        }
    }
    return 1.5 * T + Ee  + Efi;
}

double MixData::s()
{
    double vFree = vfree();
    double _xe = xe();

    double Se = 0;
    if(_xe > 0)
    {
        Se = sqrt(2.0) / (M_PI*M_PI) * pow(T, 1.5) * vFree * (5.0 / 3.0 * I15mu_d_t(T, vFree, _xe) - mu(T,vFree,_xe) / T * I05mu_d_t(T, vFree, _xe));
    }

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
    return Si + Se;
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

    result.Z = maxZ;
    result.T = T;
    result.V = V;
    result.E = e();
    result.P = p();
    result.S = s();
    result.Xe = xe();
    result.M = mu(T, vFree, result.Xe);
    result.F = result.E - T * result.S;
    result.vFactor = vFree / V;
    result.IMu = I05mu_d_t(T,vFree,result.Xe);
    result.K = pow(result.IMu,1.5);
    result.vError = 0;
    result.auxIt = 0;
    result.zd = getZd(result.Xe);
    result.rd = getRd(V);

    return result;
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

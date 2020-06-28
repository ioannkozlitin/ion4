#ifndef FINDROOT_H
#define FINDROOT_H

#include <functional>

class FindRoot
{
public:
    double operator()(double a, double b, const std::function<double(double)> &F, double eps, double T, double V);
    double operator()(double logA, double logB, double x0, const std::function<double(double)> &F, double eps, double T, double V);

private:
    double core(double a, double b, double fa, double fb, const std::function<double(double)> &F, double eps, double T, double V);
    double xroot(double x1, double y1, double x2, double y2,double x3, double y3);
    double invxroot(double x1, double y1, double x2, double y2, double x3, double y3);
    double chord(double x1, double y1, double x2, double y2);
};

#endif // FINDROOT_H

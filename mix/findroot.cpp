#include "findroot.h"
#include <cmath>

double FindRoot::operator ()(double logA, double logB, const std::function<double (double)> &F, double eps, double T, double V)
{
    double a = logA, b = logB, c, cnew;
    double fa, fb, fc;
    double root = 0;

    int cc = 0, cc2 = 0;

    cnew = b + 1;
    fa = F(exp(a));
    fb = F(exp(b));
    if(((fa >= 0) && (fb >= 0)) || ((fa <= 0) && (fb <= 0)))
    {
        if (fabs(fa) < fabs(fb)) root = exp(a);
        else root = exp(b);
    }
    else
    {
        do
        {
            if((cnew > a) && (cnew < b))
            {
                c = cnew;
                cc2++;
            }
            else
            {
                c = 0.5*(a + b);
            }

            fc = F(exp(c));

            if(cc2 < 16) cnew = xroot(a,fa,b,fb,c,fc);

            if (((fa <= 0) && (fc >= 0)) || ((fa >= 0) && (fc <= 0)))
            {
                b = c; fb = fc;
            }
            else if (((fb <= 0) && (fc >= 0)) || ((fb >= 0) && (fc <= 0)))
            {
                a = c; fa = fc;
            }
            else
            {
                printf("find root error: ln(a) = %g ln(b) = %g ln(c) = %g fa = %g fb = %g fc = %g T = %g V = %g\n", a, b, c, fa, fb, fc, T, V);
                return 0;
            }

            cc++;
        }
        while ((b - a > eps) && (cc < 60) && (fabs(fa) > eps) && (fabs(fb) > eps));
        root = exp(0.5*(a + b));
    }

    //printf("<%d>",cc);
    double Froot = F(root);

    if(fabs(Froot) > fa) {root = exp(a);Froot = fa;};
    if(fabs(Froot) > fb) {root = exp(b);Froot = fb;};

    double root2 = exp(chord(a,fa,b,fb));
    if(fabs(F(root2)) < fabs(Froot))
    {
        return root2;
    }
    else return root;
}

double FindRoot::xroot(double x1, double y1, double x2, double y2,double x3, double y3)
{
    double a=((y3-y2)/(x3-x2)-(y2-y1)/(x2-x1))/(x3-x1);
    double b=a*(x3-x2)+(y3-y2)/(x3-x2);
    double c=y3;
    double d = b*b-4*a*c;
    double r1 = (-b + sqrt(d))/(2*a);
    double r2 = (-b - sqrt(d))/(2*a);

    if(fabs(r1)<fabs(r2)) return r1+x3;
    else return r2+x3;
}

double FindRoot::chord(double x1, double y1, double x2, double y2)
{
    return x1 - (x2-x1) / (y2-y1) * y1;
}

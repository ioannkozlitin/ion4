#include "Constants.h"
#include <cmath>

// Вычисление функций Ферми-Дирака и химического потенциала
// Формулы взяты из диссертации Кузьминой стр. 41, 52-55 и препринта "Ионизационное равновесие с учетом вырождения"

/// --- I_{1/2}^{-1}(log(z)) --- ///

// Набор коэффициентов для вычисления invI05
const std::vector<fp> a =
{
    0.23960888E+0,
    0.25551970E-1,
    0.16001700E-2,
    0.62234299E-4,
    0.19789713E-5
};

const std::vector<fp> b =
{
    -0.26354443E-1
};

const std::vector<fp> alpha =
{
    0.72771980E+3,
    -0.11754360E+4,
    0.32688108E+3,
    -0.22333565E+1
};

const std::vector<fp> beta =
{
    0.33069903E+3,
    0.10575554E+1
};

inline fp F1(fp y)
{
    fp exp_y = std::exp(y);
    return std::log(4.0 / (3.0 * std::sqrt(M_PIl))) + y + std::log((1.0 + exp_y * (a[0] + exp_y * (a[1] + exp_y * (a[2] + exp_y * (a[3] + a[4] * exp_y))))) / (1.0 + b[0] * exp_y));
}


inline fp F2(fp y)
{
    fp exp_y = std::exp(y);
    return std::pow((alpha[0] + exp_y * (alpha[1] + exp_y * (alpha[2] + exp_y * (alpha[3] + exp_y)))) / (beta[0] + exp_y * (beta[1] + exp_y)), 0.25);
}

fp Constants::invI05_log(fp log_z) // Обратная функция к функции Ферми-Дирака I_{1/2}(log(z))
{
    fp result;

    if (log_z < std::log((3.02 * 3.02) / 1.5))
    {
        result = F1(std::log(1.5) + log_z);
    }
    else
    {
        result = F2(4.0 / 3.0 * (std::log(1.5) + log_z));
    }

    return result;
}

/// --- I_{3/2}(I_{1/2}^{-1}(z)) --- ///

// Набор коэффициентов для вычисления I15_invI05
const std::vector<fp> aa =
{
    1,
    0.50625059,
    0.93656560E-1,
    0.74482877E-2,
    0.21519676E-3
};

const std::vector<fp> bbb =
{
    1,
    0.10730909,
    0.33951152E-2
};

const std::vector<fp> aal =
{
    0.58432930E+5,
    0.21851397E+5,
    0.18917921E+4,
    0.57806497E+2
};

const std::vector<fp> bbt =
{
    0.11348583E+4,
    0.41356686E+2
};

fp I15low(fp y2) // Ветвь при малых значениях аргумента
{
    fp hs = aa[0] + y2 * (aa[1] + y2 * (aa[2] + y2 * (aa[3] + y2 * aa[4])));
    fp zn = bbb[0] + y2 * (bbb[1] + y2 * bbb[2]);
    return y2 * std::pow(hs / zn, 1.0 / 3.0);
}

fp I15hi(fp y) // Ветвь при больших значениях аргумента
{
    fp y8d3 = std::pow(y, 8.0 / 3.0);
    fp hs = aal[0] + y8d3 * (aal[1] + y8d3 * (aal[2] + y8d3 * (aal[3] + y8d3)));
    fp zn = bbt[0] + y8d3 * (bbt[1] + y8d3);
    return 0.4 * y * y * std::pow(hs / zn, 0.25);
}

fp Constants::I15_invI05(fp z) // Функция Ферми-Дирака I_{3/2}(I_{1/2}^{-1}(z))
{
    fp y = std::sqrt(1.5 * z);
    if (y < 3.75)
    {
        return I15low(1.5 * z);
    }
    else
    {
        return I15hi(y);
    }
}

fp Constants::mu_log(fp T, fp log_V, fp log_xe) // mu(log(xe))
{
    return T * invI05_log(std::log(M_PIl * M_PIl / (std::sqrt(2.0) * std::pow(T, 1.5))) - log_V + log_xe);
}

fp Constants::I05_mu_T(fp T, fp V, fp xe) // Величина I_{1/2}(mu/T)
{
    return M_PIl * M_PIl * xe / (std::sqrt(2.0) * std::pow(T, 1.5) * V);
}

fp Constants::I15_mu_T(fp T, fp V, fp xe) // Величина I_{3/2}(mu/T)
{
    return I15_invI05(I05_mu_T(T, V, xe));
}

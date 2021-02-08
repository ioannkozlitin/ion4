#include "calculator.h"
#include <cmath>
#include <fstream>
#include <QDebug>

Calculator::Calculator()
{

}

Calculator::~Calculator()
{
    _isRun = false;
    exit(-1);
    wait();
}

void Calculator::start(const std::vector<unsigned int> &Z, const std::vector<double> &x, double rCoeff, double lgRhoMin, double lgRhoMax, double lgRhoStep, double lgTMin, double lgTMax, double lgTStep, const std::string &filePath)
{
    _Z = Z;
    _x = x;
    _rCoeff = rCoeff;
    _lgRhoMin = lgRhoMin;
    _lgRhoMax = lgRhoMax;
    _lgRhoStep = lgRhoStep;
    _lgTMin = lgTMin;
    _lgTMax = lgTMax;
    _lgTStep = lgTStep;
    _filePath = filePath;

    _isRun = true;
    QThread::start();
}

void Calculator::stop()
{
    _isRun = false;
}

void Calculator::outputArray(std::ostream &os, const std::string dataName, const std::vector<double> &data)
{
    os << dataName << " = [";
    for(std::vector<double>::const_iterator it = data.begin(); it!= data.end(); ++it) {os << *it << " ";};
    os << "];" << std::endl;
}

void Calculator::outputTable(std::ostream &os, const std::string &tableName, const std::vector<std::vector<SahaPoint>> table, std::function<double(const SahaPoint&)> accessor)
{
    os << tableName << " = [" << std::endl;
    for (size_t i = 0; i < table.size(); ++i)
    {
        const std::vector<SahaPoint>& line = table[i];
        for (size_t j = 0; j < line.size(); ++j)
        {
            os << accessor(line[j]) << " ";
        }
        os << std::endl;
    }
    os << "];" << std::endl;
}

void Calculator::run()
{
    double iterationsNum = ((_lgTMax - _lgTMin) / _lgTStep + 1.0) * ((_lgRhoMax - _lgRhoMin) / _lgRhoStep + 1.0);
    int iterationNum = 0;

    SahaMixSolver mixSolver;

    std::vector<double> lgTPhys;
    std::vector<double> lgVa;
    std::vector<double> _lgRho;

    std::vector<std::vector<SahaPoint>> fullTable;

    MixData md(_Z, _x, _rCoeff, true, true, 1, 1);

    for (double lgT = _lgTMax; lgT > _lgTMin - _lgTStep / 2.0; lgT -= _lgTStep)
    {
        lgTPhys.push_back(lgT);
        std::vector<SahaPoint> fullLine;

        bool fillFlag = _lgRho.empty();

        for (double lgRho = _lgRhoMax; lgRho > _lgRhoMin - _lgRhoStep / 2.0; lgRho -= _lgRhoStep)
        {
            if (!_isRun)
            {
                exit(-1);
                return;
            }

            md.SetTeVRho(pow(10, lgT), pow(10, lgRho));
            mixSolver.GetFullIonizationInfo(md);
            fullLine.push_back(md.GetSahaPoint());

            if(fillFlag)
            {
                _lgRho.push_back(lgRho);
                lgVa.push_back(log10(md.GetFullV()));
            }

            iterationNum++;
            emit setProgress(iterationNum / iterationsNum * 100.0);
        }

        fullTable.push_back(fullLine);
    }

    std::fstream f(_filePath.c_str(), std::fstream::out);
    f << std::scientific;

    f << "Z=[";for(auto &z : _Z) f << z << " ";f << "];" << std::endl;
    f << "x=[";for(auto &x : _x) f << x << " ";f << "];" << std::endl;

    outputArray(f, "lgT", lgTPhys);
    outputArray(f, "lgV", lgVa);
    outputArray(f, "lgRho", _lgRho);
    outputTable(f, "xe", fullTable, std::mem_fn(&SahaPoint::Xe));
    outputTable(f, "P", fullTable, std::mem_fn(&SahaPoint::P));
    outputTable(f, "E", fullTable, std::mem_fn(&SahaPoint::E));
    outputTable(f, "S", fullTable, std::mem_fn(&SahaPoint::S));

    f.close();
    exit(0);
}

void Calculator::exit(int retcode)
{
    emit finished(retcode == 0);
    QThread::exit(retcode);
}

#ifndef CALCULATOR_H
#define CALCULATOR_H

#include <QThread>
#include "mix/sahamixsolver.h"

class Calculator : public QThread
{
    Q_OBJECT

public:
    Calculator();
    ~Calculator();

public slots:
    void start(const std::vector<unsigned int> &Z, const std::vector<double> &x, double rCoeff, double lgRhoMin, double lgRhoMax, double lgRhoStep, double lgTMin, double lgTMax, double lgTStep, const std::string &filePath);
    void stop();

signals:
    void setProgress(int percentage);
    void finished(bool isSuccess);

private:
    void outputArray(std::ostream &os, const std::string dataName, const std::vector<double> &data);
    void outputTable(std::ostream &os, const std::string &tableName, const std::vector<std::vector<SahaPoint>> table, std::function<double(const SahaPoint&)> accessor);
    void run();
    void exit(int retcode);

    bool _isRun = false;

    std::vector<unsigned int> _Z;
    std::vector<double> _x;
    double _rCoeff;
    double _lgRhoMin;
    double _lgRhoMax;
    double _lgRhoStep;
    double _lgTMin;
    double _lgTMax;
    double _lgTStep;
    std::string _filePath;
};

#endif // CALCULATOR_H
